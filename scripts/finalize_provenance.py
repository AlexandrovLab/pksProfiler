#!/usr/bin/env python3
"""Finish the run record, and write the run report a person actually reads.

F16 asks for: the production commit, parameters, enabled modules, reference digests,
sample/data identifiers, tool versions, trace, and output checksums.

The workflow writes the launch half (commit, parameters, reference digests,
environment) before any task starts, so a run that dies still says what it was
attempting. This adds what only exists afterwards -- outcome, which steps ran, the
software conda installed, the tool versions each task actually ran (and whether they
match their pins), the samples, the output checksums -- and writes both forms:

  runs/provenance/<sessionId>.json   authoritative, complete, machine-readable
  runs/RUN_REPORT.txt                the same run in reading order, appended per session

One record per session, because with `-resume` a results tree is usually the product of
several sessions, and a record describing only the last one would misdescribe it.

Invoked by the workflow.onComplete handler in nextflow.config; safe to run by hand.
"""
import argparse
import csv
import datetime
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

MAX_HASH_BYTES = 512 * 1024 * 1024

# What a reader wants to see a version for, out of the few hundred packages conda pulls
# in. Everything else stays in the JSON.
NOTABLE = (
    "python", "r-base", "pandas", "numpy", "scipy", "matplotlib-base", "pysam",
    "samtools", "bcftools", "bedtools", "minimap2", "bowtie2", "fastp", "subread",
    "hmmer", "krakenuniq", "kraken2", "bracken", "megahit", "spades", "metabat2",
    "checkm2", "gtdbtk", "genomad", "mlst", "prokka", "diamond", "seqkit", "htslib",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default=None, type=Path,
                        help="results root; RUN_REPORT.txt for the latest run is written here")
    parser.add_argument("--runs-dir", required=True, type=Path)
    parser.add_argument("--cohort-dir", required=True, type=Path)
    parser.add_argument("--session", required=True)
    parser.add_argument("--success", default="false")
    parser.add_argument("--exit-status", default="")
    parser.add_argument("--duration", default="")
    parser.add_argument("--error-message", default="")
    parser.add_argument("--trace-file", default=None, type=Path)
    parser.add_argument("--conda-cache-dir", default=None, type=Path)
    # From workflow.stats: correct even when the trace file is absent or still being
    # flushed when this runs.
    parser.add_argument("--succeeded", default=None)
    parser.add_argument("--cached", default=None)
    parser.add_argument("--failed", default=None)
    parser.add_argument("--ignored", default=None)
    # F16: strict exits non-zero when a task ran a tool that is not its pinned build from
    # its own conda env, or left no record of what it ran; warn only reports it.
    parser.add_argument("--tool-version-check", default="strict", choices=("strict", "warn"))
    parser.add_argument("--repo-dir", default=Path(__file__).resolve().parents[1], type=Path,
                        help="pipeline source, to know which processes record tool versions")
    return parser.parse_args()


# ---------------------------------------------------------------- gathering

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cohort_checksums(cohort_dir):
    """Cohort tables only: per-sample output is one file per sample per stage."""
    if not cohort_dir.is_dir():
        return []
    entries = []
    for path in sorted(p for p in cohort_dir.rglob("*") if p.is_file()):
        size = path.stat().st_size
        entries.append({
            "path": str(path.relative_to(cohort_dir)),
            "bytes": size,
            "sha256": sha256(path) if size <= MAX_HASH_BYTES
                      else f"skipped: larger than {MAX_HASH_BYTES} bytes",
        })
    return entries


# A cohort-scale trace is ~330k rows and ~45 MB. It is read once, streamed, and only
# the per-process tallies and the failures are kept -- loading it into a list of dicts
# cost ~450 MB on the head node at the end of every run, for nothing.
MAX_FAILURES_SHOWN = 20


def scan_trace(trace_file, work_dir):
    """
    One pass over Nextflow's trace: what ran, grouped by process, and what failed.

    Returns (steps, failures). `name` is "processName (tag)", where the tag is the
    sample and so varies per task; the process name is what groups.
    """
    if not trace_file or not Path(trace_file).is_file():
        return [], []

    by_process = defaultdict(Counter)
    failures = []
    with Path(trace_file).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = (row.get("name") or row.get("process") or "").strip()
            process = name.split(" (")[0].strip()
            if not process:
                continue
            status = (row.get("status") or "UNKNOWN").strip()
            by_process[process][status] += 1

            if status in ("FAILED", "ABORTED"):
                task_hash = (row.get("hash") or "").strip()
                failures.append({
                    "task": name,
                    "status": status,
                    "exit": (row.get("exit") or "").strip(),
                    "hash": task_hash,
                    "stderr": f"{work_dir}/{task_hash}*/.command.err"
                              if work_dir and task_hash else None,
                })

    steps = [{"process": process,
              "tasks": sum(statuses.values()),
              "statuses": dict(sorted(statuses.items()))}
             for process, statuses in sorted(by_process.items())]
    return steps, failures


def conda_inventory(cache_dir):
    """
    Versions conda actually installed, read from each environment's conda-meta.

    The filenames there are `<name>-<version>-<build>.json`, so no conda call and no
    import is needed -- which matters, since this runs after the pipeline, outside any
    of its environments.
    """
    if not cache_dir or not Path(cache_dir).is_dir():
        return {}
    inventory = {}
    for env_dir in sorted(p for p in Path(cache_dir).iterdir() if (p / "conda-meta").is_dir()):
        packages = {}
        for meta in sorted((env_dir / "conda-meta").glob("*.json")):
            parts = meta.stem.rsplit("-", 2)
            if len(parts) == 3:
                packages[parts[0]] = parts[1]
        if packages:
            inventory[env_dir.name] = packages
    return inventory


# ---------------------------------------------------------------- tool versions (F16)
#
# conda_inventory() above is what conda INSTALLED. What RAN is recorded by the task
# itself: scripts/record_tool_versions.sh, rendered into every process that runs an
# external tool by lib/ToolVersions.groovy, asks each binary on the task's PATH -- after
# conda activation -- for its version and location, and writes the answer to the work
# dir's .tool_versions.tsv. The two differ exactly when something earlier on PATH
# shadows the environment: on TSCC a user minimap2 2.24 ran every host-depletion pass
# while this report said 2.31 (29 Sep).

TOOL_RECORD = ".tool_versions.tsv"
# Statuses the task itself assigns that mean "not the pinned build from its own env".
BLOCKING_STATUSES = ("OUTSIDE_ENV", "VERSION_MISMATCH", "NOT_FOUND", "NO_ENV")
# A task that finished without a record, although its process records one.
MISSING = "MISSING"
# The binary that ran reports a version the environment's conda-meta does not hold.
INSTALLED_MISMATCH = "INSTALLED_MISMATCH"
RECORDED_STATUSES = ("COMPLETED", "CACHED")
PROCESS_RE = re.compile(r"^\s*process\s+(\w+)\s*\{", re.M)
ALIAS_RE = re.compile(r"\b(\w+)\s+as\s+(\w+)\b")


def version_matches(observed, expected):
    """True when `observed` is `expected` or a more specific release of it (x.y.z of x.y)."""
    return observed == expected or observed.startswith(expected + ".")


# Tool-bump (Oct 2026): conda builds whose binary prints another version than the package
# (bracken 3.1p1 prints 3.0.1). The same file scripts/record_tool_versions.sh reads.
ALIASES_FILE = Path(__file__).resolve().parent / "tool_version_aliases.tsv"


def version_aliases(path=ALIASES_FILE):
    """{(package, conda_version): reported_version} from tool_version_aliases.tsv."""
    aliases = {}
    if not Path(path).is_file():
        return aliases
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) >= 3 and fields[0] and fields[1] and fields[2]:
            aliases[(fields[0], fields[1])] = fields[2]
    return aliases


def version_satisfies(observed, expected, package, aliases):
    """version_matches, or the known self-report of `package` at conda version `expected`."""
    if version_matches(observed, expected):
        return True
    alias = aliases.get((package, expected))
    return bool(alias) and version_matches(observed, alias)


def conda_pins(env_file):
    """{package: version} for the exact `- package=version` pins in a conda env yml."""
    pins = {}
    path = Path(env_file) if env_file else None
    if not path or path.suffix not in (".yml", ".yaml") or not path.is_file():
        return pins
    for line in path.read_text().splitlines():
        spec = line.split("#", 1)[0].strip()
        if not spec.startswith("-"):
            continue
        spec = spec[1:].replace(" ", "").split("::")[-1]
        match = re.match(r"^([A-Za-z0-9_.+-]+?)(==?)([^=<>!~,|]+)", spec)
        if match:
            pins[match.group(1)] = match.group(3).rstrip("*").rstrip(".")
    return pins


def recording_processes(repo_dir):
    """
    Names, as the trace shows them, of the processes that write a tool-version record:
    the declared name of every process whose body calls ToolVersions.record, and every
    `include { declared as alias }` of one (the trace records the alias).
    """
    names = set()
    sources = sorted(Path(repo_dir).glob("Modules/*.nf")) + [Path(repo_dir) / "main.nf"]
    texts = [p.read_text() for p in sources if p.is_file()]
    for text in texts:
        starts = [m for m in PROCESS_RE.finditer(text)]
        for index, match in enumerate(starts):
            end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
            if "ToolVersions.record(" in text[match.end():end]:
                names.add(match.group(1))
    for text in texts:
        for include in re.findall(r"include\s*\{([^}]*)\}", text):
            for declared, alias in ALIAS_RE.findall(include):
                if declared in names:
                    names.add(alias)
    return names


def task_work_dir(row, work_dir):
    """The task's work directory: the trace's workdir column, else found from its hash."""
    workdir = (row.get("workdir") or "").strip()
    if workdir:
        return Path(workdir)
    task_hash = (row.get("hash") or "").strip()
    if work_dir and "/" in task_hash:
        hits = sorted(Path(work_dir).glob(f"{task_hash}*"))
        if len(hits) == 1:
            return hits[0]
    return None


def installed_version(conda_prefix, package, inventory):
    """What conda-meta says the task's own environment holds for `package`."""
    if not conda_prefix:
        return None
    name = Path(conda_prefix).name
    if name not in inventory:
        # Not under --conda-cache-dir (a hand run, or an env given as a path): read its
        # own conda-meta the same way conda_inventory() does.
        packages = {}
        for meta in sorted((Path(conda_prefix) / "conda-meta").glob("*.json")):
            parts = meta.stem.rsplit("-", 2)
            if len(parts) == 3:
                packages[parts[0]] = parts[1]
        inventory[name] = packages
    return inventory[name].get(package)


def collect_tool_versions(trace_file, work_dir, repo_dir, inventory, mode="strict"):
    """
    Read every task's record, compare it with its pin and with what conda installed, and
    say what does not agree.

    Returns {"check", "tools", "problems", "tasks_recorded", "tasks_unreadable"}.
    `tools` has one row per distinct (tool, version, path, environment); `problems` one
    per distinct problem, each with the processes and task count it covers.
    """
    summary = {"check": mode, "tools": [], "problems": [],
               "tasks_recorded": 0, "tasks_unreadable": 0}
    if not trace_file or not Path(trace_file).is_file():
        return summary
    recording = recording_processes(repo_dir)
    aliases = version_aliases()
    pins_by_env = {}
    tools, problems = {}, {}

    def problem(status, key, detail, process):
        entry = problems.setdefault((status,) + key, {"status": status, "detail": detail,
                                                      "processes": set(), "tasks": 0})
        entry["processes"].add(process)
        entry["tasks"] += 1

    with Path(trace_file).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = (row.get("name") or "").strip()
            process = name.split(" (")[0].strip()
            status = (row.get("status") or "").strip()
            simple = process.rsplit(":", 1)[-1]
            directory = task_work_dir(row, work_dir)
            record = directory / TOOL_RECORD if directory else None
            if record is None or not directory.is_dir():
                if simple in recording and status in RECORDED_STATUSES:
                    summary["tasks_unreadable"] += 1
                continue
            if not record.is_file():
                if simple in recording and status in RECORDED_STATUSES:
                    problem(MISSING, (simple,),
                            f"{simple}: finished without a tool-version record "
                            f"(e.g. {directory})", simple)
                continue
            summary["tasks_recorded"] += 1
            with record.open(newline="") as rec:
                for entry in csv.DictReader(rec, delimiter="\t"):
                    tool = entry.get("tool", "")
                    package = entry.get("package") or tool
                    version = entry.get("version", "")
                    path = entry.get("path", "")
                    prefix = entry.get("conda_prefix", "")
                    env_file = entry.get("env", "")
                    if env_file not in pins_by_env:
                        pins_by_env[env_file] = conda_pins(env_file)
                    pinned = pins_by_env[env_file].get(package, entry.get("pinned", ""))
                    installed = installed_version(prefix, package, inventory) or ""
                    task_status = entry.get("status", "")

                    # The task judged against the yml it was given; judge again here,
                    # against the same yml, and against what conda installed.
                    verdict = task_status
                    if task_status not in BLOCKING_STATUSES:
                        if pinned and version and not version_satisfies(version, pinned,
                                                                        package, aliases):
                            verdict = "VERSION_MISMATCH"
                        elif installed and version and not (
                                version_matches(installed, version)
                                or version_satisfies(version, installed, package, aliases)):
                            verdict = INSTALLED_MISMATCH

                    key = (tool, version, path, prefix)
                    row_out = tools.setdefault(key, {
                        "tool": tool, "package": package, "version": version,
                        "pinned": pinned, "installed": installed, "path": path,
                        "conda_prefix": prefix, "env": env_file, "status": verdict,
                        "raw": entry.get("raw", ""), "processes": set(), "tasks": 0})
                    row_out["processes"].add(simple)
                    row_out["tasks"] += 1

                    if verdict in BLOCKING_STATUSES:
                        where = f"from {path}" if path else "not on PATH"
                        detail = {
                            "OUTSIDE_ENV": f"{tool} {version or '(version unread)'} ran {where}, "
                                           f"outside the task's conda env {prefix}"
                                           + (f"; pinned {package}={pinned}" if pinned else ""),
                            "VERSION_MISMATCH": f"{tool} reports {version} {where}; "
                                                f"{Path(env_file).name or 'its env'} pins "
                                                f"{package}={pinned}",
                            "NOT_FOUND": f"{tool} was not on the task's PATH",
                            "NO_ENV": f"{tool} {version} ran {where} with no conda env active",
                        }[verdict]
                        problem(verdict, key, detail, simple)
                    elif verdict == INSTALLED_MISMATCH:
                        problem(verdict, key,
                                f"{tool} reports {version} {('from ' + path) if path else ''}, "
                                f"but conda installed {package} {installed} in {prefix}", simple)

    summary["tools"] = [dict(row, processes=sorted(row["processes"]))
                        for _key, row in sorted(tools.items())]
    summary["problems"] = [dict(entry, processes=sorted(entry["processes"]))
                           for _key, entry in sorted(problems.items())]
    return summary


def render_tool_versions(summary):
    """RUN_REPORT's TOOL VERSIONS section, and TOOL VERSION MISMATCH when there is one."""
    lines = section("TOOL VERSIONS (as run)")
    tools = summary.get("tools") or []
    if not tools:
        lines.append("  no task left a tool-version record (a run from before F16, or no trace)")
    else:
        lines.append("  reported by each binary from inside its task, after conda activation"
                     " (per task: <work dir>/.tool_versions.tsv)")
        lines.append(f"    {'tool':<34}{'version':<11}{'pinned':<10}{'installed':<11}"
                     f"{'tasks':>6}   ran from")
        for row in tools:
            prefix = row.get("conda_prefix") or ""
            path = row.get("path") or "-"
            where = (f"env {Path(prefix).name}" if prefix and path.startswith(prefix.rstrip("/") + "/")
                     else path)
            flag = "" if row["status"] in ("ok", "unpinned", "version_unread") else f"   {row['status']}"
            lines.append(f"    {row['tool']:<34}{row['version'] or '?':<11}{row['pinned'] or '-':<10}"
                         f"{row['installed'] or '-':<11}{row['tasks']:>6}   {where}{flag}")
    if summary.get("tasks_unreadable"):
        lines.append(f"  {summary['tasks_unreadable']} task(s) could not be checked: their work"
                     " directory was not found (cleaned up, or another machine)")

    problems = summary.get("problems") or []
    if problems:
        lines += section("TOOL VERSION MISMATCH")
        if summary.get("check") == "strict":
            lines.append("  --tool_version_check strict: a task that finds one of these stops the"
                         " run (exit 78)")
        else:
            lines.append("  --tool_version_check warn: the run carried on; these results were NOT"
                         " produced by the pinned tools")
        for entry in problems:
            lines.append(f"  {entry['status']:<20}{entry['detail']}")
            lines.append(f"  {'':<20}{', '.join(entry['processes'])}; {entry['tasks']} task(s)")
    return lines


def sample_identifiers(sheet_path):
    """
    Which samples went in. F16 asks for sample/data identifiers, and the sheet path
    alone does not survive the sheet being edited.
    """
    if not sheet_path:
        return {}
    path = Path(sheet_path)
    if not path.is_file():
        return {"sheet": str(sheet_path), "error": "sample sheet not found at finalisation"}
    try:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
    except OSError as problem:
        return {"sheet": str(sheet_path), "error": f"unreadable: {problem}"}

    ids = [str(row.get("patient", "")).strip() for row in rows]
    ids = [i for i in ids if i]
    joined = "\n".join(sorted(ids))
    return {
        "sheet": str(sheet_path),
        "count": len(ids),
        "identifiers": sorted(ids),
        "identifier_digest": hashlib.sha256(joined.encode()).hexdigest()[:12],
    }


# Parameters grouped the way someone reads them, rather than the ~90 resolved ones in
# alphabetical order. Anything not listed here stays in the JSON sidecar.
#
# Every name here must be a real params.* in main.nf (tests/test_provenance_record.py
# checks this). Ludmil, 2026-09-25 report P0-1: the tier group asked for
# tumor_*_min_genes while the parameters are tumor_*_min_clb_genes, and a name that
# does not exist was skipped exactly like a parameter left unset -- so RUN_REPORT
# silently dropped the gene thresholds from the tier definition.
TIER_NAMES = ("extensive_island", "broad_island", "multi_gene")
TIER_FIELDS = (("min_pks_reads", "pks reads"), ("min_clb_genes", "clb genes"),
               ("min_breadth", "island breadth"))

PARAMETER_GROUPS = (
    ("INPUT", ("sample", "input_data_type", "sample_type", "cram_reference", "outdir",
               # Ludmil, 26 Sep mode memo section 4: RNA library type is recorded.
               "rna_library_type",
               # Owner decision, 29 Sep: long-read input.
               "read_type",
               # E1 (Ludmil, 25 Sep): low-score mapped extraction from BAM/CRAM.
               "extract_low_score_mapped", "extract_min_as")),
    ("PROFILING", ("profiling_method", "hmm_evalue", "hmm_protein_evalue", "prefilter_mode",
                   "pks_contig", "pks_start_1based", "pks_end_1based",
                   "pks_island_len_1based",
                   # A4: the primary alignment's similarity rule.
                   "pks_align_bowtie2_mode", "pks_align_bowtie2_preset",
                   "pks_align_bowtie2_score_min", "pks_align_min_mapq",
                   # E0: host-depletion fragment mode.
                   "minimap2_frag_mode", "long_read_min_length")),
    ("EVIDENCE TIERS", tuple(f"tumor_{tier}_{field}" for tier in TIER_NAMES
                             for field, _label in TIER_FIELDS)),
    # Ludmil, 25 Sep A1/A2: the targeted-recruitment rule and where the island sits in
    # the island+flank recruitment reference.
    ("TARGETED RECRUITMENT", ("pks_recruit_min_aligned", "pks_recruit_min_identity",
                              "recruit_min_query_coverage", "recruit_boundary_tolerance_bp",
                              "recruit_min_mapq", "recruit_include_flank_only",
                              "pks_recruit_flank_bp")),
    ("STAGES", ("enable_mags", "tumor_targeted_assembly",
                "diamond_rescue", "pks_community_taxa",
                "enable_strain_typing", "pks_taxa")),
    # Ludmil, 26 Sep mode memo sections 3 and 4: the WES / RNA QC definitions.
    ("MODE QC", ("rna_transcript_min_reads", "rna_residual_human_sample_reads",
                 "rna_residual_human_args")),
    ("DATABASES", ("hg38_db", "rna_splice_db", "rna_junctions_bed", "wes_capture_bed",
                   "hg38_db_long", "rna_splice_db_long", "t2t_phix_db_long",
                   "t2t_phix_db", "pangenome_db", "kraken_db", "genomad_db",
                   "gtdbtk_db", "checkm2_db", "bracken_read_length")),
    ("RESOURCES", ("max_cpus", "max_memory", "max_time", "conda_cache_dir",
                   # F16: what a tool-version mismatch does to the run.
                   "tool_version_check")),
)


def is_set(value):
    return value not in (None, "", "null")


def grouped_parameters(params):
    """(group, [(name, value)]) for the parameters that were actually set."""
    grouped = []
    for title, names in PARAMETER_GROUPS:
        rows = [(name, params[name]) for name in names if is_set(params.get(name))]
        if rows:
            grouped.append((title, rows))
    return grouped


def unknown_parameter_names(params):
    """
    Report names this script asks for that the run's params map does not contain at all.

    A set-to-null parameter is still a key in the map; only a name that is not a
    pipeline parameter is absent. An empty map (a launch record written before params
    were recorded) cannot tell the two apart, so it is not checked.
    """
    if not params:
        return []
    return [name for _title, names in PARAMETER_GROUPS for name in names
            if name not in params]


def long_read_type(params):
    """The long read type (ont/hifi) the run used, or None for short reads."""
    value = params.get("read_type")
    return value if is_set(value) and value != "short" else None


# Owner decision, 29 Sep: which aligner made the profiling alignment. Mirrors
# lib/HostDepletion.DNA_PRESET for the long-read presets.
LONG_READ_PRESET = {"ont": "map-ont", "hifi": "map-hifi"}


def profiling_aligner(params):
    mapq = params.get("pks_align_min_mapq") or "?"
    long_read = long_read_type(params)
    if long_read:
        return (f"minimap2 -ax {LONG_READ_PRESET.get(long_read, '?')} against the IHE3034 genome "
                f"FASTA, MAPQ >= {mapq}, featureCounts -L")
    return (f"bowtie2 --{params.get('pks_align_bowtie2_mode') or '?'} "
            f"--{params.get('pks_align_bowtie2_preset') or '?'} "
            f"--score-min {params.get('pks_align_bowtie2_score_min') or '?'}, MAPQ >= {mapq}")


def island_interval(params):
    """The canonical analytical interval, as 'contig:start-end (N bp)', or None."""
    start, end = params.get("pks_start_1based"), params.get("pks_end_1based")
    if not (is_set(start) and is_set(end)):
        return None
    try:
        length = int(end) - int(start) + 1
    except ValueError:
        return None
    contig = params.get("pks_contig") if is_set(params.get("pks_contig")) else "?"
    return f"{contig}:{start}-{end}  (1-based, inclusive; {length:,} bp)"


def tier_definition(params):
    """One line per tier: every threshold it needs, or a note saying which are missing."""
    lines = []
    for tier in TIER_NAMES:
        parts, missing = [], []
        for field, label in TIER_FIELDS:
            value = params.get(f"tumor_{tier}_{field}")
            if is_set(value):
                parts.append(f"{label} >= {value}")
            else:
                missing.append(field)
        text = ", ".join(parts)
        if missing:
            text += f"  [MISSING: {', '.join(missing)}]"
        lines.append(f"    {tier:<20}{text}")
    return lines


def command_line_parameters(command_line):
    """
    The flags the operator actually typed -- the SigProfiler-style question of "what did
    I choose", as opposed to the ~90 resolved parameters, which are in the JSON.
    """
    if not command_line:
        return []
    tokens = command_line.split()
    chosen = []
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token.startswith("--") or (token.startswith("-") and len(token) > 1 and not token[1].isdigit()):
            value = ""
            if index + 1 < len(tokens) and not tokens[index + 1].startswith("-"):
                value = tokens[index + 1]
                index += 1
            chosen.append(f"{token} {value}".strip())
        index += 1
    return chosen


# ---------------------------------------------------------------- rendering

def rule(char="="):
    return char * 78


def section(title):
    return ["", title, rule("-")]


def render(record):
    """The run in reading order: did it work, what was asked, what ran, what came out."""
    run = record.get("run", {})
    code = record.get("code", {})
    invocation = record.get("invocation", {})
    environment = record.get("environment") or {}
    samples = record.get("samples") or {}
    lines = [rule(),
             f"pksProfiler {code.get('manifest_version', '?')}",
             rule()]

    # ---- JOB STATUS ----
    status = record.get("status", "unknown").upper()
    lines += section("JOB STATUS")
    # workflow.exitStatus is empty for an abort that never reached a task, which is
    # what a failed preflight is; "exit null" reads like a bug in the report.
    exit_status = run.get("exit_status")
    has_exit = exit_status not in (None, "", "null")
    lines.append(f"  {status}" + (f" -- exit {exit_status}"
                                  if run.get("success") is False and has_exit else ""))
    if run.get("success") is False and run.get("error_message"):
        lines.append(f"    {run['error_message']}")
    lines.append(f"  started    {run.get('started')}")
    lines.append(f"  finished   {run.get('completed')}  ({run.get('duration') or 'unknown'})")
    if environment.get("hostname"):
        lines.append(f"  host       {environment['hostname']}")
    lines.append(f"  results    {invocation.get('outdir') or 'unknown'}")
    # Ludmil, 26 Sep mode memo section 6: the sample mode, where nobody can miss it.
    run_params = invocation.get("params") or {}
    if is_set(run_params.get("sample_type")):
        lines.append(f"  mode       sample_type={run_params['sample_type']}")

    # ---- EXECUTION PARAMETERS ----
    lines += section("EXECUTION PARAMETERS")
    lines.append("  as given on the command line")
    lines.append(f"    {invocation.get('command_line')}")
    for item in invocation.get("parameters_given") or []:
        lines.append(f"      {item}")
    for title, rows in grouped_parameters(invocation.get("params") or {}):
        lines.append("")
        lines.append(f"  {title}")
        for name, value in rows:
            lines.append(f"    {name:<38}{value}")

    params = invocation.get("params") or {}
    if params:
        lines.append("")
        lines.append(f"  profiling aligner  {profiling_aligner(params)}")
        if params.get("profiling_method") != "hmm":
            # Ludmil, 25 Sep A5: say what the divergence columns are before anyone reads
            # them as more.
            lines.append("  divergence QC      NM / alignment columns of the clb-assigned reads: "
                         "divergence from IHE3034")
            lines.append("                     only (alignment/ihe3034_divergence.qc.tsv). Not "
                         "specificity or organism")
            lines.append("                     identity; specificity is the competitive-mapping "
                         "lane, relative to its panel.")
    interval = island_interval(params)
    long_read = long_read_type(params)
    if interval or params.get("sample_type") == "wgs":
        lines += section("PKS ISLAND AND TIER DEFINITION")
        lines.append(f"  analytical interval  {interval or 'NOT RECORDED'}")
        if params.get("sample_type") == "wgs" and long_read:
            lines.append(f"  evidence tiers       not applied: long-read WGS (read_type={long_read}) is "
                         "evidence-only; the tiers were calibrated on short-read Bowtie2")
        elif params.get("sample_type") == "wgs":
            lines.append("  evidence tiers (wgs; all thresholds must be met)")
            lines += tier_definition(params)

    unknown = unknown_parameter_names(params)
    if unknown:
        lines += section("PROVENANCE WARNING")
        lines.append("  This report asked for parameters the run does not define. The report")
        lines.append("  is incomplete; fix PARAMETER_GROUPS in scripts/finalize_provenance.py:")
        for name in unknown:
            lines.append(f"    {name}")

    # ---- SAMPLES ----
    lines += section("SAMPLES")
    sheet = invocation.get("sample_sheet") or {}
    if sheet.get("path"):
        lines.append(f"  sheet      {sheet['path']}")
        lines.append(f"  digest     {sheet.get('digest')}")
    if samples.get("count") is not None:
        lines.append(f"  samples    {samples['count']}"
                     f"   (identifiers: digest {samples.get('identifier_digest')},"
                     " full list in the JSON record)")
    elif samples.get("error"):
        lines.append(f"  {samples['error']}")
    elif not sheet.get("path"):
        lines.append("  no sample sheet recorded")

    # ---- INPUT ALIGNMENTS (E3) ----
    headers = record.get("input_headers")
    if headers:
        lines += section("INPUT ALIGNMENTS")
        files = next((r for r in headers if r.get("section") == "files"), {})
        lines.append(f"  BAM/CRAM headers read   {files.get('n_files', '?')}"
                     "   (per sample: by_sample/<sample>/provenance/input_header.tsv)")
        refs = [r for r in headers if r.get("section") == "reference"]
        lines.append(f"  distinct @SQ fingerprints  {len(refs)}")
        for r in refs:
            lines.append(f"    {r.get('key', '')[:16]}  {r.get('detail', '')}  files={r.get('n_files')}")
        chains = [r for r in headers if r.get("section") == "aligner_chain"]
        lines.append(f"  distinct @PG chains     {len(chains)}   (PN:VN in PP order)")
        for r in chains:
            lines.append(f"    files={r.get('n_files'):<6} {r.get('key', '')}")
        # E1: per-sample extraction rule, counted (on_auto_bwa_mem, off_auto_not_bwa_mem, ...).
        low_score = [r for r in headers if r.get("section") == "extract_low_score"]
        if low_score:
            lines.append("  low-score mapped extraction (AS below --extract_min_as; per sample in"
                         " input_header.tsv)")
            for r in low_score:
                lines.append(f"    files={r.get('n_files'):<6} {r.get('key', '')}")

    # ---- L4 / V2, from the cohort QC table (counts only) ----
    lines += render_qc_rollup(record.get("qc_rollup") or {})

    # ---- SOFTWARE ----
    software = record.get("software") or {}
    lines += section("SOFTWARE")
    if software:
        seen = {}
        for packages in software.values():
            for name, version in packages.items():
                if name in NOTABLE:
                    seen.setdefault(name, version)
        for name in sorted(seen):
            lines.append(f"    {name:<38}{seen[name]}")
    else:
        lines.append("  conda environments were not found, so installed versions are not recorded")
    for key in ("nextflow_version", "java_version"):
        if environment.get(key):
            lines.append(f"    {key:<38}{environment[key]}")
    if software:
        lines.append(f"  ({len(software)} conda environments; every package is in the JSON record)")
        lines.append("  what conda installed; what each task actually ran is under TOOL VERSIONS")

    # ---- TOOL VERSIONS (F16) ----
    lines += render_tool_versions(record.get("tool_versions") or {})

    # ---- WHAT RAN ----
    steps = record.get("steps") or []
    counts = record.get("task_counts") or {}
    lines += section("WHAT RAN")
    if counts:
        lines.append("  " + "   ".join(f"{count} {name}" for name, count in counts.items()))
        lines.append("")
    if steps:
        lines.append(f"    {'process':<38}{'tasks':>6}   outcome")
        for step in steps:
            outcome = ", ".join(f"{count} {status.lower()}"
                                for status, count in step["statuses"].items())
            lines.append(f"    {step['process']:<38}{step['tasks']:>6}   {outcome}")
    else:
        lines.append("  per-process detail is unavailable: runs/trace.txt was not readable"
                     " when the report was written")

    # ---- WHAT FAILED ----
    failures = record.get("failures") or []
    if failures:
        lines += section("WHAT FAILED")
        for failure in failures[:MAX_FAILURES_SHOWN]:
            lines.append(f"  {failure['task']}   {failure['status'].lower()}, exit {failure['exit']}")
            if failure.get("stderr"):
                lines.append(f"    stderr: {failure['stderr']}")
        if len(failures) > MAX_FAILURES_SHOWN:
            lines.append(f"  ... and {len(failures) - MAX_FAILURES_SHOWN} more"
                         " (all of them are in the JSON record)")

    # ---- OUTPUTS ----
    outputs = (record.get("outputs") or {}).get("cohort_files") or []
    lines += section("OUTPUTS")
    if outputs:
        lines.append(f"  cohort tables, sha256 (full digests in the JSON record)")
        for entry in outputs:
            lines.append(f"    {str(entry.get('sha256'))[:16]}  {entry.get('path')}")
    else:
        lines.append("  no cohort-level output was written")

    # ---- REPRODUCIBILITY ----
    lines += section("REPRODUCIBILITY")
    lines.append(f"    {'commit':<38}{code.get('commit')} ({code.get('commit_source')})")
    if code.get("dirty"):
        lines.append(f"    {'working tree':<38}UNCOMMITTED CHANGES -- this run is not exactly any commit")
        for entry in code.get("modified_files", []):
            lines.append(f"      {entry}")
    else:
        lines.append(f"    {'working tree':<38}clean")
    for key, value in sorted((record.get("dependencies") or {}).items()):
        lines.append(f"    {key:<38}{value}")
    lines.append(f"    {'session':<38}{run.get('session_id')}")
    lines.append(f"    {'work directory':<38}{invocation.get('work_dir')}")
    lines.append(f"    {'full record':<38}runs/provenance/{run.get('session_id')}.json")

    if record.get("notes"):
        lines += section("NOTES") + [f"  - {note}" for note in record["notes"]]

    return "\n".join(lines) + "\n"


# L4 (Ludmil, 1 Oct): the low-score rescue at each stage, and V2: how many tiered samples'
# tier would move if counted as fragments or deduplicated. Summed from the cohort QC
# table -- counts only, never a sample identifier.
RESCUE_STAGES = (("extracted_low_score_mapped_reads", "extracted (AS below --extract_min_as)"),
                 ("low_score_rescued_after_fastp", "after fastp"),
                 ("low_score_rescued_after_hg38", "after hg38 depletion"),
                 ("low_score_rescued_after_t2t_phix", "after T2T+phiX depletion"),
                 ("low_score_rescued_profiled", "in the profiling BAM"),
                 ("low_score_rescued_clb_assigned", "assigned to a clb gene"))


def qc_rollup(cohort_dir, started=None):
    """Sums over the cohort QC table's L4 columns, and the V2 tier-robustness tally."""
    path = Path(cohort_dir) / "qc" / "pks.qc.summary.tsv"
    if not path.is_file():
        return {}
    if started:
        try:
            if path.stat().st_mtime < datetime.datetime.fromisoformat(started).timestamp():
                return {}
        except (TypeError, ValueError):
            pass
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rollup = {}
    if rows and "low_score_rescued_after_fastp" in rows[0]:
        traced = [r for r in rows if (r.get("low_score_rescued_after_fastp") or "NA") != "NA"]
        totals = {}
        for column, _label in RESCUE_STAGES:
            values = [r.get(column) for r in traced]
            totals[column] = sum(int(v) for v in values if v not in (None, "", "NA"))
        rollup["low_score_rescue"] = {"samples": len(traced), "totals": totals}
    if rows and "tier_robustness" in rows[0]:
        tiered = [r.get("tier_robustness") for r in rows if r.get("tier_robustness") not in (None, "", "NA")]
        rollup["tier_robustness"] = {"tiered_samples": len(tiered),
                                     "stable": tiered.count("stable"),
                                     "unstable": tiered.count("unstable")}
    return rollup


def render_qc_rollup(rollup):
    lines = []
    rescue = rollup.get("low_score_rescue")
    if rescue:
        lines += section("LOW-SCORE RESCUE (L4)")
        lines.append(f"  reads E1 added at extraction, and how many are left after each stage;"
                     f" summed over the {rescue['samples']} sample(s)")
        lines.append("  that rescued any (per sample: cohort/qc/pks.qc.summary.tsv)")
        for column, label in RESCUE_STAGES:
            lines.append(f"    {label:<38}{rescue['totals'].get(column, 0)}")
    robustness = rollup.get("tier_robustness")
    if robustness:
        lines += section("TIER ROBUSTNESS (V2)")
        lines.append("  read_evidence is the tier on raw reads; this is whether counting fragments,"
                     " removing")
        lines.append("  duplicates, or both would give the same tier (tier_if_* in read_evidence.tsv)")
        lines.append(f"    {'tiered samples':<38}{robustness['tiered_samples']}")
        lines.append(f"    {'same tier every way (stable)':<38}{robustness['stable']}")
        lines.append(f"    {'tier would change (unstable)':<38}{robustness['unstable']}")
    return lines


def input_header_summary(cohort_dir, started=None):
    """E3: the counts-only input-header summary (header_provenance.py summarize).

    Read only if this session wrote it: a FASTQ run into an --outdir that once held a
    BAM run would otherwise report the old run's references as its own.
    """
    path = Path(cohort_dir) / "provenance" / "input_header_summary.tsv"
    if not path.is_file():
        return None
    if started:
        try:
            start = datetime.datetime.fromisoformat(started).timestamp()
            if path.stat().st_mtime < start:
                return None
        except (TypeError, ValueError):
            pass
    with path.open(newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t")]


# ---------------------------------------------------------------- driver

def main():
    args = parse_args()
    target = args.runs_dir / "provenance" / f"{args.session}.json"
    if not target.exists():
        # --help, or a failure before the launch record was written. Not an error.
        return 0

    record = json.loads(target.read_text())
    success = args.success.strip().lower() == "true"
    already_recorded = bool(record.get("finalized"))

    record["status"] = "completed" if success else "failed"
    record.setdefault("run", {}).update({
        "completed": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "duration": args.duration or None,
        "success": success,
        "exit_status": None if args.exit_status in ("", "null") else args.exit_status,
        "error_message": args.error_message or None,
    })
    invocation = record.setdefault("invocation", {})
    invocation["parameters_given"] = command_line_parameters(invocation.get("command_line"))
    record["samples"] = sample_identifiers((invocation.get("sample_sheet") or {}).get("path"))
    record["steps"], record["failures"] = scan_trace(args.trace_file, invocation.get("work_dir"))
    record["task_counts"] = {name: value for name, value in (
        ("succeeded", args.succeeded), ("cached", args.cached),
        ("failed", args.failed), ("ignored", args.ignored)) if value not in (None, "")}
    record["software"] = conda_inventory(args.conda_cache_dir)
    record["tool_versions"] = collect_tool_versions(
        args.trace_file, invocation.get("work_dir"), args.repo_dir,
        dict(record["software"]), args.tool_version_check)
    record["outputs"] = {"cohort_files": cohort_checksums(args.cohort_dir)}
    record["input_headers"] = input_header_summary(args.cohort_dir,
                                                   (record.get("run") or {}).get("started"))
    record["qc_rollup"] = qc_rollup(args.cohort_dir, (record.get("run") or {}).get("started"))
    record["finalized"] = True

    target.write_text(json.dumps(record, indent=2) + "\n")

    tool_problems = record["tool_versions"]["problems"]
    for entry in tool_problems:
        print(f"[provenance] TOOL VERSION MISMATCH {entry['status']}: {entry['detail']}",
              file=sys.stderr)
    # strict: a mismatch is a failure of the run's provenance, and says so in the exit
    # status for anyone running this by hand or from a driver script. (Nextflow ignores
    # an onComplete handler's status; there the task itself already stopped the run.)
    status = 1 if tool_problems and args.tool_version_check == "strict" else 0

    if already_recorded:
        print(f"[provenance] {target} (already reported, history not appended again)")
        return status

    report_text = render(record)
    unknown = unknown_parameter_names((record.get("invocation") or {}).get("params") or {})
    if unknown:
        print("[provenance] WARNING: RUN_REPORT asked for undefined parameters: "
              + ", ".join(unknown), file=sys.stderr)

    # The latest run, where someone handed a results directory will find it.
    if args.outdir:
        (args.outdir / "RUN_REPORT.txt").write_text(report_text)

    # Every run, kept: with -resume a results tree is the product of several sessions.
    history_dir = args.runs_dir / "reports"
    history_dir.mkdir(parents=True, exist_ok=True)
    index = len(list(history_dir.glob("*.txt"))) + 1
    (history_dir / f"{index:02d}-{args.session}.txt").write_text(report_text)

    print(f"[provenance] {(args.outdir / 'RUN_REPORT.txt') if args.outdir else history_dir}")
    return status


if __name__ == "__main__":
    sys.exit(main())
