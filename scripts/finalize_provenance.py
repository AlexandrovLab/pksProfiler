#!/usr/bin/env python3
"""Finish the run record, and write the run report a person actually reads.

F16 asks for: the production commit, parameters, enabled modules, reference digests,
sample/data identifiers, tool versions, trace, and output checksums.

The workflow writes the launch half (commit, parameters, reference digests,
environment) before any task starts, so a run that dies still says what it was
attempting. This adds what only exists afterwards -- outcome, which steps ran, the
software conda actually installed, the samples, the output checksums -- and writes
both forms:

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
PARAMETER_GROUPS = (
    ("INPUT", ("sample", "input_data_type", "sample_type", "cram_reference", "outdir")),
    ("PROFILING", ("profiling_method", "hmm_evalue", "hmm_protein_evalue", "prefilter_mode",
                   "pks_island_len", "pks_contig")),
    ("EVIDENCE TIERS", ("tumor_extensive_island_min_pks_reads", "tumor_extensive_island_min_genes",
                        "tumor_extensive_island_min_breadth", "tumor_broad_island_min_pks_reads",
                        "tumor_broad_island_min_genes", "tumor_broad_island_min_breadth",
                        "tumor_multi_gene_min_pks_reads", "tumor_multi_gene_min_genes",
                        "tumor_multi_gene_min_breadth")),
    ("STAGES", ("enable_mags", "tumor_enable_mags", "tumor_targeted_assembly",
                "tumor_full_contig_context", "diamond_rescue", "pks_community_taxa",
                "enable_strain_typing", "pks_taxa")),
    ("DATABASES", ("hg38_db", "t2t_phix_db", "pangenome_db", "kraken_db", "genomad_db",
                   "gtdbtk_db", "checkm2_db", "bracken_read_length")),
    ("RESOURCES", ("max_cpus", "max_memory", "max_time", "conda_cache_dir")),
)


def grouped_parameters(params):
    """(group, [(name, value)]) for the parameters that were actually set."""
    grouped = []
    for title, names in PARAMETER_GROUPS:
        rows = [(name, params[name]) for name in names
                if params.get(name) not in (None, "", "null")]
        if rows:
            grouped.append((title, rows))
    return grouped


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
    record["outputs"] = {"cohort_files": cohort_checksums(args.cohort_dir)}
    record["finalized"] = True

    target.write_text(json.dumps(record, indent=2) + "\n")

    if already_recorded:
        print(f"[provenance] {target} (already reported, history not appended again)")
        return 0

    report_text = render(record)

    # The latest run, where someone handed a results directory will find it.
    if args.outdir:
        (args.outdir / "RUN_REPORT.txt").write_text(report_text)

    # Every run, kept: with -resume a results tree is the product of several sessions.
    history_dir = args.runs_dir / "reports"
    history_dir.mkdir(parents=True, exist_ok=True)
    index = len(list(history_dir.glob("*.txt"))) + 1
    (history_dir / f"{index:02d}-{args.session}.txt").write_text(report_text)

    print(f"[provenance] {(args.outdir / 'RUN_REPORT.txt') if args.outdir else history_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
