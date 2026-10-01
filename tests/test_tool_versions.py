"""Every task records the tool versions it actually ran, and a mismatch fails loudly.

F16 (Ludmil) asks for tool versions. RUN_REPORT's SOFTWARE section listed what conda
installed, read from conda-meta after the run. On TSCC, ~/.bashrc put a user minimap2
2.24 ahead of the conda base on PATH, so host depletion ran 2.24 in every task while the
report said 2.31 (29 Sep). It reported what was installed, not what ran.

Three layers, tested here without Nextflow:

  scripts/record_tool_versions.sh   in each task: version + path of every tool it runs,
                                    compared with the task's own env yml; strict exits 78
  lib/ToolVersions.groovy           renders that call into every process that runs a tool
  scripts/finalize_provenance.py    reads every task's record through the trace, compares
                                    it with the yml pins and conda-meta, and reports

The four cases asked for: a match, a version mismatch, a tool that ran from outside the
task's env (the 2.24 scenario), and a task that left no record.
"""
import csv
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/record_tool_versions.sh"
FINALIZE = ROOT / "scripts/finalize_provenance.py"
LIB = (ROOT / "lib/ToolVersions.groovy").read_text()
MAIN = (ROOT / "main.nf").read_text()
CONFIG = (ROOT / "nextflow.config").read_text()
BASE_CONFIG = (ROOT / "conf/base.config").read_text()
MODULES = {p.name: p.read_text() for p in sorted((ROOT / "Modules").glob("*.nf"))}

sys.path.insert(0, str(ROOT / "scripts"))
import finalize_provenance as fp  # noqa: E402

MINIMAP2_YML = ROOT / "conda_envs/minimap2_env.yml"


def fake_tool(directory, name, output):
    """An executable that prints `output`, the way the real tool reports its version."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text("#!/bin/sh\ncat <<'EOF'\n" + output + "\nEOF\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def read_record(path):
    with Path(path).open(newline="") as handle:
        return {row["tool"]: row for row in csv.DictReader(handle, delimiter="\t")}


class Workspace(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        # A conda env as Nextflow's cache lays it out: <cache>/env-<hash>/{bin,conda-meta}.
        self.env = self.root / "conda_cache" / "env-0123456789abcdef"
        (self.env / "conda-meta").mkdir(parents=True)
        for package, version in (("minimap2", "2.31"), ("samtools", "1.21"), ("htslib", "1.21")):
            (self.env / "conda-meta" / f"{package}-{version}-h0_0.json").write_text("{}")
        fake_tool(self.env / "bin", "minimap2", "2.31-r1297")
        fake_tool(self.env / "bin", "samtools", "samtools 1.21\nUsing htslib 1.21")
        fake_tool(self.env / "bin", "bgzip", "bgzip (htslib) 1.21")
        # The TSCC trap: a user build earlier on PATH than the environment.
        self.user_bin = self.root / "home" / "tools" / "minimap2"
        fake_tool(self.user_bin, "minimap2", "2.24-r1122")

    def run_helper(self, mode="strict", path_dirs=None, env_file=MINIMAP2_YML,
                   tools=("minimap2", "samtools", "bgzip:htslib"), conda_prefix=None):
        out = self.root / "work" / ".tool_versions.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        dirs = path_dirs if path_dirs is not None else [self.env / "bin"]
        environ = {"PATH": ":".join([str(d) for d in dirs] + ["/usr/bin", "/bin"]),
                   "HOME": str(self.root)}
        prefix = self.env if conda_prefix is None else conda_prefix
        if prefix:
            environ["CONDA_PREFIX"] = str(prefix)
        result = subprocess.run(
            ["bash", str(HELPER), mode, str(env_file), "mapReads", str(out), *tools],
            capture_output=True, text=True, env=environ)
        return result, (read_record(out) if out.exists() else {})


class TheTaskRecordsWhatItRuns(Workspace):
    def test_match_the_pinned_build_from_the_tasks_own_env(self):
        result, record = self.run_helper()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record["minimap2"]["status"], "ok")
        self.assertEqual(record["minimap2"]["version"], "2.31")
        self.assertEqual(record["minimap2"]["pinned"], "2.31")
        self.assertEqual(record["minimap2"]["path"], str(self.env / "bin" / "minimap2"))
        self.assertEqual(record["minimap2"]["conda_prefix"], str(self.env))
        self.assertEqual(record["samtools"]["status"], "ok")
        self.assertEqual(record["samtools"]["version"], "1.21")
        # htslib is listed without a version in minimap2_env.yml.
        self.assertEqual(record["bgzip"]["package"], "htslib")
        self.assertEqual(record["bgzip"]["status"], "unpinned")
        self.assertEqual(result.stderr, "")

    def test_path_outside_env_the_minimap2_224_scenario_stops_the_task(self):
        result, record = self.run_helper(path_dirs=[self.user_bin, self.env / "bin"])
        self.assertEqual(result.returncode, 78)
        self.assertEqual(record["minimap2"]["status"], "OUTSIDE_ENV")
        self.assertEqual(record["minimap2"]["version"], "2.24")
        self.assertEqual(record["minimap2"]["path"], str(self.user_bin / "minimap2"))
        self.assertIn("ERROR: tool-version check, mapReads", result.stderr)
        self.assertIn(str(self.user_bin / "minimap2"), result.stderr)
        self.assertIn("pinned minimap2=2.31", result.stderr)
        # The others are still recorded, and still fine.
        self.assertEqual(record["samtools"]["status"], "ok")

    def test_warn_records_it_and_carries_on(self):
        result, record = self.run_helper(mode="warn", path_dirs=[self.user_bin, self.env / "bin"])
        self.assertEqual(result.returncode, 0)
        self.assertEqual(record["minimap2"]["status"], "OUTSIDE_ENV")
        self.assertIn("WARNING: tool-version check", result.stderr)

    def test_version_mismatch_inside_the_env_stops_the_task(self):
        fake_tool(self.env / "bin", "minimap2", "2.30-r1287")
        result, record = self.run_helper()
        self.assertEqual(result.returncode, 78)
        self.assertEqual(record["minimap2"]["status"], "VERSION_MISMATCH")
        self.assertIn("pins minimap2=2.31", result.stderr)

    def test_a_more_specific_release_of_the_pin_matches(self):
        fake_tool(self.env / "bin", "samtools", "samtools 1.21.1")
        result, record = self.run_helper()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record["samtools"]["status"], "ok")

    def test_a_pin_prefix_is_not_a_version_prefix(self):
        # 2.3 must not accept 2.31.
        env_file = self.root / "env.yml"
        env_file.write_text("dependencies:\n  - minimap2=2.3\n")
        result, record = self.run_helper(env_file=env_file, tools=("minimap2",))
        self.assertEqual(record["minimap2"]["status"], "VERSION_MISMATCH")
        self.assertEqual(result.returncode, 78)

    def test_a_tool_not_on_path(self):
        result, record = self.run_helper(tools=("minimap2", "featureCounts:subread"))
        self.assertEqual(result.returncode, 78)
        self.assertEqual(record["featureCounts"]["status"], "NOT_FOUND")
        self.assertEqual(record["featureCounts"]["package"], "subread")

    def test_no_conda_env_active(self):
        result, record = self.run_helper(conda_prefix="")
        self.assertEqual(result.returncode, 78)
        self.assertEqual(record["minimap2"]["status"], "NO_ENV")

    def test_an_env_given_as_a_directory_is_path_checked_only(self):
        result, record = self.run_helper(env_file=self.env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record["minimap2"]["status"], "unpinned")

    def test_pin_forms_conda_accepts(self):
        env_file = self.root / "env.yml"
        env_file.write_text(
            "name: x\nchannels:\n  - bioconda\ndependencies:\n"
            "  - bioconda::minimap2==2.31   # channel prefix, double equals\n"
            "  - samtools=1.21=h50ea8bc_0\n"
            "  - htslib>=1.20\n")
        result, record = self.run_helper(env_file=env_file)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record["minimap2"]["pinned"], "2.31")
        self.assertEqual(record["samtools"]["pinned"], "1.21")
        # A range is not a pin; only the path is checked.
        self.assertEqual(record["bgzip"]["pinned"], "")
        self.assertEqual(fp.conda_pins(env_file),
                         {"minimap2": "2.31", "samtools": "1.21"})

    def test_a_path_in_the_version_output_is_not_read_as_the_version(self):
        # bowtie2 prints its own absolute path first; an env dir named python3.10 must not win.
        env = self.root / "envs" / "python3.10"
        fake_tool(env / "bin", "bowtie2",
                  f"{env}/bin/bowtie2-align-s version 2.5.4\n64-bit\nBuilt on host")
        env_file = self.root / "env.yml"
        env_file.write_text("dependencies:\n  - bowtie2=2.5.4\n")
        result, record = self.run_helper(path_dirs=[env / "bin"], env_file=env_file,
                                         tools=("bowtie2",), conda_prefix=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record["bowtie2"]["version"], "2.5.4")

    def test_a_bad_mode_is_refused(self):
        result, _record = self.run_helper(mode="lenient")
        self.assertEqual(result.returncode, 2)
        self.assertIn("MODE must be strict or warn", result.stderr)


class TheReportComparesWhatRanWithThePins(Workspace):
    """scripts/finalize_provenance.py, end to end, on records the helper really wrote."""

    def make_run(self, tasks, mode="strict", with_workdir_column=True):
        """tasks: [(name, status, workdir or None)]. Returns (result, record, report)."""
        runs = self.root / "results" / "runs"
        (runs / "provenance").mkdir(parents=True)
        (self.root / "results" / "cohort").mkdir(parents=True)
        (runs / "provenance" / "S.json").write_text(json.dumps({
            "code": {}, "invocation": {"params": {}, "work_dir": str(self.root / "work")},
            "run": {"session_id": "S"}}))
        trace = self.root / "trace.txt"
        header = ["task_id", "hash", "name", "status", "exit"] + (["workdir"] if with_workdir_column else [])
        rows = ["\t".join(header)]
        for index, (name, status, workdir) in enumerate(tasks, start=1):
            fields = [str(index), f"{index:02d}/abcdef", name, status, "0"]
            if with_workdir_column:
                fields.append(str(workdir) if workdir else "")
            rows.append("\t".join(fields))
        trace.write_text("\n".join(rows) + "\n")
        result = subprocess.run(
            [sys.executable, str(FINALIZE), "--outdir", str(self.root / "results"),
             "--runs-dir", str(runs), "--cohort-dir", str(self.root / "results" / "cohort"),
             "--session", "S", "--success", "true", "--trace-file", str(trace),
             "--conda-cache-dir", str(self.env.parent), "--tool-version-check", mode],
            capture_output=True, text=True)
        record = json.loads((runs / "provenance" / "S.json").read_text())
        report = (self.root / "results" / "RUN_REPORT.txt").read_text()
        return result, record, report

    def task_dir(self, name, **helper):
        """Run the helper for one task and return its work dir."""
        result, _record = self.run_helper(**helper)
        directory = self.root / "tasks" / name
        directory.mkdir(parents=True)
        (self.root / "work" / ".tool_versions.tsv").rename(directory / fp.TOOL_RECORD)
        return directory, result

    def test_match_is_listed_and_passes(self):
        one, _ = self.task_dir("a")
        two, _ = self.task_dir("b")
        result, record, report = self.make_run([("mapReads (S1)", "COMPLETED", one),
                                                ("mapReads (S2)", "CACHED", two)])
        self.assertEqual(result.returncode, 0, result.stderr)
        versions = record["tool_versions"]
        self.assertEqual(versions["problems"], [])
        self.assertEqual(versions["tasks_recorded"], 2)
        minimap2 = next(t for t in versions["tools"] if t["tool"] == "minimap2")
        self.assertEqual((minimap2["version"], minimap2["pinned"], minimap2["installed"],
                          minimap2["tasks"], minimap2["status"]),
                         ("2.31", "2.31", "2.31", 2, "ok"))
        self.assertIn("TOOL VERSIONS (as run)", report)
        self.assertRegex(report, r"minimap2\s+2\.31\s+2\.31\s+2\.31\s+2\s+env env-0123456789abcdef")
        self.assertNotIn("TOOL VERSION MISMATCH", report)

    def test_path_outside_env_fails_the_report_the_224_scenario(self):
        # warn, so the task carried on and the report is the only place this shows.
        bad, result = self.task_dir("a", mode="warn", path_dirs=[self.user_bin, self.env / "bin"])
        self.assertEqual(result.returncode, 0)
        for mode, expected_status in (("strict", 1), ("warn", 0)):
            with self.subTest(mode=mode):
                shutil_rmtree(self.root / "results")
                result, record, report = self.make_run([("mapReads (S1)", "COMPLETED", bad)], mode=mode)
                self.assertEqual(result.returncode, expected_status, result.stderr)
                problem = record["tool_versions"]["problems"][0]
                self.assertEqual(problem["status"], "OUTSIDE_ENV")
                self.assertIn(str(self.user_bin / "minimap2"), problem["detail"])
                self.assertIn("TOOL VERSION MISMATCH", report)
                self.assertRegex(report, r"OUTSIDE_ENV\s+minimap2 2\.24 ran from .*home/tools/minimap2/minimap2")
                self.assertIn("pinned minimap2=2.31", report)
                self.assertIn("TOOL VERSION MISMATCH OUTSIDE_ENV", result.stderr)

    def test_version_mismatch_against_the_yml_is_caught_even_if_the_task_passed_it(self):
        # The task judged against the yml it was given; the report re-reads the yml.
        good, _ = self.task_dir("a")
        with (good / fp.TOOL_RECORD).open() as handle:
            text = handle.read()
        env_file = self.root / "repinned.yml"
        env_file.write_text("dependencies:\n  - minimap2=2.28\n")
        (good / fp.TOOL_RECORD).write_text(text.replace(str(MINIMAP2_YML), str(env_file)))
        result, record, report = self.make_run([("mapReads (S1)", "COMPLETED", good)])
        self.assertEqual(result.returncode, 1)
        statuses = {p["status"] for p in record["tool_versions"]["problems"]}
        self.assertEqual(statuses, {"VERSION_MISMATCH"})
        self.assertIn("repinned.yml pins minimap2=2.28", report)

    def test_a_binary_conda_did_not_install_is_caught(self):
        good, _ = self.task_dir("a")
        (self.env / "conda-meta" / "minimap2-2.31-h0_0.json").rename(
            self.env / "conda-meta" / "minimap2-2.30-h0_0.json")
        result, record, report = self.make_run([("mapReads (S1)", "COMPLETED", good)])
        self.assertEqual(result.returncode, 1)
        problem = record["tool_versions"]["problems"][0]
        self.assertEqual(problem["status"], "INSTALLED_MISMATCH")
        self.assertIn("conda installed minimap2 2.30", problem["detail"])

    def test_missing_record_from_a_process_that_writes_one(self):
        empty = self.root / "tasks" / "empty"
        empty.mkdir(parents=True)
        # pksProfilerAlign is an include alias of pksProfiler_align: the trace shows the alias.
        result, record, report = self.make_run([("mapReads (S1)", "COMPLETED", empty),
                                                ("pksProfilerAlign (S1)", "CACHED", empty)])
        self.assertEqual(result.returncode, 1)
        missing = {p["processes"][0] for p in record["tool_versions"]["problems"]
                   if p["status"] == "MISSING"}
        self.assertEqual(missing, {"mapReads", "pksProfilerAlign"})
        self.assertRegex(report, r"MISSING\s+mapReads: finished without a tool-version record")

    def test_missing_is_not_claimed_for_processes_that_run_no_tool_or_did_not_finish(self):
        empty = self.root / "tasks" / "empty"
        empty.mkdir(parents=True)
        result, record, _report = self.make_run([
            ("pksMAG:magSampleStatus (S1)", "COMPLETED", empty),   # runs no external tool
            ("mapReads (S2)", "FAILED", empty),                     # its failure is reported elsewhere
        ])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record["tool_versions"]["problems"], [])

    def test_a_cleaned_work_dir_is_counted_not_failed(self):
        result, record, report = self.make_run([("mapReads (S1)", "COMPLETED", self.root / "gone")])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record["tool_versions"]["tasks_unreadable"], 1)
        self.assertIn("1 task(s) could not be checked", report)

    def test_a_trace_without_the_workdir_column_falls_back_to_the_hash(self):
        work = self.root / "work" / "01"
        work.mkdir(parents=True)
        good, _ = self.task_dir("a")
        good.rename(work / "abcdef0123456789")
        result, record, _report = self.make_run([("mapReads (S1)", "COMPLETED", None)],
                                                with_workdir_column=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record["tool_versions"]["tasks_recorded"], 1)

    def test_no_records_at_all_says_so(self):
        result, record, report = self.make_run([])
        self.assertEqual(result.returncode, 0)
        self.assertIn("no task left a tool-version record", report)


def shutil_rmtree(path):
    import shutil
    shutil.rmtree(path, ignore_errors=True)


# ---------------------------------------------------------------- construction

def process_blocks(text):
    starts = list(re.finditer(r"^\s*process\s+(\w+)\s*\{", text, re.M))
    for index, match in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        yield match.group(1), text[match.end():end]


def script_block(body):
    """The process's own script: the last triple-quoted string after `script:`."""
    tail = body.split("script:", 1)[1]
    blocks = re.findall(r'"""(.*?)"""', tail, re.S)
    return blocks[-1] if blocks else ""


def env_files(conda_expr, module_text):
    """The conda_envs/*.yml a conda directive can resolve to."""
    files = set(re.findall(r"conda_envs/([\w.-]+\.yml)", conda_expr))
    for name in re.findall(r"params\.(\w+_env)\b", conda_expr):
        for source in (MAIN, module_text):
            for value in re.findall(rf"^params\.{name}\s*=\s*(.+)$", source, re.M):
                files.update(re.findall(r"conda_envs/([\w.-]+\.yml)", value))
    return files


# Processes that declare a conda env but run no external tool of their own (gzip, awk
# and the shell are the host's and are not pinned anywhere).
RUNS_NO_PINNED_TOOL = {
    "mergePksCandidates": "gzip/awk/cat only",
    "strainTypeSummary": "awk/cat only",
}


class EveryToolRunningProcessRecords(unittest.TestCase):
    def test_every_process_with_a_conda_env_records_its_tools(self):
        missing = []
        for fname, text in MODULES.items():
            for name, body in process_blocks(text):
                if not re.search(r'^\s*conda\s', body, re.M) or name in RUNS_NO_PINNED_TOOL:
                    continue
                if script_block(body).count("ToolVersions.record(") != 1:
                    missing.append(f"{fname}:{name}")
        self.assertEqual(missing, [], "these run tools without recording which build ran")

    def test_the_exemptions_really_run_no_pinned_tool(self):
        for fname, text in MODULES.items():
            for name, body in process_blocks(text):
                if name in RUNS_NO_PINNED_TOOL:
                    code = "\n".join(l for l in script_block(body).splitlines()
                                     if not l.strip().startswith("#"))
                    for tool in ("python", "samtools", "minimap2", "bowtie2", "Rscript", "bracken"):
                        self.assertNotRegex(code, rf"(?<![\w.-]){tool}\b", f"{name} runs {tool}")

    def test_each_recorded_package_is_in_the_env_the_process_runs_in(self):
        # Would have caught scoreCompetitiveRealignment: python3 in minimap2_env.yml, which
        # has no python, so it ran whichever python the launching shell had.
        problems = []
        for fname, text in MODULES.items():
            for name, body in process_blocks(text):
                call = re.search(r"ToolVersions\.record\(task, params, (.+?)\)\}", body)
                if not call:
                    continue
                conda = re.search(r'^\s*conda\s+(.+)$', body, re.M).group(1)
                envs = env_files(conda, text)
                self.assertTrue(envs, f"{name}: cannot resolve {conda}")
                deps = set()
                for env in envs:
                    deps |= {re.split(r"[=<>!~ ]", line.strip()[1:].strip())[0].split("::")[-1]
                             for line in (ROOT / "conda_envs" / env).read_text().splitlines()
                             if line.strip().startswith("- ")}
                tools = " ".join(re.findall(r"'([^']*)'", call.group(1))).split()
                for token in tools:
                    package = token.split(":", 1)[-1]
                    if package not in deps:
                        problems.append(f"{fname}:{name} records {token}, not in {sorted(envs)}")
        self.assertEqual(problems, [])

    def test_the_record_is_inside_the_process_script_not_a_groovy_string(self):
        # filterReads builds a nested """...""" for its awk pipe; the call belongs in the
        # task script, not in there.
        for fname, text in MODULES.items():
            for name, body in process_blocks(text):
                if "ToolVersions.record(" in body:
                    self.assertIn("ToolVersions.record(", script_block(body), f"{fname}:{name}")

    def test_the_rendered_call_writes_next_to_command_sh_and_propagates_its_status(self):
        # Under `scratch true` the current directory is node-local and discarded; $0 is
        # the work dir's .command.sh. And .command.sh runs with -u but not -e.
        self.assertIn('\\"\\$(dirname \\"\\$0\\")/.tool_versions.tsv\\"', LIB)
        self.assertIn("|| exit \\$?", LIB)
        self.assertIn("task.conda", LIB)
        self.assertIn("params.tool_version_check", LIB)


class StandaloneModulesRecordToo(unittest.TestCase):
    """pks_specificity.nf and pks_carrier.nf run on their own (tests/check_v8_lanes_execution.sh)."""

    def test_standalone_modules_find_the_lib_and_default_the_mode(self):
        link = ROOT / "Modules" / "lib"
        self.assertTrue(link.is_symlink())
        self.assertEqual(Path(os.readlink(link)), Path("../lib"))
        for name in ("pks_specificity.nf", "pks_carrier.nf"):
            self.assertIn('params.tool_version_check = params.tool_version_check ?: "strict"', MODULES[name])

    def test_the_v8_standalone_runs_use_conda(self):
        # No nextflow.config is loaded for a module run on its own, so conda was off and
        # the tasks ran the launching shell's tools.
        script = (ROOT / "tests/check_v8_lanes_execution.sh").read_text()
        self.assertIn("conda.enabled = true", script)
        self.assertIn("conda.cacheDir", script)


class StrictIsTheDefaultAndStopsTheRun(unittest.TestCase):
    def test_default_is_strict_in_both_places_it_is_read(self):
        self.assertRegex(MAIN, r'params\.tool_version_check\s*=\s*"strict"')
        self.assertRegex(CONFIG, r"tool_version_check\s*=\s*'strict'")

    def test_main_refuses_an_unknown_mode(self):
        self.assertIn("ToolVersions.validate(params.tool_version_check)", MAIN)
        self.assertIn("static final List<String> MODES = ['strict', 'warn']", LIB)

    def test_exit_78_terminates_in_every_error_strategy(self):
        strategies = re.findall(r"errorStrategy\s*=\s*\{([^}]*)\}", BASE_CONFIG)
        self.assertEqual(len(strategies), 3)
        for strategy in strategies:
            self.assertIn("task.exitStatus == 78 ? 'terminate'", strategy)
        self.assertIn("EXIT_STATUS = 78", LIB)
        self.assertIn("exit 78", HELPER.read_text())

    def test_the_trace_carries_workdir_and_the_report_gets_the_mode(self):
        fields = re.search(r"trace\.fields\s*=\s*'([^']+)'", CONFIG).group(1).split(",")
        for field in ("hash", "name", "status", "exit", "workdir"):
            self.assertIn(field, fields)
        self.assertIn("'--tool-version-check'", CONFIG)

    def test_the_helper_never_pipes_into_head(self):
        # head under pipefail: the writer can die of SIGPIPE and fail the pipe.
        code = "\n".join(l for l in HELPER.read_text().splitlines() if not l.lstrip().startswith("#"))
        self.assertNotRegex(code, r"\|\s*head\b")


if __name__ == "__main__":
    unittest.main()
