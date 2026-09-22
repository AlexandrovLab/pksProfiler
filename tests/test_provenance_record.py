"""Every run records what produced it, and a resumed tree records every session.

F16, from Ludmil's 2026-09-19 report: nothing tied the released source to the code that
actually produced a cohort. main.nf writes a JSON record at launch; this finaliser adds
the outcome and the cohort checksums and appends a readable entry to RUN_REPORT.txt.

The launch half is Groovy (lib/RunRecord.groovy) because only Nextflow knows the
workflow metadata. This half is Python so it can be tested here -- Nextflow's config
parser cannot see lib/, which is what forced the split.
"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FINALIZE = ROOT / "scripts/finalize_provenance.py"
CONFIG = (ROOT / "nextflow.config").read_text()
MAIN = (ROOT / "main.nf").read_text()

LAUNCH_RECORD = {
    "record_version": 1,
    "status": "started",
    "code": {
        "commit": "d2d813691aa3ef942b52002a6ab3766c13231c47",
        "commit_source": "git",
        "dirty": True,
        "modified_files": ["M main.nf"],
        "manifest_version": "0.0.2dev",
        "project_dir": "/repo",
    },
    "invocation": {
        "command_line": "nextflow run main.nf --sample sheet.csv",
        "sample_sheet": {"path": "sheet.csv", "digest": "e28d07fcd51e/meta"},
        "stages_requested": ["profiling:bowtie2", "sample_type:tumor_wgs", "enable_strain_typing"],
        "params": {"outdir": "results"},
    },
    "dependencies": {"hg38_db": "4bf068530323/meta", "scripts": "f4d1377c8e52/content"},
    "environment": {"nextflow_version": "26.04.6"},
    "run": {"session_id": "SESSION-A", "started": "2026-09-21T09:00:00-07:00"},
    "outputs": {},
    "notes": ["a note"],
}


class ProvenanceFinalise(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.runs = self.root / "runs"
        self.cohort = self.root / "cohort"
        (self.runs / "provenance").mkdir(parents=True)
        (self.cohort / "qc").mkdir(parents=True)
        (self.cohort / "qc" / "pks.qc.summary.tsv").write_text("Sample\tstatus\nS1\tcomplete\n")

    def write_record(self, session="SESSION-A", **overrides):
        record = json.loads(json.dumps(LAUNCH_RECORD))
        record["run"]["session_id"] = session
        record.update(overrides)
        path = self.runs / "provenance" / f"{session}.json"
        path.write_text(json.dumps(record, indent=2) + "\n")
        return path

    def finalise(self, session="SESSION-A", success="true", exit_status="0",
                 duration="1h 2m", error="", extra=()):
        return subprocess.run(
            [sys.executable, str(FINALIZE), "--outdir", str(self.root),
             "--runs-dir", str(self.runs),
             "--cohort-dir", str(self.cohort), "--session", session,
             "--success", success, "--exit-status", exit_status,
             "--duration", duration, "--error-message", error, *extra],
            capture_output=True, text=True)

    def write_trace(self):
        path = self.root / "trace.txt"
        rows = ["task_id\thash\tnative_id\tname\tstatus\texit"]
        rows += [f"{i}\tab/cd\t1\textractReads (S{i})\tCOMPLETED\t0" for i in (1, 2)]
        rows.append("3\tab/ce\t1\tmapReads (S1)\tCACHED\t0")
        rows.append("4\tab/cf\t1\tpksProfiler_align (S2)\tFAILED\t137")
        path.write_text("\n".join(rows) + "\n")
        return path

    def write_conda_cache(self):
        cache = self.root / "condacache"
        meta = cache / "env-abc123" / "conda-meta"
        meta.mkdir(parents=True)
        for stem in ("python-3.10.18-h4de0772_0", "samtools-1.22-h96c455f_0",
                     "pandas-2.2.2-py310hf9f9076_1", "some-internal-lib-0.4-habc_0"):
            (meta / f"{stem}.json").write_text("{}")
        return cache

    def history(self):
        """The latest run's report, which is what a reader opens."""
        return (self.root / "RUN_REPORT.txt").read_text()

    def kept_reports(self):
        """One per session, kept: a resumed tree is the product of several."""
        return sorted((self.runs / "reports").glob("*.txt"))

    def test_a_successful_run_is_recorded_as_completed_with_checksums(self):
        path = self.write_record()
        result = self.finalise()
        self.assertEqual(result.returncode, 0, result.stderr)
        record = json.loads(path.read_text())
        self.assertEqual(record["status"], "completed")
        self.assertTrue(record["run"]["success"])
        self.assertEqual(record["run"]["duration"], "1h 2m")
        files = record["outputs"]["cohort_files"]
        self.assertEqual([f["path"] for f in files], ["qc/pks.qc.summary.tsv"])
        self.assertRegex(files[0]["sha256"], r"^[0-9a-f]{64}$")

    def test_a_failed_run_says_so_instead_of_staying_started(self):
        path = self.write_record()
        self.finalise(success="false", exit_status="137", error="mapReads killed")
        record = json.loads(path.read_text())
        self.assertEqual(record["status"], "failed")
        self.assertFalse(record["run"]["success"])
        self.assertEqual(record["run"]["exit_status"], "137")
        self.assertIn("mapReads killed", self.history())

    def test_finalising_twice_does_not_keep_two_copies_of_one_session(self):
        self.write_record()
        self.finalise()
        self.finalise()
        self.assertEqual(len(self.kept_reports()), 1)

    def test_every_session_is_kept_and_the_latest_is_the_headline(self):
        # A resumed tree is the product of several sessions; none may be lost, and
        # RUN_REPORT.txt must be the most recent one.
        for session in ("SESSION-A", "SESSION-B"):
            self.write_record(session)
            self.finalise(session)
        self.assertEqual([p.name for p in self.kept_reports()],
                         ["01-SESSION-A.txt", "02-SESSION-B.txt"])
        self.assertIn("SESSION-B", self.history())
        self.assertNotIn("SESSION-A", self.history())

    def test_no_launch_record_is_not_an_error(self):
        # --help, or a failure before the record was written.
        result = self.finalise(session="NEVER-STARTED")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse((self.root / "RUN_REPORT.txt").exists())

    def test_the_rendered_entry_carries_what_a_reader_needs(self):
        self.write_record()
        self.finalise()
        history = self.history()
        for expected in ("d2d813691aa3ef942b52002a6ab3766c13231c47",
                         "UNCOMMITTED CHANGES",
                         "e28d07fcd51e/meta",
                         "EXECUTION PARAMETERS",
                         "hg38_db",
                         "qc/pks.qc.summary.tsv"):
            with self.subTest(expected=expected):
                self.assertIn(expected, history)


class WhatRanAndWithWhat(unittest.TestCase):
    """The three things F16 asks for that the launch record cannot know."""

    def setUp(self):
        ProvenanceFinalise.setUp(self)
        self.sheet = self.root / "sheet.csv"
        self.sheet.write_text("patient,bam\nCASE_01,/data/a.bam\nCASE_02,/data/b.bam\n")

    write_record = ProvenanceFinalise.write_record
    finalise = ProvenanceFinalise.finalise
    write_trace = ProvenanceFinalise.write_trace
    write_conda_cache = ProvenanceFinalise.write_conda_cache
    history = ProvenanceFinalise.history

    def record_with_sheet(self):
        record = json.loads(json.dumps(LAUNCH_RECORD))
        record["invocation"]["sample_sheet"] = {"path": str(self.sheet), "digest": "abc/meta"}
        path = self.runs / "provenance" / "SESSION-A.json"
        path.write_text(json.dumps(record, indent=2) + "\n")
        return path

    def test_steps_are_grouped_by_process_with_their_outcome(self):
        path = self.record_with_sheet()
        self.finalise(extra=["--trace-file", str(self.write_trace())])
        steps = {s["process"]: s for s in json.loads(path.read_text())["steps"]}
        self.assertEqual(steps["extractReads"]["tasks"], 2)
        self.assertEqual(steps["extractReads"]["statuses"], {"COMPLETED": 2})
        self.assertEqual(steps["mapReads"]["statuses"], {"CACHED": 1})
        self.assertEqual(steps["pksProfiler_align"]["statuses"], {"FAILED": 1})
        self.assertIn("extractReads", self.history())

    def test_task_counts_do_not_depend_on_the_trace_file(self):
        # onComplete can run before the trace observer has flushed.
        self.record_with_sheet()
        self.finalise(extra=["--succeeded", "7", "--cached", "2", "--failed", "1", "--ignored", "0"])
        self.assertIn("7 succeeded", self.history())
        self.assertIn("per-process detail is unavailable", self.history())

    def test_installed_versions_come_from_the_environments_conda_built(self):
        path = self.record_with_sheet()
        self.finalise(extra=["--conda-cache-dir", str(self.write_conda_cache())])
        software = json.loads(path.read_text())["software"]
        self.assertEqual(software["env-abc123"]["python"], "3.10.18")
        self.assertEqual(software["env-abc123"]["samtools"], "1.22")
        history = self.history()
        self.assertRegex(history, r"python\s+3\.10\.18")
        # the full list is in the JSON; the report shows the ones a reader recognises
        self.assertNotIn("some-internal-lib", history)

    def test_sample_identifiers_are_recorded_from_the_sheet(self):
        path = self.record_with_sheet()
        self.finalise()
        samples = json.loads(path.read_text())["samples"]
        self.assertEqual(samples["count"], 2)
        self.assertEqual(samples["identifiers"], ["CASE_01", "CASE_02"])
        self.assertRegex(samples["identifier_digest"], r"^[0-9a-f]{12}$")
        # the report gives the count and the digest, not the list
        self.assertRegex(self.history(), r"samples\s+2\s+\(identifiers: digest")
        self.assertNotIn("CASE_01", self.history())

    def test_a_missing_sheet_is_reported_not_crashed_on(self):
        record = json.loads(json.dumps(LAUNCH_RECORD))
        record["invocation"]["sample_sheet"] = {"path": str(self.root / "gone.csv")}
        (self.runs / "provenance" / "SESSION-A.json").write_text(json.dumps(record) + "\n")
        result = self.finalise()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("not found", self.history())

    def test_the_parameters_the_operator_typed_are_singled_out(self):
        path = self.record_with_sheet()
        self.finalise()
        given = json.loads(path.read_text())["invocation"]["parameters_given"]
        self.assertIn("--sample sheet.csv", given)
        self.assertIn("as given on the command line", self.history())


class ReportStructure(unittest.TestCase):
    """The report has to be readable, which means grouped and pointing somewhere."""

    setUp = ProvenanceFinalise.setUp
    write_record = ProvenanceFinalise.write_record
    finalise = ProvenanceFinalise.finalise
    write_trace = ProvenanceFinalise.write_trace
    history = ProvenanceFinalise.history

    def test_parameters_are_grouped_not_listed_flat(self):
        record = json.loads(json.dumps(LAUNCH_RECORD))
        record["invocation"]["params"] = {
            "sample": "sheet.csv", "sample_type": "tumor_wgs",
            "profiling_method": "bowtie2", "enable_mags": "false",
            "hg38_db": "/dbs/hg38.mmi", "max_cpus": "64",
            "some_internal_thing": "42",
        }
        (self.runs / "provenance" / "SESSION-A.json").write_text(json.dumps(record) + "\n")
        self.finalise()
        report = self.history()
        for heading in ("INPUT", "PROFILING", "STAGES", "DATABASES", "RESOURCES"):
            with self.subTest(heading=heading):
                self.assertIn(f"  {heading}\n", report)
        # ungrouped parameters stay in the JSON rather than padding the report
        self.assertNotIn("some_internal_thing", report)

    def test_a_failed_task_says_where_its_stderr_is(self):
        record = json.loads(json.dumps(LAUNCH_RECORD))
        record["invocation"]["work_dir"] = "/scratch/pks/work"
        (self.runs / "provenance" / "SESSION-A.json").write_text(json.dumps(record) + "\n")
        self.finalise(extra=["--trace-file", str(self.write_trace())])
        report = self.history()
        self.assertIn("WHAT FAILED", report)
        self.assertIn("pksProfiler_align (S2)", report)
        self.assertIn("/scratch/pks/work/ab/cf*/.command.err", report)

    def test_a_clean_run_has_no_failure_section(self):
        self.write_record()
        self.finalise()
        self.assertNotIn("WHAT FAILED", self.history())

    def test_the_report_points_at_the_json_sidecar(self):
        self.write_record()
        self.finalise()
        self.assertIn("runs/provenance/SESSION-A.json", self.history())


class Wiring(unittest.TestCase):
    def test_the_workflow_writes_a_launch_record(self):
        self.assertIn("RunRecord.write(run_record_path, RunRecord.start(", MAIN)
        self.assertIn("dependencies: params.dep_digest", MAIN)

    def test_the_config_finalises_on_completion(self):
        self.assertIn("workflow.onComplete", CONFIG)
        self.assertIn("finalize_provenance.py", CONFIG)

    def test_the_config_passes_what_only_the_session_knows(self):
        for flag in ("--trace-file", "--conda-cache-dir", "--succeeded", "--failed"):
            with self.subTest(flag=flag):
                self.assertIn(flag, CONFIG)

    def test_a_provenance_failure_cannot_mask_the_runs_own_outcome(self):
        handler = CONFIG.split("workflow.onComplete", 1)[1]
        self.assertIn("catch", handler)


if __name__ == "__main__":
    unittest.main()


class CondaCacheAgreement(unittest.TestCase):
    """The cache Nextflow writes to and the one the report reads must be the same.

    The 2026-09-21 smoke test recorded "conda environments were not found" because
    params.conda_cache_dir was passed to the finaliser and never to Nextflow, which
    used its own <work-dir>/conda default. Two settings, one of them ignored.
    """

    def test_nextflow_is_told_where_the_cache_goes(self):
        self.assertIn("conda.cacheDir = params.conda_cache_dir", CONFIG)

    def test_the_report_reads_the_same_location(self):
        cache_lines = [l for l in CONFIG.splitlines() if "conda_cache" in l]
        # both the cacheDir setting and the --conda-cache-dir argument, same fallback
        self.assertGreaterEqual(len(cache_lines), 2)
        fallbacks = [l for l in cache_lines if "conda_cache" in l and "?:" in l]
        self.assertEqual(len(fallbacks), 2, f"expected both to share a fallback: {cache_lines}")
        self.assertTrue(all("conda_cache" in l for l in fallbacks))

    def test_the_default_is_not_inside_the_work_directory(self):
        # <work-dir>/conda means a deleted work dir rebuilds every environment.
        cacheline = [l for l in CONFIG.splitlines() if "conda.cacheDir" in l][0]
        self.assertNotIn("workDir", cacheline)
