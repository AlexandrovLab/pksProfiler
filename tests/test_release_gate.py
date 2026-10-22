"""A tag build must run the real checks, and must not call a skipped suite green.

rc2 review section 2, items 2.2 and 2.3:

    Run real-execution tests on the actual release tag. The release CI should execute the
    end-to-end and mode-specific checks on the exact v1.0.0 tag, not only unit
    tests/graph construction.

    Fail the release if critical tests are skipped. The fast test environment can skip
    rescue, fragment-support and competitive-realignment tests when external tools are
    absent. Add an explicit skip-count/required-test assertion to the release job so
    'tests pass' cannot exclude the most consequential coverage.

Our own rc2 run reported `OK (skipped=13)`. Nothing looked at the 13, and they were the
tests that execute the AS<30 rescue, the fragment/duplicate collapse behind
tier_robustness, and the competitive-realignment score.

Three things are tested here:

  * the workflow file wires it up -- a tag triggers the execution job, the skip assertion
    runs where the tools exist, and the gate job exists and depends on the right jobs;
  * assert_required_coverage.py's logic, against synthetic results: a skipped required
    test fails, an unlisted skip fails, an over-budget skip count fails, an allowed skip
    passes;
  * verify_local_execution_record.py's logic, against synthetic records: a record from a
    different source tree is refused, a missing or failed required check is refused, a
    zero assertion count is refused, an unsigned record is refused.

The point of the last two is that these two scripts are the release gate. If they can be
satisfied by a run that did not happen, the gate is decoration.
"""
import importlib.util
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CI = (ROOT / ".github/workflows/ci.yml").read_text()


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tests" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


coverage = load("assert_required_coverage")
record = load("verify_local_execution_record")


class FakeResult:
    """Just enough of unittest's result for assert_required_coverage.check()."""

    def __init__(self, ran, skipped_by_module, failures=(), errors=()):
        self.ran_modules = set(ran)
        self.skipped_by_module = dict(skipped_by_module)
        self.failures = list(failures)
        self.errors = list(errors)
        self.skipped = [entry for skips in self.skipped_by_module.values() for entry in skips]
        self.testsRun = 1000


def all_required_ran(**overrides):
    return FakeResult(ran=set(coverage.REQUIRED) | {"test_other"},
                      skipped_by_module=overrides.get("skipped", {}),
                      failures=overrides.get("failures", ()),
                      errors=overrides.get("errors", ()))


class TheWorkflowRunsTheRealChecksOnATag(unittest.TestCase):
    def test_a_tag_triggers_the_workflow_at_all(self):
        self.assertIn('tags: ["v*"]', CI)

    def test_a_tag_triggers_the_execution_job(self):
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        condition = job.split("steps:", 1)[0]
        self.assertIn("startsWith(github.ref, 'refs/tags/v')", condition,
                      "a v* tag must run the real execution checks, not only unit tests")

    def test_the_execution_job_runs_every_execution_check_the_runner_can(self):
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        for script in ("check_e2e_execution.sh", "check_malformed_input_execution.sh",
                       "check_mode_execution.sh", "check_v8_lanes_execution.sh",
                       "check_resume_reference_change.sh"):
            with self.subTest(script=script):
                self.assertIn(script, job)

    def test_the_skip_assertion_runs_where_the_tools_exist(self):
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        self.assertIn("assert_required_coverage.py", job)
        # After the checks, so the conda environments it resolves tools from are built.
        self.assertGreater(job.index("assert_required_coverage.py"),
                           job.index("check_resume_reference_change.sh"))

    def test_the_gate_job_exists_and_is_tag_only(self):
        self.assertIn("  release-gate:", CI)
        job = CI.split("  release-gate:", 1)[1]
        self.assertIn("if: startsWith(github.ref, 'refs/tags/v')", job)

    def test_the_gate_waits_for_the_jobs_that_do_the_running(self):
        job = CI.split("  release-gate:", 1)[1]
        self.assertIn("needs: [unit-tests, workflow-validation, e2e-execution]", job)

    def test_the_gate_verifies_the_recorded_local_checks(self):
        job = CI.split("  release-gate:", 1)[1]
        self.assertIn("verify_local_execution_record.py", job)

    def test_every_script_the_workflow_names_exists(self):
        import re
        for script in sorted(set(re.findall(r"tests/[\w.]+\.(?:sh|py)", CI))):
            with self.subTest(script=script):
                self.assertTrue((ROOT / script).exists(), f"{script} is named by CI but absent")


class TheSkipAssertionCannotBeSatisfiedByASkippedSuite(unittest.TestCase):
    def test_a_clean_run_passes(self):
        self.assertEqual(coverage.check(all_required_ran(), max_skips=1), [])

    def test_the_three_areas_ludmil_named_are_all_required(self):
        joined = " ".join(coverage.REQUIRED)
        self.assertIn("test_low_score_extraction", joined)        # rescue
        self.assertIn("test_fragment_support_variants", joined)   # fragment support
        self.assertIn("test_competitive_realign_score", joined)   # competitive realignment

    def test_every_required_module_is_a_real_test_file(self):
        for module in sorted(coverage.REQUIRED):
            with self.subTest(module=module):
                self.assertTrue((ROOT / "tests" / f"{module}.py").exists())

    def test_every_allowed_skip_is_a_real_test_file_with_a_written_reason(self):
        for module, reason in sorted(coverage.ALLOWED_SKIPS.items()):
            with self.subTest(module=module):
                self.assertTrue((ROOT / "tests" / f"{module}.py").exists())
                self.assertGreater(len(reason), 40, "give a real reason, not a label")

    def test_a_required_module_that_did_not_run_fails(self):
        result = all_required_ran()
        result.ran_modules.discard("test_fragment_support_variants")
        problems = coverage.check(result, max_skips=1)
        self.assertTrue(any("test_fragment_support_variants" in p for p in problems))

    def test_a_skipped_required_test_fails_even_within_budget(self):
        result = all_required_ran(skipped={
            "test_competitive_realign_score": [("test_competitive_realign_score.K.t",
                                                "minimap2 not found")]})
        problems = coverage.check(result, max_skips=99)
        self.assertTrue(any("required test skipped" in p for p in problems))

    def test_a_skip_nobody_wrote_down_fails(self):
        result = all_required_ran(skipped={
            "test_something_new": [("test_something_new.K.t", "tool missing")]})
        self.assertTrue(any("unexpected skip" in p
                            for p in coverage.check(result, max_skips=99)))

    def test_an_allowed_skip_passes(self):
        module = sorted(coverage.ALLOWED_SKIPS)[0]
        result = all_required_ran(skipped={module: [(f"{module}.K.t", "Rscript not on PATH")]})
        self.assertEqual(coverage.check(result, max_skips=1), [])

    def test_too_many_skips_fails_even_if_each_was_allowed(self):
        module = sorted(coverage.ALLOWED_SKIPS)[0]
        result = all_required_ran(skipped={
            module: [(f"{module}.K.t{i}", "Rscript not on PATH") for i in range(3)]})
        self.assertTrue(any("budget" in p for p in coverage.check(result, max_skips=1)))

    def test_failures_and_errors_still_fail(self):
        self.assertTrue(coverage.check(all_required_ran(failures=[("t", "boom")]), 1))
        self.assertTrue(coverage.check(all_required_ran(errors=[("t", "boom")]), 1))

    def test_the_default_budget_is_the_size_of_the_written_down_list(self):
        # So widening the budget means adding a named, reasoned entry to ALLOWED_SKIPS --
        # a visible edit -- rather than quietly bumping a number in the workflow file.
        source = (ROOT / "tests/assert_required_coverage.py").read_text()
        self.assertIn('default=len(ALLOWED_SKIPS)', source)
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        self.assertNotIn("--max-skips", job,
                         "CI must take the written-down default, not pass its own budget")


HEADER = textwrap.dedent("""\
    # commit: {commit}
    # source_digest: {digest}
    # recorded_on: 2026-10-22
    # recorded_by: Someone <someone@example.org>
    # host: tscc-login
    # signed_off: Someone ran these and confirms the results
    """)
ROW = "{check}\t{status}\t{tests}\t{seconds}\t{log}\n"


class TheRecordedResultCannotBeFaked(unittest.TestCase):
    DIGEST = "a" * 64
    COMMIT = "b" * 40

    def build(self, rows=None, **header):
        fields = {"commit": self.COMMIT, "digest": self.DIGEST}
        text = HEADER.format(**fields)
        for key, value in header.items():
            text = "\n".join(line for line in text.splitlines()
                             if not line.startswith(f"# {key}:")) + "\n"
            if value is not None:
                text += f"# {key}: {value}\n"
        text += "check\tstatus\ttests\tduration_s\tlog\n"
        if rows is None:
            rows = [{"check": name, "status": "pass", "tests": "12",
                     "seconds": "900", "log": f"logs/{Path(name).name}.log"}
                    for name in record.REQUIRED_CHECKS]
        for row in rows:
            text += ROW.format(**row)
        return text

    def verify(self, text, digest=None):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "r.tsv"
            path.write_text(text)
            header, rows = record.parse(path)
        return record.verify(header, rows, digest or self.DIGEST)

    def test_a_good_record_passes(self):
        self.assertEqual(self.verify(self.build()), [])

    def test_the_required_checks_are_the_ones_the_runner_cannot_do(self):
        self.assertIn("tests/check_metagenome_mag_execution.sh", record.REQUIRED_CHECKS)
        self.assertIn("tests/assert_required_coverage.py", record.REQUIRED_CHECKS)
        for name in record.REQUIRED_CHECKS:
            with self.subTest(check=name):
                self.assertTrue((ROOT / name).exists())

    def test_a_record_from_a_different_source_tree_is_refused(self):
        problems = self.verify(self.build(), digest="c" * 64)
        self.assertTrue(any("source" in p for p in problems), problems)

    def test_a_missing_required_check_is_refused(self):
        rows = [{"check": "tests/assert_required_coverage.py", "status": "pass",
                 "tests": "1391", "seconds": "70", "log": "logs/x.log"}]
        problems = self.verify(self.build(rows=rows))
        self.assertTrue(any("check_metagenome_mag_execution" in p for p in problems))

    def test_a_failed_check_is_refused(self):
        rows = [{"check": name, "status": "fail", "tests": "12", "seconds": "9",
                 "log": "logs/x.log"} for name in record.REQUIRED_CHECKS]
        problems = self.verify(self.build(rows=rows))
        self.assertTrue(any("not releasable" in p for p in problems))

    def test_a_check_that_asserted_nothing_is_refused(self):
        rows = [{"check": name, "status": "pass", "tests": "0", "seconds": "9",
                 "log": "logs/x.log"} for name in record.REQUIRED_CHECKS]
        problems = self.verify(self.build(rows=rows))
        self.assertTrue(any("asserted nothing" in p for p in problems))

    def test_an_unsigned_record_is_refused(self):
        problems = self.verify(self.build(signed_off=None))
        self.assertTrue(any("signed_off" in p for p in problems))

    def test_a_record_with_no_source_digest_is_refused(self):
        problems = self.verify(self.build(source_digest=None))
        self.assertTrue(any("source_digest" in p for p in problems))

    def test_a_malformed_column_header_is_refused_rather_than_misread(self):
        text = self.build().replace("check\tstatus\ttests", "check\tresult\ttests")
        with self.assertRaises(SystemExit):
            self.verify(text)


class TheSourceDigestCoversWhatChangesARun(unittest.TestCase):
    def test_it_covers_the_workflow_the_modules_the_scripts_and_the_environments(self):
        for name in ("main.nf", "Modules", "scripts", "conda_envs", "conf", "lib", "tests"):
            with self.subTest(path=name):
                self.assertIn(name, record.SOURCE_PATHS)

    def test_it_excludes_the_record_itself_so_committing_it_does_not_invalidate_it(self):
        self.assertIn("tests/LOCAL_EXECUTION_CHECKS.tsv", record.SOURCE_EXCLUDE)

    def test_it_changes_when_a_source_file_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "Modules").mkdir()
            (root / "main.nf").write_text("workflow {}\n")
            (root / "Modules" / "a.nf").write_text("process A {}\n")
            before = record.source_digest(root)
            (root / "Modules" / "a.nf").write_text("process A { label 'x' }\n")
            self.assertNotEqual(before, record.source_digest(root))

    def test_it_does_not_change_when_only_the_record_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "tests").mkdir()
            (root / "main.nf").write_text("workflow {}\n")
            (root / "tests" / "LOCAL_EXECUTION_CHECKS.tsv").write_text("one\n")
            before = record.source_digest(root)
            (root / "tests" / "LOCAL_EXECUTION_CHECKS.tsv").write_text("two\n")
            self.assertEqual(before, record.source_digest(root))


if __name__ == "__main__":
    unittest.main()
