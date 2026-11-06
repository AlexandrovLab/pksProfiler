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

    def test_a_failed_run_never_saves_its_conda_environments(self):
        """actions/cache saves in a post step even when the job failed.

        That kept the v1.0.0 builds red for several rounds: a failed run saved its
        half-built environments, the next run restored them, `conda activate` on a
        partial environment prepended a bin directory with no python in it, and the
        lookup fell through to base's python -- reported as a PATH problem when the
        real fault was the environment not matching its yml. Bumping the cache key did
        not help, because the next failure re-poisoned the new key.

        So: restore and save are separate steps, and the save is guarded by success().
        """
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        self.assertNotIn("uses: actions/cache@", job,
                         "the combined cache action saves on failure too; use "
                         "cache/restore plus a cache/save guarded by success()")
        self.assertIn("actions/cache/restore@", job, "the job should restore a cache")
        save = job.split("actions/cache/save@", 1)
        self.assertEqual(len(save), 2, "the job should save the cache in its own step")
        # the guard must be on the save step, in the few lines before it
        self.assertIn("if: success()", save[0][-400:],
                      "the cache save must be guarded by success()")

    def test_a_tag_runs_the_gate_and_the_gate_verifies_the_recorded_checks(self):
        """rc2 review 2.2, under the v1.0.0 arrangement.

        2.2 asks that the real execution checks run on the exact tag. For v1.0.0 they do
        NOT run on the runner: e2e-execution has never passed on a GitHub-hosted image
        since October, for reasons specific to that image rather than to the pipeline,
        and it is dispatch-only (see its comment in ci.yml). The review's own fallback
        applies instead -- a committed, signed-off local result whose source digest this
        gate checks against the tagged tree. That substitution is stated in
        RELEASE_NOTES_v1.0.0.md and in the 6 Nov status report, not left implicit.

        What must therefore hold: a tag runs the gate, and the gate verifies the record.
        """
        self.assertIn('tags: ["v*"]', CI)
        gate = CI.split("  release-gate:", 1)[1]
        self.assertIn("startsWith(github.ref, 'refs/tags/v')", gate.split("steps:", 1)[0])
        self.assertIn("verify_local_execution_record.py", gate)

    def test_the_gate_does_not_depend_on_a_job_that_cannot_run(self):
        """A `needs:` on a skipped job skips the gate, leaving the tag verified by nothing.

        e2e-execution is dispatch-only for v1.0.0, so it must not be a dependency of the
        gate; the two jobs that do run on a tag must be.
        """
        needs = CI.split("  release-gate:", 1)[1].split("steps:", 1)[0]
        self.assertIn("needs: [unit-tests, workflow-validation]", needs)
        self.assertNotIn("e2e-execution]", needs)

    def test_the_gate_says_the_execution_checks_did_not_run_here(self):
        """Silence would read as though they had. The summary has to say otherwise."""
        gate = CI.split("  release-gate:", 1)[1]
        self.assertIn("did NOT run on this runner", gate)

    def test_the_execution_job_survives_for_the_next_release(self):
        """Dispatch-only, not deleted: v1.0.1 should fix it, which needs it to still exist."""
        self.assertIn("  e2e-execution:", CI)
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        self.assertIn("workflow_dispatch", job.split("runs-on:", 1)[0])

    def test_the_execution_job_runs_every_execution_check_the_runner_can(self):
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        for script in ("check_e2e_execution.sh", "check_malformed_input_execution.sh",
                       "check_mode_execution.sh", "check_v8_lanes_execution.sh",
                       "check_resume_reference_change.sh"):
            with self.subTest(script=script):
                self.assertIn(script, job)

    def test_the_job_installs_every_external_tool_its_checks_resolve(self):
        """The v1.0.0 tag build failed because the job installed minimap2 but not samtools.

        The execution checks need two tools that a stock runner does not carry and that
        Nextflow does not provide for them: minimap2, to build the .mmi fixtures, and
        samtools, to sort and index the BAM fixtures. Nextflow builds its own copies
        inside the pipeline's conda envs, but the fixture builders run BEFORE the
        pipeline and resolve these from the environment. Installing one and not the
        other fails only at tag time, in the job that nothing else exercises.
        """
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        install = [l for l in job.splitlines()
                   if "conda install" in l or "conda create" in l]
        self.assertTrue(install, "the job installs no fixture tools at all")
        installed = " ".join(install)
        for tool in ("minimap2", "samtools"):
            with self.subTest(tool=tool):
                self.assertIn(tool, installed,
                              f"{tool} is resolved by tests/check_*.sh before the pipeline "
                              f"runs, so the job has to install it")

    def test_the_resolved_tool_paths_do_not_depend_on_PATH_surviving_the_step(self):
        """Each `run:` is its own login shell; base-env PATH did not carry over at v1.0.0."""
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        self.assertIn("MINIMAP2_BIN=", job,
                      "resolve the binary once and export it, rather than trusting PATH")
        self.assertIn("GITHUB_ENV", job,
                      "the exported path has to reach the later steps via GITHUB_ENV")

    def test_fixture_tools_are_not_installed_into_the_base_environment(self):
        """Installing them into base is what shadowed every task's pinned environment.

        setup-miniconda auto-activates base in every `shell: bash -l` step, so anything
        in base lands on PATH, is inherited by Nextflow, and sits ahead of each task's
        own conda environment. samtools also brings bgzip via htslib, so one `-n base`
        install displaced three pinned tools at once and --tool_version_check strict
        refused the run. Two v1.0.0 tags went red before this was understood.
        """
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        for line in job.splitlines():
            if "conda install" in line or "conda create" in line:
                with self.subTest(line=line.strip()[:70]):
                    self.assertNotIn("-n base", line,
                                     "fixture tools must go in their own environment")
                    self.assertNotIn("--name base", line,
                                     "fixture tools must go in their own environment")

    def test_the_job_never_prepends_anything_to_PATH_for_later_steps(self):
        """Adding to PATH shadows every task's pinned conda environment; removing is fine.

        The second v1.0.0 tag attempt wrote `PATH=<miniconda>/bin:$PATH` to GITHUB_ENV so
        later steps could find minimap2. Nextflow tasks then resolved minimap2, samtools
        and bgzip from the base environment instead of their own pinned ones -- samtools
        1.21 against a pinned 1.24 -- and --tool_version_check strict refused the run.

        The distinction that matters is prepending versus stripping, not whether PATH is
        exported at all: the job legitimately exports a PATH with the stray interpreters
        REMOVED. A prepend keeps `$PATH` on the right-hand side; a strip does not.
        """
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        for line in job.splitlines():
            stripped = line.strip()
            if "GITHUB_ENV" not in stripped or "PATH=" not in stripped:
                continue
            if "$PATH" in stripped.split("PATH=", 1)[1]:
                self.fail(f"the job prepends to PATH for later steps, which shadows "
                          f"every task's pinned environment: {stripped}")

    def test_the_execution_job_runs_every_execution_check_the_runner_can(self):
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        for script in ("check_e2e_execution.sh", "check_malformed_input_execution.sh",
                       "check_mode_execution.sh", "check_v8_lanes_execution.sh",
                       "check_resume_reference_change.sh"):
            with self.subTest(script=script):
                self.assertIn(script, job)

    def test_the_job_installs_every_external_tool_its_checks_resolve(self):
        """The v1.0.0 tag build failed because the job installed minimap2 but not samtools.

        The execution checks need two tools that a stock runner does not carry and that
        Nextflow does not provide for them: minimap2, to build the .mmi fixtures, and
        samtools, to sort and index the BAM fixtures. Nextflow builds its own copies
        inside the pipeline's conda envs, but the fixture builders run BEFORE the
        pipeline and resolve these from the environment. Installing one and not the
        other fails only at tag time, in the job that nothing else exercises.
        """
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        install = [l for l in job.splitlines()
                   if "conda install" in l or "conda create" in l]
        self.assertTrue(install, "the job installs no fixture tools at all")
        installed = " ".join(install)
        for tool in ("minimap2", "samtools"):
            with self.subTest(tool=tool):
                self.assertIn(tool, installed,
                              f"{tool} is resolved by tests/check_*.sh before the pipeline "
                              f"runs, so the job has to install it")

    def test_the_resolved_tool_paths_do_not_depend_on_PATH_surviving_the_step(self):
        """Each `run:` is its own login shell; base-env PATH did not carry over at v1.0.0."""
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        self.assertIn("MINIMAP2_BIN=", job,
                      "resolve the binary once and export it, rather than trusting PATH")
        self.assertIn("GITHUB_ENV", job,
                      "the exported path has to reach the later steps via GITHUB_ENV")

    def test_fixture_tools_are_not_installed_into_the_base_environment(self):
        """Installing them into base is what shadowed every task's pinned environment.

        setup-miniconda auto-activates base in every `shell: bash -l` step, so anything
        in base lands on PATH, is inherited by Nextflow, and sits ahead of each task's
        own conda environment. samtools also brings bgzip via htslib, so one `-n base`
        install displaced three pinned tools at once and --tool_version_check strict
        refused the run. Two v1.0.0 tags went red before this was understood.
        """
        job = CI.split("  e2e-execution:", 1)[1].split("\n  release-gate:", 1)[0]
        for line in job.splitlines():
            if "conda install" in line or "conda create" in line:
                with self.subTest(line=line.strip()[:70]):
                    self.assertNotIn("-n base", line,
                                     "fixture tools must go in their own environment")
                    self.assertNotIn("--name base", line,
                                     "fixture tools must go in their own environment")

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

    WAIVER = ("conda cache quota blocks rebuilding one env; the check last passed "
              "2026-10-22 and nothing it exercises changed since")

    def skipped_rows(self):
        return [{"check": "tests/assert_required_coverage.py", "status": "pass",
                 "tests": "1423", "seconds": "74", "log": "logs/a.log"},
                {"check": "tests/check_metagenome_mag_execution.sh", "status": "skipped",
                 "tests": "0", "seconds": "0", "log": "-"}]

    def test_a_skip_nobody_waived_is_refused(self):
        problems = self.verify(self.build(rows=self.skipped_rows()))
        self.assertTrue(any("skip_waiver" in p for p in problems), problems)

    def test_a_waived_skip_of_the_mag_check_is_accepted(self):
        record_text = self.build(
            rows=self.skipped_rows(),
            skip_waiver=f"tests/check_metagenome_mag_execution.sh: {self.WAIVER}")
        self.assertEqual(self.verify(record_text), [])

    def test_a_flippant_waiver_is_refused(self):
        record_text = self.build(
            rows=self.skipped_rows(),
            skip_waiver="tests/check_metagenome_mag_execution.sh: because")
        problems = self.verify(record_text)
        self.assertTrue(any("substantive" in p for p in problems), problems)

    def test_the_unit_gate_can_never_be_waived(self):
        rows = [{"check": "tests/assert_required_coverage.py", "status": "skipped",
                 "tests": "0", "seconds": "0", "log": "-"},
                {"check": "tests/check_metagenome_mag_execution.sh", "status": "pass",
                 "tests": "25", "seconds": "64", "log": "logs/m.log"}]
        record_text = self.build(
            rows=rows, skip_waiver=f"tests/assert_required_coverage.py: {self.WAIVER}")
        problems = self.verify(record_text)
        self.assertTrue(any("unit gate" in p for p in problems), problems)

    def test_a_waived_skip_is_still_not_a_pass_for_a_failed_check(self):
        # a fail with a waiver header is still a fail
        rows = [{"check": name, "status": "fail", "tests": "0", "seconds": "9",
                 "log": "logs/x.log"} for name in record.REQUIRED_CHECKS]
        record_text = self.build(
            rows=rows,
            skip_waiver=f"tests/check_metagenome_mag_execution.sh: {self.WAIVER}")
        problems = self.verify(record_text)
        self.assertTrue(any("not releasable" in p for p in problems), problems)

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
