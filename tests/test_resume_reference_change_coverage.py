"""tests/check_resume_reference_change.sh: the execution proof of the F15 mechanism.

F15 is the finding that a reference, annotation or helper referenced only as a path string
inside a script block is neither a declared input nor part of the rendered script, so
`-resume` reuses tasks that no longer correspond to the current references.
lib/Provenance.groovy fixes it by interpolating digests into those scripts, and
tests/test_dependency_digests.py checks that the digests are there -- as source text.

Whether Nextflow's resume then behaves as intended is a property of a real run's task
hashing, and F21 asks for it explicitly ("resume/reference-change tests"). This file is
the structural coverage for that execution script: that it exists, runs a control, changes
content rather than a timestamp, checks both directions, and is wired into CI.

The execution script itself is not run here (it needs Nextflow, conda and a few minutes),
following the same convention as tests/test_e2e_execution_coverage.py.
"""
import stat
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "tests/check_resume_reference_change.sh"
SCRIPT = SCRIPT_PATH.read_text()
ASSERTIONS_PATH = ROOT / "tests/fixtures/assert_resume_reference_change.py"
ASSERTIONS = ASSERTIONS_PATH.read_text()
CI = (ROOT / ".github/workflows/ci.yml").read_text()


class ScriptExists(unittest.TestCase):
    def test_the_script_and_its_assertions_exist(self):
        self.assertTrue(SCRIPT_PATH.is_file())
        self.assertTrue(ASSERTIONS_PATH.is_file())

    def test_the_script_is_executable_by_bash_without_a_shebang_surprise(self):
        self.assertTrue(SCRIPT.startswith("#!/usr/bin/env bash"))

    def test_it_reuses_the_existing_tiny_fixtures_rather_than_new_ones(self):
        # F21 asks for CI to run on lightweight fixtures; a second fixture generator
        # would be a second thing to keep in step with the pipeline.
        self.assertIn("tests/fixtures/generate_e2e_fixtures.py", SCRIPT)


class ItRunsARealResume(unittest.TestCase):
    def test_it_runs_the_pipeline_for_real_not_in_preview(self):
        self.assertIn("nextflow run", SCRIPT)
        self.assertNotIn("-preview", SCRIPT)

    def test_it_resumes_into_the_same_work_directory(self):
        # A resume against a fresh work dir would cache nothing and prove nothing.
        self.assertIn('-work-dir "$work/work"', SCRIPT)
        self.assertIn("-resume", SCRIPT)

    def test_it_runs_a_control_resume_with_nothing_changed(self):
        self.assertIn("run_pipeline unchanged -resume", SCRIPT)

    def test_it_changes_content_not_a_timestamp(self):
        # Provenance.code() hashes content. A `touch` would not exercise it, and a
        # metadata-only check would pass for the wrong reason.
        self.assertNotIn("touch ", SCRIPT)
        self.assertIn(">> \"$annotation\"", SCRIPT)

    def test_it_edits_a_private_copy_not_the_repository_reference(self):
        self.assertIn('cp "$repo_dir/ref/annotations/IHE3034.clbA-clbS.gff" "$annotation"',
                       SCRIPT)
        self.assertIn('--pks_genome_annotation "$annotation"', SCRIPT)

    def test_it_captures_a_trace_per_run(self):
        self.assertIn('-with-trace "$work/trace_${label}.txt"', SCRIPT)


class TheAssertionsCheckBothDirections(unittest.TestCase):
    def test_the_control_requires_everything_cached(self):
        self.assertIn("resume with nothing changed re-ran", ASSERTIONS)

    def test_a_carrier_that_stayed_cached_is_a_failure(self):
        self.assertIn("stayed CACHED after the", ASSERTIONS)
        self.assertIn("which is exactly F15", ASSERTIONS)

    def test_an_unrelated_process_that_re_ran_is_also_a_failure(self):
        # Otherwise "re-ran everything" would pass, and the digest would be doing nothing.
        self.assertIn("digest is invalidating more than it should", ASSERTIONS)
        self.assertIn("UPSTREAM_PROCESSES", ASSERTIONS)

    def test_the_affected_processes_are_discovered_from_the_modules_not_hardcoded(self):
        self.assertIn("def processes_carrying(", ASSERTIONS)
        self.assertIn('needle = f"dep_digest?.{digest}"', ASSERTIONS)

    def test_it_resolves_include_aliases_before_matching_the_trace(self):
        # main.nf imports pksProfiler_align as pksProfilerAlign, and the trace records
        # the alias. Matching the declared name alone finds nothing and reports a
        # false "this check cannot conclude anything".
        self.assertIn("def include_aliases(", ASSERTIONS)
        self.assertIn("carriers |= aliases.get(name, {name})", ASSERTIONS)
        main = (ROOT / "main.nf").read_text()
        self.assertIn("pksProfiler_align as pksProfilerAlign", main)

    def test_a_run_with_no_carrier_at_all_is_a_failure_not_a_pass(self):
        self.assertIn("the check cannot conclude anything", ASSERTIONS)

    def test_the_digest_it_uses_is_actually_carried_by_a_module(self):
        # The default --digest must correspond to a real dependency, or the execution
        # test would fail for a reason that has nothing to do with resume.
        modules = "\n".join(p.read_text() for p in (ROOT / "Modules").glob("*.nf"))
        self.assertIn("dep_digest?.pks_annotation", modules)


class CIWiring(unittest.TestCase):
    def test_ci_runs_the_resume_check(self):
        self.assertIn("check_resume_reference_change.sh", CI)

    def test_it_runs_in_the_scheduled_execution_job(self):
        job = CI[CI.index("e2e-execution:"):]
        self.assertIn("check_resume_reference_change.sh", job)
        self.assertIn("github.event_name == 'schedule'", job)
        self.assertIn("github.event_name == 'push' && github.ref == 'refs/heads/main'", job)


if __name__ == "__main__":
    unittest.main()
