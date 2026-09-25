"""A real, executing test for malformed input exists (F21's third scenario).

F21 (Ludmil's audit): "...add a small reference-backed integration workflow covering
BAM, CRAM, single and paired FASTQ, zero, positive and malformed inputs." Every other
real-execution check in this suite proves the pipeline computes the right answer on
good input; tests/check_malformed_input_execution.sh is the one that proves it fails
visibly and informatively on bad input -- a truncated BAM and a BAM with corrupted
leading magic bytes, run alongside a healthy sample in one cohort -- rather than
hanging, crashing without explanation, or silently reporting a failed sample as though
it had succeeded.

This file is the static half of tests/check_malformed_input_execution.sh,
tests/fixtures/generate_malformed_input_fixtures.py and
tests/fixtures/assert_malformed_input_outputs.py -- the same relationship
test_e2e_execution_coverage.py has to tests/check_e2e_execution.sh. The empirical proof
is the scripts themselves, run live:

    CONDA_CACHE_DIR=/path/to/cache bash tests/check_malformed_input_execution.sh

Not run from here: real conda environments (minimap2, samtools, bowtie2, fastp).
"""
import stat
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tests/check_malformed_input_execution.sh"
GENERATOR = ROOT / "tests/fixtures/generate_malformed_input_fixtures.py"
ASSERTIONS = ROOT / "tests/fixtures/assert_malformed_input_outputs.py"
EXTRACT_READS = (ROOT / "Modules/extract_reads.nf").read_text()
CI = (ROOT / ".github/workflows/ci.yml").read_text()


def code_only(text):
    return "\n".join(line for line in text.splitlines() if not line.strip().startswith("#"))


class ScriptsExistAndAreExecutable(unittest.TestCase):
    def test_the_three_files_exist(self):
        for path in (SCRIPT, GENERATOR, ASSERTIONS):
            self.assertTrue(path.is_file(), f"missing: {path}")

    def test_they_are_executable(self):
        for path in (SCRIPT, GENERATOR, ASSERTIONS):
            mode = path.stat().st_mode
            self.assertTrue(mode & stat.S_IXUSR, f"not executable: {path}")


class RealExecutionNotAPreview(unittest.TestCase):
    def setUp(self):
        self.script = SCRIPT.read_text()
        self.code = code_only(self.script)

    def test_preview_flag_is_never_used_in_the_actual_invocation(self):
        self.assertNotIn("-preview", self.code)

    def test_it_actually_invokes_nextflow_run(self):
        self.assertIn("nextflow run", self.script)

    def test_it_is_bounded_by_a_timeout_so_a_real_hang_cannot_block_forever(self):
        self.assertIn("timeout", self.code)
        self.assertIn("NEXTFLOW_TIMEOUT", self.script)

    def test_a_timeout_is_reported_as_its_own_distinct_failure_mode(self):
        # A hang and an ordinary nonzero exit are different bugs; conflating them would
        # hide the exact failure mode this check exists to catch.
        self.assertIn("-eq 124", self.script)
        self.assertIn("HUNG", self.script)


class MalformedFixturesMatchTheRealGuardClauses(unittest.TestCase):
    """The two malformed fixtures must actually exercise extractReads.nf's two
    distinct guard clauses, not one generic "broken file" thrown at both."""

    def setUp(self):
        self.generator = GENERATOR.read_text()

    def test_truncation_keeps_the_leading_magic_bytes_intact(self):
        # Cuts the tail (`valid_bytes[:keep]`), not the head -- htsfile must still see
        # BAM so the run reaches samtools quickcheck, not the format-detection branch.
        self.assertIn("valid_bytes[:keep]", self.generator)

    def test_bad_magic_corrupts_only_the_leading_bytes(self):
        self.assertIn("corrupted[0:4]", self.generator)

    def test_the_control_sample_is_real_non_truncated_input(self):
        # Proves the corruption -- not the fixture's own construction -- is what fails;
        # a fixture that could never have worked wouldn't demonstrate that distinction.
        self.assertIn("valid.bam", self.generator)
        self.assertIn("check=True", self.generator)  # samtools view -bS must itself succeed


class ExtractReadsHasTheGuardClausesTheseFixturesAreBuiltAgainst(unittest.TestCase):
    """If either guard clause in extract_reads.nf disappears, these fixtures would stop
    testing anything real; pin their existence here too."""

    def test_htsfile_format_detection_exists(self):
        self.assertIn("htsfile", EXTRACT_READS)
        self.assertIn("neither BAM nor CRAM", EXTRACT_READS)

    def test_quickcheck_runs_before_decode(self):
        code = code_only(EXTRACT_READS)
        self.assertIn("samtools quickcheck", code)
        # quickcheck must run before samtools fastq, or a truncated file would reach
        # decode first and fail with a less specific error.
        self.assertLess(code.index("samtools quickcheck"), code.index("samtools fastq"))

    def test_the_script_runs_under_set_dash_e_so_quickcheck_failure_is_fatal(self):
        self.assertIn("set -euo pipefail", EXTRACT_READS)


class AssertionsCheckPerSampleFailureIsolationAndInformativeRecording(unittest.TestCase):
    """params.sample_failure_strategy defaults to 'ignore' (main.nf): a malformed
    sample must fail visibly, not silently, and not take the healthy sample with it --
    the live half of Ludmil's finding-8 acceptance test
    (test_sample_stage_failure_policy.py's own docstring says this needs a real run)."""

    def setUp(self):
        self.assertions = ASSERTIONS.read_text()

    def test_it_checks_the_healthy_sample_completes_normally(self):
        self.assertIn('"complete"', self.assertions)

    def test_it_checks_the_malformed_samples_are_not_reported_as_complete(self):
        self.assertIn('!= "complete"', self.assertions)

    def test_it_checks_no_fabricated_tier_output_for_a_failed_sample(self):
        self.assertIn("read_evidence.tsv", self.assertions)
        self.assertIn("not read_evidence.exists()", self.assertions)

    def test_it_reads_run_report_what_failed_section(self):
        self.assertIn("RUN_REPORT.txt", self.assertions)
        self.assertIn("WHAT FAILED", self.assertions)

    def test_it_checks_both_malformed_samples_are_recorded_as_separate_failures(self):
        self.assertIn(">= 2", self.assertions)

    def test_it_checks_the_real_specific_diagnostic_text_not_just_a_nonzero_exit(self):
        # A failure could be nonzero for the wrong reason; this must confirm the actual
        # tool message reached the log, not merely that *something* failed.
        self.assertIn("EOF block", self.assertions)
        self.assertIn("neither BAM nor CRAM", self.assertions)


class CIWiring(unittest.TestCase):
    def test_ci_references_the_malformed_input_script(self):
        self.assertIn("check_malformed_input_execution.sh", CI)

    def test_it_is_not_on_every_push_or_pull_request(self):
        # Same reasoning as check_e2e_execution.sh: no persistent conda cache on a
        # GitHub-hosted runner, gated to manual dispatch until that cost is measured.
        script_ref = CI[CI.index("check_malformed_input_execution.sh"):]
        preceding = CI[:CI.index("check_malformed_input_execution.sh")]
        self.assertIn("workflow_dispatch", preceding + script_ref)


if __name__ == "__main__":
    unittest.main()
