"""A real, executing end-to-end test exists and is not accidentally a third -preview.

Ludmil's audit, "not yet confirmed" item u-3: "The test suite is much stronger, but
several historical failures only appeared when Nextflow channel wiring and external
tools interacted. Add one tiny positive and negative reference-backed run covering
BAM, CRAM, single FASTQ, paired FASTQ, HMM, and cohort aggregation. This should also
test the mate-pair check above [u-1]."

Everything in this file is a static check on tests/check_e2e_execution.sh,
tests/fixtures/generate_e2e_fixtures.py and tests/fixtures/assert_e2e_outputs.py --
the same kind of source-pattern check every other test_*.py here already is. The
empirical proof that the real run behaves as these scripts claim is the scripts
themselves, run live: `bash tests/check_e2e_execution.sh`. That is deliberately not
run from here -- unlike tests/check_workflow_construction.sh, it needs real conda
environments (bowtie2, minimap2, fastp, samtools) and takes about two minutes even
with every environment already cached, so it is not part of the fast tier
`python -m unittest discover` runs, and not part of tests/run_checks.sh either. See
.github/workflows/ci.yml for how it is (and, so far, is not) wired into CI.

Scope note: the first pass (commit 73ed947) covered paired FASTQ and BAM input,
--sample_type tumor_wgs, --profiling_method bowtie2, positive and negative fixtures,
and the u-1 mate-suffix survival check through host depletion. This follow-up closes
the rest of Ludmil's wish list: CRAM input (extractReads' --cram_reference branch,
validate_cram_reference.py's @SQ/MD5 check), --profiling_method hmm alone
(pksProfilerHMM/hmm_best_hit.py, real nhmmscan hits and a real zero-hit run), and a
5-sample cohort spanning broad_island/extensive_island/localized_indeterminate on top
of the original multi_gene/negative pair, with exact (not just "both are present")
cohort-membership checks on the roll-up tables.
"""
import stat
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tests/check_e2e_execution.sh"
GENERATOR = ROOT / "tests/fixtures/generate_e2e_fixtures.py"
ASSERTIONS = ROOT / "tests/fixtures/assert_e2e_outputs.py"
CI = (ROOT / ".github/workflows/ci.yml").read_text()


class ScriptsExistAndAreExecutable(unittest.TestCase):
    def test_the_three_scripts_exist(self):
        for path in (SCRIPT, GENERATOR, ASSERTIONS):
            self.assertTrue(path.is_file(), f"missing: {path}")

    def test_they_are_executable(self):
        for path in (SCRIPT, GENERATOR, ASSERTIONS):
            mode = path.stat().st_mode
            self.assertTrue(mode & stat.S_IXUSR, f"not executable: {path}")


def code_only(text):
    """Strip full-line comments. The script's own header deliberately says '-preview'
    several times, contrasting itself with tests/check_workflow_construction.sh -- the
    thing that must never contain that flag is the actual invocation, not the prose."""
    return "\n".join(line for line in text.splitlines() if not line.strip().startswith("#"))


class RealExecutionNotAThirdPreview(unittest.TestCase):
    """The entire point of u-3 is that this run must not be -preview again."""

    def setUp(self):
        self.script = SCRIPT.read_text()
        self.code = code_only(self.script)

    def test_preview_flag_is_never_used_in_the_actual_invocation(self):
        self.assertNotIn("-preview", self.code)

    def test_it_actually_invokes_nextflow_run(self):
        self.assertIn("nextflow run", self.script)

    def test_it_does_not_stand_in_placeholder_databases(self):
        # tests/check_workflow_construction.sh's `need()`/stand-in-file pattern is
        # correct there (it never reads the databases), but would silently defeat the
        # point here -- a run whose --hg38_db is a one-byte placeholder never really
        # depletes anything, and a placeholder minimap2 index would fail to load
        # rather than deplete zero reads, so the checks below could not tell "the
        # index is fake" from "depletion correctly removed nothing".
        self.assertNotIn("printf 'x' >", self.script)

    def test_it_builds_real_minimap2_indices(self):
        self.assertIn("-d", self.script)
        self.assertIn(".mmi", self.script)
        self.assertIn("minimap2", self.script)

    def test_it_requires_sample_type_and_a_real_profiling_method(self):
        self.assertIn("--sample_type tumor_wgs", self.script)
        self.assertIn("--profiling_method bowtie2", self.script)

    def test_it_publishes_intermediates_so_host_depletion_output_is_checkable(self):
        self.assertIn("--save_intermediates true", self.script)


class CramInputExecutionCoverage(unittest.TestCase):
    """u-3 follow-up: CRAM was the one input form the original run left uncovered."""

    def setUp(self):
        self.script = SCRIPT.read_text()

    def test_it_converts_the_bam_fixtures_to_cram_with_real_samtools(self):
        self.assertIn("view -C -T", self.script)

    def test_it_passes_cram_reference(self):
        self.assertIn("--input_data_type cram", self.script)
        self.assertIn("--cram_reference", self.script)

    def test_cram_reference_is_the_same_fasta_the_cram_was_built_against(self):
        # extractReads.nf's validate_cram_reference.py compares the CRAM's embedded
        # @SQ M5 against `samtools dict` of whatever --cram_reference names; passing
        # anything other than the exact FASTA samtools view -C -T was given would be
        # asking that check to fail, not exercising the success path.
        self.assertIn('--cram_reference "$fixtures/synthetic_hg38.fa"', self.script)

    def test_it_checks_cram_output_with_the_same_f01_metric_as_bam(self):
        cram_section = self.script[self.script.index("Third input form: CRAM"):]
        self.assertIn("--expect-extracted-unmapped-reads 8", cram_section)


class HmmProfilingExecutionCoverage(unittest.TestCase):
    """u-3 follow-up: --profiling_method hmm had never been executed by any test."""

    def setUp(self):
        self.script = SCRIPT.read_text()
        self.assertions = ASSERTIONS.read_text()

    def test_it_runs_profiling_method_hmm(self):
        self.assertIn("--profiling_method hmm", self.script)

    def test_it_uses_sample_type_metagenome_not_tumor_wgs(self):
        # main.nf: "--sample_type tumor_wgs requires --profiling_method bowtie2 or
        # both" -- tumor_wgs rejects an hmm-alone run outright, so this must be a
        # different sample_type or the run never reaches pksProfilerHMM at all. The
        # code (not the prose above it, which also says "--profiling_method hmm")
        # is what has to pair it with metagenome.
        code = code_only(self.script)
        hmm_section = code[code.index("--profiling_method hmm"):]
        self.assertIn("--sample_type metagenome", hmm_section[:hmm_section.index("outdir")])

    def test_assertions_support_an_hmm_profiling_mode(self):
        self.assertIn("--profiling-method", self.assertions)
        self.assertIn("check_hmm_evidence", self.assertions)

    def test_assertions_check_the_real_hmm_output_files(self):
        self.assertIn("hmm_counts.tsv", self.assertions)
        self.assertIn("hmm.qc.tsv", self.assertions)

    def test_assertions_check_ambiguity_is_reported_separately_from_the_assigned_count(self):
        # hmm_best_hit.py's F07 fix: ambiguous reads (a tie across two genes) count for
        # neither gene and are reported on their own metric, not folded into either.
        self.assertIn("hmm_ambiguous_reads", self.assertions)


class CohortScalingExecutionCoverage(unittest.TestCase):
    """u-3 follow-up: a 2-sample cohort cannot show a sorting/dedup/boundary bug that
    only appears at N>2; this run's 5-sample cohort can."""

    def setUp(self):
        self.script = SCRIPT.read_text()
        self.assertions = ASSERTIONS.read_text()

    def test_it_runs_a_five_sample_cohort(self):
        for var in ("POS_SAMPLE", "NEG_SAMPLE", "BROAD_SAMPLE", "EXTENSIVE_SAMPLE", "BORDERLINE_SAMPLE"):
            self.assertIn(var, self.script)

    def test_it_passes_cohort_samples_to_the_assertions(self):
        self.assertIn("--cohort-samples", self.script)

    def test_assertions_support_naming_extra_cohort_samples_and_their_tiers(self):
        self.assertIn("parse_cohort_samples", self.assertions)
        for tier in ("broad_island", "extensive_island", "localized_indeterminate"):
            self.assertIn(tier, self.assertions)

    def test_assertions_check_exact_cohort_membership_not_just_containment(self):
        # The original 2-sample check only asked "are both of these somewhere in the
        # table" -- true even if a sample were duplicated or a stray row present. The
        # cohort-scaling assertion must ask for the exact set, once each.
        self.assertIn("check_cohort_membership", self.assertions)
        self.assertIn("strict", self.assertions)

    def test_assertions_check_the_gene_by_sample_matrix_columns_too(self):
        # masterTableAlign's pks.gene.counts.align.txt is Gene x Sample -- membership
        # there is column headers, not rows, and the original check never looked at it.
        self.assertIn("pks.gene.counts.align.txt", self.assertions)


class FixtureGeneratorBuildsPositiveAndNegativeCases(unittest.TestCase):
    def setUp(self):
        self.generator = GENERATOR.read_text()

    def test_positive_reads_come_from_the_real_shipped_reference(self):
        self.assertIn("GCF_000025745.1", self.generator)

    def test_positive_reads_span_multiple_clb_genes(self):
        # multi_gene tier requires >=3 distinct genes (classify_tumor_pks_evidence.py);
        # the fixture is built from 4 so a small alignment loss doesn't flake the tier.
        self.assertIn("LOCI", self.generator)
        genes_mentioned = sum(1 for gene in ("clbS", "clbN", "clbD", "clbA") if gene in self.generator)
        self.assertGreaterEqual(genes_mentioned, 3)

    def test_negative_reads_are_independent_pseudorandom_sequence(self):
        self.assertIn("build_negative", self.generator)
        self.assertIn("random", self.generator.lower())

    def test_synthetic_host_references_are_built_not_shipped(self):
        # --hg38_db/--t2t_phix_db have no small real equivalent in this repo (see the
        # task's own framing); the fixture must synthesize them rather than silently
        # reuse something else in indices/ that happens to exist.
        self.assertIn("build_synthetic_host", self.generator)

    def test_read_names_have_no_pre_existing_mate_suffix(self):
        # filterReads (Modules/filter_reads.nf) adds /1 and /2 itself; pre-suffixing
        # the fixture would let that step go untested, which is exactly what u-1 is
        # about (Modules/filter_reads.nf, the mate-suffix-tagging awk block).
        self.assertIn("WITHOUT a /1 or /2 suffix", self.generator)


class FixtureGeneratorBuildsTheScaledCohort(unittest.TestCase):
    """u-3 follow-up: broad_island/extensive_island/localized_indeterminate on top of
    the original multi_gene/negative pair, for a 5-sample cohort."""

    def setUp(self):
        self.generator = GENERATOR.read_text()

    def test_it_defines_a_second_validated_window_pool(self):
        self.assertIn("TIER_WINDOW_POOL", self.generator)

    def test_it_builds_broad_island_extensive_island_and_localized_indeterminate(self):
        for name in ("BROAD_WINDOWS", "EXTENSIVE_WINDOWS", "BORDERLINE_WINDOWS"):
            self.assertIn(name, self.generator)

    def test_extensive_pads_read_count_with_real_duplicate_fragments_not_new_positions(self):
        # Padding the read-count floor with an unvalidated new genomic position would
        # reintroduce exactly the multi-mapping risk clbK was dropped over; a
        # duplicate read at an already-validated position cannot.
        self.assertIn("EXTENSIVE_DUPLICATE_WINDOWS", self.generator)

    def test_it_uses_a_shared_builder_not_copy_pasted_per_tier_logic(self):
        self.assertIn("def build_tier_sample", self.generator)


class AssertionsCoverTheMatePairRegressionGuard(unittest.TestCase):
    """Ludmil: "This should also test the mate-pair check above [u-1]."."""

    def setUp(self):
        self.assertions = ASSERTIONS.read_text()

    def test_it_reads_the_real_host_depleted_fastq(self):
        self.assertIn("host_depleted.fastq.gz", self.assertions)

    def test_it_checks_both_mate_suffixes_survive(self):
        self.assertIn("mate_suffix_census", self.assertions)
        self.assertIn('"1", "2"', self.assertions.replace("'", '"'))

    def test_it_checks_the_tumour_read_evidence_tier(self):
        self.assertIn("read_evidence", self.assertions)
        self.assertIn("multi_gene", self.assertions)
        self.assertIn("negative", self.assertions)

    def test_it_checks_the_cohort_qc_summary_not_just_the_per_sample_files(self):
        self.assertIn("pks.qc.summary.tsv", self.assertions)

    def test_it_checks_cohort_aggregation_output(self):
        # Ludmil's wish list explicitly names "cohort aggregation"; this run only has
        # two samples, but the roll-up tables must still carry both of them.
        self.assertIn("pks_cohort_report.tsv", self.assertions)
        self.assertIn("pks.master_summary.tsv", self.assertions)


class CIWiring(unittest.TestCase):
    """Not wired into the fast, always-on jobs -- see the module docstring -- but
    present in ci.yml so it is at least runnable without local setup."""

    def test_ci_references_the_e2e_script(self):
        self.assertIn("check_e2e_execution.sh", CI)

    def test_the_e2e_job_is_not_on_every_push_or_pull_request(self):
        # A cold GitHub-hosted runner has no persistent conda cache (unlike
        # CONDA_CACHE_DIR on TSCC) and would solve bowtie2/minimap2/fastp/samtools/
        # subread/deeptools/R+circlize from scratch -- an unknown, possibly large,
        # multiple of the ~2 minutes this takes locally with every env pre-built.
        # Gated to manual dispatch until that cold-start cost is actually measured.
        e2e_job = CI[CI.index("check_e2e_execution.sh"):]
        preceding = CI[:CI.index("check_e2e_execution.sh")]
        self.assertIn("workflow_dispatch", preceding + e2e_job)


if __name__ == "__main__":
    unittest.main()
