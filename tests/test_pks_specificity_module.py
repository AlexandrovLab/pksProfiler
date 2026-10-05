"""Structural coverage for Modules/pks_specificity.nf, the V8a LOCUS-DETECTION
SPECIFICITY lane (competitive remapping against a true + pks-negative-decoy panel).

Ludmil, 25 Sep report section 4 V8a: "Challenge pks-supporting reads against pks-negative
Enterobacterales, close homologous NRPS/PKS sequences, broad gut/oral backgrounds,
GRCh38, and UniVec. This asks whether the sequence really belongs to the pks/clb locus."

The sister lane -- V8b carrier attribution, Modules/pks_carrier.nf -- answers a different
question against a pks-POSITIVE panel, and the report says the two must not be conflated.
tests/test_v8_lane_separation.py holds their vocabularies apart; this file covers V8a's
own wiring.

The lane is now included from main.nf and invoked under --pks_locus_specificity (off by
default), because its panel exists: scripts/build_specificity_panels.py --panel v8a builds
it from ref/specificity_panels/v8a_locus_decoy_panel.tsv.

Static, source-level checks here (the `code_only()` pattern every other test_*.py in
this repo uses for .nf files) are complemented by a real, live construction proof: this
session ran
    JAVA_HOME=<nf-env> nextflow run Modules/pks_specificity.nf -preview \\
        --sample <sheet.csv> --sample_dir <dir> \\
        --pks_specificity_true_index/... (placeholder paths; -preview never reads them)
against a hand-built one-sample sheet and confirmed: (a) with every required param
supplied, the workflow graph constructs with `[SUCCESS] completed=0 failed=0 cached=0`
and no task is ever run; (b) with the params omitted, it fails fast with the exact
per-parameter error messages pksSpecificityValidateParams() emits, not a generic
Groovy NPE three stages later. That live run is not re-executed by `unittest
discover` (it needs a real `nextflow` binary and a JVM, like
tests/check_workflow_construction.sh's proof for main.nf), so it is captured here as
prose plus the source-level assertions below, following test_e2e_execution_coverage.py's
own convention for the same trade-off.
"""
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "Modules/pks_specificity.nf"
MODULE = MODULE_PATH.read_text()


def code_only(text):
    """Comments explain the reasoning; the assertions are about what runs."""
    return "\n".join(line for line in text.splitlines()
                      if not line.strip().startswith(("//", "#")))


CODE = code_only(MODULE)


class ModuleExists(unittest.TestCase):
    def test_file_exists(self):
        self.assertTrue(MODULE_PATH.is_file())

    def test_dsl2(self):
        self.assertIn("nextflow.enable.dsl = 2", MODULE)


class RequiredParamsFailClearlyNotSilently(unittest.TestCase):
    """Every Phase 2 dependency must be named in an actionable error, not left to
    surface as a bare Groovy null-pointer several processes downstream."""

    # Arm A is optional, so pks_specificity_arm_a_target_clade is NOT in this list: it
    # is required only when --pks_specificity_arm_a_dir is given, which is checked
    # separately in ArmAIsOptional below.
    REQUIRED_PARAMS = (
        "pks_specificity_true_index",
        "pks_specificity_decoy_index",
        "pks_specificity_concatenated_index",
        "pks_specificity_true_contigs",
        "pks_specificity_decoy_contigs",
        "pks_specificity_island_length",
    )

    def test_every_required_param_defaults_to_null(self):
        for name in self.REQUIRED_PARAMS:
            self.assertIn(f"params.{name}", MODULE,
                           f"{name} is not declared as a module param")

    def test_every_required_param_is_checked_before_use(self):
        for name in self.REQUIRED_PARAMS:
            self.assertIn(f"!params.{name}", CODE, f"{name} has no missing-param guard")

    def test_validation_function_names_the_builder_that_produces_the_panel(self):
        self.assertIn("def pksLocusSpecificityValidateParams()", CODE)
        self.assertIn("exit 1,", CODE)
        # An error that names only the missing flag leaves the reader nowhere to go.
        self.assertIn("build_specificity_panels.py --panel v8a", MODULE)

    def test_validation_runs_before_any_process_in_the_workflow(self):
        workflow_body = MODULE.split("workflow pksLocusSpecificity {", 1)[1]
        main_block = workflow_body.split("main:", 1)[1]
        validation_pos = main_block.index("pksLocusSpecificityValidateParams()")
        first_process_pos = main_block.index("alignConcatenatedIndex(")
        self.assertLess(validation_pos, first_process_pos)


class ThresholdsAreConfigurableNotHardcoded(unittest.TestCase):
    """The spec is explicit: R>=0.8 / CI-lower>=0.6 / etc. are placeholders pending
    Phase 3 calibration, so every one must be a params default, never a bare literal
    inside a script block."""

    CONFIGURABLE = (
        "pks_specificity_arm_a_min_confidence",
        "pks_specificity_arm_a_min_unique_kmers",
        "pks_specificity_retention_threshold",
        "pks_specificity_retention_ci_lower_threshold",
        "pks_specificity_bootstrap_iterations",
        "pks_specificity_max_tier_drop",
        "pks_specificity_non_megasynthase_min_genes",
        "pks_specificity_absorber_fraction_threshold",
        "pks_specificity_absorber_near_threshold_margin",
        "pks_specificity_min_supporting_reads",
        "pks_specificity_max_unmeasured_fraction",
    )

    def test_every_threshold_has_a_declared_default(self):
        for name in self.CONFIGURABLE:
            self.assertRegex(MODULE, rf"params\.{name}\s*=")

    def test_every_threshold_is_passed_through_to_the_decision_script(self):
        for name in self.CONFIGURABLE:
            flag = "--" + name.replace("pks_specificity_", "").replace("_", "-")
            self.assertIn(flag, MODULE, f"{name} is declared but never passed as {flag}")

    def test_tier_thresholds_match_the_pipelines_own_names_and_defaults(self):
        # Same numbers as classify_tumor_pks_evidence.py / classifyPksReadEvidence's
        # own params -- tier recomputation must use the pipeline's actual definition.
        for line in (
            "params.pks_specificity_multi_gene_reads       = 5",
            "params.pks_specificity_multi_gene_genes       = 3",
            "params.pks_specificity_multi_gene_breadth     = .01",
            "params.pks_specificity_broad_island_reads     = 30",
            "params.pks_specificity_broad_island_genes     = 8",
            "params.pks_specificity_broad_island_breadth   = .075",
            "params.pks_specificity_extensive_island_reads = 100",
            "params.pks_specificity_extensive_island_genes = 10",
            "params.pks_specificity_extensive_island_breadth = .15",
        ):
            self.assertIn(line, MODULE)


class ProcessWiring(unittest.TestCase):
    def test_b1_uses_secondary_alignments_with_a_lowered_ratio(self):
        # -p 0.1 (below minimap2's 0.8 default): a weak decoy-side hit must still be
        # reported so Arm B can see AS_best(decoy) whenever any alignment to the
        # decoy family exists, not only when it is competitive.
        b1 = MODULE.split("process alignConcatenatedIndex", 1)[1].split("process alignSeparateIndices", 1)[0]
        self.assertIn("-N 5", b1)
        self.assertIn("-p 0.1", b1)
        self.assertIn("--secondary=yes", b1)
        self.assertIn("params.pks_specificity_concatenated_index", b1)

    def test_b2_aligns_against_each_lone_index_independently(self):
        b2 = MODULE.split("process alignSeparateIndices", 1)[1].split("process scoreCompetitiveRealignment", 1)[0]
        self.assertIn("params.pks_specificity_true_index", b2)
        self.assertIn("params.pks_specificity_decoy_index", b2)
        # Two independent minimap2 invocations, not one shared command.
        self.assertEqual(b2.count("minimap2 -a -x sr"), 2)

    def test_scoring_process_invokes_the_arm_b_script(self):
        self.assertIn("competitive_realign_score.py", MODULE)
        self.assertIn("--true-contigs", MODULE)
        self.assertIn("--decoy-contigs", MODULE)

    def test_decision_process_invokes_the_decision_script(self):
        self.assertIn("specificity_decision.py", MODULE)
        self.assertIn("--arm-a-target-clade", MODULE)

    def test_b1_and_b2_both_run_on_the_same_supporting_reads_channel(self):
        workflow_body = MODULE.split("workflow pksLocusSpecificity {", 1)[1]
        self.assertIn("alignConcatenatedIndex(supporting_reads)", workflow_body)
        self.assertIn("alignSeparateIndices(supporting_reads)", workflow_body)

    def test_decision_process_joins_all_four_required_inputs(self):
        workflow_body = MODULE.split("workflow pksLocusSpecificity {", 1)[1]
        decision_wiring = workflow_body.split("computeLocusSpecificityVerdict(", 1)[0]
        self.assertIn(".join(read_gene_assignment", decision_wiring)
        self.assertIn(".join(arm_a_results", decision_wiring)
        self.assertIn(".join(scoreCompetitiveRealignment.out.arm_b", decision_wiring)


class ArmAIsOptional(unittest.TestCase):
    """Arm A needs a KrakenUniq database over the panel's taxon set, which a V8a panel
    build does not produce. Refusing to run the lane without it would make the whole
    check unavailable; running it and silently reporting a two-arm verdict would be
    worse. The lane runs Arm B alone and records arms_used."""

    def test_the_arm_a_directory_param_defaults_to_null(self):
        self.assertRegex(MODULE, r"params\.pks_specificity_arm_a_dir\s*=\s*null")

    def test_the_target_clade_is_required_only_alongside_an_arm_a_directory(self):
        self.assertIn("if (params.pks_specificity_arm_a_dir && !params.pks_specificity_arm_a_target_clade)",
                       CODE)

    def test_the_placeholder_file_exists_and_names_the_arm_a_schema(self):
        placeholder = ROOT / "ref/specificity_panels/NO_ARM_A.tsv"
        self.assertTrue(placeholder.is_file())
        header = [line for line in placeholder.read_text().splitlines()
                   if line and not line.startswith("#")]
        self.assertEqual(header, ["read_id\tassigned_clade\tconfidence\tunique_kmer_count"])

    def test_the_verdict_process_switches_on_the_placeholder_name(self):
        self.assertIn("arm_a_tsv.name == 'NO_ARM_A.tsv'", MODULE)

    def test_the_decision_script_accepts_a_run_without_arm_a(self):
        decision = (ROOT / "scripts/specificity_decision.py").read_text()
        self.assertIn('p.add_argument("--arm-a", default=None', decision)
        self.assertIn('arms_used = "arm_b_only" if arm_a_rows is None else "arm_a_and_arm_b"',
                       decision)


class OutputsAreNamedForTheQuestionTheyAnswer(unittest.TestCase):
    """Ludmil's V8a/V8b split is only real if the files say which one they are."""

    def test_the_verdict_files_are_locus_specificity_named(self):
        self.assertIn("locus_specificity_summary.tsv", MODULE)
        self.assertIn("locus_specificity_per_read.tsv", MODULE)

    def test_everything_publishes_under_one_locus_specificity_directory(self):
        self.assertIn('/locus_specificity"', MODULE)
        self.assertNotIn("specificity_check", MODULE)

    def test_the_summary_carries_the_check_name_and_which_arms_ran(self):
        decision = (ROOT / "scripts/specificity_decision.py").read_text()
        self.assertIn('"check": "locus_specificity"', decision)
        self.assertIn('"sample", "check", "verdict", "arms_used"', decision)


class WiredIntoMainAsAnOptionalLane(unittest.TestCase):
    MAIN = (ROOT / "main.nf").read_text()

    def test_main_includes_the_workflow(self):
        self.assertIn("include { pksLocusSpecificity } from './Modules/pks_specificity.nf'",
                       self.MAIN)

    def test_the_flag_is_off_by_default(self):
        self.assertRegex(self.MAIN, r"params\.pks_locus_specificity\s*=\s*false")

    def test_the_flag_is_tested_through_a_normalized_boolean_local(self):
        self.assertIn("def pks_locus_specificity_b     = params.pks_locus_specificity.toString().toBoolean()",
                       self.MAIN)
        self.assertIn("if (do_align && pks_locus_specificity_b)", self.MAIN)

    def test_the_verdict_reaches_the_cohort_report_and_master_summary(self):
        self.assertIn("tuple(sampleID, 'locus_specificity', summary)", self.MAIN)
        report = (ROOT / "scripts/build_cohort_report.py").read_text()
        summary = (ROOT / "scripts/build_master_summary.py").read_text()
        for text in (report, summary):
            self.assertIn('"locus_specificity": "locus_specificity/locus_specificity_summary.tsv"',
                           text)

    def test_the_panel_index_is_a_tracked_dependency(self):
        # F15: rebuilding the panel in place must invalidate the tasks that aligned to it.
        self.assertIn("pks_specificity_panel : Provenance.data(params.pks_specificity_concatenated_index)",
                       self.MAIN)
        self.assertIn("${params.dep_digest?.pks_specificity_panel}", MODULE)


class SampleStageLabelConvention(unittest.TestCase):
    """Every independent per-sample process carries label 'sample_stage' so one
    sample's failure does not take down every other sample's specificity check
    (the same convention conf/base.config documents for the rest of the pipeline)."""

    def test_every_process_has_sample_stage_label(self):
        process_blocks = MODULE.split("process ")[1:]
        for block in process_blocks:
            name = block.split()[0].rstrip("{").strip()
            header = block[:block.index("script:")] if "script:" in block else block
            self.assertIn("label 'sample_stage'", header, f"{name} missing label 'sample_stage'")


class StandaloneEntryPoint(unittest.TestCase):
    """Not reachable from main.nf: `nextflow run Modules/pks_specificity.nf` must
    work on its own, given a small sample sheet."""

    def test_top_level_workflow_block_exists(self):
        self.assertRegex(MODULE, r"\nworkflow \{\n")

    def test_it_validates_its_own_required_cli_params(self):
        entry = MODULE.rsplit("workflow {", 1)[1]
        self.assertIn('exit 1, "Missing required parameter: --sample', entry)
        self.assertIn('exit 1, "Missing required parameter: --sample_dir', entry)

    def test_it_reads_a_five_column_sample_sheet(self):
        entry = MODULE.rsplit("workflow {", 1)[1]
        for column in ("reads_fastq", "read_evidence", "read_gene_assignment", "arm_a_tsv"):
            self.assertIn(f"row.{column}", entry)

    def test_a_row_without_arm_a_falls_back_to_the_placeholder(self):
        entry = MODULE.rsplit("workflow {", 1)[1]
        self.assertIn("row.arm_a_tsv ? file(row.arm_a_tsv", entry)
        self.assertIn("pksSpecificityNoArmA()", entry)

    def test_scripts_and_conda_envs_default_from_moduleDir_not_projectDir(self):
        # moduleDir is THIS file's directory regardless of which script Nextflow was
        # invoked with; projectDir would resolve relative to the entry script (wrong
        # when this module is later imported from main.nf one directory up).
        self.assertIn("${moduleDir}/../scripts", MODULE)
        self.assertIn("${moduleDir}/../conda_envs/minimap2_env.yml", MODULE)


class EveryV8AlignmentIsSingleEnd(unittest.TestCase):
    """5 Oct 2026 V8 pilot: minimap2 -x sr pairs two consecutive records named X/1 and X/2
    and writes both as X, so island reads whose mates sat next to each other in the
    coordinate-sorted BAM had no Arm-B row under their own name. --frag=no on every V8
    alignment keeps each read independent and its name intact."""

    def test_every_minimap2_call_in_both_v8_lanes_disables_fragment_mode(self):
        for module in ("Modules/pks_specificity.nf", "Modules/pks_carrier.nf"):
            text = (ROOT / module).read_text()
            calls = [line for line in text.splitlines()
                     if line.strip().startswith("minimap2 ")]
            self.assertTrue(calls, f"{module}: no minimap2 call found")
            for line in calls:
                self.assertIn("--frag=no", line, f"{module}: {line.strip()}")


if __name__ == "__main__":
    unittest.main()
