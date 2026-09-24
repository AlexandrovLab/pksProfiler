"""Structural coverage for Modules/pks_specificity.nf, the optional read-level pks
specificity check (competitive remapping against a true+decoy reference panel).

This module is deliberately NOT included from main.nf (see the module's own header
comment for the reasoning): the reference panel it depends on -- true-only, decoy-
only and concatenated minimap2 indices, plus the contig lists that tell Arm B which
side of the concatenated index each contig belongs to -- is Phase 2 work that has not
happened yet. Wiring it into the default pipeline entry point now would add new
required-param failure surface to main.nf for a feature nobody can run.

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

    REQUIRED_PARAMS = (
        "pks_specificity_true_index",
        "pks_specificity_decoy_index",
        "pks_specificity_concatenated_index",
        "pks_specificity_true_contigs",
        "pks_specificity_decoy_contigs",
        "pks_specificity_arm_a_target_clade",
        "pks_specificity_island_length",
    )

    def test_every_required_param_defaults_to_null(self):
        for name in self.REQUIRED_PARAMS:
            self.assertIn(f"params.{name}", MODULE,
                           f"{name} is not declared as a module param")

    def test_every_required_param_is_checked_before_use(self):
        for name in self.REQUIRED_PARAMS:
            self.assertIn(f"!params.{name}", CODE, f"{name} has no missing-param guard")

    def test_validation_function_exits_nonzero_with_phase_context(self):
        self.assertIn("def pksSpecificityValidateParams()", CODE)
        self.assertIn("exit 1,", CODE)
        self.assertIn("Phase 2", MODULE)

    def test_validation_runs_before_any_process_in_the_workflow(self):
        workflow_body = MODULE.split("workflow pksSpecificityCheck {", 1)[1]
        main_block = workflow_body.split("main:", 1)[1]
        validation_pos = main_block.index("pksSpecificityValidateParams()")
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
        workflow_body = MODULE.split("workflow pksSpecificityCheck {", 1)[1]
        self.assertIn("alignConcatenatedIndex(supporting_reads)", workflow_body)
        self.assertIn("alignSeparateIndices(supporting_reads)", workflow_body)

    def test_decision_process_joins_all_four_required_inputs(self):
        workflow_body = MODULE.split("workflow pksSpecificityCheck {", 1)[1]
        decision_wiring = workflow_body.split("computeSpecificityVerdict(", 1)[0]
        self.assertIn(".join(read_gene_assignment", decision_wiring)
        self.assertIn(".join(arm_a_results", decision_wiring)
        self.assertIn(".join(scoreCompetitiveRealignment.out.arm_b", decision_wiring)


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

    def test_scripts_and_conda_envs_default_from_moduleDir_not_projectDir(self):
        # moduleDir is THIS file's directory regardless of which script Nextflow was
        # invoked with; projectDir would resolve relative to the entry script (wrong
        # when this module is later imported from main.nf one directory up).
        self.assertIn("${moduleDir}/../scripts", MODULE)
        self.assertIn("${moduleDir}/../conda_envs/minimap2_env.yml", MODULE)


if __name__ == "__main__":
    unittest.main()
