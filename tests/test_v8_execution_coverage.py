"""tests/check_v8_lanes_execution.sh: what the real V8a/V8b runs must prove.

Ludmil, 25 Sep report section 4, asks for two checks with two panels and warns that they
must not be conflated. The decision rules are unit-tested on hand-written SAM records
(tests/test_specificity_decision.py, tests/test_carrier_attribution.py) and the lane
separation is enforced on source (tests/test_v8_lane_separation.py). What neither can
show is that the lanes, run end to end through real minimap2 and real Nextflow wiring,
actually separate island reads from a homologous look-alike and actually attribute a
carrier read to the right organism.

This file is the structural coverage for that execution script -- the same convention
tests/test_e2e_execution_coverage.py follows. The script itself needs Nextflow, conda and
a couple of minutes, so it is not run from `unittest discover`.
"""
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "tests/check_v8_lanes_execution.sh"
SCRIPT = SCRIPT_PATH.read_text()
GENERATOR_PATH = ROOT / "tests/fixtures/generate_v8_fixtures.py"
GENERATOR = GENERATOR_PATH.read_text()
ASSERTIONS_PATH = ROOT / "tests/fixtures/assert_v8_outputs.py"
ASSERTIONS = ASSERTIONS_PATH.read_text()
CI = (ROOT / ".github/workflows/ci.yml").read_text()


class TheScriptAndItsFixturesExist(unittest.TestCase):
    def test_all_three_files_exist(self):
        for path in (SCRIPT_PATH, GENERATOR_PATH, ASSERTIONS_PATH):
            self.assertTrue(path.is_file(), f"{path} is missing")

    def test_the_fixtures_come_from_the_shipped_reference_not_from_sample_data(self):
        self.assertIn("indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic.fna",
                       GENERATOR)
        self.assertIn("random.Random(args.seed)", GENERATOR)

    def test_the_fixture_island_coordinates_match_the_pipelines_own(self):
        main = (ROOT / "main.nf").read_text()
        self.assertIn("params.pks_start_1based    = 2193827", main)
        self.assertIn("params.pks_end_1based      = 2244594", main)
        self.assertIn("ISLAND_START_1BASED = 2193827", GENERATOR)
        self.assertIn("ISLAND_END_1BASED = 2244594", GENERATOR)


class BothLanesRunForReal(unittest.TestCase):
    def test_it_runs_each_module_through_nextflow_not_the_scripts_directly(self):
        self.assertIn("nextflow run \"$repo_dir/Modules/pks_specificity.nf\"", SCRIPT)
        self.assertIn("nextflow run \"$repo_dir/Modules/pks_carrier.nf\"", SCRIPT)
        self.assertNotIn("-preview", SCRIPT)

    def test_it_builds_real_minimap2_indices_for_both_panels(self):
        for name in ("v8a_true", "v8a_decoy", "v8a_concatenated", "v8b_carriers"):
            self.assertIn(name, SCRIPT)
        self.assertIn('-x sr -d "$fixtures/$name.mmi"', SCRIPT)

    def test_it_checks_the_v8a_lane_refuses_to_run_without_a_panel(self):
        # A null panel path would otherwise produce an empty check that looks like a pass.
        self.assertIn("did not refuse a run with no panel", SCRIPT)
        self.assertIn("--pks_specificity_concatenated_index", SCRIPT)

    def test_v8a_is_run_without_arm_a_so_the_optional_arm_is_exercised(self):
        v8a_block = SCRIPT.split("Running the V8a", 1)[1].split("Running the V8b", 1)[0]
        self.assertNotIn("--pks_specificity_arm_a_dir", v8a_block)


class TheFixtureKnowsTheCorrectAnswer(unittest.TestCase):
    def test_the_v8a_decoy_is_a_homologous_slice_not_only_random_sequence(self):
        # Random background alone would make V8a trivial: the point is the conserved
        # NRPS/PKS domain look-alike, which is the failure mode the check exists for.
        self.assertIn("homolog", GENERATOR)
        self.assertIn("--homolog-divergence", GENERATOR)
        self.assertIn("decoy__nrps_pks_homolog", GENERATOR)
        self.assertIn("decoy__random_background", GENERATOR)

    def test_the_homolog_comes_from_a_different_part_of_the_island(self):
        # Otherwise a homolog read's best home would be ambiguous by construction and the
        # expected answer would be a coin flip.
        self.assertIn("island[30000:36000]", GENERATOR)
        self.assertIn("island[:20000]", GENERATOR)

    def test_the_v8b_panel_holds_two_carriers_of_the_same_island(self):
        self.assertIn("ecoli_test__clb_island", GENERATOR)
        self.assertIn("klebsiella_test__clb_island", GENERATOR)
        self.assertIn("--carrier-divergence", GENERATOR)

    def test_the_expected_answers_are_written_beside_the_fixtures(self):
        self.assertIn("expected.tsv", GENERATOR)
        self.assertIn("expected_carrier\\tklebsiella_test", GENERATOR)
        self.assertIn("read_expected", ASSERTIONS)


class TheAssertionsCheckBothDirectionsAndBothLanes(unittest.TestCase):
    def test_island_reads_must_be_confirmed_true(self):
        self.assertIn('("island", true_rows, "confirmed_true", "true_reads")', ASSERTIONS)

    def test_homolog_reads_must_be_confirmed_decoy(self):
        self.assertIn('("homolog", decoy_rows, "confirmed_decoy", "decoy_reads")',
                       ASSERTIONS)

    def test_retention_must_land_in_between_not_at_either_extreme(self):
        # Half the reads are homolog reads, so a lane that confirmed everything (1.0) or
        # rejected everything (0.0) would be caught even if the per-read shares passed.
        self.assertIn("0.2 <= retention <= 0.8", ASSERTIONS)

    def test_a_carrier_read_must_be_attributed_to_its_own_carrier(self):
        self.assertIn('r["best_carrier"] == wanted', ASSERTIONS)
        self.assertIn("expected_genus", ASSERTIONS)

    def test_an_arm_b_only_run_must_say_so(self):
        self.assertIn('row.get("arms_used") != "arm_b_only"', ASSERTIONS)

    def test_neither_lanes_output_may_carry_the_others_vocabulary(self):
        self.assertIn('if "carrier" in ",".join(row.keys()).lower()', ASSERTIONS)
        self.assertIn('for banned in ("SURVIVES", "NOT_DISTINGUISHABLE"', ASSERTIONS)


class CIWiring(unittest.TestCase):
    def test_ci_runs_the_v8_lane_check_in_the_scheduled_execution_job(self):
        self.assertIn("check_v8_lanes_execution.sh", CI)
        job = CI[CI.index("e2e-execution:"):]
        self.assertIn("check_v8_lanes_execution.sh", job)
        self.assertIn("github.event_name == 'schedule'", job)


if __name__ == "__main__":
    unittest.main()
