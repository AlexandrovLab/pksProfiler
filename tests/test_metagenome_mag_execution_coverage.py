"""A real, executing test for the metagenome MAG "zero real bins" scenario exists.

Found this session via a real (non-mocked) execution of a since-removed tumour-MAG
combination: build_mag_summary.py's per-sample summary table
(genomes/pks_mag_summary.tsv) requires a real CheckM2 report and a real GTDB-Tk summary
before it reports anything for a sample. checkm2Predict/gtdbtkClassify only ever run on
real bins (Modules/pks_mag.nf, metabat2Bin.out.bins, optional: true) -- so a sample with
zero real bins never had either report, and the summary table was silently never
written for it at all, even carrying a real, positive call from U1's own "unbinned
pool" detection (metabat2Bin's --unbinned output, added this same session for exactly
the "MetaBAT2 recovered nothing" case). This is the metagenome lane, still fully active
and unaffected by the (separately removed) tumour-MAG combination that first surfaced
the general "process/channel wiring only tested via real execution" lesson.

This file is the static half of tests/check_metagenome_mag_execution.sh,
tests/fixtures/generate_metagenome_mag_fixtures.py and
tests/fixtures/assert_metagenome_mag_outputs.py -- the same relationship
test_e2e_execution_coverage.py has to tests/check_e2e_execution.sh. The empirical proof
is the scripts themselves, run live:

    DB_ROOT=/path/to/dbs bash tests/check_metagenome_mag_execution.sh

Not run from here: real conda environments (megahit, metabat2, checkm2, gtdbtk, prokka,
hmmer, geNomad) and the real production GTDB-Tk/CheckM2/geNomad databases. See
.github/workflows/ci.yml.
"""
import stat
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tests/check_metagenome_mag_execution.sh"
GENERATOR = ROOT / "tests/fixtures/generate_metagenome_mag_fixtures.py"
ASSERTIONS = ROOT / "tests/fixtures/assert_metagenome_mag_outputs.py"
CI = (ROOT / ".github/workflows/ci.yml").read_text()
PKS_MAG = (ROOT / "Modules/pks_mag.nf").read_text()
BUILD_MASTER_SUMMARY = (ROOT / "scripts/build_master_summary.py").read_text()


class ScriptsExistAndAreExecutable(unittest.TestCase):
    def test_the_three_files_exist(self):
        for path in (SCRIPT, GENERATOR, ASSERTIONS):
            self.assertTrue(path.is_file(), f"missing: {path}")

    def test_they_are_executable(self):
        for path in (SCRIPT, GENERATOR, ASSERTIONS):
            mode = path.stat().st_mode
            self.assertTrue(mode & stat.S_IXUSR, f"not executable: {path}")


def code_only(text):
    return "\n".join(line for line in text.splitlines() if not line.strip().startswith("#"))


class RealExecutionNotAPreview(unittest.TestCase):
    def setUp(self):
        self.script = SCRIPT.read_text()
        self.code = code_only(self.script)

    def test_preview_flag_is_never_used_in_the_actual_invocation(self):
        self.assertNotIn("-preview", self.code)

    def test_it_actually_invokes_nextflow_run(self):
        self.assertIn("nextflow run", self.script)

    def test_it_requires_the_real_mag_databases_rather_than_skipping(self):
        self.assertIn("GTDBTK_DB", self.script)
        self.assertIn("CHECKM2_DB", self.script)
        self.assertIn("GENOMAD_DB", self.script)

    def test_it_sets_enable_mags_true_for_metagenome(self):
        self.assertIn("--sample_type metagenome", self.script)
        self.assertIn("--enable_mags true", self.script)


class FixtureIsSizedToStayUnbinned(unittest.TestCase):
    def setUp(self):
        self.generator = GENERATOR.read_text()

    def test_it_reuses_the_base_generators_helpers_rather_than_duplicating_them(self):
        self.assertIn("from generate_e2e_fixtures import", self.generator)

    def test_it_covers_the_complete_pks_island(self):
        self.assertIn("2193827", self.generator)
        self.assertIn("2244594", self.generator)

    def test_it_documents_the_200kb_minclssize_floor(self):
        self.assertIn("200 kb", self.generator)
        self.assertIn("minClsSize", self.generator)

    def test_tiling_overlaps_so_megahit_can_actually_bridge_reads(self):
        self.assertIn("TILE_STEP", self.generator)
        step_line = next(ln for ln in self.generator.splitlines() if ln.startswith("TILE_STEP"))
        step = int(step_line.split("=", 1)[1].split("#", 1)[0].strip())
        self.assertLess(step, 300)


class AssertionsCoverTheRegression(unittest.TestCase):
    def setUp(self):
        self.assertions = ASSERTIONS.read_text()

    def test_it_confirms_zero_real_bins_rather_than_assuming_it(self):
        self.assertIn("bin_count", self.assertions)

    def test_it_checks_pks_mag_summary_exists_despite_zero_bins(self):
        self.assertIn("pks_mag_summary.tsv", self.assertions)
        self.assertIn("unbinned", self.assertions)

    def test_it_checks_no_real_checkm2_or_gtdbtk_output(self):
        self.assertIn("checkm2", self.assertions.lower())
        self.assertIn("gtdbtk", self.assertions.lower())

    def test_it_checks_the_u1_master_summary_columns(self):
        self.assertIn("mag_unbinned_locus_tier", self.assertions)
        self.assertIn("mag_unbinned_pks_positive", self.assertions)


class TheFixStaysInPlace(unittest.TestCase):
    """Static half of the regression guard: the fix exists in the source, wired the way
    the execution assertions above describe."""

    def test_stub_process_exists(self):
        self.assertIn("process stubBinlessMagQuality {", PKS_MAG)
        self.assertIn("val sampleID", PKS_MAG)

    def test_stub_process_has_no_publish_directive(self):
        # No real tool ran; publishing a fake checkm2/gtdbtk directory for a sample
        # that had zero bins would misrepresent what happened.
        block = PKS_MAG[PKS_MAG.index("process stubBinlessMagQuality {"):]
        block = block[:block.index("\nprocess ") if "\nprocess " in block[1:] else len(block)]
        self.assertNotIn("publishDir", block)

    def test_summary_input_ch_uses_the_merged_channels(self):
        workflow = PKS_MAG[PKS_MAG.index("workflow pksMAG {"):]
        self.assertIn("checkm2_report_ch = checkm2Predict.out.report.mix(stubBinlessMagQuality.out.checkm2)", workflow)
        self.assertIn("gtdbtk_summary_ch = gtdbtkClassify.out.summary.mix(stubBinlessMagQuality.out.gtdbtk)", workflow)
        self.assertIn("summary_input_ch = checkm2_report_ch", workflow)

    def test_tblouts_gained_a_remainder_true_fallback(self):
        # Before this fix, tblouts_per_sample_ch was also a strict join -- a sample
        # with zero real bins has no tblouts either, so it needed the same
        # remainder:true + `?: []` treatment contexts_per_sample_ch already had.
        workflow = PKS_MAG[PKS_MAG.index("workflow pksMAG {"):]
        self.assertIn(".join(tblouts_per_sample_ch, by: 0, remainder: true)", workflow)
        self.assertIn("tblouts ?: []", workflow)


class ConstructionStillPasses(unittest.TestCase):
    """This fix must not change what tests/check_workflow_construction.sh's
    "metagenome, MAGs" case builds -- it only changes what happens once metabat2
    actually runs, which -preview never does."""

    def test_metagenome_mags_construction_case_still_present(self):
        construction = (ROOT / "tests/check_workflow_construction.sh").read_text()
        self.assertIn('construct "metagenome, MAGs"', construction)
        self.assertIn("--enable_mags true", construction)


class CIWiring(unittest.TestCase):
    def test_ci_references_the_metagenome_mag_script(self):
        self.assertIn("check_metagenome_mag_execution.sh", CI)

    def test_the_job_is_not_on_every_push_or_pull_request(self):
        job = CI[CI.index("check_metagenome_mag_execution.sh"):]
        preceding = CI[:CI.index("check_metagenome_mag_execution.sh")]
        self.assertIn("workflow_dispatch", preceding + job)


if __name__ == "__main__":
    unittest.main()
