"""A real execution exists for every mode (Ludmil, 26 Sep mode memo section 6).

The static half of tests/check_mode_execution.sh, tests/fixtures/generate_mode_fixtures.py
and tests/fixtures/assert_mode_outputs.py -- the same relationship
test_e2e_execution_coverage.py has to check_e2e_execution.sh. The proof is the script,
run live (real conda environments, real minimap2/bowtie2/fastp/fastplong):

    CONDA_CACHE_DIR=/path/to/cache bash tests/check_mode_execution.sh
"""
import stat
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tests/check_mode_execution.sh"
GENERATOR = ROOT / "tests/fixtures/generate_mode_fixtures.py"
ASSERTIONS = ROOT / "tests/fixtures/assert_mode_outputs.py"


def code_only(text):
    return "\n".join(line for line in text.splitlines() if not line.strip().startswith("#"))


class ModeExecutionCoverage(unittest.TestCase):
    def setUp(self):
        self.script = SCRIPT.read_text()
        self.code = code_only(self.script)
        self.assertions = ASSERTIONS.read_text()

    def test_files_exist_and_are_executable(self):
        for path in (SCRIPT, GENERATOR, ASSERTIONS):
            self.assertTrue(path.is_file(), path)
            self.assertTrue(path.stat().st_mode & stat.S_IXUSR, f"not executable: {path}")

    def test_real_runs_not_previews(self):
        self.assertIn("nextflow run", self.code)
        self.assertNotIn("-preview", self.code)

    def test_every_mode_and_long_read_lane_is_run(self):
        self.assertIn('MODES=${MODES:-"wgs wgs_assembly wes rna metagenome wgs_ont rna_ont normal_parity"}',
                      self.script)
        for mode, flags in (("wgs", "--sample_type wgs"), ("wes", "--sample_type wes"),
                            ("rna", "--sample_type rna"), ("metagenome", "--sample_type metagenome"),
                            ("wgs_ont", "--sample_type wgs --read_type ont"),
                            ("rna_ont", "--sample_type rna --read_type ont")):
            with self.subTest(mode=mode):
                self.assertIn(f"run_mode {mode} {flags}", self.code)

    def test_targeted_recruitment_runs_for_real(self):
        # Ludmil, 25 Sep A1/A2: an extensive_island sample with flank-only, boundary and
        # chimeric pairs runs recruitment + both assemblers for real.
        self.assertIn("run_mode wgs_assembly --sample_type wgs --profiling_method bowtie2", self.code)
        self.assertIn("recruitment_R1.fastq.gz", self.code)
        generator = GENERATOR.read_text()
        for needle in ("sim_recruit_flank_", "sim_recruit_boundary_left", "sim_recruit_boundary_right",
                       "sim_recruit_chimera_"):
            self.assertIn(needle, generator)
        for needle in ("island_supporting_fragments", "fragments_flank_only", "fragments_boundary_spanning",
                       "reads_rejected_low_query_coverage", "recruitment Bowtie2 has no -k 1",
                       "assembly_input_diagnostic"):
            self.assertIn(needle, self.assertions)

    def test_rna_uses_a_real_splice_index_and_junctions(self):
        self.assertIn('-x splice:sr -d "$fixtures/human_splice_sr.mmi"', self.code)
        self.assertIn('--rna_junctions_bed "$fixtures/junctions.bed"', self.code)
        generator = GENERATOR.read_text()
        self.assertIn('"GT"', generator)
        self.assertIn('"AG"', generator)

    def test_wes_and_rna_qc_run_for_real(self):
        # Ludmil, 26 Sep mode memo sections 3 and 4: WES runs from a BAM with 0x400 flags
        # and a capture BED; RNA checks clb_transcripts_detected and residual human.
        self.assertIn('--wes_capture_bed "$fixtures/wes_capture.bed"', self.code)
        self.assertIn('"$fixtures/wes_input.bam"', self.code)
        self.assertIn('--input_data_type "$input_type"', self.code)
        generator = GENERATOR.read_text()
        self.assertIn("1024 if i < WES_HOST_DUPLICATES", generator)
        self.assertIn("for flag in (256, 2048)", generator)
        for needle in ("wes_input_duplicate_fraction", "wes_capture_off_target_fraction",
                       "wes_pks_island_duplicate_fraction", "clb_transcripts_detected",
                       "rna_residual_human_fraction", "Residual-human check", "-n 2 -m 20 -s 30"):
            with self.subTest(needle=needle):
                self.assertIn(needle, self.assertions)

    def test_assertions_cover_branch_disabled_features_preset_and_state(self):
        for needle in ("EVIDENCE_FILES", "no targeted assembly", "no MAG lane", "-ax splice:sr ",
                       "-ax map-ont ", "-ax splice ", "terminal_state", "pks_reads_per_million",
                       "targeted recruitment and assembly switched off", "master summary mode"):
            with self.subTest(needle=needle):
                self.assertIn(needle, self.assertions)


if __name__ == "__main__":
    unittest.main()
