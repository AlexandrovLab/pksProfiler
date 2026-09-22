"""One preflight contract, run before any task, reporting everything at once.

F20, from Ludmil's 2026-09-19 report. Validation existed but asked only whether a flag
was set, never whether what it pointed at was usable, and it aborted at the first
problem. `--hg38_db /typo/human.mmi` passed every check and failed in the first
mapReads task, after the whole cohort had been extracted.
"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT = ROOT / "scripts/preflight.py"
BOWTIE2_SUFFIXES = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
HMM_PRESS = (".h3f", ".h3i", ".h3m", ".h3p")
CLB_GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]


class PreflightTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.d = Path(self.tmp.name)

        # a run that should pass cleanly, which every test then breaks in one way
        (self.d / "a.bam").write_text("x")
        (self.d / "sheet.csv").write_text(
            f"patient,bam\nCASE_01,{self.d}/a.bam\nCASE_02,{self.d}/a.bam\n")
        for name in ("hg38.mmi", "t2t.mmi", "ref.fna", "clb.gff", "st_lookup.tsv"):
            (self.d / name).write_text("x")
        for suffix in BOWTIE2_SUFFIXES:
            (self.d / f"index{suffix}").write_text("x")
            (self.d / f"recruit{suffix}").write_text("x")
        model = self.d / "clb.hmm"
        model.write_text("".join(f"NAME  {gene}\n" for gene in CLB_GENES))
        for suffix in HMM_PRESS:
            (self.d / f"clb.hmm{suffix}").write_text("x")

        kraken = self.d / "kraken"
        kraken.mkdir()
        for name in ("database.kdb", "database.idx",
                     "database150mers.kmer_distrib", "database100mers.kmer_distrib"):
            (kraken / name).write_text("x")

        self.config = {
            "sample": str(self.d / "sheet.csv"),
            "input_data_type": "bam",
            "profiling_method": "both",
            "hg38_db": str(self.d / "hg38.mmi"),
            "t2t_phix_db": str(self.d / "t2t.mmi"),
            "pks_reference_fasta": str(self.d / "ref.fna"),
            "pks_genome_annotation": str(self.d / "clb.gff"),
            "pks_genome": str(self.d / "index"),
            "pks_recruit_index": str(self.d / "recruit"),
            "hmm_model": str(model),
            "kraken_db": str(kraken),
            "bracken_read_length": 150,
            "pks_taxa": True,
            "enable_strain_typing": True,
            "st_phylogroup_lookup": str(self.d / "st_lookup.tsv"),
        }

    def run_preflight(self, **overrides):
        config = dict(self.config)
        config.update(overrides)
        path = self.d / "config.json"
        path.write_text(json.dumps(config))
        return subprocess.run([sys.executable, str(PREFLIGHT), "--config", str(path)],
                              capture_output=True, text=True)

    def test_a_complete_configuration_passes(self):
        result = self.run_preflight()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("preflight passed", result.stdout)

    def test_a_reference_that_is_set_but_absent_is_caught(self):
        # The case that motivated F20: set, so every old check passed.
        result = self.run_preflight(hg38_db=str(self.d / "typo.mmi"))
        self.assertEqual(result.returncode, 1)
        self.assertIn("--hg38_db does not exist", result.stdout)

    def test_every_problem_is_reported_in_one_pass(self):
        result = self.run_preflight(hg38_db=str(self.d / "typo.mmi"),
                                    t2t_phix_db=str(self.d / "also_typo.mmi"),
                                    bracken_read_length=151)
        self.assertEqual(result.returncode, 1)
        self.assertIn("found 3 problems", result.stdout)
        for expected in ("--hg38_db", "--t2t_phix_db", "--bracken_read_length"):
            with self.subTest(expected=expected):
                self.assertIn(expected, result.stdout)

    def test_a_bracken_length_the_database_lacks_is_caught_with_the_alternatives(self):
        result = self.run_preflight(bracken_read_length=151)
        self.assertEqual(result.returncode, 1)
        self.assertIn("--bracken_read_length 151 has no distribution", result.stdout)
        self.assertIn("100, 150", result.stdout)

    def test_an_incomplete_bowtie2_index_is_caught(self):
        (self.d / "index.rev.2.bt2").unlink()
        result = self.run_preflight()
        self.assertEqual(result.returncode, 1)
        self.assertIn(".rev.2.bt2", result.stdout)

    def test_an_unpressed_hmm_model_is_caught(self):
        (self.d / "clb.hmm.h3i").unlink()
        result = self.run_preflight()
        self.assertEqual(result.returncode, 1)
        self.assertIn("hmmpress", result.stdout)

    def test_an_hmm_model_missing_clb_genes_is_caught(self):
        (self.d / "clb.hmm").write_text("NAME  clbA\nNAME  clbB\n")
        result = self.run_preflight()
        self.assertEqual(result.returncode, 1)
        self.assertIn("does not cover", result.stdout)

    def test_hmm_checks_are_skipped_when_hmm_profiling_is_off(self):
        (self.d / "clb.hmm").unlink()
        result = self.run_preflight(profiling_method="bowtie2")
        self.assertEqual(result.returncode, 0, result.stdout)

    def test_sheet_identifier_rules_are_enforced_before_anything_runs(self):
        (self.d / "sheet.csv").write_text(
            f"patient,bam\nCASE_01,{self.d}/a.bam\nCASE_01,{self.d}/a.bam\nbad id!,{self.d}/a.bam\n")
        result = self.run_preflight()
        self.assertEqual(result.returncode, 1)
        self.assertIn("duplicate patient value", result.stdout)
        self.assertIn("invalid patient value", result.stdout)

    def test_an_input_file_named_in_the_sheet_but_absent_is_caught(self):
        (self.d / "sheet.csv").write_text(f"patient,bam\nCASE_01,{self.d}/gone.bam\n")
        result = self.run_preflight()
        self.assertEqual(result.returncode, 1)
        self.assertIn("bam does not exist", result.stdout)

    def test_cram_input_without_a_reference_is_caught(self):
        (self.d / "s.cram").write_text("x")
        (self.d / "sheet.csv").write_text(f"patient,cram\nCASE_01,{self.d}/s.cram\n")
        result = self.run_preflight(input_data_type="cram")
        self.assertEqual(result.returncode, 1)
        self.assertIn("--cram_reference", result.stdout)

    def test_mag_stages_require_their_databases_to_exist(self):
        result = self.run_preflight(enable_mags=True, gtdbtk_db=str(self.d / "gtdb"),
                                    checkm2_db=str(self.d / "checkm2.dmnd"),
                                    genomad_db=str(self.d / "genomad"))
        self.assertEqual(result.returncode, 1)
        for expected in ("--gtdbtk_db", "--checkm2_db", "--genomad_db"):
            with self.subTest(expected=expected):
                self.assertIn(expected, result.stdout)

    def test_stage_databases_are_not_demanded_when_the_stage_is_off(self):
        result = self.run_preflight(enable_mags=False)
        self.assertEqual(result.returncode, 0, result.stdout)


class WiringTests(unittest.TestCase):
    MAIN = (ROOT / "main.nf").read_text()

    def test_the_workflow_runs_preflight_and_stops_on_failure(self):
        self.assertIn("Preflight.run(", self.MAIN)
        self.assertIn("preflight_status != 0", self.MAIN)

    def test_preflight_runs_before_any_process_is_invoked(self):
        body = self.MAIN.split("workflow {", 1)[1]
        preflight_at = body.index("Preflight.run(")
        first_process = body.index("extractReads(")
        self.assertLess(preflight_at, first_process,
                        "preflight must gate the run, not follow it")


if __name__ == "__main__":
    unittest.main()


class AdaptersAreOptional(unittest.TestCase):
    """fastp detects adapters itself; the 234-sequence list is opt-in.

    Matching every read against that file was a measurable share of the Hartwig
    extraction runtime for no change in what survived, so it is no longer the default.
    """

    FILTER = (ROOT / "Modules/filter_reads.nf").read_text()
    MAIN = (ROOT / "main.nf").read_text()

    def test_the_default_is_unset(self):
        self.assertIn("params.adapters     = null", self.MAIN)

    def test_the_flag_is_omitted_when_unset(self):
        self.assertIn('params.adapters ? "--adapter_fasta', self.FILTER)
        # never interpolated unconditionally into the command
        command = self.FILTER.split("fastp \\", 1)[1]
        self.assertNotIn('--adapter_fasta "${params.adapters}"', command)

    def test_preflight_does_not_demand_one(self):
        source = (ROOT / "scripts/preflight.py").read_text()
        self.assertIn('if config.get("adapters"):', source)

    def test_but_checks_it_when_given(self):
        import json, subprocess, sys, tempfile
        from pathlib import Path as P
        with tempfile.TemporaryDirectory() as directory:
            d = P(directory)
            (d / "sheet.csv").write_text("patient,bam\n")
            config = {"sample": str(d / "sheet.csv"), "input_data_type": "bam",
                      "profiling_method": "bowtie2",
                      "adapters": str(d / "not_here.fna")}
            (d / "c.json").write_text(json.dumps(config))
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/preflight.py"),
                 "--config", str(d / "c.json")], capture_output=True, text=True)
            self.assertIn("--adapters does not exist", result.stdout)
