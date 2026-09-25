"""Strain typing joins an MLST call to the 45,761-genome reference panel.

Reads are never typed. A partial allele profile yields an explicit "insufficient" status
rather than a guessed ST, and a pks+ assignment outside phylogroup B2 is flagged, since
B2 accounts for 4386/4760 (92%) of pks+ genomes in the panel.
"""
import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/assign_strain_type.py"
LOOKUP = ROOT / "ref/typing/st_phylogroup.tsv"

spec = importlib.util.spec_from_file_location("assign_strain_type", SCRIPT)
assign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(assign)

ALLELES_FULL = ["adk(6)", "fumC(4)", "gyrB(12)", "icd(1)", "mdh(20)", "purA(13)", "recA(7)"]


def run(tmp, st, alleles, label="bin.1"):
    mlst = Path(tmp) / "in.tsv"
    mlst.write_text("\t".join(["assembly.fa", "ecoli_achtman_4", st] + alleles) + "\n")
    out = Path(tmp) / "out.tsv"
    subprocess.run([sys.executable, str(SCRIPT), "--mlst", str(mlst), "--lookup", str(LOOKUP),
                    "--sample", "S1", "--label", label, "--output", str(out)], check=True)
    with out.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))[0]


class LookupTests(unittest.TestCase):
    def test_lookup_is_present_and_covers_the_dominant_pks_lineages(self):
        table = {r["ST"]: r for r in csv.DictReader(LOOKUP.open(), delimiter="\t")}
        self.assertGreater(len(table), 2000)
        # ST73 and ST95 are the two dominant pks+ lineages in the panel.
        for st in ("73", "95"):
            self.assertIn(st, table)
            self.assertEqual(table[st]["phylogroup"], "B2")
            self.assertGreater(float(table[st]["pks_fraction"]), 0.0)


class ClassificationTests(unittest.TestCase):
    def test_full_profile_with_known_st_is_typed(self):
        with tempfile.TemporaryDirectory() as tmp:
            row = run(tmp, "73", ALLELES_FULL)
        self.assertEqual(row["status"], "typed")
        self.assertEqual(row["ST"], "73")
        self.assertEqual(row["phylogroup"], "B2")
        self.assertEqual(row["loci_called"], "7")
        self.assertEqual(row["note"], "")

    def test_partial_profile_is_refused_not_guessed(self):
        with tempfile.TemporaryDirectory() as tmp:
            row = run(tmp, "-", ALLELES_FULL[:4] + ["-", "-", "-"])
        self.assertEqual(row["status"], "insufficient_loci")
        self.assertEqual(row["ST"], "NA")
        self.assertEqual(row["phylogroup"], "NA")
        self.assertIn("too little sequence", row["note"])

    def test_full_profile_without_an_st_is_a_novel_combination(self):
        with tempfile.TemporaryDirectory() as tmp:
            row = run(tmp, "-", ALLELES_FULL)
        self.assertEqual(row["status"], "novel_allele_combination")
        self.assertIn("not in the scheme", row["note"])

    def test_empty_mlst_output_is_insufficient(self):
        with tempfile.TemporaryDirectory() as tmp:
            mlst = Path(tmp) / "empty.tsv"; mlst.write_text("")
            out = Path(tmp) / "out.tsv"
            subprocess.run([sys.executable, str(SCRIPT), "--mlst", str(mlst), "--lookup", str(LOOKUP),
                            "--sample", "S1", "--label", "bin.1", "--output", str(out)], check=True)
            row = list(csv.DictReader(out.open(), delimiter="\t"))[0]
        self.assertEqual(row["status"], "insufficient_loci")

    def test_non_b2_phylogroup_is_flagged(self):
        table = {r["ST"]: r for r in csv.DictReader(LOOKUP.open(), delimiter="\t")}
        st = next(s for s, r in table.items() if r["phylogroup"] == "B1" and int(r["n_genomes"]) > 50)
        with tempfile.TemporaryDirectory() as tmp:
            row = run(tmp, st, ALLELES_FULL)
        self.assertEqual(row["status"], "typed")
        self.assertIn("atypical for pks+", row["note"])

    def test_st_absent_from_the_panel_is_reported_as_such(self):
        with tempfile.TemporaryDirectory() as tmp:
            row = run(tmp, "9999999", ALLELES_FULL)
        self.assertEqual(row["status"], "typed")
        self.assertIn("not represented in the reference panel", row["note"])


class WiringTests(unittest.TestCase):
    MAIN = (ROOT / "main.nf").read_text()
    MAG = (ROOT / "Modules/pks_mag.nf").read_text()

    def test_typing_runs_on_assemblies_not_reads(self):
        # The call site takes an assembled unit, not a reads channel.
        self.assertIn("magStrainTyping(bins_flat_ch", self.MAG)
        self.assertNotIn("StrainTyping(reads", self.MAG)

    def test_assembly_floor_is_enforced_in_the_module(self):
        typing = (ROOT / "Modules/pks_typing.nf").read_text()
        self.assertIn("typing_min_assembly_bp", typing)

    def test_lookup_default_points_into_the_repo(self):
        self.assertIn('params.st_phylogroup_lookup   = "${projectDir}/ref/typing/st_phylogroup.tsv"', self.MAIN)


if __name__ == "__main__":
    unittest.main()
