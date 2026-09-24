"""A provirus call needs length and hallmarks, not just a geNomad row.

T3, from the internal v0.0.2 audit. geNomad reported 134 "viral contigs" on AA-3850,
the top hits 369 bp with one gene and one hallmark. The count tracked how fragmented
the assembly was rather than any biology, and mobility claims were built on top of it
-- `island_in_prophage` is only meaningful if the prophage is real.

The metagenome MAG lane applies these floors: params.provirus_min_length_bp and
params.provirus_min_hallmarks. (A second, tumour-contig lane applied the same floors
until it was removed entirely -- no evidence MAG/contig binning could work on a
tumour sample's tiny bacterial fraction; see the removal commit for the rationale.)
"""
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# these scripts import mag_utils as a sibling, as they do when Nextflow runs them
sys.path.insert(0, str(ROOT / "scripts"))
MAG = (ROOT / "Modules/pks_mag.nf").read_text()
MAIN = (ROOT / "main.nf").read_text()


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"scripts/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


context = load("extract_genomic_context")
community = load("build_community_prophage")

HEADER = ["seq_name", "length", "topology", "coordinates", "n_genes",
          "genetic_code", "virus_score", "fdr", "n_hallmarks", "taxonomy"]


def summary(rows):
    lines = ["\t".join(HEADER)]
    for name, length, hallmarks in rows:
        lines.append("\t".join([f"{name}|provirus_1", str(length), "provirus",
                                "100-500", "4", "11", "0.95", "0.01",
                                str(hallmarks), "Caudoviricetes"]))
    return "\n".join(lines) + "\n"


class NoiseIsNotAProphage(unittest.TestCase):
    def write(self, rows):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        handle.write(summary(rows))
        handle.close()
        self.addCleanup(lambda: Path(handle.name).unlink(missing_ok=True))
        return handle.name

    def test_the_aa3850_pattern_is_rejected(self):
        # 369 bp, one hallmark -- the top hit among 134 "viral contigs".
        path = self.write([("k141_1", 369, 1)])
        self.assertEqual(context.parse_genomad_proviruses(path), {})
        self.assertEqual(community.read_prophages(path), [])

    def test_a_real_prophage_is_kept(self):
        path = self.write([("k141_2", 42000, 6)])
        self.assertEqual(len(context.parse_genomad_proviruses(path)), 1)
        self.assertEqual(len(community.read_prophages(path)), 1)

    def test_long_but_one_hallmark_is_rejected(self):
        path = self.write([("k141_3", 40000, 1)])
        self.assertEqual(context.parse_genomad_proviruses(path), {})

    def test_many_hallmarks_but_short_is_rejected(self):
        path = self.write([("k141_4", 500, 8)])
        self.assertEqual(context.parse_genomad_proviruses(path), {})

    def test_the_floors_are_adjustable(self):
        path = self.write([("k141_5", 369, 1)])
        self.assertEqual(len(context.parse_genomad_proviruses(path, 100, 1)), 1)

    def test_both_lanes_agree_on_the_same_input(self):
        rows = [("k141_1", 369, 1), ("k141_2", 42000, 6), ("k141_3", 40000, 1)]
        path = self.write(rows)
        self.assertEqual(len(context.parse_genomad_proviruses(path)),
                         len(community.read_prophages(path)))


class Wiring(unittest.TestCase):
    def test_both_thresholds_are_declared_once(self):
        self.assertIn("params.provirus_min_length_bp = 3000", MAIN)
        self.assertIn("params.provirus_min_hallmarks = 2", MAIN)

    def test_the_invocation_passes_them(self):
        self.assertEqual(MAG.count("--min-provirus-length"), 1)
        self.assertEqual(MAG.count("--min-provirus-hallmarks"), 1)


if __name__ == "__main__":
    unittest.main()
