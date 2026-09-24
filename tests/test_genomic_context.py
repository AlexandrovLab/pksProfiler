"""Nearby-gene context: mobility markers and geNomad provirus overlap.

A clb gene inside a predicted provirus, or flanked by an integrase or a tRNA, is
evidence the island itself is mobile. The canonical pks island integrates at a tRNA
locus, so nearby_trna is a specific signal rather than a generic annotation.
"""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/extract_genomic_context.py"

GFF = """\
##gff-version 3
ctg1\tProkka\tCDS\t1000\t2000\t.\t+\t0\tID=PK_00001;gene=clbB;product=colibactin polyketide synthase
ctg1\tProkka\tCDS\t3000\t4000\t.\t+\t0\tID=PK_00002;gene=intA;product=Integrase core domain protein
ctg1\tProkka\ttRNA\t5000\t5080\t.\t+\t0\tID=PK_t0001;product=tRNA-Asn(gtt)
ctg1\tProkka\tCDS\t900000\t901000\t.\t+\t0\tID=PK_00003;gene=xyz;product=hypothetical protein
ctg2\tProkka\tCDS\t1000\t2000\t.\t+\t0\tID=PK_00010;gene=clbS;product=colibactin self-resistance protein
"""

# hmmsearch --tblout needs >=19 whitespace fields; col1 target, col3 query, col5 E-value.
def tblout(rows):
    head = "# target name        accession  query name  accession    E-value\n"
    body = "".join(
        f"{tag} - {gene} - {ev} 100.0 0.0 {ev} 100.0 0.0 1 1 0 0 1 1 1 1 -\n"
        for tag, gene, ev in rows
    )
    return head + body


GENOMAD = """\
seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy
ctg1|provirus_800_6000\t5200\tprovirus\t800-6000\t9\t11\t0.9710\tNA\t3\t12.5\tCaudoviricetes
ctg2\t3000\tNo terminal repeats\tNA\t4\t11\t0.5100\tNA\t0\t1.2\tUnclassified
"""


def run(tmp, genomad=None, window=50000):
    d = Path(tmp)
    (d / "in.gff").write_text(GFF)
    (d / "in.tblout").write_text(tblout([("PK_00001", "clbB", "1e-40"), ("PK_00010", "clbS", "1e-30")]))
    cmd = [sys.executable, str(SCRIPT), "--gff", str(d / "in.gff"), "--tblout", str(d / "in.tblout"),
           "--evalue", "1e-5", "--window", str(window), "--out", str(d / "out.tsv")]
    if genomad:
        (d / "genomad.tsv").write_text(GENOMAD)
        cmd += ["--genomad", str(d / "genomad.tsv")]
    subprocess.run(cmd, check=True, capture_output=True)
    with (d / "out.tsv").open() as fh:
        return {r["clb_gene"]: r for r in csv.DictReader(fh, delimiter="\t")}


class ContextColumnTests(unittest.TestCase):
    def test_new_columns_are_emitted(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = run(tmp, genomad=True)
        for column in ("nearby_trna", "in_prophage", "prophage_id",
                       "prophage_virus_score", "prophage_taxonomy"):
            self.assertIn(column, rows["clbB"])

    def test_integrase_and_trna_in_the_flank_are_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = run(tmp, genomad=True)
        self.assertEqual(rows["clbB"]["has_integrase"], "True")
        self.assertEqual(rows["clbB"]["nearby_trna"], "True")

    def test_clb_gene_inside_a_provirus_is_flagged(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = run(tmp, genomad=True)
        self.assertEqual(rows["clbB"]["in_prophage"], "True")
        self.assertEqual(rows["clbB"]["prophage_id"], "ctg1|provirus_800_6000")
        self.assertEqual(rows["clbB"]["prophage_virus_score"], "0.9710")
        self.assertEqual(rows["clbB"]["prophage_taxonomy"], "Caudoviricetes")

    def test_non_provirus_topology_is_not_counted(self):
        # ctg2 is a whole-contig virus call, not a provirus, so clbS gets no interval.
        with tempfile.TemporaryDirectory() as tmp:
            rows = run(tmp, genomad=True)
        self.assertEqual(rows["clbS"]["in_prophage"], "False")
        self.assertEqual(rows["clbS"]["prophage_id"], "")

    def test_without_genomad_the_columns_are_empty_not_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = run(tmp, genomad=None)
        self.assertEqual(rows["clbB"]["in_prophage"], "False")
        self.assertEqual(rows["clbB"]["prophage_id"], "")
        # Prokka-derived markers still work with no geNomad input.
        self.assertEqual(rows["clbB"]["has_integrase"], "True")

    def test_window_bounds_the_flank(self):
        with tempfile.TemporaryDirectory() as tmp:
            narrow = run(tmp, genomad=True, window=500)
        # At 500 bp neither the integrase (3 kb away) nor the tRNA (5 kb) is in range.
        self.assertEqual(narrow["clbB"]["has_integrase"], "False")
        self.assertEqual(narrow["clbB"]["nearby_trna"], "False")


class WiringTests(unittest.TestCase):
    MAG = (ROOT / "Modules/pks_mag.nf").read_text()
    MAIN = (ROOT / "main.nf").read_text()

    def test_window_is_a_parameter(self):
        self.assertIn("params.context_window_bp = 50000", self.MAIN)
        self.assertEqual(self.MAG.count("--window ${params.context_window_bp}"), 1)


if __name__ == "__main__":
    unittest.main()
