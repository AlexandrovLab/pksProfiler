"""One row per sample, joining every stage that ran.

Stages that were not enabled must contribute no columns at all, rather than a wall of NA —
so the width of the table tells you what the run actually did.
"""
import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/build_master_summary.py"


def write(path, header, rows, delim="\t"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter=delim)
        w.writerow(header)
        w.writerows(rows)


def minimal(root):
    """Only read profiling — the smallest real run."""
    write(root / "cohort/qc/pks.qc.summary.tsv",
          ["Sample", "input_alignment_records", "reads_after_fastp", "reads_after_hg38",
           "reads_after_t2t_phix", "num_clb_genes_align", "reads_clb_genes_align",
           "num_clb_genes_hmm", "reads_clb_genes_hmm"],
          [["S1", "1000", "990", "500", "480", "15", "1497", "NA", "NA"],
           ["S2", "1000", "980", "400", "390", "0", "0", "NA", "NA"]])


def add_everything(root):
    write(root / "by_sample/S1/read_evidence.tsv",
          ["sample", "read_evidence", "pks_reads", "clb_genes_detected", "island_breadth_1x"],
          [["S1", "extensive_island", "1497", "15", "0.731800"]])
    write(root / "by_sample/S1/contigs/final_evidence/final_pks_evidence.tsv",
          ["sample", "final_structural_evidence", "assembler_agreement"],
          [["S1", "near_complete_island", "both"]])
    write(root / "by_sample/S1/genomes/pks_mag_summary.tsv",
          ["sample", "bin_id", "classification", "completeness", "distinct_clb_genes"],
          [["S1", "bin.1", "d__Bacteria;s__Escherichia coli", "96.4", "17"],
           ["S1", "bin.2", "d__Bacteria;s__Bacteroides fragilis", "88.1", "0"]])
    write(root / "by_sample/S1/community/community_prophage_inventory.tsv",
          ["sample", "bin_id", "prophage_id"],
          [["S1", "bin.2", "p1"], ["S1", "bin.2", "p2"]])
    write(root / "by_sample/S1/community/pks_community_interactions.tsv",
          ["sample", "pks_producer_bin", "recipient_bin", "recipient_prophage_count"],
          [["S1", "bin.1", "bin.2", "2"], ["S1", "bin.1", "bin.3", "0"]])
    write(root / "by_sample/S1/community/pks_island_mobility.tsv",
          ["sample", "bin_id", "nearby_integrase", "nearby_trna"],
          [["S1", "bin.1", "True", "False"]])
    write(root / "by_sample/S1/strain/strain_types.tsv",
          ["sample", "unit", "status", "ST", "clonal_complex", "phylogroup"],
          [["S1", "bin.1", "typed", "73", "ST73_Cplx", "B2"]])


def run(root):
    out = root / "master.tsv"
    r = subprocess.run([sys.executable, str(SCRIPT), "--results", str(root), "--output", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with out.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows, r.stderr


class MasterSummaryTests(unittest.TestCase):
    def test_minimal_run_has_only_profiling_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); minimal(root)
            rows, log = run(root)
        self.assertEqual([r["sample"] for r in rows], ["S1", "S2"])
        self.assertIn("reads_clb_genes_align", rows[0])
        for absent in ("read_evidence", "phylogroup", "mag_bins_total", "community_prophages_total"):
            self.assertNotIn(absent, rows[0], f"{absent} present although that stage did not run")
        self.assertIn("not run", log)

    def test_full_run_joins_every_stage_onto_one_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); minimal(root); add_everything(root)
            rows, _ = run(root)
        s1 = next(r for r in rows if r["sample"] == "S1")
        self.assertEqual(s1["reads_clb_genes_align"], "1497")
        self.assertEqual(s1["read_evidence"], "extensive_island")
        self.assertEqual(s1["final_structural_evidence"], "near_complete_island")
        self.assertEqual(s1["mag_bins_total"], "2")
        self.assertEqual(s1["mag_bins_pks_positive"], "1")
        self.assertIn("Escherichia coli", s1["pks_mag_taxonomy"])
        self.assertEqual(s1["pks_mag_clb_genes"], "17")
        self.assertEqual(s1["community_prophages_total"], "2")
        self.assertEqual(s1["community_neighbours_assessed"], "2")
        self.assertEqual(s1["community_neighbours_with_prophage"], "1")
        self.assertEqual(s1["island_nearby_integrase"], "yes")
        self.assertEqual(s1["island_nearby_trna"], "no")
        self.assertEqual(s1["ST"], "73")
        self.assertEqual(s1["phylogroup"], "B2")

    def test_sample_missing_from_later_stages_is_kept_with_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); minimal(root); add_everything(root)
            rows, _ = run(root)
        s2 = next(r for r in rows if r["sample"] == "S2")
        # S2 was profiled and found negative, so it has no MAGs or typing — but must not vanish.
        self.assertEqual(s2["reads_clb_genes_align"], "0")
        self.assertEqual(s2["read_evidence"], "NA")
        self.assertEqual(s2["phylogroup"], "NA")

    def test_the_pks_positive_genome_is_the_one_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); minimal(root); add_everything(root)
            rows, _ = run(root)
        s1 = next(r for r in rows if r["sample"] == "S1")
        # bin.2 is more numerous in the file but carries no clb genes; bin.1 must win.
        self.assertNotIn("Bacteroides", s1["pks_mag_taxonomy"])

    def test_empty_results_directory_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run([sys.executable, str(SCRIPT), "--results", tmp,
                                "--output", os.path.join(tmp, "o.tsv")],
                               capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("No per-sample results", r.stderr)


if __name__ == "__main__":
    unittest.main()
