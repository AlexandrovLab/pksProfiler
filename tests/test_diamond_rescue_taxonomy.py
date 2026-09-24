#!/usr/bin/env python3
"""summarize_diamond_rescue_taxonomy.py: item 2 of the metagenome-report-visibility work.

select_diamond_reads.py keeps only the read ID from diamond.tsv and krakenuniq.output.txt
is never rejoined to it, so a rescued read's clb-gene hit and its community classifier
taxid were both computed and then discarded. This is the join that was missing.
"""
import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/summarize_diamond_rescue_taxonomy.py"

spec = importlib.util.spec_from_file_location("summarize_diamond_rescue_taxonomy", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

TAXDB = """\
1\t1\troot\tno rank
131567\t1\tcellular organisms\tno rank
2\t131567\tBacteria\tsuperkingdom
1224\t2\tProteobacteria\tphylum
1236\t1224\tGammaproteobacteria\tclass
91347\t1236\tEnterobacterales\torder
543\t91347\tEnterobacteriaceae\tfamily
561\t543\tEscherichia\tgenus
562\t561\tEscherichia coli\tspecies
570\t543\tKlebsiella\tgenus
573\t570\tKlebsiella pneumoniae\tspecies
"""


def write_taxdb(directory):
    (directory / "taxDB").write_text(TAXDB)
    return directory


class ClbGeneFromSseqid(unittest.TestCase):
    def test_extracts_the_gene_token(self):
        sseqid = "SAMN01|SAMN01.contig00001_14|clbB|score=7820.10|hmmcov=0.9991"
        self.assertEqual(mod.clb_gene_from_sseqid(sseqid), "clbB")

    def test_unknown_when_no_clb_token_is_present(self):
        self.assertEqual(mod.clb_gene_from_sseqid("contig1|gene2|score=1.0"), "unknown")


class LoadDiamondHits(unittest.TestCase):
    def _write(self, tmp, lines):
        path = Path(tmp) / "diamond.tsv"
        path.write_text("\n".join(lines) + "\n" if lines else "")
        return path

    def test_a_passing_hit_is_kept_with_its_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(tmp, [
                "read1\tid|id|clbB|score=1|hmmcov=1\t90\t30\t100\t170\t1e-20\t80"])
            hits = mod.load_diamond_hits(path, 1e-5, 25, 0.5)
        self.assertEqual(hits, {"read1": "clbB"})

    def test_a_failing_hit_is_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            # length 10 aa, below the min-aa floor of 25.
            path = self._write(tmp, [
                "read1\tid|id|clbB|score=1|hmmcov=1\t90\t10\t100\t170\t1e-20\t80"])
            hits = mod.load_diamond_hits(path, 1e-5, 25, 0.5)
        self.assertEqual(hits, {})

    def test_the_lowest_evalue_hit_wins_when_a_read_has_two(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(tmp, [
                "read1\tid|id|clbA|score=1|hmmcov=1\t90\t30\t100\t170\t1e-8\t80",
                "read1\tid|id|clbB|score=1|hmmcov=1\t90\t30\t100\t170\t1e-20\t80",
            ])
            hits = mod.load_diamond_hits(path, 1e-5, 25, 0.5)
        self.assertEqual(hits, {"read1": "clbB"})

    def test_an_empty_file_yields_no_hits(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(tmp, [])
            hits = mod.load_diamond_hits(path, 1e-5, 25, 0.5)
        self.assertEqual(hits, {})

    def test_a_malformed_row_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(tmp, ["not\tenough\tcolumns"])
            with self.assertRaises(ValueError):
                mod.load_diamond_hits(path, 1e-5, 25, 0.5)


class LoadReadTaxids(unittest.TestCase):
    def test_only_classified_reads_are_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "kraken.tsv"
            path.write_text("C\tread1\t562\t100\tA:1\nU\tread2\t0\t100\tA:0\n")
            taxids = mod.load_read_taxids(path)
        self.assertEqual(taxids, {"read1": "562"})

    def test_taxid_zero_is_not_a_classification(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "kraken.tsv"
            path.write_text("C\tread1\t0\t100\tA:0\n")
            taxids = mod.load_read_taxids(path)
        self.assertEqual(taxids, {})


class BuildTable(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        write_taxdb(Path(self.tmp.name))
        self.parents, self.ranks, self.names = mod.load_taxonomy(Path(self.tmp.name))

    def test_a_read_resolves_to_its_species_name(self):
        counts = mod.build_table({"r1": "clbB"}, {"r1": "573"},
                                 self.parents, self.ranks, self.names)
        self.assertEqual(counts, {("Klebsiella pneumoniae", "573", "clbB"): 1})

    def test_two_reads_same_species_same_gene_are_summed(self):
        counts = mod.build_table({"r1": "clbB", "r2": "clbB"}, {"r1": "573", "r2": "573"},
                                 self.parents, self.ranks, self.names)
        self.assertEqual(counts[("Klebsiella pneumoniae", "573", "clbB")], 2)

    def test_an_unclassified_read_is_labelled_not_dropped(self):
        counts = mod.build_table({"r1": "clbA"}, {}, self.parents, self.ranks, self.names)
        self.assertEqual(counts, {("Unclassified", "0", "clbA"): 1})

    def test_a_taxid_with_no_species_ancestor_is_unresolved_not_dropped(self):
        counts = mod.build_table({"r1": "clbA"}, {"r1": "131567"},
                                 self.parents, self.ranks, self.names)
        self.assertEqual(counts, {("Unresolved_at_species", "-1", "clbA"): 1})


class EndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        write_taxdb(self.dir)
        (self.dir / "diamond.tsv").write_text(
            "readA\tid|id|clbB|score=1|hmmcov=1\t90\t30\t100\t170\t1e-20\t80\n"
            "readB\tid|id|clbA|score=1|hmmcov=1\t90\t30\t100\t170\t1e-20\t80\n")
        (self.dir / "kraken.tsv").write_text(
            "C\treadA\t573\t100\tA:1\n"
            "C\treadB\t562\t100\tA:1\n")
        self.out = self.dir / "out.tsv"

    def _run(self):
        return subprocess.run(
            [sys.executable, str(SCRIPT), "--sample", "S1",
             "--diamond", str(self.dir / "diamond.tsv"),
             "--kraken-output", str(self.dir / "kraken.tsv"),
             "--kraken-db", str(self.dir),
             "--output", str(self.out)], capture_output=True, text=True)

    def test_it_runs_and_writes_the_declared_columns(self):
        run = self._run()
        self.assertEqual(run.returncode, 0, run.stderr)
        with self.out.open(newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(list(rows[0].keys()),
                         ["sample", "organism", "taxid", "clb_gene", "read_count"])
        by_organism = {r["organism"]: r for r in rows}
        self.assertEqual(by_organism["Klebsiella pneumoniae"]["clb_gene"], "clbB")
        self.assertEqual(by_organism["Escherichia coli"]["clb_gene"], "clbA")
        self.assertEqual(by_organism["Klebsiella pneumoniae"]["read_count"], "1")

    def test_every_row_is_specifically_a_rescued_read(self):
        # Nothing in this table is anything but a DIAMOND rescue -- select_diamond_reads.py
        # only ever sees reads krakenPrefilter's target-taxon routing did not keep.
        self._run()
        with self.out.open(newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(len(rows), 2)

    def test_no_diamond_hits_skips_loading_the_taxonomy_and_still_writes_a_header(self):
        (self.dir / "empty.tsv").write_text("")
        out = self.dir / "empty_out.tsv"
        run = subprocess.run(
            [sys.executable, str(SCRIPT), "--sample", "S1",
             "--diamond", str(self.dir / "empty.tsv"),
             "--kraken-output", str(self.dir / "kraken.tsv"),
             "--kraken-db", str(self.dir / "does_not_exist"),
             "--output", str(out)], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        with out.open(newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(rows, [])


if __name__ == "__main__":
    unittest.main()
