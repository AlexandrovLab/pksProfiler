#!/usr/bin/env python3
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

STAGE = Path(__file__).resolve().parents[1]

def load(name):
    path = STAGE / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

route = load("route_kraken_reads")
import sys
sys.modules["route_kraken_reads"] = route
diamond = load("select_diamond_reads")

class PrefilterTests(unittest.TestCase):
    def test_kraken_routes_mates_independently(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d / "nodes.dmp").write_text("1\t|\t1\t|\n91347\t|\t1\t|\n543\t|\t91347\t|\n562\t|\t543\t|\n")
            (d / "kraken.tsv").write_text("C\tpair/1\t562\t100\tA:1\nU\tpair/2\t0\t100\tA:0\n")
            (d / "reads.fastq").write_text("@pair/1\nACGT\n+\nIIII\n@pair/2\nTGCA\n+\nIIII\n")
            counts = route.route(d / "reads.fastq", d / "kraken.tsv", d / "nodes.dmp", 91347, d / "primary.gz", d / "other.gz")
            self.assertEqual(counts["primary"], 1)
            self.assertEqual(counts["non_target"], 1)

    def test_diamond_thresholds_and_ids_are_read_specific(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d / "hits.tsv").write_text(
                "pair/1\tclbS\t90\t30\t100\t170\t1e-20\t80\n"
                "pair/2\tclbS\t90\t10\t100\t170\t1e-20\t30\n"
            )
            selected = diamond.load_hits(d / "hits.tsv", 1e-5, 25, 0.5)
            self.assertEqual(selected, {"pair/1"})

class ReferenceTests(unittest.TestCase):
    def test_clb_protein_reference_has_exactly_a_through_s(self):
        fasta = STAGE / "ref" / "clb_reference_proteins.faa"
        labels = [line[1:].split("|", 1)[0] for line in fasta.read_text().splitlines() if line.startswith(">")]
        self.assertEqual(sorted(labels), [f"clb{x}" for x in "ABCDEFGHIJKLMNOPQRS"])

if __name__ == "__main__":
    unittest.main()
