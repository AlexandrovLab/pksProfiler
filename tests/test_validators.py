#!/usr/bin/env python3

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPO = Path(__file__).resolve().parents[1]
GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]


def load_script(name):
    path = REPO / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ValidatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bedgraph = load_script("validate_bedgraph")
        cls.featurecounts = load_script("validate_featurecounts")
        cls.cram = load_script("validate_cram_reference")

    def test_bedgraph_accepts_scientific_notation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "coverage.bedgraph"
            path.write_text("contig\t0\t10\t1.2e-07\ncontig\t10\t20\t3E+2\n")
            self.bedgraph.validate(path)

    def test_bedgraph_rejects_invalid_values_and_empty_data(self):
        for content in ("", "contig\t0\t10\tNaN\n", "contig\t0\t10\tInf\n", "contig\t0\t10\t-1e-3\n"):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "coverage.bedgraph"
                path.write_text(content)
                with self.assertRaises(ValueError):
                    self.bedgraph.validate(path)

    def test_annotation_requires_exact_unique_19_gene_set(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "clb.gff"
            path.write_text("".join(
                f"ctg\ttest\tgene\t{i}\t{i + 9}\t.\t+\t.\tID={gene};Name={gene}\n"
                for i, gene in enumerate(reversed(GENES), 1)
            ))
            self.featurecounts.annotation_genes(path)
            path.write_text(path.read_text().replace("Name=clbA", "Name=clbB"))
            with self.assertRaises(ValueError):
                self.featurecounts.annotation_genes(path)

    def test_cram_reference_match_and_mismatch(self):
        cram_header = "@SQ\tSN:ctg\tLN:100\tM5:abc\n"
        matching_dict = "@SQ\tSN:ctg\tLN:100\tM5:abc\n"
        mismatch_dict = "@SQ\tSN:ctg\tLN:101\tM5:def\n"

        def result(text):
            return subprocess.CompletedProcess([], 0, stdout=text, stderr="")

        with patch.object(self.cram.subprocess, "run", side_effect=[result(cram_header), result(matching_dict)]):
            self.cram.validate("sample.cram", "ref.fa")
        with patch.object(self.cram.subprocess, "run", side_effect=[result(cram_header), result(mismatch_dict)]):
            with self.assertRaises(ValueError):
                self.cram.validate("sample.cram", "ref.fa")

    def test_merge_rejects_a_missing_gene(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            counts = directory / "sample.counts.txt"
            output = directory / "merged.tsv"
            rows = ["Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample"]
            rows.extend(f"{gene}\tctg\t1\t10\t+\t10\t0" for gene in GENES[:-1])
            counts.write_text("\n".join(rows) + "\n")
            result = subprocess.run(
                [sys.executable, REPO / "scripts" / "mergeGeneCounts.py", counts, output],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
