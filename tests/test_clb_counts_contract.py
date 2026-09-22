"""One contract for a per-sample clb count table, held by both methods.

F06, from Ludmil's 2026-09-19 report: "The HMM merger accepts an unknown gene (`BAD`),
negative counts, a missing count converted to zero, a header-only sample, and
fractional counts. Duplicate genes are grouped and summed. Alignment validators are
stronger, but still accept fractional values even though this pipeline does not
request fractional feature assignment."

His acceptance test: "All five malformed HMM fixtures and fractional alignment counts
must fail. Valid zero samples must pass. Duplicate genes, duplicate sample IDs, extra
columns, and incomplete files need explicit tests."
"""
import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

spec = importlib.util.spec_from_file_location(
    "clb_counts_schema", ROOT / "scripts/clb_counts_schema.py")
schema = importlib.util.module_from_spec(spec)
spec.loader.exec_module(schema)

GENES = schema.CLB_GENES
HMM = ROOT / "scripts/build_hmm_matrix.py"
MERGE = ROOT / "scripts/mergeGeneCounts.py"


class TheContract(unittest.TestCase):
    def check(self, rows):
        return schema.validate_clb_counts(rows, "sample.tsv")

    def test_nineteen_genes_with_counts_pass(self):
        counts = self.check([(g, i) for i, g in enumerate(GENES)])
        self.assertEqual(list(counts), GENES)          # canonical order
        self.assertEqual(counts["clbA"], 0)

    def test_a_sample_of_all_zeros_is_valid_and_not_empty(self):
        # The distinction the finding asks for: a negative sample is nineteen zeros.
        self.assertEqual(set(self.check([(g, 0) for g in GENES]).values()), {0})

    def test_an_empty_file_is_not_a_zero_sample(self):
        with self.assertRaises(schema.CountSchemaError) as caught:
            self.check([])
        self.assertIn("truncated", str(caught.exception))

    def test_an_unknown_gene_is_refused(self):
        with self.assertRaises(schema.CountSchemaError):
            self.check([(g, 1) for g in GENES] + [("BAD", 1)])

    def test_a_missing_gene_is_refused_rather_than_filled_with_zero(self):
        with self.assertRaises(schema.CountSchemaError) as caught:
            self.check([(g, 1) for g in GENES[:-1]])
        self.assertIn("clbS", str(caught.exception))

    def test_a_negative_count_is_refused(self):
        with self.assertRaises(schema.CountSchemaError):
            self.check([(g, -1 if g == "clbA" else 1) for g in GENES])

    def test_a_fractional_count_is_refused(self):
        with self.assertRaises(schema.CountSchemaError):
            self.check([(g, 3.7 if g == "clbA" else 1) for g in GENES])

    def test_a_duplicate_gene_is_refused_not_summed(self):
        with self.assertRaises(schema.CountSchemaError) as caught:
            self.check([(g, 1) for g in GENES] + [("clbA", 9)])
        self.assertIn("not a count to be summed", str(caught.exception))

    def test_an_empty_field_is_not_a_zero(self):
        with self.assertRaises(schema.CountSchemaError):
            self.check([(g, "" if g == "clbA" else 1) for g in GENES])

    def test_a_non_numeric_count_is_refused(self):
        with self.assertRaises(schema.CountSchemaError):
            self.check([(g, "many" if g == "clbA" else 1) for g in GENES])


class BothMergersHoldIt(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.d = Path(self.tmp.name)

    def hmm_file(self, name, rows):
        path = self.d / f"{name}.hmm_counts.tsv"
        path.write_text("Gene\tCount\n" + "".join(f"{g}\t{v}\n" for g, v in rows))
        return path

    def featurecounts_file(self, name, values):
        path = self.d / f"{name}.counts.txt"
        rows = ["Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample"]
        rows += [f"{g}\tctg\t1\t10\t+\t10\t{v}" for g, v in zip(GENES, values)]
        path.write_text("\n".join(rows) + "\n")
        return path

    def run_hmm(self, path):
        return subprocess.run(
            [sys.executable, str(HMM), "--inputs", str(path), "--out",
             str(self.d / "m.tsv"), "--strip-suffix", ".hmm_counts.tsv"],
            capture_output=True, text=True)

    def run_merge(self, path):
        return subprocess.run(
            [sys.executable, str(MERGE), str(path), str(self.d / "a.tsv")],
            capture_output=True, text=True)

    def test_hmm_accepts_a_valid_zero_sample(self):
        result = self.run_hmm(self.hmm_file("zero", [(g, 0) for g in GENES]))
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_hmm_refuses_each_malformed_fixture(self):
        fixtures = {
            "unknown": [(g, 1) for g in GENES] + [("BAD", 1)],
            "negative": [(g, -1 if g == "clbA" else 1) for g in GENES],
            "missing": [(g, 1) for g in GENES[:-1]],
            "header_only": [],
            "fractional": [(g, 3.7 if g == "clbA" else 1) for g in GENES],
            "duplicate": [(g, 1) for g in GENES] + [("clbA", 9)],
        }
        for name, rows in fixtures.items():
            with self.subTest(fixture=name):
                result = self.run_hmm(self.hmm_file(name, rows))
                self.assertNotEqual(result.returncode, 0,
                                    f"{name} was accepted: {result.stdout}")

    def test_hmm_rows_come_out_in_canonical_gene_order(self):
        path = self.hmm_file("ordered", [(g, 1) for g in reversed(GENES)])
        self.assertEqual(self.run_hmm(path).returncode, 0)
        with (self.d / "m.tsv").open(newline="") as handle:
            rows = [row[0] for row in csv.reader(handle, delimiter="\t")][1:]
        self.assertEqual(rows, GENES)

    def test_the_alignment_merger_now_refuses_fractional_counts(self):
        path = self.featurecounts_file("frac", [2.5] + [1] * 18)
        self.assertNotEqual(self.run_merge(path).returncode, 0)

    def test_the_alignment_merger_still_accepts_a_real_table(self):
        path = self.featurecounts_file("ok", list(range(19)))
        self.assertEqual(self.run_merge(path).returncode, 0)

    def test_both_mergers_import_the_same_validator(self):
        for script in (HMM, MERGE):
            with self.subTest(script=script.name):
                self.assertIn("from clb_counts_schema import", script.read_text())


if __name__ == "__main__":
    unittest.main()
