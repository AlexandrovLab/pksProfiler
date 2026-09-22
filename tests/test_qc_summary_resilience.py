#!/usr/bin/env python3
"""The cohort QC table must survive a sample that failed or never ran.

Losing 96 samples' QC because one hit a time cap is the failure these guard
against. Corruption stays fatal; incompleteness gets recorded.
"""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_qc_summary", ROOT / "scripts/build_qc_summary.py")
build_qc_summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build_qc_summary)

COMPLETE = [
    ("input_alignment_records", 100),
    ("extracted_unmapped_reads", 90),
    ("filter_input_reads", 90),
    ("reads_after_fastp", 80),
    ("reads_after_hg38", 70),
    ("reads_after_t2t_phix", 60),
]


def fragment(directory, sample, metrics):
    path = Path(directory) / f"{sample}.qc.tsv"
    lines = ["Sample\tMetric\tValue"]
    lines += [f"{sample}\t{metric}\t{value}" for metric, value in metrics]
    path.write_text("\n".join(lines) + "\n")
    return path


def summarise(directory, fragments, expected=None):
    out = Path(directory) / "summary.tsv"
    exp = None
    if expected is not None:
        exp = Path(directory) / "expected.txt"
        exp.write_text("\n".join(expected) + "\n")
    metrics = build_qc_summary.load_fragments(fragments)
    build_qc_summary.write_summary(metrics, out, build_qc_summary.read_expected(exp))
    with out.open(newline="") as handle:
        return {row["Sample"]: row for row in csv.DictReader(handle, delimiter="\t")}


class PartialCohortTests(unittest.TestCase):
    def test_one_failed_sample_does_not_lose_the_others(self):
        with tempfile.TemporaryDirectory() as d:
            good = fragment(d, "good", COMPLETE)
            rows = summarise(d, [good], expected=["good", "failed"])
            self.assertEqual(sorted(rows), ["failed", "good"])
            self.assertEqual(rows["good"]["status"], "complete")
            self.assertEqual(rows["good"]["reads_after_t2t_phix"], "60")

    def test_a_sample_that_produced_nothing_is_still_reported(self):
        with tempfile.TemporaryDirectory() as d:
            good = fragment(d, "good", COMPLETE)
            rows = summarise(d, [good], expected=["good", "failed"])
            self.assertEqual(rows["failed"]["status"], "no_qc_produced")
            self.assertEqual(rows["failed"]["input_alignment_records"], "NA")

    def test_a_partially_run_sample_is_marked_incomplete_not_fatal(self):
        with tempfile.TemporaryDirectory() as d:
            partial = fragment(d, "partial", COMPLETE[:2])
            rows = summarise(d, [partial], expected=["partial"])
            self.assertEqual(rows["partial"]["status"], "incomplete")
            self.assertEqual(rows["partial"]["reads_after_t2t_phix"], "NA")

    def test_every_expected_sample_gets_exactly_one_row(self):
        with tempfile.TemporaryDirectory() as d:
            good = fragment(d, "good", COMPLETE)
            rows = summarise(d, [good], expected=["good", "a", "b", "c"])
            self.assertEqual(sorted(rows), ["a", "b", "c", "good"])

    def test_no_fragments_at_all_still_produces_the_table(self):
        with tempfile.TemporaryDirectory() as d:
            rows = summarise(d, [], expected=["a", "b"])
            self.assertEqual(sorted(rows), ["a", "b"])
            self.assertTrue(all(r["status"] == "no_qc_produced" for r in rows.values()))


class CorruptionIsStillFatalTests(unittest.TestCase):
    def test_conflicting_values_for_one_metric_still_raise(self):
        with tempfile.TemporaryDirectory() as d:
            a = fragment(d, "s", [("reads_after_fastp", 10)])
            b = fragment(d, "s2", [("reads_after_fastp", 20)])
            b.write_text(b.read_text().replace("s2", "s"))
            with self.assertRaises(ValueError):
                build_qc_summary.load_fragments([a, b])

    def test_counts_increasing_down_the_chain_still_raise(self):
        with tempfile.TemporaryDirectory() as d:
            bad = fragment(d, "s", [
                ("input_alignment_records", 10),
                ("extracted_unmapped_reads", 10),
                ("filter_input_reads", 10),
                ("reads_after_fastp", 10),
                ("reads_after_hg38", 10),
                ("reads_after_t2t_phix", 999),
            ])
            metrics = build_qc_summary.load_fragments([bad])
            with self.assertRaises(ValueError):
                build_qc_summary.output_row("s", metrics["s"])

    def test_an_out_of_range_gene_count_still_raises(self):
        with tempfile.TemporaryDirectory() as d:
            bad = fragment(d, "s", COMPLETE + [("num_clb_genes_align", 42)])
            metrics = build_qc_summary.load_fragments([bad])
            with self.assertRaises(ValueError):
                build_qc_summary.output_row("s", metrics["s"])

    def test_nothing_supplied_and_nothing_expected_still_raises(self):
        self.assertEqual(build_qc_summary.read_expected(None), [])


if __name__ == "__main__":
    unittest.main()
