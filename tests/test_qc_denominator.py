"""QC reports a real library denominator, and every column has exactly one meaning.

F09, from Ludmil's 2026-09-19 report. The summary read `bam_input_primary_records`
and **nothing ever wrote it**, so for BAM/CRAM the `input_reads` column silently fell
back to the post-extraction count. A column that is sometimes the library and
sometimes what survived extraction cannot show how much of the library was dropped --
which is why F01 went unnoticed in QC for as long as it did.

The first test here is the one that matters: it compares what the summary consumes
against what the modules emit. The old test suite had a fixture that fed
`bam_input_primary_records` in by hand and asserted the summary handled it, which
proved the reader worked and never asked whether anything wrote it.
"""
import csv
import importlib.util
import re
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_qc_summary", ROOT / "scripts/build_qc_summary.py")
build_qc_summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build_qc_summary)

QC_SUMMARY_SOURCE = (ROOT / "scripts/build_qc_summary.py").read_text()
MODULE_SOURCES = {p.name: p.read_text() for p in (ROOT / "Modules").glob("*.nf")}


def emitted_metrics():
    """Metric names any module writes into a QC fragment."""
    names = set()
    for text in MODULE_SOURCES.values():
        # The value may be a %s or a literal -- the zero-hit branch of the HMM module
        # writes `printf "%s\treads_clb_genes_hmm\t0\n"`. Match on the metric name only.
        names.update(re.findall(r'printf\s+"%s\\*t(\w+)', text))
    return names


def consumed_metrics():
    """Metric names the summary reads out of the fragments."""
    return set(re.findall(r'values\.get\("(\w+)"', QC_SUMMARY_SOURCE))


class MetricsExist(unittest.TestCase):
    def test_every_metric_the_summary_reads_is_written_by_some_module(self):
        missing = sorted(consumed_metrics() - emitted_metrics())
        self.assertEqual(missing, [], f"consumed but never emitted: {missing}")

    def test_extraction_reports_a_denominator_not_derived_from_itself(self):
        extract = MODULE_SOURCES["extract_reads.nf"]
        self.assertIn("input_alignment_records", extract)
        self.assertIn("idxstats", extract)

    def test_the_exact_count_is_opt_in_because_it_costs_a_second_decode(self):
        extract = MODULE_SOURCES["extract_reads.nf"]
        self.assertIn("total_primary_reads", extract)
        self.assertIn("params.exact_input_counts", extract)
        self.assertIn("params.exact_input_counts = false", (ROOT / "main.nf").read_text())

    def test_the_reference_mapped_count_is_exported_not_only_tested(self):
        # F09 asked for it to be exported; the align module writes it into the
        # per-sample QC fragment.
        align = MODULE_SOURCES["pksProfiler_align.nf"]
        self.assertIn("reads_mapped_ihe3034", align)

    def test_it_is_not_a_column_in_the_cohort_table(self):
        # Reads mapping anywhere on the 5.1 Mb genome, of which the island is 1%.
        # Beside the clb columns it reads as though it were the same measurement.
        self.assertNotIn('"reads_mapped_ihe3034",', QC_SUMMARY_SOURCE.split("OUTPUT_COLUMNS = [", 1)[1].split("]", 1)[0])
        self.assertIn("reads_mapped_ihe3034", QC_SUMMARY_SOURCE)  # still bounds the clb count


class NoSubstitution(unittest.TestCase):
    """A missing metric is NA. It is never filled in from a different stage."""

    def summarise(self, metrics, sample="S1"):
        with tempfile.TemporaryDirectory() as directory:
            fragment = Path(directory) / f"{sample}.qc.tsv"
            lines = ["Sample\tMetric\tValue"]
            lines += [f"{sample}\t{metric}\t{value}" for metric, value in metrics]
            fragment.write_text("\n".join(lines) + "\n")
            output = Path(directory) / "summary.tsv"
            build_qc_summary.write_summary(
                build_qc_summary.load_fragments([fragment]), output, [sample])
            with output.open(newline="") as handle:
                return list(csv.DictReader(handle, delimiter="\t"))[0]

    def test_a_missing_alignment_denominator_stays_na(self):
        row = self.summarise([
            ("filter_input_reads", 400), ("reads_after_fastp", 360),
            ("reads_after_hg38", 300), ("reads_after_t2t_phix", 250),
        ])
        self.assertEqual(row["input_alignment_records"], "NA")
        self.assertEqual(row["total_primary_reads"], "NA")
        self.assertEqual(row["extracted_unmapped_reads"], "NA")
        self.assertEqual(row["filter_input_reads"], "400")

    def test_the_denominator_and_the_extracted_count_are_both_kept(self):
        row = self.summarise([
            ("input_alignment_records", 1000), ("extracted_unmapped_reads", 200),
            ("filter_input_reads", 200), ("reads_after_fastp", 180),
            ("reads_after_hg38", 150), ("reads_after_t2t_phix", 120),
        ])
        self.assertEqual(row["input_alignment_records"], "1000")
        self.assertEqual(row["extracted_unmapped_reads"], "200")
        self.assertEqual(row["status"], "complete")

    def test_the_chain_still_refuses_counts_that_grow(self):
        with self.assertRaises(ValueError):
            self.summarise([
                ("input_alignment_records", 100), ("extracted_unmapped_reads", 200),
                ("filter_input_reads", 200), ("reads_after_fastp", 180),
                ("reads_after_hg38", 150), ("reads_after_t2t_phix", 120),
            ])

    def test_the_mapped_count_bounds_the_clb_count(self):
        with self.assertRaises(ValueError):
            self.summarise([
                ("filter_input_reads", 400), ("reads_after_fastp", 360),
                ("reads_after_hg38", 300), ("reads_after_t2t_phix", 250),
                ("reads_mapped_ihe3034", 10), ("reads_clb_genes_align", 99),
            ])


if __name__ == "__main__":
    unittest.main()
