"""A missing metric must not crash the cohort QC summary.

`output_row` compared reads_clb_genes_hmm against the depletion metric while guarding only
the first operand for None. A sample missing its depletion metric raised

    TypeError: '>' not supported between instances of 'int' and 'NoneType'

and took masterQCSummary down with it -- and every cohort table is a single task, so the
whole cohort's QC summary is lost with it.

That is reachable in production, not just in a hand-assembled run: `sample_failure_strategy`
defaults to `ignore` for the per-sample lane, so a failed filterReads leaves the sample with
no reads_after_t2t_phix and the run continues to the summary, which then dies. The strategy
exists to stop one bad sample killing the cohort tables; this crash defeated it.

Found on 2026-09-22 building the cohort tables by hand from a partially complete run.
"""
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("bqs", ROOT / "scripts/build_qc_summary.py")
bqs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bqs)


def row(**overrides):
    values = {c: None for c in bqs.OUTPUT_COLUMNS}
    values.update(overrides)
    return values


class MissingMetricsDoNotCrash(unittest.TestCase):
    def test_hmm_reads_without_a_depletion_metric(self):
        out = bqs.output_row("S1", row(reads_clb_genes_hmm=12))
        self.assertEqual(out["reads_clb_genes_hmm"], 12)
        self.assertEqual(out["reads_after_t2t_phix"], "NA")

    def test_the_impossible_count_check_still_fires_when_both_exist(self):
        with self.assertRaises(ValueError):
            bqs.output_row("S1", row(reads_clb_genes_hmm=500,
                                     reads_after_t2t_phix=100))

    def test_a_consistent_pair_passes(self):
        out = bqs.output_row("S1", row(reads_clb_genes_hmm=10,
                                       reads_after_t2t_phix=1000))
        self.assertEqual(out["reads_clb_genes_hmm"], 10)

    def test_a_row_with_nothing_in_it_still_produces_every_column(self):
        out = bqs.output_row("S1", row())
        self.assertEqual(set(out), set(bqs.OUTPUT_COLUMNS))
