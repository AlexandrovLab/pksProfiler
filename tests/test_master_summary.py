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
          ["sample", "final_structural_evidence", "assembler_agreement",
           "megahit_reference_covered_bp", "megahit_reference_coverage",
           "megahit_supporting_contigs", "metaspades_reference_covered_bp",
           "metaspades_reference_coverage", "metaspades_supporting_contigs",
           "recruited_fragment_ids", "paired_fragments"],
          [["S1", "near_complete_island", "both", "412", "0.812000", "3",
            "398", "0.789000", "2", "9", "3"]])
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
        self.assertEqual(s1["megahit_reference_covered_bp"], "412")
        self.assertEqual(s1["megahit_reference_coverage"], "0.812000")
        self.assertEqual(s1["megahit_supporting_contigs"], "3")
        self.assertEqual(s1["metaspades_reference_covered_bp"], "398")
        self.assertEqual(s1["metaspades_reference_coverage"], "0.789000")
        self.assertEqual(s1["metaspades_supporting_contigs"], "2")
        self.assertEqual(s1["recruited_fragment_ids"], "9")
        self.assertEqual(s1["paired_fragments"], "3")
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

    def test_unbinned_row_is_reported_but_never_counted_as_a_bin(self):
        """U1: pks_mag_summary.tsv can carry an extra unit_type=unbinned row for the
        per-sample pool of contigs MetaBAT2 never placed in any bin -- a positive call
        there must reach the table (a real detection blind spot otherwise) but must
        never inflate mag_bins_total/mag_bins_pks_positive, which describe recovered
        genomes specifically.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); minimal(root); add_everything(root)
            write(root / "by_sample/S1/genomes/pks_mag_summary.tsv",
                  ["sample", "bin_id", "unit_type", "classification", "completeness",
                   "distinct_clb_genes", "locus_tier"],
                  [["S1", "bin.1", "bin", "d__Bacteria;s__Escherichia coli", "96.4", "17", "extensive_island"],
                   ["S1", "bin.2", "bin", "d__Bacteria;s__Bacteroides fragilis", "88.1", "0", "negative"],
                   ["S1", "unbinned", "unbinned", "NA", "NA", "NA", "broad_island"]])
            rows, _ = run(root)
        s1 = next(r for r in rows if r["sample"] == "S1")
        # Still 2, not 3: the unbinned row is not a bin.
        self.assertEqual(s1["mag_bins_total"], "2")
        self.assertEqual(s1["mag_bins_pks_positive"], "1")
        self.assertEqual(s1["mag_unbinned_locus_tier"], "broad_island")
        self.assertEqual(s1["mag_unbinned_pks_positive"], "yes")

    def test_unbinned_columns_are_na_when_there_is_no_unbinned_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); minimal(root); add_everything(root)
            rows, _ = run(root)
        s1 = next(r for r in rows if r["sample"] == "S1")
        self.assertEqual(s1["mag_unbinned_locus_tier"], "NA")
        self.assertEqual(s1["mag_unbinned_pks_positive"], "NA")

    def test_a_negative_unbinned_tier_reads_as_not_positive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); minimal(root); add_everything(root)
            write(root / "by_sample/S1/genomes/pks_mag_summary.tsv",
                  ["sample", "bin_id", "unit_type", "locus_tier"],
                  [["S1", "unbinned", "unbinned", "negative"]])
            rows, _ = run(root)
        s1 = next(r for r in rows if r["sample"] == "S1")
        self.assertEqual(s1["mag_bins_total"], "0")
        self.assertEqual(s1["mag_unbinned_locus_tier"], "negative")
        self.assertEqual(s1["mag_unbinned_pks_positive"], "no")

    def test_empty_results_directory_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run([sys.executable, str(SCRIPT), "--results", tmp,
                                "--output", os.path.join(tmp, "o.tsv")],
                               capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("No per-sample results", r.stderr)


class TheSummaryNeverShowsAStaleSample(unittest.TestCase):
    """Automating this script (dead-code item d-2) reintroduces finding 5's exact
    bug class -- glob-by-path admits a directory left over from an older run at
    the same --results -- unless it gets the same --expected-samples filter
    build_cohort_report.py already has.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        minimal(self.root)
        stale = self.root / "by_sample/S_stale_old_run"
        write(stale / "read_evidence.tsv",
              ["sample", "read_evidence", "pks_reads", "clb_genes_detected", "island_breadth_1x"],
              [["S_stale_old_run", "extensive_island", "9999", "19", "0.99"]])
        self.expected = self.root / "expected_samples.txt"
        self.expected.write_text("S1\nS2\n")

    def run_with_expected(self):
        out = self.root / "master.tsv"
        r = subprocess.run(
            [sys.executable, str(SCRIPT), "--results", str(self.root), "--output", str(out),
             "--expected-samples", str(self.expected)],
            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with out.open() as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_without_the_flag_the_stale_sample_leaks_through(self):
        out = self.root / "master.tsv"
        r = subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.root),
                            "--output", str(out)], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with out.open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertIn("S_stale_old_run", {r["sample"] for r in rows},
                      "fixture assumption: an unfiltered scan does pick it up")

    def test_with_the_flag_the_stale_sample_is_gone(self):
        rows = self.run_with_expected()
        self.assertEqual({r["sample"] for r in rows}, {"S1", "S2"})

    def test_the_flag_is_optional(self):
        # Still safe to run by hand against a finished tree, exactly as before.
        out = self.root / "master.tsv"
        r = subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.root),
                            "--output", str(out)], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)


class Wiring(unittest.TestCase):
    """The process is now automatic, gated on the lanes it actually reads from."""

    PLOTTING = (ROOT / "Modules/plotting.nf").read_text()
    MAIN = (ROOT / "main.nf").read_text()
    PKS_MAG = (ROOT / "Modules/pks_mag.nf").read_text()

    def test_the_process_exists_and_is_invoked(self):
        self.assertIn("process masterSummary", self.PLOTTING)
        self.assertIn("masterSummary(", self.MAIN)

    def test_it_takes_the_same_three_inputs_as_cohortReport(self):
        process = self.PLOTTING.split("process masterSummary {", 1)[1]
        self.assertIn("path(results_marker)", process)
        self.assertIn("path(report_script)", process)
        self.assertIn("path(expected_samples)", process)
        self.assertIn('--expected-samples "${expected_samples}"', process)

    def test_it_runs_after_the_cohort_report(self):
        body = self.MAIN.split("workflow {", 1)[1]
        self.assertLess(body.index("cohortReport("), body.index("masterSummary("))

    def test_its_gate_covers_cohort_report_gate_plus_the_optional_lanes(self):
        body = self.MAIN.split("workflow {", 1)[1]
        call = body[body.index("masterSummary("):body.index("masterSummary(") + 300]
        self.assertIn("masterQCSummary.out", call)
        self.assertIn("cohort_report_gate", call)
        self.assertIn("optional_lane_gate", call)
        self.assertIn("EXPECTED_SAMPLE_IDS", call)

    def test_optional_lane_gate_defaults_empty_outside_if_do_align(self):
        # Declared where cohort_report_gate is, for the identical reason: an
        # HMM-only run enables none of these lanes and genuinely has nothing to
        # wait on.
        body = self.MAIN.split("workflow {", 1)[1]
        decl = body.index("def optional_lane_gate = channel.empty()")
        first_if_do_align = body.index("if (do_align) {")
        self.assertLess(decl, first_if_do_align)

    def test_the_metagenome_mag_lane_feeds_the_gate(self):
        body = self.MAIN.split("workflow {", 1)[1]

        mag_block = body[body.index("if (enable_mags_b) {"):]
        mag_block = mag_block[:mag_block.index("\n    }")]
        self.assertIn("pksMAG.out.mag_summary", mag_block)
        self.assertIn("pksMAG.out.community_summary", mag_block)
        self.assertIn("pksMAG.out.strain_summary", mag_block)

    def test_pksMAG_emits_what_the_gate_and_the_script_both_need(self):
        workflow = self.PKS_MAG[self.PKS_MAG.index("workflow pksMAG {"):]
        emit_block = workflow[workflow.index("\n    emit:"):]
        self.assertIn("mag_summary       = magSummaryTable.out.summary", emit_block)
        self.assertIn("community_summary = communityProphageSummary.out.inventory", emit_block)
        self.assertIn("strain_summary    = strain_typing_summary_ch", emit_block)

    def test_pksMAG_strain_summary_is_empty_when_typing_is_off(self):
        workflow = self.PKS_MAG[self.PKS_MAG.index("workflow pksMAG {"):
                               self.PKS_MAG.index("\n    emit:", self.PKS_MAG.index("workflow pksMAG {"))]
        # Lowercase channel.empty(), not Channel.empty(): Nextflow 24.10.0 (the
        # manifest's declared floor, and the version CI pins) fails to compile
        # `def <newVar> = Channel.something()` -- see m-1's fix in this same commit.
        self.assertIn("def strain_typing_summary_ch = channel.empty()", workflow)
        # Assigned only inside the flag check, matching cohort_report_gate's own pattern.
        flag_block = workflow[workflow.index("if (params.enable_strain_typing"):]
        self.assertIn("strain_typing_summary_ch = magStrainTyping.out.summary", flag_block)


if __name__ == "__main__":
    unittest.main()
