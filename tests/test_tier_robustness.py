"""V2 robustness columns (Ammal's decision A, 1 Oct): the tier stays on raw reads.

scripts/compute_fragment_support_variants.py counts a sample's clb evidence four ways --
raw reads, duplicates removed, DNA fragments, both -- and the tier each implies. It was an
on-demand reviewer check; it now runs inside classifyPksReadEvidence for every tiered
(short-read wgs) sample, and classify_tumor_pks_evidence.py writes the three alternative
counts and tiers beside read_evidence:

    pks_fragments        tier_if_fragments
    pks_reads_dedup      tier_if_dedup
    pks_fragments_dedup  tier_if_fragments_dedup
    tier_robustness      stable | unstable

The requirement this file guards: read_evidence (and pks_reads) are byte-for-byte what they
were without the variants, whatever the variants say; the new columns are present, filled
from the variants file, and carried into the QC summary, master summary and cohort report.
The variants script's own arithmetic is tests/test_fragment_support_variants.py.
"""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLASSIFY = ROOT / "scripts/classify_tumor_pks_evidence.py"
TARGETED = (ROOT / "Modules/pks_targeted.nf").read_text()
sys.path.insert(0, str(ROOT / "scripts"))
import build_qc_summary as bqs  # noqa: E402
import classify_tumor_pks_evidence as classify  # noqa: E402

ISLAND = 1000
NEW = classify.ROBUSTNESS_FIELDS


def inputs(directory, reads=120, genes=12, covered=250):
    qc = directory / "qc.tsv"
    qc.write_text(f"Sample\tMetric\tValue\nS1\treads_clb_genes_align\t{reads}\n")
    counts = directory / "counts.tsv"
    counts.write_text("Geneid\tCount\n" + "".join(
        f"clb{x}\t{1 if i < genes else 0}\n" for i, x in enumerate("ABCDEFGHIJKLMNOPQRS")))
    depth = directory / "depth.tsv"
    depth.write_text("".join(f"NC_017628.1\t{p + 1}\t{1 if p < covered else 0}\n" for p in range(ISLAND)))
    return qc, counts, depth


def variants_file(directory, rows, verdict="PASS"):
    """compute_fragment_support_variants.py's report layout (write_report)."""
    path = directory / "variants.tsv"
    lines = ["sample\tvariant\treads_or_fragments\tgenes_detected\tisland_breadth_1x\ttier"]
    lines += [f"S1\t{name}\t{count}\t12\t0.250000\t{tier}" for name, count, tier in rows]
    path.write_text("\n".join(lines) + f"\n\n# verdict\t{verdict}\n")
    return path


def run(directory, *extra, mode="wgs", reads=120):
    qc, counts, depth = inputs(directory, reads=reads)
    out, fragment = directory / "evidence.tsv", directory / "evidence.qc.tsv"
    result = subprocess.run(
        [sys.executable, str(CLASSIFY), "--sample", "S1", "--qc", str(qc), "--counts", str(counts),
         "--depth", str(depth), "--island-length", str(ISLAND), "--mode", mode,
         "--depletion-qc", str(write_depletion(directory)), "--output", str(out),
         "--qc-fragment", str(fragment), *extra], capture_output=True, text=True)
    if result.returncode:
        return result, None, None
    with out.open() as handle:
        row = list(csv.DictReader(handle, delimiter="\t"))[0]
    with fragment.open() as handle:
        metrics = {r["Metric"]: r["Value"] for r in csv.DictReader(handle, delimiter="\t")}
    return result, row, metrics


def write_depletion(directory):
    path = directory / "depletion.tsv"
    path.write_text("Sample\tMetric\tValue\nS1\treads_after_t2t_phix\t1000\n")
    return path


STABLE = [("raw", 120, "extensive_island"), ("dup_marked", 110, "extensive_island"),
          ("fragment", 101, "extensive_island"), ("dup_collapsed", 100, "extensive_island")]
UNSTABLE = [("raw", 120, "extensive_island"), ("dup_marked", 60, "broad_island"),
            ("fragment", 101, "extensive_island"), ("dup_collapsed", 55, "broad_island")]


class TheTierStaysOnRawReads(unittest.TestCase):
    def test_tier_and_counts_are_unchanged_by_any_variant(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _, without, _ = run(d)
            for rows in (STABLE, UNSTABLE):
                with self.subTest(rows=rows[1]):
                    _, row, _ = run(d, "--fragment-variants", str(variants_file(d, rows)))
                    for column in classify.READ_EVIDENCE_FIELDS:
                        self.assertEqual(row[column], without[column], column)
                    self.assertEqual(row["read_evidence"], "extensive_island")

    def test_the_new_columns_are_present_and_filled(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _, row, metrics = run(d, "--fragment-variants", str(variants_file(d, UNSTABLE)))
        self.assertEqual(list(row)[-len(NEW):], list(NEW))
        self.assertEqual((row["pks_fragments"], row["tier_if_fragments"]), ("101", "extensive_island"))
        self.assertEqual((row["pks_reads_dedup"], row["tier_if_dedup"]), ("60", "broad_island"))
        self.assertEqual((row["pks_fragments_dedup"], row["tier_if_fragments_dedup"]), ("55", "broad_island"))
        self.assertEqual(row["tier_robustness"], "unstable")
        self.assertEqual(metrics["tier_robustness"], "unstable")
        self.assertEqual(metrics["pks_reads_dedup"], "60")
        self.assertEqual(metrics["evidence_outcome"], "positive")

    def test_stable_when_every_count_gives_the_same_tier(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _, row, _ = run(d, "--fragment-variants", str(variants_file(d, STABLE)))
        self.assertEqual(row["tier_robustness"], "stable")

    def test_a_raw_row_that_disagrees_with_read_evidence_is_refused(self):
        drifted = [("raw", 119, "extensive_island")] + STABLE[1:]
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            result, _, _ = run(d, "--fragment-variants", str(variants_file(d, drifted)))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("disagrees with read_evidence", result.stderr)

    def test_no_reads_is_zero_negative_stable_without_running_anything(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, row, metrics = run(Path(tmp), reads=0)
        self.assertEqual(row["read_evidence"], "negative")
        self.assertEqual((row["pks_fragments"], row["tier_if_fragments"], row["tier_robustness"]),
                         ("0", "negative", "stable"))
        self.assertEqual(metrics["tier_if_dedup"], "negative")

    def test_untiered_rows_say_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, row, metrics = run(Path(tmp), "--no-tiers")
        self.assertEqual(row["read_evidence"], "untiered")
        self.assertEqual({row[c] for c in NEW}, {"NA"})
        self.assertNotIn("tier_robustness", metrics)

    def test_other_modes_have_no_robustness_columns(self):
        for mode in ("metagenome", "wes", "rna"):
            with self.subTest(mode=mode):
                self.assertFalse(set(NEW) & set(classify.SCHEMAS[mode]))


class TheProcessRunsTheVariants(unittest.TestCase):
    BLOCK = TARGETED.split("process classifyPksReadEvidence", 1)[1].split("\nprocess ", 1)[0]

    def test_only_for_tiered_samples_and_only_with_reads(self):
        self.assertIn("def tiered = params.sample_type.toString() == 'wgs' && !HostDepletion.isLong(params.read_type)",
                      self.BLOCK)
        self.assertIn('if [[ "${tiered}" == "true" ]] && awk', self.BLOCK)
        self.assertIn('\\$2 == "reads_clb_genes_align" && \\$3 + 0 > 0', self.BLOCK)

    def test_same_thresholds_island_depth_and_overlap_as_the_tier(self):
        call = self.BLOCK.split('python "${params.scripts}/compute_fragment_support_variants.py"', 1)[1]
        call = call.split("\n", 1)[0]
        for needle in ('--coverage "${sampleID}.island.raw_depth.tsv"', "--min-overlap ${params.min_gene_overlap_bp}",
                       '--pks-contig "${params.pks_contig}"', "--pks-start ${params.pks_start_1based}",
                       "--pks-end ${params.pks_end_1based}", "${thresholds}"):
            self.assertIn(needle, call)
        for param in ("tumor_multi_gene_min_pks_reads", "tumor_broad_island_min_breadth",
                      "tumor_extensive_island_min_clb_genes"):
            self.assertIn(param, self.BLOCK.split("def thresholds", 1)[1].split("\n", 1)[0])
        self.assertIn('"\\${VARIANT_ARGS[@]}"', self.BLOCK)


class TheColumnsReachTheTables(unittest.TestCase):
    def fragment(self, tmp, metrics):
        path = Path(tmp) / f"f{len(list(Path(tmp).iterdir()))}.tsv"
        path.write_text("Sample\tMetric\tValue\n" + "".join(f"S1\t{k}\t{v}\n" for k, v in metrics.items()))
        return path

    BASE = {"filter_input_reads": 10, "reads_after_fastp": 10, "reads_after_hg38": 10,
            "reads_after_t2t_phix": 10, "reads_clb_genes_align": 5, "num_clb_genes_align": 3}

    def test_qc_summary_gets_the_columns_only_when_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            plain = bqs.load_fragments([self.fragment(tmp, self.BASE)])
            self.assertFalse(set(NEW) & set(bqs.table_columns([], plain)))
            with_v2 = bqs.load_fragments([self.fragment(tmp, dict(
                self.BASE, pks_fragments=4, tier_if_fragments="multi_gene", pks_reads_dedup=5,
                tier_if_dedup="multi_gene", pks_fragments_dedup=4,
                tier_if_fragments_dedup="localized_indeterminate", tier_robustness="unstable"))])
            columns = bqs.table_columns([], with_v2)
            self.assertEqual([c for c in columns if c in NEW], list(NEW))
            self.assertEqual(bqs.output_row("S1", with_v2["S1"])["tier_if_fragments_dedup"],
                             "localized_indeterminate")

    def test_a_tier_value_is_checked(self):
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(ValueError):
            bqs.load_fragments([self.fragment(tmp, dict(self.BASE, tier_if_dedup="probably"))])

    def test_master_summary_and_cohort_report_carry_them(self):
        master = (ROOT / "scripts/build_master_summary.py").read_text()
        self.assertIn('("tier robustness", TIER_ROBUSTNESS_COLUMNS, read_tier_robustness)', master)
        report = (ROOT / "scripts/build_cohort_report.py").read_text()
        self.assertIn("tier_robustness_section(records)", report)
        self.assertIn("fields += list(ROBUSTNESS_FIELDS)", report)
        self.assertEqual(tuple(bqs.TIER_ROBUSTNESS_COLUMNS), tuple(NEW))
        import build_cohort_report
        self.assertEqual(tuple(build_cohort_report.ROBUSTNESS_FIELDS), tuple(NEW))

    def test_documented(self):
        doc = (ROOT / "docs/output.md").read_text()
        for column in NEW:
            self.assertIn(column, doc)


if __name__ == "__main__":
    unittest.main()
