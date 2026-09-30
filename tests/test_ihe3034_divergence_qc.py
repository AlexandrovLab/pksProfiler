"""Ludmil, 25 Sep A5: divergence from IHE3034 of the clb-supporting reads.

MAPQ >= 40 within IHE3034 is not uniqueness across microbes. NM / alignment length
quantifies divergence from IHE3034 and nothing more, so it is reported as a QC
distribution -- never as specificity, uniqueness or carrier identity -- and it must not
alter any primary call. Acceptance: "Per-sample divergence QC is generated without
altering primary calls."

Everything here is synthetic SAM against a synthetic 19-gene annotation.
"""
import csv
import importlib.util
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/ihe3034_divergence_qc.py"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


div = load("ihe3034_divergence_qc", SCRIPT)
bqs = load("bqs_a5", ROOT / "scripts/build_qc_summary.py")

GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]
CONTIG = "NC_017628.1"


def write_gff(path):
    """clbA at 1001-2000, clbB at 2001-3000, ... each 1 kb, back to back."""
    with open(path, "w") as out:
        for index, gene in enumerate(GENES):
            start = 1001 + 1000 * index
            out.write(f"{CONTIG}\tsynthetic\tgene\t{start}\t{start + 999}\t.\t+\t.\t"
                      f"ID=gene-{gene};Name={gene}\n")


def sam(name, pos, cigar, nm, flag=0, mapq=42, contig=CONTIG):
    tags = "" if nm is None else f"\tNM:i:{nm}"
    return f"{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t*\t*{tags}\n"


class Fixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.gff = self.dir / "clb_genes.gff"
        write_gff(self.gff)
        self.genes = div.load_genes(self.gff)

    def tearDown(self):
        self.tmp.cleanup()

    def collect(self, lines, min_mapq=40):
        return div.collect(lines, self.genes, min_mapq)


class ReadSelection(Fixture):
    def test_only_primary_confident_clb_assigned_reads_are_measured(self):
        lines = ["@HD\tVN:1.6\n",
                 sam("keep", 1101, "100M", 1),
                 sam("secondary", 1101, "100M", 1, flag=256),
                 sam("supplementary", 1101, "100M", 1, flag=2048),
                 sam("unmapped", 1101, "*", None, flag=4, mapq=0),
                 sam("low_mapq", 1101, "100M", 1, mapq=39),
                 sam("outside_island", 50, "100M", 1),
                 sam("other_contig", 1101, "100M", 1, contig="plasmid")]
        reads = self.collect(lines)
        self.assertEqual([r[0] for r in reads], ["clbA"])

    def test_largest_overlap_assigns_like_featurecounts_and_ties_are_dropped(self):
        # 1951-2050 spans clbA/clbB 50/50: a tie, unassigned. 1961-2060: 40 clbA, 60 clbB.
        reads = self.collect([sam("tie", 1951, "100M", 0), sam("mostly_b", 1961, "100M", 0)])
        self.assertEqual([r[0] for r in reads], ["clbB"])

    def test_a_missing_nm_tag_is_an_error_not_a_zero(self):
        with self.assertRaisesRegex(ValueError, "no NM tag"):
            self.collect([sam("r", 1101, "100M", None)])

    def test_an_nm_larger_than_the_alignment_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "alignment columns"):
            self.collect([sam("r", 1101, "10M", 11)])


class PerReadArithmetic(unittest.TestCase):
    def test_columns_length_and_span(self):
        # 5S 50M 2I 40M 3D 8M 4H: columns 50+2+40+3+8, read 5+50+2+40+8+4, span 50+40+3+8
        self.assertEqual(div.read_metrics("5S50M2I40M3D8M4H"), (103, 109, 101))

    def test_a_malformed_cigar_is_rejected(self):
        with self.assertRaises(ValueError):
            div.read_metrics("10M5")

    def test_quantiles_are_linear_interpolation(self):
        values = [1, 2, 3, 4]
        self.assertEqual(div.quantile(values, 0.5), 2.5)
        self.assertEqual(div.quantile(values, 0.25), 1.75)
        self.assertEqual(div.quantile(values, 0.75), 3.25)
        self.assertEqual(div.quantile([7], 0.25), 7)


class Summary(Fixture):
    def test_hand_computed_distribution(self):
        # NM over 100 alignment columns: 0, 2 (exactly 2%: not above), 3, 6 (above 5%).
        lines = [sam("a", 1101, "100M", 0), sam("b", 2101, "100M", 2),
                 sam("c", 3101, "100M", 3), sam("d", 4101, "10S100M", 6)]
        metrics = div.summarise(self.collect(lines))
        self.assertEqual(metrics["divergence_ihe3034_reads"], "4")
        self.assertEqual(metrics["nm_ihe3034_median"], "2.500000")
        self.assertEqual(metrics["divergence_ihe3034_median"], "0.025000")
        self.assertEqual(metrics["divergence_ihe3034_q1"], "0.015000")
        self.assertEqual(metrics["divergence_ihe3034_q3"], "0.037500")
        self.assertEqual(metrics["divergence_ihe3034_frac_gt_2pct"], "0.500000")
        self.assertEqual(metrics["divergence_ihe3034_frac_gt_5pct"], "0.250000")
        # read length counts the soft clip: 100, 100, 100, 110
        self.assertEqual(metrics["clb_read_length_median"], "100")
        self.assertEqual(metrics["clb_read_length_min"], "100")
        self.assertEqual(metrics["clb_read_length_max"], "110")
        self.assertEqual(metrics["clb_read_length_q3"], "102.500000")

    def test_indels_count_in_the_denominator(self):
        # 50M 2D 50M, NM 2: 2 / 102 columns -- under 2%, whereas 2 / 100 read bases is not
        metrics = div.summarise(self.collect([sam("a", 1101, "50M2D50M", 2)]))
        self.assertEqual(metrics["divergence_ihe3034_frac_gt_2pct"], "0")
        self.assertEqual(metrics["divergence_ihe3034_median"], f"{2 / 102:.6f}")

    def test_no_reads_writes_only_the_count(self):
        self.assertEqual(div.summarise([]), {"divergence_ihe3034_reads": "0"})

    def test_per_gene_rows_cover_all_19_genes(self):
        rows = div.per_gene(self.collect([sam("a", 1101, "100M", 1), sam("b", 1201, "100M", 3)]))
        self.assertEqual([r[0] for r in rows], GENES)
        self.assertEqual(rows[0], ("clbA", 2, "2", "0.020000"))
        self.assertEqual(rows[1], ("clbB", 0, "NA", "NA"))

    def test_an_annotation_missing_a_gene_is_rejected(self):
        broken = self.dir / "broken.gff"
        broken.write_text("".join(l for l in self.gff.read_text().splitlines(True)
                                  if "Name=clbS" not in l))
        with self.assertRaisesRegex(ValueError, "clbS"):
            div.load_genes(broken)


class CommandLineIntoCohortTable(Fixture):
    """The fragment the script writes is one build_qc_summary.py accepts and tabulates."""

    def run_script(self, lines):
        sam_path = self.dir / "in.sam"
        sam_path.write_text("".join(lines))
        qc, genes = self.dir / "S1.div.qc.tsv", self.dir / "S1.div.genes.tsv"
        subprocess.run([sys.executable, str(SCRIPT), "--sample", "S1",
                        "--annotation", str(self.gff), "--sam", str(sam_path),
                        "--qc-out", str(qc), "--genes-out", str(genes)], check=True)
        return qc, genes

    def cohort_row(self, fragments):
        out = self.dir / "qc.tsv"
        metrics = bqs.load_fragments(fragments)
        bqs.write_summary(metrics, out)
        with open(out) as handle:
            return next(csv.DictReader(handle, delimiter="\t"))

    def alignment_fragment(self, mapped):
        path = self.dir / "S1.alignment.qc.tsv"
        path.write_text("Sample\tMetric\tValue\n"
                        f"S1\treads_mapped_ihe3034_mapq40\t{mapped}\n"
                        f"S1\treads_clb_genes_align\t{mapped}\n"
                        "S1\tnum_clb_genes_align\t2\n")
        return path

    def test_fragment_reaches_the_cohort_columns(self):
        qc, genes = self.run_script([sam("a", 1101, "100M", 1), sam("b", 2101, "100M", 7)])
        row = self.cohort_row([qc, self.alignment_fragment(2)])
        self.assertEqual(row["divergence_ihe3034_reads"], "2")
        self.assertEqual(row["divergence_ihe3034_median"], "0.040000")
        self.assertEqual(row["divergence_ihe3034_frac_gt_5pct"], "0.500000")
        self.assertEqual(row["clb_read_length_median"], "100")
        with open(genes) as handle:
            self.assertEqual(handle.readline().rstrip("\n").split("\t"), div.GENE_FIELDS)

    def test_zero_reads_is_zero_and_na_not_a_crash(self):
        qc, _genes = self.run_script(["@HD\tVN:1.6\n"])
        self.assertEqual(qc.read_text(), "Sample\tMetric\tValue\nS1\tdivergence_ihe3034_reads\t0\n")
        row = self.cohort_row([qc, self.alignment_fragment(0)])
        self.assertEqual(row["divergence_ihe3034_reads"], "0")
        self.assertEqual(row["divergence_ihe3034_median"], "NA")


class CohortSchema(unittest.TestCase):
    def row(self, **values):
        base = {c: None for c in bqs.OUTPUT_COLUMNS}
        base.update(values)
        return bqs.output_row("S1", base)

    def test_columns_sit_beside_the_alignment_counts_in_order(self):
        cols = bqs.OUTPUT_COLUMNS
        start = cols.index("reads_clb_genes_align") + 1
        self.assertEqual(cols[start:start + len(div.QC_METRICS)], div.QC_METRICS)
        self.assertEqual(bqs.DIVERGENCE_COLUMNS, div.QC_METRICS)

    def test_absent_is_na(self):
        out = self.row()
        for column in div.QC_METRICS:
            self.assertEqual(out[column], "NA")

    def test_more_measured_reads_than_mapped_reads_is_impossible(self):
        with self.assertRaisesRegex(ValueError, "divergence_ihe3034_reads"):
            self.row(divergence_ihe3034_reads=5, reads_mapped_ihe3034_mapq40=4)

    def test_quantiles_out_of_order_are_impossible(self):
        with self.assertRaisesRegex(ValueError, "divergence_ihe3034_q1"):
            self.row(divergence_ihe3034_reads=3, divergence_ihe3034_q1=0.03,
                     divergence_ihe3034_median=0.02)
        with self.assertRaisesRegex(ValueError, "frac_gt_5pct"):
            self.row(divergence_ihe3034_reads=3, divergence_ihe3034_frac_gt_2pct=0.1,
                     divergence_ihe3034_frac_gt_5pct=0.2)

    def test_decimals_without_measured_reads_are_impossible(self):
        with self.assertRaisesRegex(ValueError, "divergence_ihe3034_median"):
            self.row(divergence_ihe3034_reads=0, divergence_ihe3034_median=0.01)

    def test_fragment_values_are_range_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            for metric, value in (("divergence_ihe3034_frac_gt_2pct", "1.5"),
                                  ("divergence_ihe3034_median", "-0.1"),
                                  ("nm_ihe3034_median", "nan"),
                                  ("clb_read_length_median", "abc"),
                                  ("divergence_ihe3034_reads", "1.5")):
                path = Path(tmp) / "f.tsv"
                path.write_text(f"Sample\tMetric\tValue\nS1\t{metric}\t{value}\n")
                with self.subTest(metric=metric), self.assertRaises(ValueError):
                    bqs.load_fragments([path])

    def test_names_say_divergence_not_specificity(self):
        for column in bqs.DIVERGENCE_COLUMNS:
            self.assertNotRegex(column, r"specific|unique|carrier|identity")


class Wiring(unittest.TestCase):
    ALIGN = (ROOT / "Modules/pksProfiler_align.nf").read_text()
    MAIN = (ROOT / "main.nf").read_text()

    def test_the_alignment_task_runs_it_on_the_counts_bam(self):
        self.assertIn('ihe3034_divergence_qc.py', self.ALIGN)
        self.assertRegex(self.ALIGN, r'samtools view "\$\{bam\}" \\\s*\n\s*\| python')
        self.assertIn("--min-mapq ${params.pks_align_min_mapq}", self.ALIGN)
        self.assertIn("--annotation clb_genes.gff", self.ALIGN)
        self.assertIn("emit: divergence", self.ALIGN)

    def test_it_feeds_the_qc_table_and_nothing_else(self):
        uses = re.findall(r"ALIGN_OUT\.divergence[^\n]*", self.MAIN)
        self.assertEqual(len(uses), 1, uses)
        block = self.MAIN[self.MAIN.index("ALIGN_OUT.divergence") - 200:
                          self.MAIN.index("ALIGN_OUT.divergence")]
        self.assertIn("QC_FRAGMENTS = QC_FRAGMENTS.mix(", block)

    def test_it_runs_after_the_counts_are_fixed(self):
        # The clb counts come from featureCounts' table; the divergence step reads the
        # BAM afterwards and writes only its own two files.
        self.assertLess(self.ALIGN.index("CLB_READS="), self.ALIGN.index("ihe3034_divergence_qc.py"))
        self.assertNotRegex(self.ALIGN, r'divergence_qc[^\n]*>>\s*"\$\{qc\}"')


if __name__ == "__main__":
    unittest.main()
