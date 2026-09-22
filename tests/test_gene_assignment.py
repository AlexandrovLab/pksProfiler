"""One definition of which clb gene a read belongs to, used everywhere.

F14, from Ludmil's 2026-09-19 report. Three places decided this differently:
featureCounts with --largestOverlap; the alignment QC with `bedtools intersect -u`,
which counted every read *touching* any clb gene rather than assigning it to one; and
the taxonomy lane with its own awk over a BED built from every named feature, not
only clbA-clbS. The three numbers sat in the same tables and were not comparable.
"""
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BED = ROOT / "scripts/clb_gene_bed.py"
ASSIGN = ROOT / "scripts/assign_reads_to_genes.py"
CLB_GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]

spec = importlib.util.spec_from_file_location("assign_reads_to_genes", ASSIGN)
assign_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(assign_module)


def gff(genes, extra=()):
    lines = ["##gff-version 3"]
    start = 100
    for gene in genes:
        lines.append(f"ctg\tsrc\tgene\t{start}\t{start + 99}\t.\t+\t.\tID=g{gene};Name={gene}")
        start += 200
    for name in extra:
        lines.append(f"ctg\tsrc\tgene\t{start}\t{start + 99}\t.\t+\t.\tID=g{name};Name={name}")
        start += 200
    return "\n".join(lines) + "\n"


class CanonicalIntervals(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.d = Path(self.tmp.name)

    def build(self, text, *flags):
        (self.d / "clb.gff").write_text(text)
        result = subprocess.run(
            [sys.executable, str(BED), "--annotation", str(self.d / "clb.gff"),
             "--output", str(self.d / "out.bed"), *flags],
            capture_output=True, text=True)
        return result, (self.d / "out.bed")

    def test_only_clb_genes_are_included(self):
        # The taxonomy lane used to take every named feature, so a read over a
        # neighbouring gene was assigned to it there and ignored in the QC.
        result, bed = self.build(gff(CLB_GENES, extra=["neighbourX", "tRNA-Asn"]))
        self.assertEqual(result.returncode, 0, result.stderr)
        names = [line.split("\t")[3] for line in bed.read_text().splitlines()]
        self.assertEqual(names, CLB_GENES)

    def test_gff_coordinates_become_half_open_bed(self):
        result, bed = self.build(gff(CLB_GENES))
        self.assertEqual(result.returncode, 0, result.stderr)
        contig, start, end, name = bed.read_text().splitlines()[0].split("\t")[:4]
        self.assertEqual((start, end, name), ("99", "199", "clbA"))

    def test_an_annotation_missing_genes_is_refused(self):
        result, _ = self.build(gff(CLB_GENES[:5]))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing 14 clb genes", result.stdout + result.stderr)

    def test_a_duplicated_gene_is_refused(self):
        result, _ = self.build(gff(CLB_GENES) + gff(["clbA"]).splitlines()[1] + "\n")
        self.assertNotEqual(result.returncode, 0)


class AssignmentRule(unittest.TestCase):
    def test_a_read_goes_to_the_gene_it_shares_most_bases_with(self):
        best = assign_module.assign([("r1", "clbA", 30), ("r1", "clbB", 70)])
        self.assertEqual(best["r1"], ("clbB", 70))

    def test_overlaps_from_several_blocks_of_one_read_are_summed(self):
        best = assign_module.assign([("r1", "clbA", 40), ("r1", "clbA", 40), ("r1", "clbB", 70)])
        self.assertEqual(best["r1"], ("clbA", 80))

    def test_ties_break_deterministically_not_by_input_order(self):
        forward = assign_module.assign([("r1", "clbQ", 50), ("r1", "clbB", 50)])
        reverse = assign_module.assign([("r1", "clbB", 50), ("r1", "clbQ", 50)])
        self.assertEqual(forward, reverse)
        self.assertEqual(forward["r1"][0], "clbB")

    def test_an_overlap_below_the_minimum_does_not_count(self):
        self.assertEqual(assign_module.assign([("r1", "clbA", 5)], min_overlap=10), {})
        self.assertIn("r1", assign_module.assign([("r1", "clbA", 15)], min_overlap=10))

    def test_each_read_is_assigned_at_most_once_so_counts_sum_to_reads(self):
        best = assign_module.assign([("r1", "clbA", 10), ("r1", "clbB", 20),
                                     ("r2", "clbC", 30)])
        self.assertEqual(len(best), 2)


class OneRulePerLane(unittest.TestCase):
    """featureCounts is the rule. Nothing recomputes it; the taxonomy lane, which has
    no featureCounts, implements the same documented rule in the shared script."""

    ALIGN = (ROOT / "Modules/pksProfiler_align.nf").read_text()
    TAXA = (ROOT / "Modules/pks_taxa.nf").read_text()

    def test_the_alignment_qc_count_is_the_featurecounts_total(self):
        # Not a second opinion computed with bedtools: the sum of the matrix column,
        # so the QC number and the published counts cannot disagree.
        self.assertNotIn("assign_reads_to_genes.py", self.ALIGN)
        self.assertNotIn("clb_genes.qc.bed", self.ALIGN)
        self.assertIn("--largestOverlap", self.ALIGN)
        body = self.ALIGN.split("CLB_READS=", 1)[1].split(")", 1)[0]
        self.assertIn("counts", body)

    def test_the_alignment_lane_no_longer_counts_bare_overlaps(self):
        # The comment explaining what was removed mentions it; the script must not run it.
        code = [line for line in self.ALIGN.splitlines()
                if "bedtools intersect" in line and not line.strip().startswith("#")]
        self.assertEqual(code, [])

    def test_both_read_counts_come_from_the_same_file_in_the_alignment_lane(self):
        # num_clb_genes_align was always derived from the counts table; the read count
        # now is too.
        for name in ("CLB_READS", "CLB_GENES_DETECTED"):
            section = self.ALIGN.split(f"{name}=", 1)[1][:400]
            self.assertIn("${counts}", section, f"{name} does not read the counts table")

    def test_the_taxonomy_lane_uses_the_shared_definition(self):
        self.assertIn("clb_gene_bed.py", self.TAXA)
        self.assertIn("assign_reads_to_genes.py", self.TAXA)
        self.assertNotIn('sub(/^Name=/, "", attributes[i])', self.TAXA)

    def test_the_overlap_threshold_is_a_declared_parameter(self):
        self.assertIn("params.min_gene_overlap_bp", self.TAXA)
        self.assertIn("params.min_gene_overlap_bp = 1", (ROOT / "main.nf").read_text())


if __name__ == "__main__":
    unittest.main()
