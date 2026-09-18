"""Coverage for the frozen three-tier tumour pks evidence classifier.

Tiers come from tier_threshold_development_20260909: all three criteria must hold, and
classification proceeds from the highest tier downward.
"""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

spec = importlib.util.spec_from_file_location("classify_tumor", ROOT / "scripts/classify_tumor_pks_evidence.py")
classify_tumor = importlib.util.module_from_spec(spec)
spec.loader.exec_module(classify_tumor)

FROZEN = {
    "multi_gene":       (5, 3, 0.01),
    "broad_island":     (30, 8, 0.075),
    "extensive_island": (100, 10, 0.15),
}

def tier(reads, genes, breadth):
    return classify_tumor.classify(reads, genes, breadth, FROZEN)


class TierCascadeTests(unittest.TestCase):
    def test_zero_reads_is_negative(self):
        self.assertEqual(tier(0, 0, 0.0), "negative")
        # Genes/breadth cannot rescue a sample with no assigned reads.
        self.assertEqual(tier(0, 19, 1.0), "negative")

    def test_below_multi_gene_is_localized_indeterminate(self):
        self.assertEqual(tier(1, 1, 0.0001), "localized_indeterminate")
        self.assertEqual(tier(4, 3, 0.05), "localized_indeterminate")   # reads short
        self.assertEqual(tier(50, 2, 0.05), "localized_indeterminate")  # genes short
        self.assertEqual(tier(50, 5, 0.005), "localized_indeterminate") # breadth short

    def test_multi_gene_boundary_is_five_reads_not_eleven(self):
        # 11 reads was candidate multi_gene_B, evaluated and rejected; the frozen
        # boundary is 5. Guards against reintroducing the pipeline's old default.
        self.assertEqual(tier(5, 3, 0.01), "multi_gene")
        self.assertEqual(tier(4, 3, 0.01), "localized_indeterminate")

    def test_each_tier_at_its_exact_boundary(self):
        self.assertEqual(tier(5, 3, 0.01), "multi_gene")
        self.assertEqual(tier(30, 8, 0.075), "broad_island")
        self.assertEqual(tier(100, 10, 0.15), "extensive_island")

    def test_all_three_criteria_are_required(self):
        # Meets broad_island reads and genes but not breadth -> falls back to multi_gene.
        self.assertEqual(tier(30, 8, 0.05), "multi_gene")
        # Meets extensive reads and breadth but not genes -> falls back to broad_island.
        self.assertEqual(tier(100, 8, 0.15), "broad_island")

    def test_classification_takes_the_highest_qualifying_tier(self):
        self.assertEqual(tier(1000, 19, 0.99), "extensive_island")

    def test_tier_order_is_highest_first(self):
        # The cascade returns on the first match, so TIERS must run highest to lowest.
        self.assertEqual(classify_tumor.TIERS, ("extensive_island", "broad_island", "multi_gene"))


class ContigGateTests(unittest.TestCase):
    """The default gate assembles only broad_island and extensive_island."""

    DEFAULT_GATE = ("broad_island", "extensive_island")

    def test_default_gate_excludes_multi_gene(self):
        main = (ROOT / "main.nf").read_text()
        self.assertIn('params.tumor_contig_tiers = "broad_island,extensive_island"', main)
        self.assertNotIn("multi_gene,broad_island,extensive_island", main.split("params.tumor_contig_tiers")[1].split("\n")[0])

    def test_multi_gene_classifies_but_does_not_assemble(self):
        self.assertEqual(tier(5, 3, 0.01), "multi_gene")
        self.assertNotIn(tier(5, 3, 0.01), self.DEFAULT_GATE)

    def test_higher_tiers_do_assemble(self):
        self.assertIn(tier(30, 8, 0.075), self.DEFAULT_GATE)
        self.assertIn(tier(100, 10, 0.15), self.DEFAULT_GATE)

    def test_indeterminate_and_negative_never_assemble(self):
        self.assertNotIn(tier(0, 0, 0.0), self.DEFAULT_GATE)
        self.assertNotIn(tier(4, 3, 0.01), self.DEFAULT_GATE)


class BooleanFlagTests(unittest.TestCase):
    """`--flag false` reaches Nextflow as the String "false", which is truthy in Groovy.

    Every v0.0.2 boolean must therefore be read through the asBool coercion, never
    tested as params.<flag> directly, or `--enable_mags false` silently enables MAGs.
    """

    MAIN = (ROOT / "main.nf").read_text()
    BODY = MAIN.split("workflow {", 1)[1]
    # Everything after the coercion block. The block itself reads params.<flag> once
    # inside asBool(), which is the only legitimate direct read.
    AFTER = BODY.split("Required input validation", 1)[1]
    FLAGS = ("enable_mags", "tumor_enable_mags", "tumor_full_contig_context",
             "tumor_targeted_assembly", "diamond_rescue", "pks_community_taxa",
             "enable_strain_typing")

    def test_coercion_uses_string_to_boolean(self):
        # Nextflow 26 strict DSL rejects a local closure called as a function,
        # so the coercion is inlined per flag rather than shared via a helper.
        self.assertIn("toString().toBoolean()", self.BODY)

    def test_every_boolean_flag_is_coerced(self):
        for flag in self.FLAGS:
            self.assertIn(f"def {flag}_b", self.BODY, f"{flag} has no coerced local")

    def test_no_flag_is_tested_directly_after_coercion(self):
        for flag in self.FLAGS:
            self.assertNotIn(f"params.{flag}", self.AFTER,
                             f"params.{flag} tested directly; use {flag}_b")

    def test_coercion_precedes_the_mag_branch(self):
        # The MAG branch is the earliest consumer, and sits before the validation block.
        coercion = self.BODY.index("def enable_mags_b")
        mag_branch = self.BODY.index("MAG_ASSEMBLY_READS")
        self.assertLess(coercion, mag_branch)


class ClassifierOutputTests(unittest.TestCase):
    def _write_inputs(self, directory, reads, gene_count, covered_positions, island_length):
        qc = directory / "qc.tsv"
        with qc.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, ["Metric", "Value"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"Metric": "reads_clb_genes_align", "Value": str(reads)})
        counts = directory / "counts.tsv"
        genes = [f"clb{x}" for x in "ABCDEFGHIJKLMNOPQRS"]
        with counts.open("w") as handle:
            handle.write("Geneid\tCount\n")
            for index, gene in enumerate(genes):
                handle.write(f"{gene}\t{1 if index < gene_count else 0}\n")
        depth = directory / "depth.tsv"
        with depth.open("w") as handle:
            for position in range(island_length):
                handle.write(f"NC_017628.1\t{position + 1}\t{1 if position < covered_positions else 0}\n")
        return qc, counts, depth

    def test_end_to_end_emits_expected_tier_and_breadth(self):
        island_length = 1000
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            qc, counts, depth = self._write_inputs(directory, reads=120, gene_count=12,
                                                   covered_positions=250, island_length=island_length)
            output = directory / "evidence.tsv"
            import sys
            argv = sys.argv
            sys.argv = ["classify", "--sample", "S1", "--qc", str(qc), "--counts", str(counts),
                        "--depth", str(depth), "--island-length", str(island_length),
                        "--output", str(output)]
            try:
                classify_tumor.main()
            finally:
                sys.argv = argv
            with output.open() as handle:
                row = list(csv.DictReader(handle, delimiter="\t"))[0]
            self.assertEqual(row["read_evidence"], "extensive_island")
            self.assertEqual(row["pks_reads"], "120")
            self.assertEqual(row["clb_genes_detected"], "12")
            self.assertEqual(row["island_breadth_1x"], "0.250000")

    def test_non_monotonic_thresholds_are_rejected_by_main(self):
        island_length = 100
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            qc, counts, depth = self._write_inputs(directory, reads=50, gene_count=8,
                                                   covered_positions=50, island_length=island_length)
            import sys
            argv = sys.argv
            sys.argv = ["classify", "--sample", "S1", "--qc", str(qc), "--counts", str(counts),
                        "--depth", str(depth), "--island-length", str(island_length),
                        "--output", str(directory / "out.tsv"),
                        "--broad-island-reads", "1"]
            try:
                with self.assertRaises(SystemExit):
                    classify_tumor.main()
            finally:
                sys.argv = argv

    def test_depth_row_count_mismatch_is_fatal(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            qc, counts, depth = self._write_inputs(directory, reads=10, gene_count=5,
                                                   covered_positions=5, island_length=100)
            with self.assertRaises(ValueError):
                classify_tumor.breadths(depth, 999)


if __name__ == "__main__":
    unittest.main()
