"""One authoritative definition of the pks island's coordinates.

Ludmil, revised report finding 6: the annotation spans 2,193,827-2,244,594
inclusive (50,768 bp). pks_shift=2193827 and pks_island_len=50767 got a
different +/-1 adjustment in every module that turned them into a region:
correct in pks_taxa.nf (pks_shift-1 for a 0-based BED start), off by one bp at
the start everywhere else -- pksProfiler_align.nf, pks_targeted.nf's two call
sites, and this session's own pks_mag.nf locus alignment -- silently dropping
the true first base of the island from a samtools region or a PAF/BED-style
0-based start.

pks_start_1based/pks_end_1based are now the single source of truth; every
correctness-sensitive call site derives its own coordinate representation from
them directly. The legacy pks_shift/pks_island_len params, last read by
plotting.nf's coverage-plot window before it moved to the canonical pair, were
removed in 1.0.0-rc2.
"""
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()
ALIGN = (ROOT / "Modules/pksProfiler_align.nf").read_text()
TARGETED = (ROOT / "Modules/pks_targeted.nf").read_text()
TAXA = (ROOT / "Modules/pks_taxa.nf").read_text()
MAG = (ROOT / "Modules/pks_mag.nf").read_text()


def code_only(text):
    return "\n".join(line for line in text.splitlines()
                      if not line.strip().startswith(("#", "//")))


class TheCanonicalValuesAreDeclaredOnce(unittest.TestCase):
    def test_start_and_end_match_the_annotation(self):
        self.assertIn("params.pks_start_1based    = 2193827", MAIN)
        self.assertIn("params.pks_end_1based      = 2244594", MAIN)

    def test_the_corrected_length_is_derived_not_hardcoded(self):
        self.assertIn(
            "params.pks_island_len_1based = params.pks_end_1based - params.pks_start_1based + 1",
            MAIN)

    def test_the_corrected_length_is_50768_not_50767(self):
        # The exact bug: pks_island_len (50767) is one short of the true,
        # inclusive-on-both-ends length.
        start, end = 2193827, 2244594
        self.assertEqual(end - start + 1, 50768)
        self.assertNotEqual(end - start + 1, 50767)

    def test_the_legacy_params_are_gone(self):
        # 1.0.0-rc2: nothing read pks_shift/pks_island_len (50,767) any more, so the
        # definitions were removed rather than left to be picked up again by mistake.
        self.assertNotIn("params.pks_shift =", MAIN)
        self.assertNotIn("params.pks_island_len =", MAIN)
        for path in sorted((ROOT / "Modules").glob("*.nf")) + sorted((ROOT / "lib").glob("*.groovy")):
            with self.subTest(file=path.name):
                text = path.read_text()
                self.assertNotRegex(text, r"params\.pks_shift\b")
                self.assertNotRegex(text, r"params\.pks_island_len\b(?!_1based)")


class EveryCorrectnessSensitiveCallSiteUsesTheCanonicalValues(unittest.TestCase):
    def test_pksProfiler_align_samtools_region_needs_no_arithmetic(self):
        code = code_only(ALIGN)
        self.assertIn('pks_start_1based}-${params.pks_end_1based', code)
        self.assertNotIn("pks_shift.toString().toInteger() + 1", code)

    def test_classifyPksReadEvidence_samtools_region_and_island_length_agree(self):
        code = code_only(TARGETED)
        self.assertIn('pks_start_1based}-${params.pks_end_1based', code)
        self.assertIn("--island-length \"${params.pks_island_len_1based}\"", code)
        # The old mismatched pair must both be gone, not just the region.
        self.assertNotIn("--island-length \"${params.pks_island_len}\"", code)

    def test_summarizeTargetedPksEvidence_converts_to_zero_based_half_open(self):
        code = code_only(TARGETED)
        self.assertIn('--island-start "${params.pks_start_1based.toString().toInteger() - 1}"', code)
        self.assertIn('--island-end "${params.pks_end_1based}"', code)

    def test_pks_taxa_bed_start_still_subtracts_one_from_the_canonical_value(self):
        code = code_only(TAXA)
        self.assertIn("pks_start_1based} - 1", code)
        self.assertIn("END=${params.pks_end_1based}", code)

    def test_mag_bin_locus_evidence_converts_to_zero_based_half_open(self):
        code = code_only(MAG)
        self.assertIn("pks_start_1based} - 1", code)
        self.assertIn("--island-end ${params.pks_end_1based}", code)

    def test_no_fixed_site_still_does_the_old_pks_shift_arithmetic(self):
        # The specific pattern each of the four sites his report named used:
        # pks_shift(+1) combined with pks_island_len to build a region/length.
        # pks_taxa.nf is exempt -- it already used pks_shift-1 correctly and is
        # migrated to the same canonical params, not to an arithmetic-free form.
        for name, code in (("pksProfiler_align.nf", code_only(ALIGN)),
                           ("pks_targeted.nf", code_only(TARGETED)),
                           ("pks_mag.nf", code_only(MAG)),
                           ("plotting.nf", code_only((ROOT / "Modules" / "plotting.nf").read_text()))):
            with self.subTest(module=name):
                self.assertNotIn("params.pks_shift.toString().toInteger() + params.pks_island_len",
                                 code, f"{name} still derives a bound from pks_shift+pks_island_len")


if __name__ == "__main__":
    unittest.main()


class GeneOverlapRuleReachesEveryLane(unittest.TestCase):
    """params.min_gene_overlap_bp must reach the primary counts, not just the side lanes.

    Ludmil's v3 review, item 12. Modules/pksProfiler_align.nf passed --largestOverlap to
    featureCounts and never --minOverlap, so the primary clb counts used featureCounts'
    own 1 bp default while Modules/pks_taxa.nf and Modules/pks_targeted.nf both honoured
    params.min_gene_overlap_bp. At the shipped default of 1 the two agree by coincidence,
    which is why this went unnoticed; at any other value the lanes disagree and the
    primary counts answer to a rule the run does not report.

    C2 confirmed every production run used 1, so no existing result is affected.
    """

    ALIGN = (ROOT / "Modules" / "pksProfiler_align.nf").read_text()

    def test_the_primary_count_passes_the_configured_overlap(self):
        self.assertIn("--minOverlap", self.ALIGN,
                      "featureCounts must be told the overlap rule, not left on its default")
        self.assertIn("--minOverlap ${params.min_gene_overlap_bp}", self.ALIGN,
                      "and it must be the pipeline's parameter, not a literal")

    def test_it_is_not_hardcoded_to_the_default(self):
        # A literal 1 would pass today and silently ignore the parameter tomorrow.
        self.assertNotIn("--minOverlap 1", self.ALIGN)

    def test_every_lane_that_assigns_reads_to_genes_uses_the_same_parameter(self):
        lanes = {
            "Modules/pksProfiler_align.nf": "--minOverlap ${params.min_gene_overlap_bp}",
            "Modules/pks_taxa.nf": "--min-overlap",
            "Modules/pks_targeted.nf": "--min-overlap",
        }
        for path, needle in lanes.items():
            with self.subTest(lane=path):
                text = (ROOT / path).read_text()
                self.assertIn(needle, text)
                self.assertIn("params.min_gene_overlap_bp", text,
                              f"{path} assigns reads to clb genes, so it has to use the "
                              f"one configured overlap rule")

    def test_a_non_default_value_would_actually_change_the_command(self):
        """The regression the review asked for: substitute a non-default and check it lands."""
        rendered = self.ALIGN.replace("${params.min_gene_overlap_bp}", "25")
        self.assertIn("--minOverlap 25", rendered)
        self.assertNotIn("--minOverlap ${params.min_gene_overlap_bp}", rendered)
