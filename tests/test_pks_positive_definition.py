"""One definition of a pks-positive bin, counted the same way everywhere.

M3, from the internal v0.0.2 audit. `grep -vc '^#'` on the hmmsearch tblout counted
protein-to-model *hit lines*: one gene matched by three predicted proteins counted
three times, no threshold beyond -E, and a single line anywhere made a bin positive.
That count decided the status file's pks_positive_bin_count and which bins got
genomic-context extraction, while build_mag_summary.py reported distinct clb genes
and community producer selection used a third rule with its own threshold.

The count is now distinct clb genes, produced where hmmsearch runs, and one
parameter -- params.mag_min_clb_genes -- is applied by every consumer.
"""
import importlib.util
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAG = (ROOT / "Modules/pks_mag.nf").read_text()
MAIN = (ROOT / "main.nf").read_text()


def code_only(text):
    """Comments explain what was removed; the assertions are about what runs."""
    return "\n".join(line for line in text.splitlines()
                      if not line.strip().startswith(("#", "//")))


MAG_CODE = code_only(MAG)
MAIN_CODE = code_only(MAIN)

spec = importlib.util.spec_from_file_location("mag_utils", ROOT / "scripts/mag_utils.py")
mag_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mag_utils)

# The awk exactly as the module ships it, with Nextflow's escaping resolved.
AWK = r"""!/^#/ && NF>=19 && ($5+0) <= 1e-5 { genes[$3]=1 }
     END { print length(genes)+0 }"""


def tblout(rows):
    """rows: (locus_tag, clb_gene, evalue). 19+ fields, as hmmsearch writes."""
    lines = ["# target name        accession  query name"]
    for locus, gene, evalue in rows:
        lines.append(" ".join([locus, "-", gene, "-", evalue, "50.0", "0.0",
                               evalue, "49.0", "0.0"] + ["1.0"] * 9 + ["description"]))
    return "\n".join(lines) + "\n"


def count_with_awk(text):
    with tempfile.NamedTemporaryFile("w", suffix=".tblout", delete=False) as handle:
        handle.write(text)
        path = handle.name
    try:
        result = subprocess.run(["awk", AWK, path], capture_output=True, text=True, check=True)
        return int(result.stdout.strip())
    finally:
        Path(path).unlink()


class TheCountIsGenesNotLines(unittest.TestCase):
    def test_one_gene_matched_by_several_proteins_counts_once(self):
        text = tblout([("p1", "clbB", "1e-30"), ("p2", "clbB", "1e-25"),
                       ("p3", "clbB", "1e-20")])
        self.assertEqual(text.count("clbB"), 3)      # three hit lines
        self.assertEqual(count_with_awk(text), 1)    # one gene

    def test_hits_above_the_evalue_cut_do_not_count(self):
        text = tblout([("p1", "clbB", "1e-30"), ("p2", "clbA", "1e-3")])
        self.assertEqual(count_with_awk(text), 1)

    def test_the_megasynthase_only_pattern_is_two_genes(self):
        # M1's false positive: clbB and clbK hit every bin by domain homology.
        text = tblout([("p1", "clbB", "1e-30"), ("p2", "clbK", "1e-28"),
                       ("p3", "clbB", "1e-22")])
        self.assertEqual(count_with_awk(text), 2)

    def test_an_empty_tblout_is_zero_not_an_error(self):
        self.assertEqual(count_with_awk("# nothing matched\n"), 0)

    def test_it_matches_what_the_summary_reports(self):
        """The invariant: the same quantity, by construction, not by coincidence."""
        rows = [("p1", "clbB", "1e-30"), ("p2", "clbB", "1e-25"),
                ("p3", "clbK", "1e-28"), ("p4", "clbA", "1e-3")]
        text = tblout(rows)
        with tempfile.NamedTemporaryFile("w", suffix=".tblout", delete=False) as handle:
            handle.write(text)
            path = handle.name
        try:
            hits = mag_utils.parse_hmmsearch_tblout(path)
            distinct = len({hit["clb_gene"] for hit in hits})
        finally:
            Path(path).unlink()
        self.assertEqual(count_with_awk(text), distinct)


class OneThresholdEverywhere(unittest.TestCase):
    def test_the_hit_line_count_is_gone(self):
        self.assertEqual([l for l in MAG_CODE.splitlines() if "grep -vc" in l], [])
        self.assertEqual([l for l in MAG_CODE.splitlines() if "hit_count" in l], [])

    def test_every_consumer_reads_the_same_parameter(self):
        # which bins get genomic context (this HMM pre-filter), and community producer
        # selection's HMM leg. mag_status's positivity count no longer reads this
        # parameter at all -- see LocusEvidenceDecidesPositivity below -- so the
        # invariant this test protects is now "the two remaining HMM-count consumers
        # agree", not "all three positivity consumers agree".
        self.assertEqual(MAG_CODE.count("params.mag_min_clb_genes"), 2)

    def test_the_parameter_is_declared_once_and_the_old_one_is_gone(self):
        self.assertIn("params.mag_min_clb_genes = 3", MAIN_CODE)
        self.assertEqual([l for l in MAIN_CODE.splitlines() if "community_min_clb_genes" in l], [])
        self.assertEqual([l for l in MAG_CODE.splitlines() if "community_min_clb_genes" in l], [])

    def test_the_emitted_count_is_named_for_what_it_holds(self):
        self.assertIn("emit: clb_gene_count", MAG)


class LocusEvidenceDecidesPositivity(unittest.TestCase):
    """M1: the HMM gene count alone still admits ERR525841's Bifidobacterium bins --
    clbB/clbH/clbK domain homology clears mag_min_clb_genes with no island-specific
    evidence. Positivity is decided by aligning each bin's own assembly to the
    canonical locus (alignMagBinToCanonicalReference/magBinLocusEvidence) and scoring
    it the same way read-level evidence already is, not by picking specific genes out
    of the HMM hits.
    """

    def test_the_locus_tier_params_are_declared_with_the_read_level_values(self):
        # Same numbers as classify_tumor_pks_evidence.py's tier thresholds (minus the
        # read-count criterion, which has no bin-level analogue) -- one evidentiary
        # standard whether the evidence is reads or an assembled bin.
        for line in (
            "params.mag_locus_multi_gene_min_genes         = 3",
            "params.mag_locus_multi_gene_min_breadth       = .01",
            "params.mag_locus_broad_island_min_genes       = 8",
            "params.mag_locus_broad_island_min_breadth     = .075",
            "params.mag_locus_extensive_island_min_genes   = 10",
            "params.mag_locus_extensive_island_min_breadth = .15",
        ):
            self.assertIn(line, MAIN_CODE)

    def test_pks_reference_fasta_feeds_the_bin_alignment(self):
        # U1: bins_flat_ch is mixed with the per-sample unbinned pool into
        # locus_alignment_units_ch before this call, so the same alignment covers
        # both -- Nextflow forbids invoking a process twice in one workflow, and
        # every bin is still in there (bins_flat_ch.mix(unbinned_flat_ch)).
        self.assertIn("alignMagBinToCanonicalReference(locus_alignment_units_ch, pks_reference_ch)", MAG_CODE)
        self.assertIn("locus_alignment_units_ch = bins_flat_ch.mix(unbinned_flat_ch)", MAG_CODE)

    def test_every_bin_is_aligned_not_only_hmm_candidates(self):
        # bins_flat_ch, not pks_pos_tblout_ch or any HMM-filtered channel.
        self.assertIn("Channel.value(file(params.pks_reference_fasta", MAG_CODE)

    def test_mag_status_count_reads_the_locus_tier_not_the_hmm_count(self):
        # U1: bin_locus_evidence_ch is magBinLocusEvidence.out.evidence with the
        # unbinned pseudo-unit's row filtered back out -- see
        # test_unbinned_pool_never_inflates_bin_positivity_counts below for why.
        self.assertIn("bin_locus_evidence_ch = magBinLocusEvidence.out.evidence", MAG_CODE)
        self.assertIn("pks_positive_bin_counts_ch = bin_locus_evidence_ch", MAG_CODE)
        self.assertNotIn("pks_positive_bin_counts_ch = hmmsearchClb.out.clb_gene_count", MAG_CODE)

    def test_positive_tiers_exclude_negative_and_indeterminate(self):
        self.assertIn(
            '(tier in ["multi_gene", "broad_island", "extensive_island"]) ? 1 : 0', MAG_CODE,
        )

    def test_community_and_summary_both_receive_locus_evidence(self):
        # communityProphageSummary reads locus_evidence_per_sample_ch (real bins only);
        # magSummaryTable reads mag_summary_locus_evidence_ch (bins + the U1 unbinned
        # pool, so build_mag_summary.py can report the unbinned call too). Two
        # different channels now, not the same one joined twice, because the
        # unbinned pseudo-unit has no taxonomy/annotation for communityProphageSummary
        # to read but does belong in the per-sample MAG summary table.
        self.assertEqual(MAG_CODE.count(".join(locus_evidence_per_sample_ch, by: 0, remainder: true)"), 1)
        self.assertEqual(MAG_CODE.count(".join(mag_summary_locus_evidence_ch, by: 0, remainder: true)"), 1)

    def test_unbinned_pool_never_inflates_bin_positivity_counts(self):
        # U1: the per-sample pool of contigs MetaBAT2 never binned is not a genome --
        # a positive call there must never count toward pks_positive_bin_count
        # (mag_status.tsv) or communityProphageSummary's producer selection.
        self.assertIn(
            'bin_locus_evidence_ch = magBinLocusEvidence.out.evidence\n'
            '        .filter { _sampleID, unitID, _evidence -> unitID != "unbinned" }',
            MAG,
        )

    def test_the_specific_gene_mechanism_is_gone(self):
        self.assertEqual([l for l in MAG_CODE.splitlines() if "has_specific_clb" in l], [])
        self.assertEqual([l for l in MAIN_CODE.splitlines() if "mag_specific_clb_genes" in l], [])


if __name__ == "__main__":
    unittest.main()
