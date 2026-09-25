"""Coverage counts every confident match; the length floor only counts contigs.

T1, from the internal v0.0.2 audit. One `--min-aligned-bp` was applied to every
alignment before anything was computed, so it decided both how much of the island
was covered and how many contigs supported the call. metaSPAdes writes shorter
contigs than MEGAHIT, so the floor removed 36-77% of its alignments against 10-38%
of MEGAHIT's, and `assembler_agreement` -- which is computed from coverage alone --
became partly an artefact of the cutoff. AF-6655 and DC-5337 were reported
`no_contig_support` with every alignment discarded.

Two further details the fix depends on: the floor was applied per alignment block,
so a contig matching in three 150 bp pieces was dropped despite carrying 450 bp of
island; and clbR is 213 bp, so a 300 bp floor can reject a contig covering a whole
gene. The floor is now 200, on each contig's total.
"""
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "summarize_tumor_pks_contigs", ROOT / "scripts/summarize_tumor_pks_contigs.py")
summarize = importlib.util.module_from_spec(spec)
spec.loader.exec_module(summarize)

ISLAND_START, ISLAND_END = 2193826, 2244593


def hit(contig, start, length, identity=0.99, mapq=60):
    return dict(contig=contig, qlen=length, strand="+", start=start, end=start + length,
                aligned=length, identity=identity, mapq=mapq)


class CoverageHasNoLengthFloor(unittest.TestCase):
    def test_short_exact_matches_still_cover_the_island(self):
        # Nine 250 bp matches: none reaches 300, all are real island sequence.
        hits = [hit(f"NODE_{i}", ISLAND_START + i * 1000, 250) for i in range(9)]
        self.assertEqual(summarize.union_length(hits, ISLAND_START, ISLAND_END), 9 * 250)

    def test_overlapping_matches_are_not_double_counted(self):
        hits = [hit("a", ISLAND_START, 500), hit("b", ISLAND_START + 200, 500)]
        self.assertEqual(summarize.union_length(hits, ISLAND_START, ISLAND_END), 700)

    def test_paf_defaults_to_no_minimum(self):
        import inspect
        self.assertEqual(inspect.signature(summarize.paf).parameters["min_aligned"].default, 0)


class SupportingContigsUseTotals(unittest.TestCase):
    def test_a_contig_is_judged_on_its_total_not_its_largest_block(self):
        # three 150 bp blocks of one contig: 450 bp of island between them
        hits = [hit("NODE_1", ISLAND_START + n * 400, 150) for n in range(3)]
        self.assertEqual(summarize.supporting_contigs(hits, ISLAND_START, ISLAND_END, 200),
                         ["NODE_1"])

    def test_a_contig_below_the_floor_in_total_does_not_count(self):
        hits = [hit("NODE_1", ISLAND_START, 90), hit("NODE_1", ISLAND_START + 500, 90)]
        self.assertEqual(summarize.supporting_contigs(hits, ISLAND_START, ISLAND_END, 200), [])

    def test_a_contig_covering_the_shortest_clb_gene_counts(self):
        # clbR is 213 bp; a 300 bp floor would have rejected this
        self.assertEqual(
            summarize.supporting_contigs([hit("NODE_1", ISLAND_START, 213)],
                                         ISLAND_START, ISLAND_END, 200),
            ["NODE_1"])

    def test_contigs_are_counted_once_however_many_blocks_they_have(self):
        hits = [hit("NODE_1", ISLAND_START + n * 400, 300) for n in range(4)]
        self.assertEqual(summarize.supporting_contigs(hits, ISLAND_START, ISLAND_END, 200),
                         ["NODE_1"])


class SupportingContigsRequireIslandOverlapNotFlankAlignment(unittest.TestCase):
    """Ludmil, revised report finding 3, his three fixtures verbatim.

    supporting_contigs() used to sum each hit's full aligned length within the
    wider +/-10kb recruitment window paf() searches, not the portion actually
    inside the canonical island -- so a contig aligning mostly to a flank could
    still clear the floor and be reported as island-supporting.
    """

    # A point comfortably inside the 5kb flank paf() searches around the island,
    # i.e. before ISLAND_START -- not part of the island itself.
    FLANK_START = ISLAND_START - 3000

    def test_fixture_a_wholly_in_flank_yields_zero_supporting_contigs(self):
        # 500 bp alignment entirely in the 5kb flank: no island overlap at all.
        hits = [hit("NODE_1", self.FLANK_START, 500)]
        self.assertEqual(summarize.supporting_contigs(hits, ISLAND_START, ISLAND_END, 200), [])

    def test_fixture_b_small_island_overlap_plus_flank_fails_the_floor(self):
        # One 500 bp block straddling the boundary: 450 bp in the flank, 50 bp
        # inside the island. Old code counted the full 500 bp; only 50 counts now.
        hits = [hit("NODE_1", ISLAND_START - 450, 500)]
        self.assertEqual(summarize.supporting_contigs(hits, ISLAND_START, ISLAND_END, 200), [])

    def test_fixture_c_disjoint_island_blocks_pass_once_totalled(self):
        # Two blocks wholly inside the island, summing to >=200 bp of real overlap.
        hits = [hit("NODE_1", ISLAND_START, 120), hit("NODE_1", ISLAND_START + 1000, 100)]
        self.assertEqual(summarize.supporting_contigs(hits, ISLAND_START, ISLAND_END, 200),
                         ["NODE_1"])

    def test_a_contig_with_plenty_of_flank_but_no_island_hits_is_never_reported(self):
        # The exact failure mode from his report: lots of "aligned" bases, none in
        # the island the count is supposed to be about.
        hits = [hit("NODE_1", self.FLANK_START, 2000)]
        self.assertEqual(summarize.supporting_contigs(hits, ISLAND_START, ISLAND_END, 200), [])
        # Coverage/breadth is unaffected -- it was never scoped to the flank.
        self.assertEqual(summarize.union_length(hits, ISLAND_START, ISLAND_END), 0)


class AgreementNoLongerFollowsTheFilter(unittest.TestCase):
    """The case from the audit, in miniature."""

    def setUp(self):
        # MEGAHIT: one long block. metaSPAdes: the same island span in short pieces.
        self.megahit = [hit("M1", ISLAND_START, 1200)]
        self.metaspades = [hit(f"S{i}", ISLAND_START + i * 260, 250) for i in range(5)]

    def coverage(self, hits):
        return summarize.union_length(hits, ISLAND_START, ISLAND_END)

    def test_before_the_fix_the_two_assemblers_looked_discordant(self):
        # Simulate the old behaviour: drop every block under 300 bp first.
        old = lambda hits: self.coverage([h for h in hits if h["aligned"] >= 300])
        self.assertEqual(old(self.megahit), 1200)
        self.assertEqual(old(self.metaspades), 0)          # every piece discarded
        call, agreement = summarize.structural_call(
            old(self.megahit) / 50767, old(self.metaspades) / 50767)
        self.assertEqual(agreement, "single_assembler_only")

    def test_after_the_fix_they_agree(self):
        mh = self.coverage(self.megahit) / 50767
        ms = self.coverage(self.metaspades) / 50767
        self.assertGreater(ms, 0)
        _call, agreement = summarize.structural_call(mh, ms)
        self.assertEqual(agreement, "concordant")

    def test_the_contig_count_does_not_affect_the_call(self):
        # It is reported, not decided on: structural_call takes coverage only.
        import inspect
        self.assertEqual(list(inspect.signature(summarize.structural_call).parameters), ["mh", "ms"])


if __name__ == "__main__":
    unittest.main()
