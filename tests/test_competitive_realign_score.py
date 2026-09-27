"""Arm B (dual-index competitive reassignment) of the read-level pks specificity check.

Ludmil's concern: MAPQ >= 40 against IHE3034 alone shows confident mapping to *that*
reference, not uniqueness against the broader microbial universe. clbB/clbC/clbH/clbI/
clbJ/clbK/clbN/clbO carry NRPS/PKS megasynthase domains conserved across unrelated
bacterial gene clusters -- the same domain-homology failure M1 already caught once at
the assembled-MAG-bin level (test_pks_positive_definition.py: "clbB and clbK hit every
bin by domain homology"). This is that failure mode at the read level.

`ExactAgreementOnRealAlignments` builds a tiny synthetic true reference and decoy
genome that share one "conserved domain" block (a stand-in for a megasynthase domain),
runs the *real* pinned minimap2 (never a hand-written SAM fixture) against a
concatenated index (B1) and each reference alone (B2), and asserts
competitive_realign_score.py reaches the correct verdict for four read origins:
a true-only region, a decoy-only region, a read mostly in the shared block but
extending into a true-side flank (small positive margin -- still true), and a read
entirely inside the shared block with no flank on either side at all (an exact tie --
this must come out "ambiguous", never force-assigned to either side; that is the
whole point of running two independent index architectures and requiring both to
agree, per the AstraZeneca-NGS/disambiguate convention this arm follows).

`ClassificationEdgeCases` covers the pure decision functions directly (B1's margin
gate, B2's AS-then-NM tiebreak, and the "both arms must agree" rule) with numbers a
real aligner would be unlikely to hand us on demand -- an exact NM tiebreak, and a
read absent from one side of the panel entirely -- without needing a second real
alignment run for each one.
"""
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import competitive_realign_score as arm_b  # noqa: E402

MINIMAP2_CANDIDATES = (
    "/tscc/projects/ps-lalexandrov/amabbasi/anaconda3/envs/minimap2_env/bin/minimap2",
    "/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/tools/minimap2/minimap2",
)


def find_minimap2():
    for candidate in MINIMAP2_CANDIDATES:
        if Path(candidate).is_file():
            return candidate
    return shutil.which("minimap2")


MINIMAP2 = find_minimap2()


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_fasta(path, header, sequence):
    with Path(path).open("w") as handle:
        handle.write(f">{header}\n")
        for i in range(0, len(sequence), 70):
            handle.write(sequence[i:i + 70] + "\n")


def write_fastq(path, reads):
    with Path(path).open("w") as handle:
        for name, seq in reads.items():
            handle.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def run_minimap2(reference, reads_fastq, output_sam, extra_args=()):
    cmd = [MINIMAP2, "-a", "-x", "sr", *extra_args, str(reference), str(reads_fastq)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"minimap2 failed: {result.stderr}"
    Path(output_sam).write_text(result.stdout)


@unittest.skipUnless(MINIMAP2, "minimap2 not found on any known path")
class ExactAgreementOnRealAlignments(unittest.TestCase):
    """Real minimap2 output, real SAM parsing -- not a hand-written fixture."""

    @classmethod
    def setUpClass(cls):
        rng = random.Random(20260923)
        # A 300bp block identical in both references: the stand-in for a conserved
        # megasynthase domain that clbB/clbK-like genes share with an unrelated cluster.
        shared_block = rand_seq(300, rng)
        true_seq = rand_seq(800, rng) + shared_block + rand_seq(800, rng)
        decoy_seq = rand_seq(800, rng) + shared_block + rand_seq(800, rng)

        cls.tmp = tempfile.TemporaryDirectory()
        d = Path(cls.tmp.name)
        write_fasta(d / "true_ref.fasta", "true_ihe3034_locus", true_seq)
        write_fasta(d / "decoy_ref.fasta", "decoy_species_X", decoy_seq)
        with (d / "concatenated.fasta").open("w") as handle:
            handle.write((d / "true_ref.fasta").read_text())
            handle.write((d / "decoy_ref.fasta").read_text())
        (d / "true_contigs.txt").write_text("true_ihe3034_locus\n")
        (d / "decoy_contigs.txt").write_text("decoy_species_X\n")

        reads = {
            # Entirely within the true-only flank: no decoy alignment should exist at all.
            "read_true_specific": true_seq[100:250],
            # Entirely within the decoy-only flank: mirror case.
            "read_decoy_specific": decoy_seq[100:250],
            # Mostly the shared block (conserved-domain reads exactly like clbB/clbK),
            # but extended 20bp into the true-side flank -- true should still win, by
            # a real, measurable margin, not by construction.
            "read_shared_leaning_true": true_seq[800 + 170:800 + 300] + true_seq[1100:1120],
            # Entirely inside the shared block, no flank on either side: identical
            # sequence exists in both references, so no aligner can prefer one over
            # the other. This is the case that must come out ambiguous.
            "read_pure_shared": true_seq[800 + 75:800 + 225],
        }
        write_fastq(d / "reads.fastq", reads)

        run_minimap2(d / "concatenated.fasta", d / "reads.fastq", d / "b1.sam",
                     ("-N", "5", "-p", "0.1", "--secondary=yes"))
        run_minimap2(d / "true_ref.fasta", d / "reads.fastq", d / "b2_true.sam")
        run_minimap2(d / "decoy_ref.fasta", d / "reads.fastq", d / "b2_decoy.sam")

        cls.d = d
        cls.true_contigs = arm_b.read_contig_list(d / "true_contigs.txt")
        cls.decoy_contigs = arm_b.read_contig_list(d / "decoy_contigs.txt")
        cls.rows = {
            row["read_id"]: row
            for row in arm_b.score_reads(
                list(arm_b.parse_sam_records(d / "b1.sam")),
                list(arm_b.parse_sam_records(d / "b2_true.sam")),
                list(arm_b.parse_sam_records(d / "b2_decoy.sam")),
                cls.true_contigs, cls.decoy_contigs,
            )
        }

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_true_only_region_confirms_true(self):
        self.assertEqual(self.rows["read_true_specific"]["classification"], "confirmed_true")

    def test_decoy_only_region_confirms_decoy(self):
        self.assertEqual(self.rows["read_decoy_specific"]["classification"], "confirmed_decoy")

    def test_conserved_domain_read_with_a_true_side_flank_still_confirms_true(self):
        row = self.rows["read_shared_leaning_true"]
        self.assertEqual(row["classification"], "confirmed_true")
        # The margin is real, not zero -- the true-side flank is what breaks the tie.
        self.assertGreater(row["b1_delta_as"], 0)

    def test_pure_conserved_domain_read_is_ambiguous_not_force_assigned(self):
        """The central case: identical sequence in both references. Neither B1 nor B2
        has any basis to prefer one side, and the decision rule must say so rather
        than defaulting to whichever reference happened to be listed first (both
        SAM records here put this read's primary alignment on the true reference --
        if this test passed by checking "is it on the true contig" instead of the
        classification field, it would hide exactly the failure mode this arm exists
        to catch)."""
        row = self.rows["read_pure_shared"]
        self.assertEqual(row["b1_delta_as"], 0)
        self.assertEqual(row["classification"], "ambiguous")

    def test_decoy_contig_is_named_for_the_reads_that_do_not_survive(self):
        # Full transparency requirement: every non-surviving read's absorbing decoy
        # must be identifiable from this table, not just a pass/fail bit.
        self.assertEqual(self.rows["read_decoy_specific"]["b1_decoy_contig"], "decoy_species_X")
        self.assertEqual(self.rows["read_pure_shared"]["b1_decoy_contig"], "decoy_species_X")


class ClassificationEdgeCases(unittest.TestCase):
    """Pure decision-function coverage for cases a real aligner will not hand us
    on demand in a four-read fixture."""

    def test_b1_only_true_family_present_is_true(self):
        self.assertEqual(arm_b.classify_b1(120, None), "true")

    def test_b1_only_decoy_family_present_is_decoy(self):
        self.assertEqual(arm_b.classify_b1(None, 90), "decoy")

    def test_b1_neither_family_present_is_ambiguous(self):
        self.assertEqual(arm_b.classify_b1(None, None), "ambiguous")

    def test_b1_margin_below_min_delta_is_ambiguous(self):
        # A real margin exists (3 AS units) but is below the configured floor.
        self.assertEqual(arm_b.classify_b1(103, 100, min_delta_as=5), "ambiguous")
        self.assertEqual(arm_b.classify_b1(103, 100, min_delta_as=2), "true")

    def test_b2_as_tiebreak_by_nm_lower_wins(self):
        # Disambiguate convention: equal AS, edit distance decides.
        self.assertEqual(arm_b.classify_b2(as_true=100, nm_true=1, as_decoy=100, nm_decoy=3), "true")
        self.assertEqual(arm_b.classify_b2(as_true=100, nm_true=3, as_decoy=100, nm_decoy=1), "decoy")

    def test_b2_exact_tie_on_as_and_nm_is_ambiguous(self):
        self.assertEqual(arm_b.classify_b2(as_true=100, nm_true=2, as_decoy=100, nm_decoy=2), "ambiguous")

    def test_b2_absent_from_one_side_uses_the_other(self):
        self.assertEqual(arm_b.classify_b2(as_true=100, nm_true=0, as_decoy=None, nm_decoy=None), "true")
        self.assertEqual(arm_b.classify_b2(as_true=None, nm_true=None, as_decoy=80, nm_decoy=1), "decoy")

    def test_final_classification_requires_both_arms_to_agree(self):
        self.assertEqual(arm_b.final_classification("true", "true"), "confirmed_true")
        self.assertEqual(arm_b.final_classification("decoy", "decoy"), "confirmed_decoy")
        # Disagreement between arms -- exactly the "single concatenated index can
        # distort MAPQ/AS by panel composition alone" risk the dual-index design
        # exists to catch -- must not be silently resolved either way.
        self.assertEqual(arm_b.final_classification("true", "decoy"), "ambiguous")
        self.assertEqual(arm_b.final_classification("true", "ambiguous"), "ambiguous")
        self.assertEqual(arm_b.final_classification("ambiguous", "ambiguous"), "ambiguous")

    def test_true_and_decoy_contigs_cannot_overlap(self):
        # A contig listed on both sides would make every comparison against it a
        # reference-against-itself tie; caught before any read is even scored.
        with self.assertRaises(ValueError):
            arm_b.validate_contig_lists({"shared_contig"}, {"shared_contig"})


if __name__ == "__main__":
    unittest.main()
