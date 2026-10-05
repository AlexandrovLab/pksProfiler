"""The read-level pks specificity decision rule: retention + CI, recomputed tier
(classify_tumor_pks_evidence's own classify(), not a reinvented parallel definition),
non-megasynthase breadth, and decoy-absorption transparency.

Every case here is built from hand-written Arm A / Arm B TSVs in the documented
schema (see scripts/specificity_decision.py's module docstring) -- exercising the
decision LOGIC, not the aligner. Real-tool-driven Arm B classification is covered in
test_competitive_realign_score.py; this file assumes that arm's output is correct and
tests what happens next.
"""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("specificity_decision", ROOT / "scripts/specificity_decision.py")
sd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sd)


def write_tsv(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


class Args:
    """A stand-in for argparse.Namespace with this script's defaults, so each test
    only overrides what it cares about."""
    arm_a_target_clade = "pks_plus_clade"
    arm_a_min_confidence = 0.65
    arm_a_min_unique_kmers = 3
    island_length = 1000
    survivor_depth = None
    multi_gene_reads = 5
    multi_gene_genes = 3
    multi_gene_breadth = .01
    broad_island_reads = 30
    broad_island_genes = 8
    broad_island_breadth = .075
    extensive_island_reads = 100
    extensive_island_genes = 10
    extensive_island_breadth = .15
    retention_threshold = 0.8
    retention_ci_lower_threshold = 0.6
    bootstrap_iterations = 1000
    bootstrap_seed = 0
    max_tier_drop = 1
    non_megasynthase_min_genes = 2
    absorber_fraction_threshold = 0.2
    absorber_near_threshold_margin = 5.0
    min_supporting_reads = 10
    max_unmeasured_fraction = 0.05


def make_reads(n, gene_cycle, overlap_bp=50):
    """n read ids r0..r(n-1), cycled across gene_cycle, each with overlap_bp bases."""
    reads = [f"r{i}" for i in range(n)]
    gene_assignment = {r: (gene_cycle[i % len(gene_cycle)], overlap_bp) for i, r in enumerate(reads)}
    return reads, gene_assignment


def arm_a_all_confirm(reads, clade="pks_plus_clade", confidence=0.9, kmers=10):
    return {r: {"read_id": r, "assigned_clade": clade, "confidence": str(confidence),
                "unique_kmer_count": str(kmers)} for r in reads}


def arm_b_row(read_id, classification, decoy_contig=None, delta_as=None):
    return {"read_id": read_id, "b1_as_true": "", "b1_as_decoy": "",
            "b1_decoy_contig": decoy_contig or "", "b1_delta_as": "" if delta_as is None else delta_as,
            "b1_mapq": "", "b1_call": "", "b2_as_true": "", "b2_nm_true": "",
            "b2_as_decoy": "", "b2_nm_decoy": "", "b2_decoy_contig": decoy_contig or "",
            "b2_call": "", "classification": classification}


class RetentionAndVerdictBoundaries(unittest.TestCase):
    """R and its CI lower bound gate SURVIVES; the spec's two named failure/underpower
    conditions are tested at their exact numeric boundaries, not just broadly."""

    def _run(self, n_supporting, n_confirmed_true, original_tier="multi_gene", **overrides):
        reads, gene_assignment = make_reads(
            n_supporting, ["clbA", "clbD", "clbG", "clbM", "clbP", "clbQ", "clbS"])
        arm_a = arm_a_all_confirm(reads)
        arm_b = {r: arm_b_row(r, "confirmed_true" if i < n_confirmed_true else "confirmed_decoy",
                              decoy_contig="decoyX")
                 for i, r in enumerate(reads)}
        evidence = {"sample": "S1", "read_evidence": original_tier, "pks_reads": str(n_supporting),
                    "clb_genes_detected": "7", "island_breadth_1x": "0.5",
                    "island_breadth_2x": "0.3", "island_breadth_3x": "0.1"}
        args = Args()
        for key, value in overrides.items():
            setattr(args, key, value)
        return sd.decide("S1", evidence, gene_assignment, arm_a, arm_b, args)

    def test_underpowered_regardless_of_a_clean_point_estimate(self):
        # 8 supporting reads, all 8 confirmed -- R=1.0, a perfect point estimate --
        # must still be UNDERPOWERED because N_supporting < 10.
        result = self._run(n_supporting=8, n_confirmed_true=8)
        self.assertEqual(result["verdict"], sd.UNDERPOWERED)

    def test_exactly_ten_supporting_reads_is_not_underpowered(self):
        result = self._run(n_supporting=10, n_confirmed_true=10)
        self.assertNotEqual(result["verdict"], sd.UNDERPOWERED)

    def test_retention_below_one_half_is_not_distinguishable(self):
        result = self._run(n_supporting=20, n_confirmed_true=9)  # R = 0.45
        self.assertEqual(result["verdict"], sd.NOT_DISTINGUISHABLE)

    def test_retention_exactly_at_point_eight_with_ci_lower_at_point_six_survives(self):
        # Same fixture proven against specificity_decision.py's real bootstrap in
        # manual smoke testing: 20 reads, 16 confirmed -> R=0.8, CI lower lands at
        # exactly 0.6 with seed 0 / 1000 iterations.
        result = self._run(n_supporting=20, n_confirmed_true=16)
        self.assertEqual(result["retention"], 0.8)
        self.assertGreaterEqual(result["retention_ci_lower"], 0.6)
        self.assertEqual(result["verdict"], sd.SURVIVES)

    def test_retention_just_under_threshold_is_not_distinguishable(self):
        result = self._run(n_supporting=20, n_confirmed_true=15)  # R = 0.75 < 0.8
        self.assertEqual(result["verdict"], sd.NOT_DISTINGUISHABLE)

    def test_high_retention_with_a_wide_ci_can_still_fail_on_ci_lower(self):
        # A small sample can hit R>=0.8 while the bootstrap CI lower bound still
        # sits below 0.6 -- the point estimate alone is not enough.
        result = self._run(n_supporting=10, n_confirmed_true=8)  # R = 0.8, small n
        if result["retention_ci_lower"] < 0.6:
            self.assertEqual(result["verdict"], sd.NOT_DISTINGUISHABLE)

    def test_verdict_has_only_three_named_outcomes(self):
        for n_supp, n_conf in ((8, 8), (20, 9), (20, 16), (20, 15)):
            result = self._run(n_supporting=n_supp, n_confirmed_true=n_conf)
            self.assertIn(result["verdict"], (sd.SURVIVES, sd.NOT_DISTINGUISHABLE, sd.UNDERPOWERED))


class TierRecomputationUsesThePipelinesOwnDefinition(unittest.TestCase):
    def test_tier_drop_of_exactly_one_still_passes(self):
        self.assertEqual(sd.tier_drop("broad_island", "multi_gene"), 1)

    def test_tier_drop_of_two_is_computed_correctly(self):
        self.assertEqual(sd.tier_drop("extensive_island", "multi_gene"), 2)

    def test_tier_unchanged_is_zero_drop(self):
        self.assertEqual(sd.tier_drop("multi_gene", "multi_gene"), 0)

    def test_recompute_tier_calls_the_pipelines_own_classify_function(self):
        # Not a reinvented parallel definition: this is a direct pass-through to
        # classify_tumor_pks_evidence.classify(), imported by file path.
        thresholds = {"multi_gene": (5, 3, 0.01), "broad_island": (30, 8, 0.075),
                      "extensive_island": (100, 10, 0.15)}
        self.assertEqual(sd.classify_tumor.classify(5, 3, 0.01, thresholds), "multi_gene")
        self.assertEqual(sd.classify_tumor.classify(0, 0, 0.0, thresholds), "negative")

    def test_two_of_seven_non_megasynthase_genes_is_the_default_floor(self):
        gene_assignment = {"r0": ("clbA", 10), "r1": ("clbD", 10), "r2": ("clbB", 10)}
        covered = sd.recompute_genes({"r0", "r1", "r2"}, gene_assignment, set(sd.NON_MEGASYNTHASE_GENES))
        # clbB is a megasynthase gene and must not count toward this floor.
        self.assertEqual(covered, {"clbA", "clbD"})


class BreadthApproximationIsMonotonicAndLabelled(unittest.TestCase):
    def test_approximate_breadth_uses_overlap_bp_sum_over_island_length(self):
        gene_assignment = {"r0": ("clbA", 200), "r1": ("clbD", 300)}
        breadth, source = sd.recompute_breadth({"r0", "r1"}, gene_assignment, 1000, None)
        self.assertAlmostEqual(breadth, 0.5)
        self.assertEqual(source, "approximate")

    def test_approximate_breadth_is_capped_at_one(self):
        gene_assignment = {"r0": ("clbA", 2000)}
        breadth, source = sd.recompute_breadth({"r0"}, gene_assignment, 1000, None)
        self.assertEqual(breadth, 1.0)

    def test_exact_breadth_uses_survivor_depth_file_when_given(self):
        with tempfile.TemporaryDirectory() as tmp:
            depth_path = Path(tmp) / "depth.tsv"
            island_length = 100
            with depth_path.open("w") as handle:
                for position in range(island_length):
                    handle.write(f"chr\t{position + 1}\t{1 if position < 40 else 0}\n")
            breadth, source = sd.recompute_breadth({"r0"}, {}, island_length, depth_path)
            self.assertEqual(source, "exact")
            self.assertAlmostEqual(breadth, 0.4)


class AbsorberReportingNeverSilentlyFoldsIntoPassFail(unittest.TestCase):
    """Criterion 4: a decoy absorbing >20% of the dropped reads at a near-threshold
    margin must be named and flagged, not folded into a bare pass/fail number."""

    def test_a_concentrated_near_threshold_absorber_is_flagged(self):
        reads, gene_assignment = make_reads(20, ["clbA", "clbD", "clbG", "clbM"])
        arm_a = arm_a_all_confirm(reads)
        arm_b = {}
        for i, r in enumerate(reads):
            if i < 16:
                arm_b[r] = arm_b_row(r, "confirmed_true")
            else:
                # 4 dropped reads, all absorbed by the same decoy, all near threshold.
                arm_b[r] = arm_b_row(r, "confirmed_decoy", decoy_contig="decoyX", delta_as=-2)
        evidence = {"sample": "S1", "read_evidence": "multi_gene", "pks_reads": "20",
                    "clb_genes_detected": "4", "island_breadth_1x": "0.1",
                    "island_breadth_2x": "0.05", "island_breadth_3x": "0.01"}
        result = sd.decide("S1", evidence, gene_assignment, arm_a, arm_b, Args())
        self.assertTrue(result["manual_followup_required"])
        self.assertIn("decoyX", result["flagged_absorbers"])

    def test_a_landslide_absorber_far_from_threshold_is_not_flagged(self):
        reads, gene_assignment = make_reads(20, ["clbA", "clbD", "clbG", "clbM"])
        arm_a = arm_a_all_confirm(reads)
        arm_b = {}
        for i, r in enumerate(reads):
            if i < 16:
                arm_b[r] = arm_b_row(r, "confirmed_true")
            else:
                # Same decoy, same fraction, but a landslide margin (delta_AS well
                # beyond the near-threshold window) -- not a close call, so not
                # flagged for manual follow-up.
                arm_b[r] = arm_b_row(r, "confirmed_decoy", decoy_contig="decoyX", delta_as=-500)
        evidence = {"sample": "S1", "read_evidence": "multi_gene", "pks_reads": "20",
                    "clb_genes_detected": "4", "island_breadth_1x": "0.1",
                    "island_breadth_2x": "0.05", "island_breadth_3x": "0.01"}
        result = sd.decide("S1", evidence, gene_assignment, arm_a, arm_b, Args())
        self.assertFalse(result["manual_followup_required"])

    def test_dropped_reads_scattered_across_many_decoys_are_not_flagged(self):
        # 30 reads, 20 confirmed_true, 10 dropped spread across 6 different decoys
        # (fractions of the 10 non-survivors: 2,2,2,2,1,1 -> 20%,20%,20%,20%,10%,10%,
        # all at or below the 20% floor, none strictly exceeding it).
        reads, gene_assignment = make_reads(30, ["clbA", "clbD", "clbG", "clbM"])
        arm_a = arm_a_all_confirm(reads)
        decoy_names = ["decoy_1", "decoy_1", "decoy_2", "decoy_2", "decoy_3", "decoy_3",
                       "decoy_4", "decoy_4", "decoy_5", "decoy_6"]
        arm_b = {}
        for i, r in enumerate(reads):
            if i < 20:
                arm_b[r] = arm_b_row(r, "confirmed_true")
            else:
                arm_b[r] = arm_b_row(r, "confirmed_decoy", decoy_contig=decoy_names[i - 20], delta_as=-1)
        evidence = {"sample": "S1", "read_evidence": "multi_gene", "pks_reads": "30",
                    "clb_genes_detected": "4", "island_breadth_1x": "0.1",
                    "island_breadth_2x": "0.05", "island_breadth_3x": "0.01"}
        result = sd.decide("S1", evidence, gene_assignment, arm_a, arm_b, Args())
        self.assertFalse(result["manual_followup_required"])

    def test_flagging_an_absorber_does_not_by_itself_block_survives(self):
        # Naming the absorber satisfies criterion 4; it does not veto an otherwise
        # passing sample (a concentrated absorber may be a real unpanelled relative,
        # which is surfaced for manual follow-up, not treated as automatic failure).
        reads, gene_assignment = make_reads(
            20, ["clbA", "clbD", "clbG", "clbM", "clbP", "clbQ", "clbS"])
        arm_a = arm_a_all_confirm(reads)
        arm_b = {}
        for i, r in enumerate(reads):
            if i < 16:
                arm_b[r] = arm_b_row(r, "confirmed_true")
            else:
                arm_b[r] = arm_b_row(r, "confirmed_decoy", decoy_contig="decoyX", delta_as=-1)
        evidence = {"sample": "S1", "read_evidence": "multi_gene", "pks_reads": "20",
                    "clb_genes_detected": "7", "island_breadth_1x": "0.5",
                    "island_breadth_2x": "0.3", "island_breadth_3x": "0.1"}
        result = sd.decide("S1", evidence, gene_assignment, arm_a, arm_b, Args())
        self.assertTrue(result["manual_followup_required"])
        self.assertEqual(result["verdict"], sd.SURVIVES)


class ArmAGating(unittest.TestCase):
    def test_confidence_below_floor_fails_arm_a(self):
        row = {"read_id": "r0", "assigned_clade": "pks_plus_clade", "confidence": "0.5",
               "unique_kmer_count": "10"}
        self.assertFalse(sd.arm_a_confirms(row, "pks_plus_clade", 0.65, 3))

    def test_unique_kmer_count_below_floor_fails_arm_a(self):
        row = {"read_id": "r0", "assigned_clade": "pks_plus_clade", "confidence": "0.9",
               "unique_kmer_count": "1"}
        self.assertFalse(sd.arm_a_confirms(row, "pks_plus_clade", 0.65, 3))

    def test_collapsed_clade_assignment_fails_arm_a(self):
        # M1's failure mode again: a genus-level LCA is not a specific pks+ call.
        row = {"read_id": "r0", "assigned_clade": "genus_only", "confidence": "0.99",
               "unique_kmer_count": "50"}
        self.assertFalse(sd.arm_a_confirms(row, "pks_plus_clade", 0.65, 3))

    def test_a_read_missing_from_arm_a_entirely_fails_closed(self):
        self.assertFalse(sd.arm_a_confirms(None, "pks_plus_clade", 0.65, 3))

    def test_exact_threshold_values_pass(self):
        row = {"read_id": "r0", "assigned_clade": "pks_plus_clade", "confidence": "0.65",
               "unique_kmer_count": "3"}
        self.assertTrue(sd.arm_a_confirms(row, "pks_plus_clade", 0.65, 3))


class BootstrapCI(unittest.TestCase):
    def test_all_survive_gives_a_degenerate_ci_at_one(self):
        lower, upper = sd.bootstrap_ci([1] * 20, iterations=200, seed=1)
        self.assertEqual((lower, upper), (1.0, 1.0))

    def test_none_survive_gives_a_degenerate_ci_at_zero(self):
        lower, upper = sd.bootstrap_ci([0] * 20, iterations=200, seed=1)
        self.assertEqual((lower, upper), (0.0, 0.0))

    def test_is_deterministic_for_a_fixed_seed(self):
        outcomes = [1] * 16 + [0] * 4
        first = sd.bootstrap_ci(outcomes, iterations=500, seed=42)
        second = sd.bootstrap_ci(outcomes, iterations=500, seed=42)
        self.assertEqual(first, second)

    def test_empty_outcomes_does_not_divide_by_zero(self):
        self.assertEqual(sd.bootstrap_ci([], iterations=100, seed=0), (0.0, 0.0))


def unaligned_row(read_id):
    """What competitive_realign_score.py writes for a read minimap2 left unmapped in all
    three alignments: no score anywhere, classification ambiguous."""
    return arm_b_row(read_id, "ambiguous")


class UnmeasuredReadsAreNotDecoyLosses(unittest.TestCase):
    """5 Oct 2026 V8 pilot: a read with no Arm-B row was a non-survivor whose absorbing
    decoy was 'unknown_missing_from_arm_b', and retention was divided by the tier's read
    count. A read the challenge never measured is not a read lost to a decoy."""

    EVIDENCE = {"sample": "S1", "read_evidence": "multi_gene", "pks_reads": "100",
                "clb_genes_detected": "7", "island_breadth_1x": "0.5",
                "island_breadth_2x": "0.3", "island_breadth_3x": "0.1"}
    GENES = ["clbA", "clbD", "clbG", "clbM", "clbP", "clbQ", "clbS"]

    def _decide(self, arm_b, gene_assignment, **overrides):
        args = Args()
        for key, value in overrides.items():
            setattr(args, key, value)
        return sd.decide("S1", self.EVIDENCE, gene_assignment, None, arm_b, args)

    def test_reads_missing_from_arm_b_are_unmeasured_not_absorbed(self):
        reads, genes = make_reads(100, self.GENES)
        # 70 confirmed, 30 never reached Arm B under their own name.
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads[:70]}
        result = self._decide(arm_b, genes)
        self.assertEqual(result["n_unmeasured_missing_from_arm_b"], 30)
        self.assertEqual(result["n_measured"], 70)
        self.assertEqual(result["n_lost"], 0)
        self.assertEqual(result["retention"], 1.0)
        self.assertAlmostEqual(result["retention_if_unmeasured_lost"], 0.7)
        self.assertNotIn("unknown_missing_from_arm_b", result["absorbers"])
        self.assertEqual(result["flagged_absorbers"], "")
        self.assertFalse(result["manual_followup_required"])

    def test_too_many_unmeasured_reads_is_incomplete_not_a_verdict(self):
        reads, genes = make_reads(100, self.GENES)
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads[:70]}
        self.assertEqual(self._decide(arm_b, genes)["verdict"], sd.INCOMPLETE)
        # ...and it is never reported as NOT_DISTINGUISHABLE, which would claim the reads
        # were challenged and lost.
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads[:40]}
        self.assertEqual(self._decide(arm_b, genes)["verdict"], sd.INCOMPLETE)

    def test_a_few_unmeasured_reads_under_the_ceiling_still_decide(self):
        reads, genes = make_reads(100, self.GENES)
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads[:96]}
        result = self._decide(arm_b, genes)
        self.assertEqual(result["n_unmeasured"], 4)
        self.assertEqual(result["verdict"], sd.SURVIVES)

    def test_a_read_aligning_to_neither_panel_side_is_unmeasured_not_absorbed(self):
        reads, genes = make_reads(100, self.GENES)
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads[:70]}
        arm_b.update({r: unaligned_row(r) for r in reads[70:]})
        result = self._decide(arm_b, genes)
        self.assertEqual(result["n_unmeasured_unaligned"], 30)
        self.assertEqual(result["n_lost"], 0)
        self.assertEqual(result["absorbers"], {})
        self.assertEqual(result["verdict"], sd.INCOMPLETE)

    def test_real_decoy_losses_are_still_counted_and_flagged(self):
        reads, genes = make_reads(100, self.GENES)
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads[:60]}
        arm_b.update({r: arm_b_row(r, "confirmed_decoy", decoy_contig="decoyX", delta_as=-2)
                      for r in reads[60:98]})
        # two missing reads alongside 38 genuine losses
        result = self._decide(arm_b, genes)
        self.assertEqual(result["n_lost"], 38)
        self.assertEqual(result["n_unmeasured"], 2)
        self.assertAlmostEqual(result["retention"], 60 / 98)
        self.assertIn("decoyX", result["flagged_absorbers"])
        self.assertNotIn("missing", result["flagged_absorbers"])
        self.assertEqual(result["verdict"], sd.NOT_DISTINGUISHABLE)

    def test_retention_and_its_ci_use_the_same_denominator(self):
        # pks_reads (100) is larger than the evaluated read set (50): the point estimate
        # must not be divided by one and the CI computed over the other.
        reads, genes = make_reads(50, self.GENES)
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads}
        result = self._decide(arm_b, genes)
        self.assertEqual(result["retention"], 1.0)
        self.assertEqual((result["retention_ci_lower"], result["retention_ci_upper"]), (1.0, 1.0))
        self.assertEqual(result["n_supporting"], 100)
        self.assertEqual(result["n_evaluated"], 50)

    def test_too_few_measured_reads_is_underpowered(self):
        reads, genes = make_reads(100, self.GENES)
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads[:5]}
        self.assertEqual(self._decide(arm_b, genes)["verdict"], sd.UNDERPOWERED)

    def test_the_per_read_table_names_each_unmeasured_kind(self):
        reads, genes = make_reads(20, self.GENES)
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads[:16]}
        arm_b[reads[16]] = unaligned_row(reads[16])
        arm_b[reads[17]] = arm_b_row(reads[17], "confirmed_decoy", decoy_contig="decoyX")
        result = self._decide(arm_b, genes)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "per_read.tsv"
            sd.write_per_read(result, {}, arm_b, path)
            with path.open() as handle:
                rows = {r["read_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(rows[reads[16]]["status"], sd.UNALIGNED)
        self.assertEqual(rows[reads[16]]["absorbing_decoy"], "")
        self.assertEqual(rows[reads[17]]["status"], sd.LOST)
        self.assertEqual(rows[reads[17]]["absorbing_decoy"], "decoyX")
        self.assertEqual(rows[reads[18]]["status"], sd.MISSING_FROM_ARM_B)
        self.assertEqual(rows[reads[18]]["absorbing_decoy"], "")
        self.assertEqual(rows[reads[0]]["status"], sd.SURVIVED)

    def test_with_arm_a_a_read_missing_from_arm_a_is_unmeasured(self):
        reads, genes = make_reads(100, self.GENES)
        arm_a = arm_a_all_confirm(reads[:90])
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads}
        result = sd.decide("S1", self.EVIDENCE, genes, arm_a, arm_b, Args())
        self.assertEqual(result["n_unmeasured_missing_from_arm_a"], 10)
        self.assertEqual(result["n_lost"], 0)
        self.assertEqual(result["verdict"], sd.INCOMPLETE)

    def test_an_arm_a_rejection_is_a_loss_named_as_such(self):
        reads, genes = make_reads(20, self.GENES)
        arm_a = arm_a_all_confirm(reads)
        arm_a[reads[0]]["confidence"] = "0.1"
        arm_b = {r: arm_b_row(r, "confirmed_true") for r in reads}
        result = sd.decide("S1", self.EVIDENCE, genes, arm_a, arm_b, Args())
        self.assertEqual(result["n_lost"], 1)
        self.assertIn("arm_a_not_confirmed", result["absorbers"])


if __name__ == "__main__":
    unittest.main()
