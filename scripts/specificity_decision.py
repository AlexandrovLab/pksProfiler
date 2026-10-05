#!/usr/bin/env python3
"""V8a decision rule: LOCUS-DETECTION SPECIFICITY for the read-level pks call.

Ludmil, 25 Sep report, section 4, V8a: "Challenge pks-supporting reads against
pks-negative Enterobacterales, close homologous NRPS/PKS sequences, broad gut/oral
backgrounds, GRCh38, and UniVec. This asks whether the sequence really belongs to the
pks/clb locus." Every member of the panel behind this decision is pks-NEGATIVE.

The sister question -- which pks-positive organism carries the island -- is V8b, and the
report says the two "should not be conflated". It has its own script
(scripts/carrier_attribution.py), its own pks-POSITIVE panel and its own vocabulary
(CARRIER_ATTRIBUTED / CARRIER_GENUS_ONLY / CARRIER_AMBIGUOUS / CARRIER_UNDETERMINED).
Nothing here reads or reports carrier identity, and nothing there reports specificity.

Combines Arm A (KrakenUniq reclassification against a true+decoy panel) and Arm B
(competitive_realign_score.py's dual-index alignment-score-margin verdict) with the
sample's existing read_evidence.tsv (classify_tumor_pks_evidence.py's tier call) to
answer one question: does the pipeline's read-level pks-positive call survive being
re-checked against a panel of pks-negative decoys, or is it indistinguishable from
cross-mapping onto conserved megasynthase domains that clbB/clbC/clbH/clbI/clbJ/clbK/
clbN/clbO share with unrelated bacterial secondary-metabolite clusters (the same
domain-homology failure M1 already caught once at the assembled-MAG-bin level --
see tests/test_pks_positive_definition.py)?

A read survives only if BOTH arms confirm it independently:
  * Arm A: its KrakenUniq LCA assignment resolves specifically to the pks+ clade
    (not collapsed to genus/family) at or above --arm-a-min-confidence, AND its
    unique-k-mer count is at or above --arm-a-min-unique-kmers (KrakenUniq's own
    statistic for "how much of this call rests on broad unique coverage versus a
    handful of shared conserved k-mers").
  * Arm B: competitive_realign_score.py classified it "confirmed_true".

From the surviving set:
  N_supporting  = pks_reads from the original read_evidence.tsv (the tier's own read
                  count; reported, and the UNDERPOWERED floor, but not the denominator)
  N_evaluated   = the reads this check is asked about: every read in
                  --read-gene-assignment (reads already passing the pipeline's own
                  MAPQ>=40-vs-IHE3034 filter; this check never re-runs primary detection)
  N_unmeasured  = evaluated reads the challenge never measured, in two named kinds:
                    missing_from_arm_b   no Arm-B row at all (the read never reached the
                                         aligner under the name it is evaluated under --
                                         see the 5 Oct note below)
                    unaligned            an Arm-B row, but no alignment to EITHER side of
                                         the panel in either architecture
                  (and, with Arm A, missing_from_arm_a: no Arm-A row). An unmeasured read
                  is not a read lost to a decoy, and is never counted as one.
  N_measured    = N_evaluated - N_unmeasured
  N_survive     = reads confirmed by both arms
  R             = N_survive / N_measured, with a bootstrap 95% CI over the SAME measured
                  reads (resampling with replacement, default 1000 iterations, stdlib
                  `random` only). retention_if_unmeasured_lost = N_survive / N_evaluated
                  is reported beside it as the fail-closed lower bound.

  5 Oct 2026 (V8 real-data pilot). Before this, R was N_survive / N_supporting while its
  CI was over the evaluated reads, and every read without an Arm-B row was a non-survivor
  whose "absorbing decoy" was unknown_missing_from_arm_b. minimap2 -x sr pairs any two
  consecutive FASTQ records whose names differ only by /1 and /2 and strips the suffix,
  so on a v0.0.1 island FASTQ (names end in /1 and /2; mates are often adjacent in a
  coordinate-sorted BAM) a third of the reads lost their name in Arm B and were scored as
  lost to a decoy. The lanes now run minimap2 with --frag=no; this script no longer lets
  any such gap pass as biology.
  recomputed tier   = classify_tumor_pks_evidence.classify() on survivor reads only,
                      so tier recomputation is the pipeline's own definition, not a
                      reinvented parallel one. Needs survivor read counts, survivor
                      clb-gene coverage (from --read-gene-assignment, the same
                      read-to-gene file assign_reads_to_genes.py already produces),
                      and survivor island breadth. Breadth is EXACT when
                      --survivor-depth is given (a samtools-depth-shaped file over
                      the survivor-only read set, same schema classify_tumor_pks_
                      evidence.breadths() already reads); otherwise it is
                      APPROXIMATED from summed read-to-gene overlap_bp among
                      survivors divided by island length, which is monotonically
                      non-increasing as reads are dropped (same direction as the
                      true quantity) but double-counts overlapping reads, so it is
                      reported with breadth_source=approximate rather than silently
                      passed off as the real thing.
  non-megasynthase breadth = how many of clbA/clbD/clbG/clbM/clbP/clbQ/clbS (the pks-
                      specific genes that do NOT share the megasynthase genes' broadly
                      conserved domains) still have >=1 survivor read assigned.

Verdict:
  UNDERPOWERED                          N_supporting < --min-supporting-reads, or
                                         N_measured below it, regardless of the point
                                         estimate.
  INCOMPLETE_MEASUREMENT                more than --max-unmeasured-fraction of the
                                         evaluated reads were never measured: the check
                                         did not run on the read set it reports on, so
                                         it says that rather than either verdict below.
  SURVIVES                              R >= --retention-threshold AND its CI lower
                                         bound >= --retention-ci-lower-threshold AND
                                         the recomputed tier drops by at most
                                         --max-tier-drop AND non-megasynthase breadth
                                         covers >= --non-megasynthase-min-genes genes.
  NOT_DISTINGUISHABLE_FROM_CROSS_MAPPING otherwise (this includes R < 0.5 as named in
                                         the spec, the tier-drop or breadth criteria
                                         failing, and the region strictly between --
                                         there is no "maybe": failing SURVIVES on reads
                                         that WERE measured, without being UNDERPOWERED,
                                         is reported as not distinguishable from cross-
                                         mapping, fail-closed rather than silently
                                         optimistic).

A single decoy genome/species is never allowed to be silently folded into a bare
pass/fail count: every measured non-surviving read's absorbing decoy is in the per-read
report (unmeasured reads are reported by kind, never as an absorber),
and any decoy absorbing more than --absorber-fraction-threshold of the non-surviving
reads at a near-threshold margin (|delta_AS| <= --absorber-near-threshold-margin, or
an Arm-B "ambiguous" call, i.e. genuinely close calls rather than a landslide) is
named explicitly in the summary with manual_followup_required=True -- it may be a
real, currently-unpanelled close relative, not generic noise.
"""
import argparse
import csv
import importlib.util
import random
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "classify_tumor_pks_evidence", SCRIPTS_DIR / "classify_tumor_pks_evidence.py")
classify_tumor = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(classify_tumor)

MEGASYNTHASE_GENES = {"clbB", "clbC", "clbH", "clbI", "clbJ", "clbK", "clbN", "clbO"}
NON_MEGASYNTHASE_GENES = ("clbA", "clbD", "clbG", "clbM", "clbP", "clbQ", "clbS")

TIER_ORDER = ("negative", "localized_indeterminate", "multi_gene", "broad_island", "extensive_island")

SURVIVES = "SURVIVES"
NOT_DISTINGUISHABLE = "NOT_DISTINGUISHABLE_FROM_CROSS_MAPPING"
UNDERPOWERED = "UNDERPOWERED"
INCOMPLETE = "INCOMPLETE_MEASUREMENT"

# Per-read status. Only LOST reads are evidence against the locus.
SURVIVED = "survived"
LOST = "lost"
MISSING_FROM_ARM_B = "unmeasured_missing_from_arm_b"
MISSING_FROM_ARM_A = "unmeasured_missing_from_arm_a"
UNALIGNED = "unmeasured_unaligned"
UNMEASURED = (MISSING_FROM_ARM_B, MISSING_FROM_ARM_A, UNALIGNED)
ARM_B_SCORE_FIELDS = ("b1_as_true", "b1_as_decoy", "b2_as_true", "b2_as_decoy")


def read_single_row_tsv(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one data row in {path}, found {len(rows)}")
    return rows[0]


def read_gene_assignment(path):
    """assign_reads_to_genes.py's read_id/Gene/overlap_bp TSV -> {read_id: (gene, overlap_bp)}."""
    assignment = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            assignment[row["read_id"]] = (row["Gene"], int(row["overlap_bp"]))
    return assignment


def read_arm_a(path):
    """One row per read: read_id, assigned_clade, confidence, unique_kmer_count."""
    rows = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows[row["read_id"]] = row
    return rows


def read_arm_b(path):
    with Path(path).open() as handle:
        rows = {row["read_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    return rows


def arm_a_confirms(row, target_clade, min_confidence, min_unique_kmers):
    if row is None:
        return False
    if row["assigned_clade"] != target_clade:
        return False
    if float(row["confidence"]) < min_confidence:
        return False
    if int(row["unique_kmer_count"]) < min_unique_kmers:
        return False
    return True


def bootstrap_ci(outcomes, iterations, seed, alpha=0.05):
    """outcomes: list of 0/1, one per original supporting read. Resample with
    replacement `iterations` times and return (lower, upper) percentiles of the
    retention fraction. stdlib-only (no numpy dependency), matching the rest of
    this repo's scripts.
    """
    n = len(outcomes)
    if n == 0:
        return 0.0, 0.0
    rng = random.Random(seed)
    fractions = []
    for _ in range(iterations):
        resample = rng.choices(outcomes, k=n)
        fractions.append(sum(resample) / n)
    fractions.sort()
    lower_index = int((alpha / 2) * iterations)
    upper_index = min(iterations - 1, int((1 - alpha / 2) * iterations))
    return fractions[lower_index], fractions[upper_index]


def recompute_breadth(survivor_ids, gene_assignment, island_length, survivor_depth_path):
    """Returns (breadth_1x, source) where source is "exact" or "approximate"."""
    if survivor_depth_path:
        b1, _b2, _b3 = classify_tumor.breadths(survivor_depth_path, island_length)
        return b1, "exact"
    if not island_length:
        return 0.0, "approximate"
    total_overlap = sum(overlap for read_id, (_gene, overlap) in gene_assignment.items()
                        if read_id in survivor_ids)
    return min(1.0, total_overlap / island_length), "approximate"


def recompute_genes(survivor_ids, gene_assignment, gene_set=None):
    genes = {gene for read_id, (gene, _overlap) in gene_assignment.items()
             if read_id in survivor_ids and (gene_set is None or gene in gene_set)}
    return genes


def arm_b_unaligned(row):
    """An Arm-B row with no alignment score against either side of the panel in either
    architecture: minimap2 wrote the read as unmapped every time. Nothing was measured,
    so nothing can be said about a decoy."""
    return all((row.get(field) or "") == "" for field in ARM_B_SCORE_FIELDS)


def read_status(read_id, arm_a_rows, arm_b_rows, survives):
    """survived / lost / one of the UNMEASURED kinds, for one evaluated read."""
    if survives:
        return SURVIVED
    row = arm_b_rows.get(read_id)
    if row is None:
        return MISSING_FROM_ARM_B
    # competitive_realign_score.py calls a read with no alignment anywhere "ambiguous";
    # a confirmed_* call always rests on an alignment score.
    if row.get("classification") == "ambiguous" and arm_b_unaligned(row):
        return UNALIGNED
    if arm_a_rows is not None and read_id not in arm_a_rows:
        return MISSING_FROM_ARM_A
    return LOST


def absorbing_decoy(row):
    """The decoy a LOST read went to. A lost read with no decoy alignment at all lost on
    the arm-A gate or on arm disagreement, and is labelled so, not 'unknown'."""
    if row is None:
        return ""
    decoy = row.get("b1_decoy_contig") or row.get("b2_decoy_contig") or ""
    if decoy:
        return decoy
    if row.get("classification") == "confirmed_true":
        return "arm_a_not_confirmed"
    return "no_decoy_alignment"


def absorber_report(lost_reads, arm_b_rows, near_threshold_margin):
    """{decoy_name: {"count": n, "near_threshold_count": n}} over LOST reads only.
    Unmeasured reads never enter it: a read the challenge did not measure was not
    absorbed by anything."""
    absorbers = {}
    for read_id in lost_reads:
        row = arm_b_rows.get(read_id)
        decoy = absorbing_decoy(row)
        delta = (row or {}).get("b1_delta_as", "")
        is_ambiguous = (row or {}).get("classification") == "ambiguous"
        try:
            near = is_ambiguous or (delta != "" and abs(float(delta)) <= near_threshold_margin)
        except ValueError:
            near = True
        entry = absorbers.setdefault(decoy, {"count": 0, "near_threshold_count": 0})
        entry["count"] += 1
        if near:
            entry["near_threshold_count"] += 1
    return absorbers


def tier_drop(original_tier, recomputed_tier):
    """How many rungs recomputed_tier falls below original_tier on TIER_ORDER.
    Negative means recomputed_tier is HIGHER (cannot normally happen, since survivors
    are a subset of the original reads, but is not treated as a failure if it does)."""
    return TIER_ORDER.index(original_tier) - TIER_ORDER.index(recomputed_tier)


def decide(sample, read_evidence, gene_assignment, arm_a_rows, arm_b_rows, args):
    """arm_a_rows is None when Arm A was not run: the survival rule is then Arm B alone,
    and `arms_used` records that. Treating a missing Arm A as "every read fails Arm A"
    would report NOT_DISTINGUISHABLE for a check that was never performed."""
    arms_used = "arm_b_only" if arm_a_rows is None else "arm_a_and_arm_b"
    n_supporting = int(read_evidence["pks_reads"])
    original_tier = read_evidence["read_evidence"]

    all_read_ids = sorted(gene_assignment) if gene_assignment else sorted(
        set(arm_a_rows or {}) | set(arm_b_rows))

    survivors = set()
    status = {}
    for read_id in all_read_ids:
        arm_a_ok = True if arm_a_rows is None else arm_a_confirms(
            arm_a_rows.get(read_id), args.arm_a_target_clade,
            args.arm_a_min_confidence, args.arm_a_min_unique_kmers)
        arm_b_row = arm_b_rows.get(read_id)
        arm_b_ok = arm_b_row is not None and arm_b_row.get("classification") == "confirmed_true"
        survives = arm_a_ok and arm_b_ok
        status[read_id] = read_status(read_id, arm_a_rows, arm_b_rows, survives)
        if survives:
            survivors.add(read_id)

    n_evaluated = len(all_read_ids)
    unmeasured_counts = {kind: sum(1 for v in status.values() if v == kind) for kind in UNMEASURED}
    n_unmeasured = sum(unmeasured_counts.values())
    measured = [r for r in all_read_ids if status[r] not in UNMEASURED]
    n_measured = len(measured)
    unmeasured_fraction = (n_unmeasured / n_evaluated) if n_evaluated else 0.0

    n_survive = len(survivors)
    # Retention and its CI over the same reads: the ones the challenge measured.
    outcomes = [1 if r in survivors else 0 for r in measured]
    retention = (n_survive / n_measured) if n_measured else 0.0
    retention_if_unmeasured_lost = (n_survive / n_evaluated) if n_evaluated else 0.0
    ci_lower, ci_upper = bootstrap_ci(outcomes, args.bootstrap_iterations, args.bootstrap_seed)

    recomputed_genes = recompute_genes(survivors, gene_assignment)
    recomputed_breadth, breadth_source = recompute_breadth(
        survivors, gene_assignment, args.island_length, args.survivor_depth)
    thresholds = {
        "multi_gene": (args.multi_gene_reads, args.multi_gene_genes, args.multi_gene_breadth),
        "broad_island": (args.broad_island_reads, args.broad_island_genes, args.broad_island_breadth),
        "extensive_island": (args.extensive_island_reads, args.extensive_island_genes, args.extensive_island_breadth),
    }
    recomputed_tier = classify_tumor.classify(n_survive, len(recomputed_genes), recomputed_breadth, thresholds)
    drop = tier_drop(original_tier, recomputed_tier)
    tier_ok = drop <= args.max_tier_drop

    non_mega_covered = recompute_genes(survivors, gene_assignment, set(NON_MEGASYNTHASE_GENES))
    breadth_ok = len(non_mega_covered) >= args.non_megasynthase_min_genes

    non_survivors = [r for r in all_read_ids if r not in survivors]
    lost = [r for r in non_survivors if status[r] == LOST]
    absorbers = absorber_report(lost, arm_b_rows, args.absorber_near_threshold_margin)
    n_lost = len(lost)
    flagged = {}
    for decoy, stats in absorbers.items():
        fraction = stats["count"] / n_lost if n_lost else 0.0
        if fraction > args.absorber_fraction_threshold and stats["near_threshold_count"] > 0:
            flagged[decoy] = {"count": stats["count"], "fraction": fraction}
    manual_followup_required = bool(flagged)

    max_unmeasured = getattr(args, "max_unmeasured_fraction", 0.05)
    if n_supporting < args.min_supporting_reads or n_measured < args.min_supporting_reads:
        verdict = UNDERPOWERED
    elif unmeasured_fraction > max_unmeasured:
        verdict = INCOMPLETE
    elif (retention >= args.retention_threshold and ci_lower >= args.retention_ci_lower_threshold
          and tier_ok and breadth_ok):
        verdict = SURVIVES
    else:
        verdict = NOT_DISTINGUISHABLE

    return {
        "sample": sample,
        # Named so a row of this table can never be mistaken for a V8b carrier row.
        "check": "locus_specificity",
        "verdict": verdict,
        "arms_used": arms_used,
        "n_supporting": n_supporting,
        "n_evaluated": n_evaluated,
        "n_measured": n_measured,
        "n_unmeasured": n_unmeasured,
        "n_unmeasured_missing_from_arm_b": unmeasured_counts[MISSING_FROM_ARM_B],
        "n_unmeasured_unaligned": unmeasured_counts[UNALIGNED],
        "n_unmeasured_missing_from_arm_a": unmeasured_counts[MISSING_FROM_ARM_A],
        "unmeasured_fraction": unmeasured_fraction,
        "max_unmeasured_fraction": max_unmeasured,
        "n_survive": n_survive,
        "n_lost": n_lost,
        "retention": retention,
        "retention_if_unmeasured_lost": retention_if_unmeasured_lost,
        "retention_ci_lower": ci_lower,
        "retention_ci_upper": ci_upper,
        "bootstrap_iterations": args.bootstrap_iterations,
        "original_tier": original_tier,
        "recomputed_tier": recomputed_tier,
        "tier_drop": drop,
        "recomputed_clb_genes_detected": len(recomputed_genes),
        "recomputed_island_breadth_1x": recomputed_breadth,
        "breadth_source": breadth_source,
        "non_megasynthase_genes_covered": len(non_mega_covered),
        "non_megasynthase_genes": ",".join(sorted(non_mega_covered)),
        "non_megasynthase_min_genes": args.non_megasynthase_min_genes,
        "manual_followup_required": manual_followup_required,
        "flagged_absorbers": ";".join(
            f"{decoy}:{stats['count']}:{stats['fraction']:.3f}" for decoy, stats in
            sorted(flagged.items(), key=lambda kv: -kv[1]["fraction"])),
        "survivors": survivors,
        "non_survivors": non_survivors,
        "status": status,
        "absorbers": absorbers,
    }


SUMMARY_FIELDS = [
    "sample", "check", "verdict", "arms_used", "n_supporting", "n_evaluated", "n_measured",
    "n_unmeasured", "n_unmeasured_missing_from_arm_b", "n_unmeasured_unaligned",
    "n_unmeasured_missing_from_arm_a", "unmeasured_fraction", "max_unmeasured_fraction",
    "n_survive", "n_lost", "retention", "retention_if_unmeasured_lost", "retention_ci_lower",
    "retention_ci_upper", "bootstrap_iterations", "original_tier", "recomputed_tier",
    "tier_drop", "recomputed_clb_genes_detected", "recomputed_island_breadth_1x",
    "breadth_source", "non_megasynthase_genes_covered", "non_megasynthase_genes",
    "non_megasynthase_min_genes", "manual_followup_required", "flagged_absorbers",
]


def write_summary(result, path):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, SUMMARY_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        row = dict(result)
        row["retention"] = f"{result['retention']:.6f}"
        row["retention_if_unmeasured_lost"] = f"{result['retention_if_unmeasured_lost']:.6f}"
        row["unmeasured_fraction"] = f"{result['unmeasured_fraction']:.6f}"
        row["retention_ci_lower"] = f"{result['retention_ci_lower']:.6f}"
        row["retention_ci_upper"] = f"{result['retention_ci_upper']:.6f}"
        row["recomputed_island_breadth_1x"] = f"{result['recomputed_island_breadth_1x']:.6f}"
        writer.writerow(row)


def write_per_read(result, arm_a_rows, arm_b_rows, path):
    with Path(path).open("w", newline="") as handle:
        fields = ["read_id", "survived", "status", "arm_a_present", "arm_b_classification",
                  "absorbing_decoy"]
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        for read_id in sorted(result["survivors"]) + sorted(result["non_survivors"]):
            survived = read_id in result["survivors"]
            status = result["status"][read_id]
            arm_b_row = arm_b_rows.get(read_id)
            writer.writerow({
                "read_id": read_id,
                "survived": survived,
                "status": status,
                "arm_a_present": read_id in arm_a_rows,
                "arm_b_classification": (arm_b_row or {}).get("classification", "missing"),
                # Only a LOST read has an absorbing decoy; an unmeasured one has none.
                "absorbing_decoy": absorbing_decoy(arm_b_row) if status == LOST else "",
            })


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--read-evidence", required=True,
                   help="the pipeline's existing per-sample read_evidence.tsv "
                        "(classify_tumor_pks_evidence.py's output)")
    p.add_argument("--read-gene-assignment", required=True,
                   help="assign_reads_to_genes.py's read_id/Gene/overlap_bp TSV for "
                        "the original supporting read set")
    p.add_argument("--arm-a", default=None,
                   help="OPTIONAL KrakenUniq reclassification TSV: read_id, "
                        "assigned_clade, confidence, unique_kmer_count. Omitted, the "
                        "decision is Arm B alone and the summary says arms_used="
                        "arm_b_only rather than reporting a two-arm verdict that was "
                        "never computed.")
    p.add_argument("--arm-b", required=True,
                   help="competitive_realign_score.py's output TSV")
    p.add_argument("--arm-a-target-clade", default=None,
                   help="required with --arm-a: the exact assigned_clade value that "
                        "counts as a specific pks+ resolution (not collapsed to "
                        "genus/family)")
    p.add_argument("--arm-a-min-confidence", type=float, default=0.65)
    p.add_argument("--arm-a-min-unique-kmers", type=int, default=3)
    p.add_argument("--island-length", type=int, default=0,
                   help="island length in bp, for breadth recomputation "
                        "(matches classify_tumor_pks_evidence.py's --island-length)")
    p.add_argument("--survivor-depth", default=None,
                   help="optional: samtools-depth-shaped file (chrom, pos, depth) over "
                        "the island computed from survivor reads only, for an exact "
                        "breadth recomputation; approximated from read-gene overlap "
                        "when omitted (see module docstring)")
    p.add_argument("--multi-gene-reads", type=int, default=5)
    p.add_argument("--multi-gene-genes", type=int, default=3)
    p.add_argument("--multi-gene-breadth", type=float, default=.01)
    p.add_argument("--broad-island-reads", type=int, default=30)
    p.add_argument("--broad-island-genes", type=int, default=8)
    p.add_argument("--broad-island-breadth", type=float, default=.075)
    p.add_argument("--extensive-island-reads", type=int, default=100)
    p.add_argument("--extensive-island-genes", type=int, default=10)
    p.add_argument("--extensive-island-breadth", type=float, default=.15)
    p.add_argument("--retention-threshold", type=float, default=0.8,
                   help="placeholder pending Phase 3 calibration against real controls")
    p.add_argument("--retention-ci-lower-threshold", type=float, default=0.6,
                   help="placeholder pending Phase 3 calibration against real controls")
    p.add_argument("--bootstrap-iterations", type=int, default=1000)
    p.add_argument("--bootstrap-seed", type=int, default=0)
    p.add_argument("--max-tier-drop", type=int, default=1)
    p.add_argument("--non-megasynthase-min-genes", type=int, default=2)
    p.add_argument("--absorber-fraction-threshold", type=float, default=0.2)
    p.add_argument("--absorber-near-threshold-margin", type=float, default=5.0,
                   help="|delta_AS| (AS units) at or below which a dropped read's "
                        "call is treated as a close call rather than a landslide")
    p.add_argument("--min-supporting-reads", type=int, default=10)
    p.add_argument("--max-unmeasured-fraction", type=float, default=0.05,
                   help="above this share of evaluated reads with no Arm-B row, no "
                        "alignment to either panel side, or (with Arm A) no Arm-A row, the "
                        "verdict is INCOMPLETE_MEASUREMENT rather than a pass or fail "
                        "computed on a different read set (default: 0.05)")
    p.add_argument("--output-summary", required=True)
    p.add_argument("--output-per-read", required=True)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.arm_a and not args.arm_a_target_clade:
        raise SystemExit("--arm-a needs --arm-a-target-clade: without the exact clade "
                          "value that counts as a specific pks+ resolution, no Arm-A row "
                          "can confirm anything and every read would fail silently.")

    read_evidence = read_single_row_tsv(args.read_evidence)
    gene_assignment = read_gene_assignment(args.read_gene_assignment)
    arm_a_rows = read_arm_a(args.arm_a) if args.arm_a else None
    arm_b_rows = read_arm_b(args.arm_b)

    result = decide(args.sample, read_evidence, gene_assignment, arm_a_rows, arm_b_rows, args)
    write_summary(result, args.output_summary)
    write_per_read(result, arm_a_rows or {}, arm_b_rows, args.output_per_read)

    print(f"{args.sample}\t{result['verdict']} ({result['arms_used']})"
          f"\tR={result['retention']:.3f} over {result['n_measured']} measured"
          f" ({result['n_unmeasured']} unmeasured)"
          f" [{result['retention_ci_lower']:.3f},{result['retention_ci_upper']:.3f}]"
          f"\ttier={result['recomputed_tier']} (was {result['original_tier']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
