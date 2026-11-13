#!/usr/bin/env python3
"""A2: how much of the cohort's positivity rests on genes seen once?

rc2 review item A2 (Ludmil, 15 Oct):

    One-read gene floor: full-cohort sensitivity at >=1, >=2, >=3 reads per gene; for every
    positive report how many qualifying genes rest on exactly one read; focus on the +20
    newly positive.

Two outcomes are reported separately at floor 2, not merged: a sample that stops being
positive at all, and one that drops a tier but stays positive. Merging them overstates how
much a cohort rests on single-read genes.

`clb_genes_detected` counts a clb gene as present at >= 1 assigned read
(scripts/classify_tumor_pks_evidence.py), and the gene count is one of the three criteria
every tier requires. So a gene carried by a single read can be what lifts a sample over a
tier boundary. This re-derives each sample's gene count and tier at floors of 1, 2 and 3
reads per gene, using the SAME tier rules, and reports how the cohort moves.

PRIVACY: this reads a per-sample count matrix, which is derived from protected data, and
writes ONLY pooled aggregates -- counts of samples, never a sample identifier, never a row.
Run it yourself; it is written so its output can be pasted into a report.

    python3 scripts/a2_gene_floor_sensitivity.py \\
        --gene-counts results/cohort/gene_counts/pks.gene.counts.align.txt \\
        --master      results/cohort/pks.master_summary.tsv \\
        --out         a2_gene_floor.tsv

With --baseline-master <an older run's master summary>, the samples positive now but not
in that run are reported as their own group: the "newly positive" he asked to focus on.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# scripts/classify_tumor_pks_evidence.py: all three criteria must hold.
# Breadth is a FRACTION in the evidence files, not a percentage: classify_tumor_pks_evidence.py
# defaults are .15 / .075 / .01 and island_breadth_1x is written with "%.6f". Using 15/7.5/1
# here silently made every sample fail the breadth criterion and reported 0 positives on a
# cohort with ~120; sanity_check() below now refuses that rather than printing it.
TIERS = [("extensive_island", 100, 10, 0.15),
         ("broad_island", 30, 8, 0.075),
         ("multi_gene", 5, 3, 0.01)]
POSITIVE = {t[0] for t in TIERS}
# TIERS is strongest-first, so a larger index is a weaker call.
TIER_RANK = {name: i for i, (name, *_ ) in enumerate(TIERS)}
FLOORS = (1, 2, 3)


def tier_for(reads, genes, breadth):
    for name, min_reads, min_genes, min_breadth in TIERS:
        if reads >= min_reads and genes >= min_genes and breadth >= min_breadth:
            return name
    return "localized_indeterminate" if reads >= 1 else "negative"


def num(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def read_matrix(path):
    """{sample: [per-gene counts]} from a gene x sample matrix (genes in rows)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip() and not l.startswith("#")]
    if not rows:
        sys.exit(f"{path}: empty")
    samples = rows[0][1:]
    per_sample = {s: [] for s in samples}
    for row in rows[1:]:
        for sample, value in zip(samples, row[1:]):
            per_sample[sample].append(num(value))
    return per_sample


def read_master(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-counts", required=True,
                    help="cohort/gene_counts/pks.gene.counts.align.txt")
    ap.add_argument("--master", required=True, help="cohort/pks.master_summary.tsv")
    ap.add_argument("--baseline-master", default=None,
                    help="an earlier run's master summary; samples positive now and not "
                         "then are reported separately as the newly positive")
    ap.add_argument("--out", default="a2_gene_floor.tsv")
    a = ap.parse_args()

    counts = read_matrix(a.gene_counts)
    master = read_master(a.master)
    baseline = read_master(a.baseline_master) if a.baseline_master else None

    shared = [s for s in master if s in counts]
    if not shared:
        sys.exit("no samples in common between the matrix and the master summary")

    newly = set()
    if baseline is not None:
        for s in shared:
            was = (baseline.get(s, {}).get("read_evidence") or "").strip()
            now = (master[s].get("read_evidence") or "").strip()
            if now in POSITIVE and was not in POSITIVE:
                newly.add(s)

    # The tiers this run already recorded, recomputed independently below. If the two
    # disagree wildly the re-derivation is wrong (a units or column mismatch), and a wrong
    # answer that looks plausible is worse than no answer -- so stop and say so.
    recorded = sum(1 for s in shared
                   if (master[s].get("read_evidence") or "").strip() in POSITIVE)
    rederived = sum(1 for s in shared
                    if tier_for(num(master[s].get("pks_reads")),
                                sum(1 for v in counts[s] if v >= 1),
                                num(master[s].get("island_breadth_1x"))) in POSITIVE)
    if recorded and abs(rederived - recorded) > max(2, 0.05 * recorded):
        sys.exit(f"re-derivation disagrees with the recorded tiers: {rederived} positives at "
                 f"floor 1 versus {recorded} recorded. Check column names and breadth units "
                 f"(breadth is a fraction, not a percent) before trusting any number here.")
    if not recorded:
        print("NOTE: no recorded positives in this cohort; nothing to cross-check against.",
              file=sys.stderr)

    groups = [("all samples", shared)]
    if baseline is not None:
        groups.append((f"newly positive (vs baseline, n={len(newly)})", sorted(newly)))

    out_rows, notes = [], []
    for label, members in groups:
        if not members:
            notes.append(f"{label}: no samples; nothing to report")
            continue
        for floor in FLOORS:
            tally, lost_positivity, demoted, genes_on_one = {}, 0, 0, 0
            for s in members:
                row = master[s]
                reads = num(row.get("pks_reads"))
                breadth = num(row.get("island_breadth_1x"))
                vals = counts[s]
                genes_at_floor = sum(1 for v in vals if v >= floor)
                tier = tier_for(reads, genes_at_floor, breadth)
                tally[tier] = tally.get(tier, 0) + 1
                if floor == 1:
                    ones = sum(1 for v in vals if v == 1)
                    genes_on_one += ones
                    # What does dropping its single-read genes actually cost this sample?
                    #
                    # These are two different things and were previously counted as one
                    # (Ludmil's v3 review, item 8). A sample moving
                    # extensive_island -> broad_island is still positive and still a
                    # confident call; a sample moving multi_gene -> localized_indeterminate
                    # has stopped being positive at all. Reporting them together overstates
                    # how much the cohort rests on single-read genes. A3 already keeps this
                    # distinction; A2 now does too.
                    if tier in POSITIVE and ones:
                        at_two = tier_for(reads, sum(1 for v in vals if v >= 2), breadth)
                        if at_two not in POSITIVE:
                            lost_positivity += 1
                        elif TIER_RANK[at_two] > TIER_RANK[tier]:
                            demoted += 1
            positives = sum(n for t, n in tally.items() if t in POSITIVE)
            for tier in ["extensive_island", "broad_island", "multi_gene",
                         "localized_indeterminate", "negative"]:
                out_rows.append([label, floor, tier, tally.get(tier, 0),
                                 f"{100*tally.get(tier,0)/len(members):.1f}"])
            out_rows.append([label, floor, "ANY_POSITIVE", positives,
                             f"{100*positives/len(members):.1f}"])
            if floor == 1:
                out_rows.append([label, 1, "genes_resting_on_exactly_one_read", genes_on_one, ""])
                out_rows.append([label, 1, "positives_that_LOSE_POSITIVITY_at_floor_2",
                                 lost_positivity, f"{100*lost_positivity/max(1,positives):.1f}"])
                out_rows.append([label, 1, "positives_DEMOTED_but_still_positive_at_floor_2",
                                 demoted, f"{100*demoted/max(1,positives):.1f}"])

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["# A2 gene-floor sensitivity. Pooled counts only, no sample identifiers."])
        w.writerow(["# Tier rules are classify_tumor_pks_evidence.py's; only the gene floor varies."])
        w.writerow(["group", "gene_floor_reads", "metric", "samples", "percent_of_group"])
        w.writerows(out_rows)

    print(f"samples analysed   {len(shared)}")
    if baseline is not None:
        print(f"newly positive     {len(newly)}")
    print()
    for label, floor, metric, n, pct in out_rows:
        if metric in ("ANY_POSITIVE", "genes_resting_on_exactly_one_read",
                      "positives_that_LOSE_POSITIVITY_at_floor_2",
                      "positives_DEMOTED_but_still_positive_at_floor_2"):
            print(f"  {label:<34} floor>={floor}  {metric:<42} {n:>6} {pct:>6}")
    for n in notes:
        print(f"  NOTE {n}")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
