#!/usr/bin/env python3
"""A3: how many positives survive fragment collapse, and duplicate collapse?

rc2 review item A3 (Ludmil, 15 Oct):

    Raw vs fragment/duplicate tiers: reconcile the denominators (100 raw positives vs 43 of
    168), report among the current positives how many change after fragment collapse and
    after duplicate collapse, and whether instability is enriched in the newly positive;
    then choose explicitly.

Every WGS run already writes, per sample, the tier under each counting rule:
`read_evidence` (raw reads, the primary tier), `tier_if_fragments`, `tier_if_dedup` and
`tier_if_fragments_dedup`, plus `tier_robustness` (stable / unstable). Nothing new has to be
computed from the data -- this reconciles what is already there and reports the movement as
pooled counts.

THE DENOMINATOR POINT, which is the part he actually flagged
------------------------------------------------------------
"100 raw positives" and "43 of 168" are not two answers to one question: they count over
different sets. This prints, explicitly, the denominator of every figure it reports -- how
many samples carry each column at all -- so the two numbers can be compared on the same
footing instead of being quoted side by side.

PRIVACY: reads the cohort master summary, which is derived from protected data, and writes
ONLY pooled counts -- no sample identifier, no row. Run it yourself.

    python3 scripts/a3_tier_collapse_reconciliation.py \\
        --master results/cohort/pks.master_summary.tsv --out a3_tier_collapse.tsv

With --baseline-master <an earlier run's summary>, samples positive now but not then are
reported as their own group, which answers "is instability enriched in the newly positive".
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

POSITIVE = {"multi_gene", "broad_island", "extensive_island"}
ORDER = ["extensive_island", "broad_island", "multi_gene",
         "localized_indeterminate", "negative"]
RULES = [("raw reads (primary)", "read_evidence"),
         ("fragment collapse", "tier_if_fragments"),
         ("duplicate collapse", "tier_if_dedup"),
         ("fragment + duplicate", "tier_if_fragments_dedup")]


def val(row, key):
    v = (row.get(key) or "").strip()
    return v if v and v != "NA" else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--baseline-master", default=None)
    ap.add_argument("--out", default="a3_tier_collapse.tsv")
    a = ap.parse_args()

    with open(a.master, encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"{a.master}: no rows")

    present = {key for _, key in RULES if any(val(r, key) for r in rows)}
    if "read_evidence" not in present:
        sys.exit("no read_evidence column: this is not a tiered (WGS) cohort summary")
    missing = [label for label, key in RULES if key not in present]

    newly = set()
    if a.baseline_master:
        with open(a.baseline_master, encoding="utf-8", errors="replace") as fh:
            base = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
        for r in rows:
            if (val(r, "read_evidence") in POSITIVE
                    and val(base.get(r["sample"], {}), "read_evidence") not in POSITIVE):
                newly.add(r["sample"])

    groups = [("all samples", rows)]
    if a.baseline_master:
        groups.append((f"newly positive (n={len(newly)})",
                       [r for r in rows if r["sample"] in newly]))

    out = []
    for label, members in groups:
        if not members:
            continue
        # Each rule reported over ITS OWN denominator: the samples that carry that column.
        for rule, key in RULES:
            if key not in present:
                continue
            scored = [r for r in members if val(r, key)]
            pos = sum(1 for r in scored if val(r, key) in POSITIVE)
            out.append([label, rule, "samples_scored", len(scored), ""])
            out.append([label, rule, "positive", pos,
                        f"{100*pos/len(scored):.1f}" if scored else ""])
            for tier in ORDER:
                n = sum(1 for r in scored if val(r, key) == tier)
                out.append([label, rule, f"tier_{tier}", n,
                            f"{100*n/len(scored):.1f}" if scored else ""])

        # Movement, counted only over samples scored under BOTH rules being compared.
        raw_pos = [r for r in members if val(r, "read_evidence") in POSITIVE]
        for rule, key in RULES[1:]:
            if key not in present:
                continue
            comparable = [r for r in raw_pos if val(r, key)]
            lost = sum(1 for r in comparable if val(r, key) not in POSITIVE)
            demoted = sum(1 for r in comparable
                          if val(r, key) in POSITIVE
                          and ORDER.index(val(r, key)) > ORDER.index(val(r, "read_evidence")))
            out.append([label, rule, "raw_positives_comparable", len(comparable), ""])
            out.append([label, rule, "raw_positives_that_become_negative", lost,
                        f"{100*lost/len(comparable):.1f}" if comparable else ""])
            out.append([label, rule, "raw_positives_demoted_a_tier", demoted,
                        f"{100*demoted/len(comparable):.1f}" if comparable else ""])

        flagged = [r for r in members if val(r, "tier_robustness")]
        if flagged:
            unstable = sum(1 for r in flagged if val(r, "tier_robustness") == "unstable")
            out.append([label, "pipeline's own flag", "samples_with_tier_robustness", len(flagged), ""])
            out.append([label, "pipeline's own flag", "unstable", unstable,
                        f"{100*unstable/len(flagged):.1f}"])

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["# A3 tier reconciliation across counting rules. Pooled counts only."])
        w.writerow(["# Every figure is reported over its own denominator (samples_scored),"])
        w.writerow(["# which is the reconciliation the review item asked for."])
        w.writerow(["group", "counting_rule", "metric", "samples", "percent_of_denominator"])
        w.writerows(out)

    for label, rule, metric, n, pct in out:
        if metric in ("samples_scored", "positive", "raw_positives_comparable",
                      "raw_positives_that_become_negative", "raw_positives_demoted_a_tier",
                      "unstable"):
            print(f"  {label:<26} {rule:<22} {metric:<36} {n:>6} {pct:>6}")
    for m in missing:
        print(f"  NOTE  '{m}' not present in this cohort; its rows are omitted rather than zeroed")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
