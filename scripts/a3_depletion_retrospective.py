#!/usr/bin/env python3
"""Tier A 3, retrospective arm: how much did each depletion stage actually remove?

The simulation benchmark (a3_depletion_floor_benchmark.py) answers what the removal
floor does to reads of KNOWN origin. This answers the companion question on real runs:
across a finished cohort, what fraction of reads did each depletion stage discard, and
is the SECOND stage still removing enough for the compounding to matter?

It needs no new computation. Every finished run already writes the per-stage survivor
counts into results/cohort/pks.master_summary.tsv, one row per sample:

    filter_input_reads -> reads_after_fastp -> reads_after_hg38 -> reads_after_t2t_phix

so each stage's removal rate is already on disk for every sample ever run.

WHY IT MATTERS. The floor is applied in series -- the human pass, then the T2T+phiX
pass, then once per pangenome index -- and each pass re-exposes every surviving read to
the same unbenchmarked rule. If the T2T+phiX pass is still removing a non-trivial
fraction of what the human pass left, the compounding is real rather than theoretical,
and the simulation's projection can be checked against it.

A KNOWN LIMIT. reads_after_pangenome is written to the per-sample depletion QC but is
not carried into the master summary, so the optional pangenome stages cannot be
measured from this file. With --pangenome_db unset -- the usual case -- there are only
the two stages below, and this covers both.

PRIVACY. This reads a table DERIVED FROM PROTECTED SAMPLES, so it is written for Ammal
to run, not for Claude. It emits only aggregate statistics -- counts, means and
quantiles of a RATIO -- and never a per-sample row, a sample identifier, or a read.
The --per-sample flag exists for local inspection and its output must stay on the
cluster: it carries identifiers and belongs in no report, deck, commit or email.

    python3 scripts/a3_depletion_retrospective.py \
        /path/to/results/cohort/pks.master_summary.tsv [more.tsv ...] \
        --out validation/a3_depletion_retrospective.tsv
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

# Consecutive counters in pipeline order. Each pair is (count reaching the stage,
# count surviving it, the stage's name).
INTERVALS = [
    ("reads_after_fastp", "reads_after_hg38", "hg38_human_pass"),
    ("reads_after_hg38", "reads_after_t2t_phix", "t2t_phix_pass"),
]
STAGE_NAMES = [stage for _, _, stage in INTERVALS]
COUNTER_COLUMNS = {key for interval in INTERVALS for key in interval[:2]}
SAMPLE_COLUMNS = ("sample", "sampleID", "Sample", "sample_id")


def as_count(value):
    """A non-negative integer, or None for NA / blank / anything unparseable."""
    if value is None:
        return None
    text = value.strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", "."}:
        return None
    try:
        number = int(float(text))
    except ValueError:
        return None
    return number if number >= 0 else None


def stage_rates(row):
    """Fraction removed BY each stage, relative to what reached it."""
    rates = {}
    for before_key, after_key, stage in INTERVALS:
        before = as_count(row.get(before_key))
        after = as_count(row.get(after_key))
        if before is None or after is None or before <= 0:
            continue
        if after > before:
            # A counter that grew across a depletion stage means the two columns did
            # not come from the same run. Skip rather than report a negative rate.
            continue
        rates[stage] = (before - after) / before
    return rates


def read_summaries(paths):
    """Every row of every master summary given, plus how many were unusable."""
    rows, unreadable = [], 0
    for path in paths:
        target = Path(path)
        if target.is_dir():
            found = sorted(target.rglob("pks.master_summary.tsv"))
            if not found:
                print(f"  no pks.master_summary.tsv under {path}", file=sys.stderr)
            rows.extend(r for f in found for r in read_summaries([f])[0])
            continue
        try:
            with target.open(encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if not reader.fieldnames:
                    unreadable += 1
                    continue
                if not (COUNTER_COLUMNS & set(reader.fieldnames)):
                    print(f"  no depletion counters in {target.name}", file=sys.stderr)
                    unreadable += 1
                    continue
                rows.extend(dict(r) for r in reader)
        except OSError:
            unreadable += 1
    return rows, unreadable


def summarise(stage, rates):
    ordered = sorted(rates)
    def at(q):
        return ordered[min(len(ordered) - 1, int(q * (len(ordered) - 1)))]
    return {
        "stage": stage,
        "samples": len(ordered),
        "mean_fraction_removed": f"{statistics.fmean(ordered):.6f}",
        "median_fraction_removed": f"{statistics.median(ordered):.6f}",
        "min_fraction_removed": f"{ordered[0]:.6f}",
        "p90_fraction_removed": f"{at(0.90):.6f}",
        "max_fraction_removed": f"{ordered[-1]:.6f}",
    }


def sample_name(row):
    for key in SAMPLE_COLUMNS:
        if row.get(key):
            return row[key]
    return "unnamed"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("summaries", nargs="+",
                        help="pks.master_summary.tsv files, or directories to search")
    parser.add_argument("--out", default="a3_depletion_retrospective.tsv")
    parser.add_argument("--per-sample", default=None,
                        help="also write a per-sample table here. It carries sample "
                             "identifiers, so it stays on the cluster and goes into no "
                             "report, deck, commit or email.")
    args = parser.parse_args(argv)

    rows, unreadable = read_summaries(args.summaries)
    if not rows:
        sys.exit("no usable master summary rows found")

    gathered = {stage: [] for stage in STAGE_NAMES}
    per_sample_rows = []
    for row in rows:
        rates = stage_rates(row)
        if not rates:
            continue
        for stage, rate in rates.items():
            gathered[stage].append(rate)
        if args.per_sample:
            entry = {"sample": sample_name(row)}
            entry.update({s: f"{r:.6f}" for s, r in rates.items()})
            per_sample_rows.append(entry)

    results = [summarise(stage, rates) for stage, rates in gathered.items() if rates]
    if not results:
        sys.exit("no sample had a usable pair of consecutive stage counts")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        handle.write("# Tier A 3 retrospective: fraction of reads removed BY each host-\n")
        handle.write("# depletion stage, relative to what reached that stage.\n")
        handle.write("# Aggregate only -- no per-sample row, no identifier.\n")
        handle.write(f"# summary rows read: {len(rows)}\n")
        writer = csv.DictWriter(handle, list(results[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    if args.per_sample and per_sample_rows:
        with open(args.per_sample, "w", newline="", encoding="utf-8") as handle:
            handle.write("# PER-SAMPLE, carries identifiers. Keep on the cluster.\n")
            writer = csv.DictWriter(handle, ["sample"] + STAGE_NAMES,
                                    delimiter="\t", restval="")
            writer.writeheader()
            writer.writerows(per_sample_rows)

    print(f"  summary rows read: {len(rows)}"
          + (f"  ({unreadable} file(s) unusable)" if unreadable else ""))
    print()
    print(f"  {'stage':<22}{'n':>6}{'mean':>10}{'median':>10}{'p90':>10}{'max':>10}")
    print("  " + "-" * 68)
    for row in results:
        print(f"  {row['stage']:<22}{row['samples']:>6}"
              f"{float(row['mean_fraction_removed']):>10.4f}"
              f"{float(row['median_fraction_removed']):>10.4f}"
              f"{float(row['p90_fraction_removed']):>10.4f}"
              f"{float(row['max_fraction_removed']):>10.4f}")
    print()

    second = next((r for r in results if r["stage"] == "t2t_phix_pass"), None)
    if second and float(second["mean_fraction_removed"]) > 0.001:
        print("  The T2T+phiX pass is still removing a measurable fraction of what the")
        print("  human pass left, so the unbenchmarked floor is applied to the same read")
        print("  more than once and the false-removal rate compounds. This is the")
        print("  real-data counterpart to the simulation benchmark.")
    elif second:
        print("  The T2T+phiX pass removes almost nothing beyond the human pass, so in")
        print("  practice the floor is effectively applied once and the compounding is")
        print("  negligible for these cohorts.")
    print(f"\n  wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
