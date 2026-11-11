#!/usr/bin/env python3
"""v3 item 18: compare the metagenome rate between --prefilter_mode off and balanced.

Reads both runs' read_evidence.tsv and reports, side by side:

    host_depleted_reads     the DENOMINATOR, taken at host depletion upstream of the
                            prefilter -- must be identical, and this asserts it
    pks_reads               the NUMERATOR, counted on PROFILING_READS, which `balanced`
                            narrows to CANDIDATE_OUT.reads
    pks_reads_per_million   the rate that follows from the two

A difference in the numerator is the finding item 18 asks for, not a failure: it is what
decides whether the old "the rate does not move with --prefilter_mode" comment was
merely imprecise or actually wrong on real behaviour. This exits non-zero only when the
comparison cannot support a conclusion:

  * a run is missing its evidence file, or the sample is absent from it
  * the denominator moved, which would contradict the pipeline's own design
  * neither mode found any pks read -- 0 vs 0 is not parity, it is no measurement

Inputs are synthetic reads over public references; the output is counts.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

EVIDENCE = "read_evidence.tsv"


def find_evidence(results, sample):
    """The sample's read_evidence.tsv under a results tree."""
    root = Path(results)
    if not root.is_dir():
        sys.exit(f"not a results directory: {results}")
    for candidate in sorted(root.rglob(EVIDENCE)):
        with candidate.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("sample") == sample:
                    return candidate, row
    sys.exit(f"no {EVIDENCE} row for sample {sample} under {results}")


def as_int(value, label):
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        sys.exit(f"{label} is not an integer: {value!r}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--off", required=True, help="results dir of the off run")
    parser.add_argument("--balanced", required=True, help="results dir of the balanced run")
    parser.add_argument("--out", default=None, help="write the comparison here as TSV")
    args = parser.parse_args(argv)

    off_path, off = find_evidence(args.off, args.sample)
    bal_path, bal = find_evidence(args.balanced, args.sample)

    off_denom = as_int(off.get("host_depleted_reads"), "off host_depleted_reads")
    bal_denom = as_int(bal.get("host_depleted_reads"), "balanced host_depleted_reads")
    off_reads = as_int(off.get("pks_reads"), "off pks_reads")
    bal_reads = as_int(bal.get("pks_reads"), "balanced pks_reads")
    off_rate = off.get("pks_reads_per_million", "NA")
    bal_rate = bal.get("pks_reads_per_million", "NA")

    rows = [
        ("host_depleted_reads", off_denom, bal_denom),
        ("pks_reads", off_reads, bal_reads),
        ("pks_reads_per_million", off_rate, bal_rate),
        ("host_depleted_stage", off.get("host_depleted_stage", "NA"),
         bal.get("host_depleted_stage", "NA")),
    ]
    width = max(len(name) for name, _, _ in rows)
    print(f"  {'metric':<{width}}  {'off':>16}  {'balanced':>16}")
    print("  " + "-" * (width + 38))
    for name, a, b in rows:
        flag = "" if str(a) == str(b) else "   <-- differs"
        print(f"  {name:<{width}}  {str(a):>16}  {str(b):>16}{flag}")
    print()

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="", encoding="utf-8") as handle:
            handle.write("# v3 item 18: metagenome rate, --prefilter_mode off vs balanced.\n")
            handle.write("# One synthetic sample, same fixture, one flag different.\n")
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["metric", "prefilter_off", "prefilter_balanced"])
            writer.writerows(rows)
        print(f"  wrote {out}")

    problems = []
    if off_denom != bal_denom:
        problems.append(
            f"the DENOMINATOR moved: host_depleted_reads {off_denom} (off) vs "
            f"{bal_denom} (balanced). It is taken at host depletion, upstream of the "
            "prefilter, so this contradicts the pipeline's design and is a real defect, "
            "not a measurement.")
    if off_reads == 0 and bal_reads == 0:
        problems.append(
            "neither mode found a pks read, so this run measures nothing. 0 vs 0 must "
            "not be reported as parity: use a fixture whose reads carry clb genes.")
    if problems:
        print()
        for problem in problems:
            print(f"  FAIL: {problem}", file=sys.stderr)
        return 1

    print(f"  Denominator identical ({off_denom}), as the design requires: it is taken")
    print("  at host depletion, upstream of the prefilter.")
    if off_reads == bal_reads:
        print(f"  Numerator identical too ({off_reads}): on this fixture the balanced")
        print("  prefilter, with the DIAMOND rescue main.nf requires, loses no clb read,")
        print("  so the rate is unchanged. That is one fixture, not a proof for every")
        print("  library -- the comment should still not claim mode-invariance outright.")
    else:
        lost = off_reads - bal_reads
        print(f"  Numerator MOVED: {off_reads} (off) vs {bal_reads} (balanced), a "
              f"difference of {lost}.")
        print("  The prefilter changes what reaches profiling, so the rate is not")
        print("  comparable across modes. This is exactly what v3 item 18 describes, and")
        print("  it confirms the old main.nf comment was wrong rather than merely loose.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
