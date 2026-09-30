#!/usr/bin/env python3
"""Assertions for tests/check_v8_lanes_execution.sh.

Checks what the two lanes actually produced from real minimap2 alignments of reads whose
correct answer the fixture knows by construction:

  V8a  a read from the true island slice must be confirmed_true; a read from the
       NRPS/PKS look-alike must be confirmed_decoy. Both directions matter -- a lane that
       called everything decoy would satisfy the second on its own.
  V8b  a read from the Klebsiella carrier copy must be attributed to that carrier, and the
       sample-level call must name Klebsiella. This is the question V8a does not answer.

It also checks that neither lane's vocabulary appears in the other's output files, which
is Ludmil's "should not be conflated" made checkable on real output rather than on source.
"""
import argparse
import csv
import sys
from pathlib import Path


def read_rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_expected(path):
    return {row["key"]: row["value"] for row in read_rows(path)}


def check_v8a(arm_b_path, summary_path, per_read_path, expected, failures,
              min_fraction):
    arm_b = read_rows(arm_b_path)
    if not arm_b:
        failures.append(f"{arm_b_path}: no rows")
        return

    true_prefix = expected["true_read_prefix"]
    decoy_prefix = expected["decoy_read_prefix"]
    true_rows = [r for r in arm_b if r["read_id"].startswith(true_prefix)]
    decoy_rows = [r for r in arm_b if r["read_id"].startswith(decoy_prefix)]

    for label, rows, wanted, count_key in (
            ("island", true_rows, "confirmed_true", "true_reads"),
            ("homolog", decoy_rows, "confirmed_decoy", "decoy_reads")):
        if len(rows) != int(expected[count_key]):
            failures.append(f"V8a: {len(rows)} {label} reads in arm_b.tsv, "
                            f"expected {expected[count_key]}")
        correct = sum(1 for r in rows if r["classification"] == wanted)
        fraction = correct / len(rows) if rows else 0.0
        print(f"V8a {label:<8} {correct}/{len(rows)} {wanted} ({fraction:.0%})")
        if fraction < min_fraction:
            got = {}
            for row in rows:
                got[row["classification"]] = got.get(row["classification"], 0) + 1
            failures.append(
                f"V8a: only {fraction:.0%} of {label} reads were {wanted} "
                f"(needed {min_fraction:.0%}); actual: {got}")

    summary = read_rows(summary_path)
    if len(summary) != 1:
        failures.append(f"{summary_path}: expected one row, found {len(summary)}")
        return
    row = summary[0]
    print(f"V8a summary: verdict={row.get('verdict')} arms={row.get('arms_used')} "
          f"retention={row.get('retention')} tier {row.get('original_tier')} -> "
          f"{row.get('recomputed_tier')}")
    if row.get("check") != "locus_specificity":
        failures.append(f"V8a summary check column is {row.get('check')!r}, "
                        "so a reader cannot tell which of the two questions it answers")
    if row.get("arms_used") != "arm_b_only":
        failures.append(f"V8a ran without Arm A but reports arms_used="
                        f"{row.get('arms_used')!r}: a one-arm verdict must say so")
    # Half the supporting reads are homolog reads by construction, so retention must land
    # near 0.5 -- not 1.0 (nothing challenged) and not 0.0 (everything rejected).
    retention = float(row.get("retention") or 0)
    if not 0.2 <= retention <= 0.8:
        failures.append(f"V8a retention {retention:.3f} is outside the range this "
                        "half-true/half-homolog fixture must produce (0.2-0.8): the "
                        "challenge either did nothing or rejected everything")

    per_read = read_rows(per_read_path)
    if len(per_read) != len(arm_b):
        failures.append(f"V8a per-read table has {len(per_read)} rows, arm_b has "
                        f"{len(arm_b)}")
    for row in per_read:
        if "carrier" in ",".join(row.keys()).lower():
            failures.append("V8a per-read output names a carrier: the two checks are "
                            "being conflated in the output, not only in the prose")
            break


def check_v8b(summary_path, per_read_path, expected, failures, min_fraction):
    per_read = read_rows(per_read_path)
    if not per_read:
        failures.append(f"{per_read_path}: no rows")
        return
    wanted = expected["expected_carrier"]
    attributed = [r for r in per_read if r["attribution"] == "attributed"]
    correct = [r for r in attributed if r["best_carrier"] == wanted]
    fraction = len(correct) / len(per_read)
    print(f"V8b per-read: {len(correct)}/{len(per_read)} attributed to {wanted} "
          f"({fraction:.0%}); {len(attributed)} attributed in total")
    if fraction < min_fraction:
        got = {}
        for row in per_read:
            key = row["best_carrier"] or row["attribution"]
            got[key] = got.get(key, 0) + 1
        failures.append(f"V8b: only {fraction:.0%} of carrier reads were attributed to "
                        f"{wanted} (needed {min_fraction:.0%}); actual: {got}")

    summary = read_rows(summary_path)
    if len(summary) != 1:
        failures.append(f"{summary_path}: expected one row, found {len(summary)}")
        return
    row = summary[0]
    print(f"V8b summary: {row.get('carrier_attribution')} top_carrier="
          f"{row.get('top_carrier')} top_genus={row.get('top_genus')} "
          f"share={row.get('top_carrier_fraction')}")
    if row.get("top_carrier") != wanted:
        failures.append(f"V8b top carrier is {row.get('top_carrier')!r}, expected "
                        f"{wanted!r}: reads from one carrier were attributed to another")
    if row.get("top_genus") != expected["expected_genus"]:
        failures.append(f"V8b top genus is {row.get('top_genus')!r}, expected "
                        f"{expected['expected_genus']!r}")
    if row.get("carrier_attribution") not in ("CARRIER_ATTRIBUTED", "CARRIER_GENUS_ONLY"):
        failures.append(f"V8b call is {row.get('carrier_attribution')!r}: reads taken "
                        "straight from one carrier's own sequence should resolve at "
                        "least to its genus")
    for banned in ("SURVIVES", "NOT_DISTINGUISHABLE", "UNDERPOWERED", "confirmed_true"):
        if banned in Path(summary_path).read_text():
            failures.append(f"V8b output contains {banned!r}, a V8a verdict: the two "
                            "checks are being conflated in the output")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixtures", required=True)
    parser.add_argument("--v8a-arm-b", required=True)
    parser.add_argument("--v8a-summary", required=True)
    parser.add_argument("--v8a-per-read", required=True)
    parser.add_argument("--v8b-summary", required=True)
    parser.add_argument("--v8b-per-read", required=True)
    parser.add_argument("--min-fraction", type=float, default=0.9,
                        help="share of reads that must land on the known-correct answer")
    args = parser.parse_args(argv)

    expected = read_expected(Path(args.fixtures) / "expected.tsv")
    failures = []
    check_v8a(args.v8a_arm_b, args.v8a_summary, args.v8a_per_read, expected, failures,
              args.min_fraction)
    print()
    check_v8b(args.v8b_summary, args.v8b_per_read, expected, failures, args.min_fraction)

    if failures:
        print("\nFAILURES:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    print("\nOK: V8a separated island reads from homolog reads, and V8b attributed the "
          "carrier reads to the right organism.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
