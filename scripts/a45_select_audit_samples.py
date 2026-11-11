#!/usr/bin/env python3
"""Pick the samples worth spending mate-loss audit decode passes on.

The audit (a45_mate_loss_audit.py) costs roughly two decode passes per alignment, so it
should start where a lost mate could actually change something:

    every positive            extensive_island / broad_island / multi_gene -- the calls
                              Ludmil's items 4-5 could have inflated-by-absence
    every near-miss           localized_indeterminate: >=1 pks read but below the
                              multi_gene floor. Recovery only ever ADDS reads, so these
                              are the samples a recovered mate could flip to positive --
                              the sharpest test of whether the defect matters
    a few clean negatives     0 pks reads, randomly drawn: the noise floor

This script joins the cohort's master summary to the run's sample sheet
(patient,alignment) and writes a ready-to-run command file of ~/mateaudit calls. It
reads a table derived from protected samples and its OUTPUT contains sample identifiers
and alignment paths, so it is for Ammal to run and the command file stays on the
cluster.

    python3 scripts/a45_select_audit_samples.py \
        --master  /path/results/cohort/pks.master_summary.tsv \
        --sheet   /path/sample_sheet.csv \
        --low-score-mode on \
        --commands validation/a45_audit_commands.sh
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

POSITIVE_TIERS = {"extensive_island", "broad_island", "multi_gene"}
NEAR_MISS = "localized_indeterminate"
SAMPLE_COLUMNS = ("sample", "sampleID", "Sample", "patient")


def sample_of(row):
    for key in SAMPLE_COLUMNS:
        if row.get(key):
            return row[key].strip()
    return None


def read_master(path):
    with open(path, encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("#")), delimiter="\t")
        if not reader.fieldnames or "read_evidence" not in reader.fieldnames:
            sys.exit(f"{path} has no read_evidence column; is it a master summary?")
        return [dict(row) for row in reader]


def read_sheet(path):
    """patient -> alignment from the v0.0.1-style sample sheet (CSV)."""
    mapping = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "alignment" not in reader.fieldnames:
            sys.exit(f"{path} has no alignment column; expected patient,alignment")
        for row in reader:
            name = sample_of(row)
            if name and row.get("alignment"):
                mapping[name] = row["alignment"].strip()
    return mapping


def select(rows, clean_negatives, rng):
    chosen, reasons = [], {"positive": 0, "near_miss": 0, "clean_negative": 0}
    negatives = []
    for row in rows:
        name = sample_of(row)
        if not name:
            continue
        evidence = (row.get("read_evidence") or "").strip()
        if evidence in POSITIVE_TIERS:
            chosen.append((name, "positive"))
            reasons["positive"] += 1
        elif evidence == NEAR_MISS:
            chosen.append((name, "near_miss"))
            reasons["near_miss"] += 1
        elif evidence == "negative":
            negatives.append(name)
    for name in rng.sample(negatives, min(clean_negatives, len(negatives))):
        chosen.append((name, "clean_negative"))
        reasons["clean_negative"] += 1
    return chosen, reasons, len(negatives)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--master", required=True, help="pks.master_summary.tsv")
    parser.add_argument("--sheet", required=True,
                        help="the run's sample sheet (patient,alignment CSV)")
    parser.add_argument("--low-score-mode", choices=["on", "off"], required=True,
                        help="the cohort's extraction mode, as C1 / provenance report it")
    parser.add_argument("--cram-reference", default=None,
                        help="appended to every command; needed for CRAM cohorts")
    parser.add_argument("--clean-negatives", type=int, default=10,
                        help="how many 0-read samples to include as controls")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--commands", required=True,
                        help="output: a shell file of ~/mateaudit calls. Contains "
                             "identifiers and paths -- it stays on the cluster.")
    args = parser.parse_args(argv)

    rows = read_master(args.master)
    alignments = read_sheet(args.sheet)
    chosen, reasons, negative_pool = select(
        rows, args.clean_negatives, random.Random(args.seed))
    if not chosen:
        sys.exit("nothing selected: no positive, near-miss or negative rows found")

    missing = [name for name, _ in chosen if name not in alignments]
    reference = f" '{args.cram_reference}'" if args.cram_reference else ""

    out = Path(args.commands)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        handle.write("#!/bin/bash\n")
        handle.write("# Mate-loss audit set. CONTAINS SAMPLE IDENTIFIERS AND PATHS --\n")
        handle.write("# this file stays on the cluster. Run lines directly, or submit\n")
        handle.write("# each as a job (see the wrapper in ~/mateaudit).\n")
        handle.write("set -euo pipefail\n\n")
        for name, reason in chosen:
            alignment = alignments.get(name)
            if not alignment:
                handle.write(f"# SKIPPED ({reason}): {name} not in the sample sheet\n")
                continue
            handle.write(f"# {reason}\n")
            handle.write(f"bash ~/mateaudit '{alignment}' "
                         f"{args.low_score_mode}{reference}\n")

    # Console output is aggregate-only: counts, never names.
    print(f"  selected {len(chosen) - len(missing)} samples"
          f" ({reasons['positive']} positive, {reasons['near_miss']} near-miss,"
          f" {reasons['clean_negative']} of {negative_pool} clean negatives)")
    if missing:
        print(f"  {len(missing)} selected sample(s) had no alignment in the sheet;"
              " marked SKIPPED in the command file")
    print(f"  wrote {out} -- identifiers inside, keep it on the cluster")
    return 0


if __name__ == "__main__":
    sys.exit(main())
