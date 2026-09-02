#!/usr/bin/env python3
"""Validate the exact clbA-clbS annotation and featureCounts outputs."""

import argparse
import csv
import math
import sys


EXPECTED = {f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"}


def validate_genes(genes, source):
    if len(genes) != len(set(genes)):
        raise ValueError(f"{source}: duplicate clb gene identifiers")
    observed = set(genes)
    if observed != EXPECTED or len(genes) != 19:
        raise ValueError(
            f"{source}: expected exactly clbA-clbS; "
            f"missing={sorted(EXPECTED - observed)}, extra={sorted(observed - EXPECTED)}"
        )


def annotation_genes(path):
    genes = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                raise ValueError(f"{path}: malformed or non-gene annotation row")
            attributes = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            genes.append(attributes.get("Name", ""))
    validate_genes(genes, path)


def counts_genes(path):
    with open(path, newline="") as handle:
        rows = list(csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    if len(rows) != 20 or not rows or rows[0][:6] != ["Geneid", "Chr", "Start", "End", "Strand", "Length"]:
        raise ValueError(f"{path}: expected a featureCounts header and exactly 19 data rows")
    genes = [row[0] for row in rows[1:]]
    validate_genes(genes, path)
    for row in rows[1:]:
        try:
            value = float(row[-1])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{path}: invalid count for {row[0]}") from exc
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{path}: count for {row[0]} must be finite and non-negative")


def summary(path):
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if len(rows) < 2 or not rows[0] or rows[0][0] != "Status":
        raise ValueError(f"{path}: missing or malformed featureCounts summary")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation")
    parser.add_argument("--counts")
    parser.add_argument("--summary")
    args = parser.parse_args()
    if args.annotation:
        annotation_genes(args.annotation)
    if args.counts:
        counts_genes(args.counts)
    if args.summary:
        summary(args.summary)
    if not any((args.annotation, args.counts, args.summary)):
        parser.error("provide at least one file to validate")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
