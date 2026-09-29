#!/usr/bin/env python3
"""Collapse `samtools idxstats` into a handful of input-reference categories.

E2 (Ludmil, 25 Sep): extraction ran idxstats and kept only the grand total, which
throws away where the input's reads actually sit -- on the primary assembly, on alt
or unplaced contigs, on decoys, on viral sequences, or unplaced altogether. Those
counts are the context the FLAG4-only extraction question (E1) needs, for a few
bytes per sample. They are QC/provenance only: they do not measure pks-read loss.

Categories, by sequence name (with or without a `chr` prefix):

    canonical   1-22, X, Y
    mito        M / MT
    alt         *_alt
    unplaced    *_random, Un_*, GL*/KI*/JH* scaffolds
    decoy       *_decoy, hs37d5, HLA-*
    viral       EBV, NC_007605*, HPV*/HBV*/HCV*/HIV*/KSHV* and similar names
    other       anything else

For every category the mapped and placed-unmapped counts (idxstats columns 3 and 4)
are summed. The `*` row -- unmapped reads with no position -- is reported alone as
`input_unmapped_unplaced`. `input_alignment_records` stays the grand total.
"""

import argparse
import re
import sys

CATEGORIES = ("canonical", "mito", "alt", "unplaced", "decoy", "viral", "other")

_CANONICAL = re.compile(r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
_MITO = re.compile(r"^(?:chr)?(?:M|MT)$")
_VIRAL = re.compile(r"^(?:chr)?(?:EBV|NC_007605|HPV|HBV|HCV|HIV|KSHV|HTLV|CMV|HHV|MCPyV|SV40)", re.I)
_UNPLACED = re.compile(r"(?:_random$|^(?:chr)?Un_|^(?:GL|KI|JH)\d)")


def category(name):
    if _CANONICAL.match(name):
        return "canonical"
    if _MITO.match(name):
        return "mito"
    if name.endswith("_alt"):
        return "alt"
    if name.endswith("_decoy") or name == "hs37d5" or name.startswith("HLA-"):
        return "decoy"
    if _VIRAL.match(name):
        return "viral"
    if _UNPLACED.search(name):
        return "unplaced"
    return "other"


def summarize(lines):
    mapped = dict.fromkeys(CATEGORIES, 0)
    placed_unmapped = dict.fromkeys(CATEGORIES, 0)
    unplaced_unmapped = 0
    for number, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 4 or not fields[2].isdigit() or not fields[3].isdigit():
            raise ValueError(f"idxstats line {number} is not name/length/mapped/unmapped")
        name, m, u = fields[0], int(fields[2]), int(fields[3])
        if name == "*":
            if m:
                raise ValueError("idxstats '*' row reports mapped reads")
            unplaced_unmapped += u
            continue
        c = category(name)
        mapped[c] += m
        placed_unmapped[c] += u
    metrics = {"input_alignment_records": sum(mapped.values()) + sum(placed_unmapped.values())
               + unplaced_unmapped}
    for c in CATEGORIES:
        metrics[f"input_mapped_{c}"] = mapped[c]
    metrics["input_unmapped_placed"] = sum(placed_unmapped.values())
    metrics["input_unmapped_unplaced"] = unplaced_unmapped
    return metrics


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--idxstats", required=True)
    p.add_argument("--sample", required=True)
    a = p.parse_args()
    with open(a.idxstats) as handle:
        metrics = summarize(handle)
    for metric, value in metrics.items():
        sys.stdout.write(f"{a.sample}\t{metric}\t{value}\n")


if __name__ == "__main__":
    main()
