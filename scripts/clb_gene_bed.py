#!/usr/bin/env python3
"""The canonical clbA-clbS interval set, in BED, from the annotation.

F14. Two modules built this themselves with near-identical awk, and they did not
agree: the alignment QC filtered to `Name=clb[A-S]`, while the taxonomy lane took
every feature with a Name -- so a read over a neighbouring gene was assigned to it
there and ignored here. One definition, one implementation.

    clb_gene_bed.py --annotation clb.gff --output clb_genes.bed

BED is 0-based half-open; GFF is 1-based inclusive, so start is start-1.
"""
import argparse
import re
import sys
from pathlib import Path

CLB_GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]
NAME = re.compile(r"(?:^|;)Name=([^;]+)")


def intervals(annotation):
    found = {}
    for line in Path(annotation).read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        match = NAME.search(fields[8])
        if not match:
            continue
        name = match.group(1).strip()
        if name not in CLB_GENES:
            continue
        if name in found:
            raise SystemExit(f"[ERROR] {annotation}: {name} appears more than once")
        found[name] = (fields[0], int(fields[3]) - 1, int(fields[4]), name, fields[6])
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--allow-partial", action="store_true",
                        help="do not fail when the annotation lacks some clb genes")
    args = parser.parse_args()

    found = intervals(args.annotation)
    missing = [gene for gene in CLB_GENES if gene not in found]
    if missing and not args.allow_partial:
        raise SystemExit(f"[ERROR] {args.annotation} is missing {len(missing)} clb genes: "
                         f"{', '.join(missing)}")

    with open(args.output, "w") as handle:
        for gene in CLB_GENES:
            if gene in found:
                contig, start, end, name, strand = found[gene]
                handle.write(f"{contig}\t{start}\t{end}\t{name}\t0\t{strand}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
