#!/usr/bin/env python3
"""The one definition of which clb gene a read belongs to.

F14. Three places decided this differently: featureCounts with --largestOverlap, the
alignment QC with a bare `bedtools intersect -u` that counted any read *touching* any
clb gene, and the taxonomy lane with its own awk that summed overlaps and broke ties
by gene name. The three numbers were not comparable, and nothing said so.

The rule, stated once:

  * a read is assigned to the gene it shares the most aligned bases with
  * overlaps from several alignment blocks of the same read are summed per gene
  * ties go to the alphabetically first gene, so the answer does not depend on input
    order
  * an overlap below --min-overlap does not count
  * a read is assigned to at most one gene, so counts sum to the number of reads
  * strand is ignored: the island is profiled unstranded, as featureCounts is run

Input is `bedtools intersect -a <reads.bed> -b <genes.bed> -wo` on stdin. Output is
read_id, Gene, overlap_bp.

    ... | assign_reads_to_genes.py --output read_clb_gene.tsv [--min-overlap 1]
"""
import argparse
import sys
from collections import defaultdict


def assign(rows, min_overlap=1):
    """rows: iterable of (read_id, gene, overlap_bp). Returns {read_id: (gene, bp)}."""
    totals = defaultdict(int)
    for read_id, gene, overlap in rows:
        totals[(read_id, gene)] += overlap

    best = {}
    for (read_id, gene), overlap in totals.items():
        if overlap < min_overlap:
            continue
        current = best.get(read_id)
        if current is None or overlap > current[1] or (overlap == current[1] and gene < current[0]):
            best[read_id] = (gene, overlap)
    return best


def parse(stream, read_name_field, gene_field, overlap_field):
    for line in stream:
        fields = line.rstrip("\n").split("\t")
        if len(fields) <= max(read_name_field, gene_field, overlap_field):
            continue
        try:
            overlap = int(fields[overlap_field])
        except ValueError:
            continue
        yield fields[read_name_field], fields[gene_field], overlap


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-overlap", type=int, default=1,
                        help="aligned bases a read must share with a gene to count")
    # bedtools intersect -wo on a 6-column read BED and a 6-column gene BED puts the
    # read name at 3, the gene name at 9 and the overlap last.
    parser.add_argument("--read-field", type=int, default=3)
    parser.add_argument("--gene-field", type=int, default=9)
    parser.add_argument("--overlap-field", type=int, default=-1)
    args = parser.parse_args()

    rows = []
    for line in sys.stdin:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        try:
            overlap = int(fields[args.overlap_field])
        except (ValueError, IndexError):
            continue
        if len(fields) <= max(args.read_field, args.gene_field):
            continue
        rows.append((fields[args.read_field], fields[args.gene_field], overlap))

    best = assign(rows, min_overlap=args.min_overlap)
    with open(args.output, "w") as handle:
        handle.write("read_id\tGene\toverlap_bp\n")
        for read_id in sorted(best):
            gene, overlap = best[read_id]
            handle.write(f"{read_id}\t{gene}\t{overlap}\n")

    print(len(best))
    return 0


if __name__ == "__main__":
    sys.exit(main())
