#!/usr/bin/env python3
"""Classify a MAG bin's colibactin-island evidence from its own assembly, aligned
against the canonical IHE3034 reference -- not from HMM hits on its proteins.

M1: HMM E-value hits to the megasynthase domains (clbB/clbC/clbH/clbI/clbJ/clbK/clbN/
clbO) cross-react with unrelated NRPS/PKS gene clusters, so a bin can clear a
gene-count floor on domain homology alone with no island carriage
(v0.0.2_functional_test_20260910, ERR525841: 7 of 8 bins, including two
Bifidobacterium bins, called pks-positive this way). Rather than requiring specific
genes by name -- an arbitrary call -- the bin's own contigs are aligned to the
reference island and scored the same way read-level evidence already is: genes an
alignment actually covers, and breadth of the island those alignments span. Tier
names and default thresholds are deliberately the same as
classify_tumor_pks_evidence.py's read-level tiers (minus the read-count criterion,
which has no bin-level analogue), so "positive" means the same thing whether the
evidence is reads or an assembled bin.

    negative                  no alignment to the island at all
    localized_indeterminate  some alignment, fails multi_gene's criteria
    multi_gene                >=3  genes, >=1%   breadth
    broad_island               >=8  genes, >=7.5% breadth
    extensive_island           >=10 genes, >=15%  breadth
"""
import argparse
import csv
from pathlib import Path

TIERS = ("extensive_island", "broad_island", "multi_gene")


def paf_hits(path, start, end, min_identity, min_mapq):
    """Alignment blocks overlapping [start, end), identity- and MAPQ-filtered.

    Mirrors summarize_tumor_pks_contigs.py's paf(): every alignment counts, not one
    per contig -- the island is repeat-rich and a bin's assembly is as fragmented as
    a targeted one, so collapsing to a single best hit per contig would discard real
    reference coverage.
    """
    hits = []
    text = Path(path).read_text() if Path(path).stat().st_size else ""
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        alen = int(fields[10])
        ident = int(fields[9]) / alen if alen else 0
        ts, te, mapq = int(fields[7]), int(fields[8]), int(fields[11])
        if te <= start or ts >= end or ident < min_identity or mapq < min_mapq:
            continue
        hits.append((max(start, ts), min(end, te)))
    return hits


def union_length(spans):
    total, right = 0, None
    for left, stop in sorted(spans):
        if right is None or left > right:
            total += stop - left
            right = stop
        elif stop > right:
            total += stop - right
            right = stop
    return total


def genes_covered(gff_path, contig, start, end, hits):
    """Distinct clb genes (GFF `Name=`/`gene=`) whose interval overlaps an
    alignment block. A gene within the window but with no aligned base over it does
    not count -- unlike summarize_tumor_pks_contigs.py's genes(), which is a display
    helper and does not filter on alignment overlap.
    """
    covered = set()
    for line in Path(gff_path).read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[0] != contig or fields[2] != "gene":
            continue
        left, right = int(fields[3]) - 1, int(fields[4])
        if right <= start or left >= end:
            continue
        attrs = dict(v.split("=", 1) for v in fields[8].split(";") if "=" in v)
        name = attrs.get("Name", attrs.get("gene", ""))
        if not name.startswith("clb"):
            continue
        if any(hit_start < right and hit_end > left for hit_start, hit_end in hits):
            covered.add(name)
    return covered


def classify(gene_count, breadth, thresholds):
    """Return the highest tier whose gene and breadth criteria both hold."""
    for tier in TIERS:
        min_genes, min_breadth = thresholds[tier]
        if gene_count >= min_genes and breadth >= min_breadth:
            return tier
    return "localized_indeterminate" if gene_count > 0 or breadth > 0 else "negative"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--bin-id", required=True)
    p.add_argument("--paf", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--contig", default="NC_017628.1")
    p.add_argument("--island-start", type=int, required=True)
    p.add_argument("--island-end", type=int, required=True)
    p.add_argument("--min-identity", type=float, default=0.90)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--multi-gene-genes", type=int, default=3)
    p.add_argument("--multi-gene-breadth", type=float, default=.01)
    p.add_argument("--broad-island-genes", type=int, default=8)
    p.add_argument("--broad-island-breadth", type=float, default=.075)
    p.add_argument("--extensive-island-genes", type=int, default=10)
    p.add_argument("--extensive-island-breadth", type=float, default=.15)
    p.add_argument("--output", required=True)
    a = p.parse_args()

    thresholds = {
        "multi_gene": (a.multi_gene_genes, a.multi_gene_breadth),
        "broad_island": (a.broad_island_genes, a.broad_island_breadth),
        "extensive_island": (a.extensive_island_genes, a.extensive_island_breadth),
    }
    for lower, higher in (("multi_gene", "broad_island"), ("broad_island", "extensive_island")):
        if any(h < l for l, h in zip(thresholds[lower], thresholds[higher])):
            raise SystemExit(f"locus tier thresholds must be non-decreasing: {lower} exceeds {higher}")

    length = a.island_end - a.island_start
    hits = paf_hits(a.paf, a.island_start, a.island_end, a.min_identity, a.min_mapq)
    breadth = union_length(hits) / length
    genes = genes_covered(a.gff, a.contig, a.island_start, a.island_end, hits)
    tier = classify(len(genes), breadth, thresholds)

    fields = ["sample", "bin_id", "locus_tier", "locus_genes_detected", "locus_genes", "locus_breadth"]
    with Path(a.output).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(dict(
            sample=a.sample, bin_id=a.bin_id, locus_tier=tier,
            locus_genes_detected=len(genes), locus_genes=",".join(sorted(genes)),
            locus_breadth=f"{breadth:.6f}",
        ))


if __name__ == "__main__":
    main()
