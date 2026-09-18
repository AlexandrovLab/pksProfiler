#!/usr/bin/env python3
"""Classify a tumour sample's read-level pks evidence into the frozen tier system.

Tiers come from tier_threshold_development_20260909. All three criteria (reads, clb
genes at >=1 assigned read, canonical-island breadth at >=1x) must hold for a tier;
classification proceeds from the highest tier downward.

    negative                 0 pks reads
    localized_indeterminate  >=1 read, fails one or more multi_gene criteria
    multi_gene               >=5   reads, >=3  genes, >=1%  breadth
    broad_island             >=30  reads, >=8  genes, >=7.5% breadth
    extensive_island         >=100 reads, >=10 genes, >=15%  breadth

localized_indeterminate is NOT a formal positive tier. These data are positive-only and
cannot establish the boundary between a true positive and background mapping.
"""
import argparse, csv
from pathlib import Path

GENES = {f"clb{x}" for x in "ABCDEFGHIJKLMNOPQRS"}

# Highest tier first; the cascade returns the first tier whose criteria all hold.
TIERS = ("extensive_island", "broad_island", "multi_gene")

def qc_metrics(path):
    with Path(path).open() as handle:
        return {row["Metric"]: row["Value"] for row in csv.DictReader(handle, delimiter="\t")}

def detected_genes(path):
    detected = 0
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("Geneid"): continue
            fields = line.rstrip().split("\t")
            if fields and fields[0] in GENES and float(fields[-1] or 0) >= 1: detected += 1
    return detected

def breadths(path, expected_length):
    covered = [0, 0, 0]; observed = 0
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) != 3: continue
            depth = int(fields[2]); observed += 1
            covered[0] += depth >= 1; covered[1] += depth >= 2; covered[2] += depth >= 3
    if observed != expected_length: raise ValueError(f"expected {expected_length} depth positions, observed {observed}")
    return [value / expected_length for value in covered]

def classify(reads, genes, breadth, thresholds):
    """Return the highest tier whose read, gene and breadth criteria all hold."""
    if reads == 0: return "negative"
    for tier in TIERS:
        min_reads, min_genes, min_breadth = thresholds[tier]
        if reads >= min_reads and genes >= min_genes and breadth >= min_breadth:
            return tier
    return "localized_indeterminate"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True); p.add_argument("--qc", required=True)
    p.add_argument("--counts", required=True); p.add_argument("--depth", required=True)
    p.add_argument("--island-length", required=True, type=int); p.add_argument("--output", required=True)
    p.add_argument("--multi-gene-reads", type=int, default=5)
    p.add_argument("--multi-gene-genes", type=int, default=3)
    p.add_argument("--multi-gene-breadth", type=float, default=.01)
    p.add_argument("--broad-island-reads", type=int, default=30)
    p.add_argument("--broad-island-genes", type=int, default=8)
    p.add_argument("--broad-island-breadth", type=float, default=.075)
    p.add_argument("--extensive-island-reads", type=int, default=100)
    p.add_argument("--extensive-island-genes", type=int, default=10)
    p.add_argument("--extensive-island-breadth", type=float, default=.15)
    a = p.parse_args()

    thresholds = {
        "multi_gene":       (a.multi_gene_reads, a.multi_gene_genes, a.multi_gene_breadth),
        "broad_island":     (a.broad_island_reads, a.broad_island_genes, a.broad_island_breadth),
        "extensive_island": (a.extensive_island_reads, a.extensive_island_genes, a.extensive_island_breadth),
    }
    for lower, higher in (("multi_gene", "broad_island"), ("broad_island", "extensive_island")):
        if any(h < l for l, h in zip(thresholds[lower], thresholds[higher])):
            raise SystemExit(f"tier thresholds must be non-decreasing: {lower} exceeds {higher}")

    metrics = qc_metrics(a.qc)
    reads = int(float(metrics.get("reads_clb_genes_align", 0))); genes = detected_genes(a.counts)
    b1, b2, b3 = breadths(a.depth, a.island_length)
    evidence = classify(reads, genes, b1, thresholds)

    fields = ["sample", "read_evidence", "pks_reads", "clb_genes_detected", "island_breadth_1x", "island_breadth_2x", "island_breadth_3x"]
    with Path(a.output).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader()
        writer.writerow(dict(sample=a.sample, read_evidence=evidence, pks_reads=reads, clb_genes_detected=genes,
            island_breadth_1x=f"{b1:.6f}", island_breadth_2x=f"{b2:.6f}", island_breadth_3x=f"{b3:.6f}"))

if __name__ == "__main__": main()
