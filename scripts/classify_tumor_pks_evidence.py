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

Ludmil, 26 Sep mode memo (sections 3, 4 and 6): the tiers are a WGS call only, and every
mode gets its own schema so that no field name implies genomic presence for WES off-target
DNA or for RNA. --mode picks the schema; every file carries `mode`:

    wgs         read_evidence.tsv  tiered (untiered with --no-tiers, e.g. long-read WGS)
    metagenome  read_evidence.tsv  same columns, read_evidence=untiered
    wes         wes_evidence.tsv   sample, mode, clb_wes_reads, wes_clb_genes_detected,
                                   wes_island_breadth_{1,2,3}x          -- never tiered
    rna         rna_evidence.tsv   sample, mode, clb_rna_reads, rna_clb_genes_detected,
                                   rna_expression_breadth_{1,2,3}x, rna_library_type
                                                                        -- never tiered

The measurements are the same three in every mode (clb-assigned reads, clb genes with >=1
read, island breadth at >=1/2/3x); only what they may be read as differs. RNA gene breadth
is expression breadth, not island completeness; a WES or RNA negative is not absence.
"""
import argparse, csv
from pathlib import Path

GENES = {f"clb{x}" for x in "ABCDEFGHIJKLMNOPQRS"}

# Highest tier first; the cascade returns the first tier whose criteria all hold.
TIERS = ("extensive_island", "broad_island", "multi_gene")

MODES = ("wgs", "metagenome", "wes", "rna")
RNA_LIBRARY_TYPES = ("polyA", "total", "unknown")

# Column order per mode. wgs and metagenome share read_evidence.tsv's schema (metagenome
# rows read untiered); wes and rna have their own, with no tier column at all.
READ_EVIDENCE_FIELDS = ["sample", "mode", "read_evidence", "pks_reads", "clb_genes_detected",
                        "island_breadth_1x", "island_breadth_2x", "island_breadth_3x"]
# Ludmil, 26 Sep mode memo section 5: for metagenomes the primary quantity is pks reads per
# million host-depleted reads (stool cohorts differ greatly in depth); the raw count stays.
METAGENOME_FIELDS = READ_EVIDENCE_FIELDS + ["host_depleted_reads", "pks_reads_per_million"]
SCHEMAS = {
    "wgs": READ_EVIDENCE_FIELDS,
    "metagenome": METAGENOME_FIELDS,
    "wes": ["sample", "mode", "clb_wes_reads", "wes_clb_genes_detected",
            "wes_island_breadth_1x", "wes_island_breadth_2x", "wes_island_breadth_3x"],
    "rna": ["sample", "mode", "clb_rna_reads", "rna_clb_genes_detected",
            "rna_expression_breadth_1x", "rna_expression_breadth_2x",
            "rna_expression_breadth_3x", "rna_library_type"],
}
# The file each mode writes, as main.nf's report manifest kind and the published name.
EVIDENCE_KIND = {"wgs": "read_evidence", "metagenome": "read_evidence",
                 "wes": "wes_evidence", "rna": "rna_evidence"}

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
    # Ludmil, 26 Sep mode memo: the tiers were derived from TCGA tumour WGS (short-read,
    # Bowtie2) and are not a metagenome, WES or RNA call. --no-tiers keeps the evidence and
    # writes read_evidence=untiered for wgs too (long-read WGS, whose profiling aligner the
    # thresholds were never calibrated on).
    p.add_argument("--mode", choices=MODES, default="wgs")
    p.add_argument("--no-tiers", action="store_true")
    p.add_argument("--rna-library-type", choices=RNA_LIBRARY_TYPES, default="unknown",
                   help="recorded in rna_evidence.tsv; --mode rna only")
    p.add_argument("--depletion-qc", default=None,
                   help="mapReads' depletion QC fragment; the metagenome reads-per-million "
                        "denominator (reads_after_pangenome, else reads_after_t2t_phix)")
    p.add_argument("--qc-fragment", default=None,
                   help="also write a Sample/Metric/Value QC fragment carrying "
                        "evidence_outcome (positive|negative|untiered) for build_qc_summary.py")
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
    if a.mode == "metagenome" and not a.depletion_qc:
        raise SystemExit("--mode metagenome needs --depletion-qc for pks_reads_per_million")
    depleted = depleted_reads(a.depletion_qc) if a.depletion_qc else None
    row = evidence_row(a.mode, a.sample, reads, genes, (b1, b2, b3), thresholds,
                       no_tiers=a.no_tiers, rna_library_type=a.rna_library_type,
                       depleted=depleted)

    with Path(a.output).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, SCHEMAS[a.mode], delimiter="\t"); writer.writeheader()
        writer.writerow(row)
    if a.qc_fragment:
        with Path(a.qc_fragment).open("w", newline="") as handle:
            handle.write("Sample\tMetric\tValue\n")
            handle.write(f"{a.sample}\tevidence_outcome\t{evidence_outcome(row)}\n")


def depleted_reads(path):
    """The final host-depleted read count from a depletion QC fragment: after the
    pangenome pass when it ran, otherwise after T2T + phiX."""
    metrics = qc_metrics(path)
    for metric in ("reads_after_pangenome", "reads_after_t2t_phix"):
        value = metrics.get(metric)
        if value not in (None, "", "NA"):
            return int(value)
    raise SystemExit(f"{path} has neither reads_after_pangenome nor reads_after_t2t_phix")


def reads_per_million(reads, depleted):
    """pks reads per million host-depleted reads; NA when nothing survived depletion."""
    if depleted is None or depleted <= 0:
        return "NA"
    return f"{reads * 1e6 / depleted:.4f}"


def evidence_outcome(row):
    """positive | negative | untiered, for the QC terminal state (GAMEPLAN rc1).

    Only a tiered wgs row is ever positive or negative: positive is a formal tier
    (multi_gene, broad_island, extensive_island); negative is 0 reads or
    localized_indeterminate, which is explicitly not a positive tier. Every other row --
    metagenome, wes, rna, long-read wgs -- is evidence without a call.
    """
    tier = row.get("read_evidence")
    if row.get("mode") != "wgs" or tier in (None, "untiered"):
        return "untiered"
    return "positive" if tier in TIERS else "negative"


def evidence_row(mode, sample, reads, genes, breadth, thresholds, no_tiers=False,
                 rna_library_type="unknown", depleted=None):
    """One output row in `mode`'s schema. Only wgs (without no_tiers) is ever tiered."""
    b1, b2, b3 = (f"{value:.6f}" for value in breadth)
    if mode in ("wgs", "metagenome"):
        tiered = mode == "wgs" and not no_tiers
        row = dict(sample=sample, mode=mode,
                   read_evidence=classify(reads, genes, breadth[0], thresholds) if tiered else "untiered",
                   pks_reads=reads, clb_genes_detected=genes,
                   island_breadth_1x=b1, island_breadth_2x=b2, island_breadth_3x=b3)
        if mode == "metagenome":
            row.update(host_depleted_reads=depleted if depleted is not None else "NA",
                       pks_reads_per_million=reads_per_million(reads, depleted))
        return row
    if mode == "wes":
        return dict(sample=sample, mode=mode, clb_wes_reads=reads, wes_clb_genes_detected=genes,
                    wes_island_breadth_1x=b1, wes_island_breadth_2x=b2, wes_island_breadth_3x=b3)
    if mode == "rna":
        return dict(sample=sample, mode=mode, clb_rna_reads=reads, rna_clb_genes_detected=genes,
                    rna_expression_breadth_1x=b1, rna_expression_breadth_2x=b2,
                    rna_expression_breadth_3x=b3, rna_library_type=rna_library_type)
    raise ValueError(f"unknown mode {mode!r}; expected one of {', '.join(MODES)}")

if __name__ == "__main__": main()
