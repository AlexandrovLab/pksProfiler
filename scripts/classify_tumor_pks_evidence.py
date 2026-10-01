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
                                   rna_expression_breadth_{1,2,3}x, rna_library_type,
                                   clb_transcripts_detected, rna_residual_human_fraction,
                                   rna_residual_checked_reads           -- never tiered

RNA QC (memo section 4; docs/running/wes_rna.md):

    clb_transcripts_detected     clb genes (clbA-clbS) with >= --rna-transcript-min-reads
                                 (default 3) featureCounts-assigned RNA reads. Gene-level:
                                 the IHE3034 annotation defines no operon or transcription
                                 unit, so each gene is its own unit. rna_clb_genes_detected
                                 is the same count at >= 1 read.
    rna_residual_human_fraction  rna_residual_human_reads / rna_residual_checked_reads from
                                 --residual-human-qc (rna_residual_human.py); NA without it
                                 or when nothing was left to check.

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
# V2 robustness (decision A, 1 Oct): the tier stays on raw reads (read_evidence, pks_reads).
# Beside it, the same evidence counted three other ways by
# compute_fragment_support_variants.py, and the tier each count would imply with the same
# genes and breadth rule: fragments (mates of one DNA fragment counted once), reads with
# duplicates (0x400 from samtools markdup) removed, and both. tier_robustness is stable when
# all three imply read_evidence, unstable otherwise. NA where there is no tier.
ROBUSTNESS_FIELDS = ["pks_fragments", "tier_if_fragments", "pks_reads_dedup", "tier_if_dedup",
                     "pks_fragments_dedup", "tier_if_fragments_dedup", "tier_robustness"]
# variant name in compute_fragment_support_variants.py -> (count column, tier column)
ROBUSTNESS_VARIANTS = {"fragment": ("pks_fragments", "tier_if_fragments"),
                       "dup_marked": ("pks_reads_dedup", "tier_if_dedup"),
                       "dup_collapsed": ("pks_fragments_dedup", "tier_if_fragments_dedup")}
SCHEMAS = {
    "wgs": READ_EVIDENCE_FIELDS + ROBUSTNESS_FIELDS,
    "metagenome": METAGENOME_FIELDS,
    "wes": ["sample", "mode", "clb_wes_reads", "wes_clb_genes_detected",
            "wes_island_breadth_1x", "wes_island_breadth_2x", "wes_island_breadth_3x"],
    "rna": ["sample", "mode", "clb_rna_reads", "rna_clb_genes_detected",
            "rna_expression_breadth_1x", "rna_expression_breadth_2x",
            "rna_expression_breadth_3x", "rna_library_type",
            "clb_transcripts_detected", "rna_residual_human_fraction",
            "rna_residual_checked_reads"],
}
# The file each mode writes, as main.nf's report manifest kind and the published name.
EVIDENCE_KIND = {"wgs": "read_evidence", "metagenome": "read_evidence",
                 "wes": "wes_evidence", "rna": "rna_evidence"}

def qc_metrics(path):
    with Path(path).open() as handle:
        return {row["Metric"]: row["Value"] for row in csv.DictReader(handle, delimiter="\t")}

def detected_genes(path, min_reads=1):
    """clb genes with >= min_reads assigned reads in a featureCounts table."""
    detected = 0
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("Geneid"): continue
            fields = line.rstrip().split("\t")
            if fields and fields[0] in GENES and float(fields[-1] or 0) >= min_reads: detected += 1
    return detected


def residual_human(path):
    """(fraction text, checked reads) from rnaResidualHuman's fragment; NA when absent."""
    if not path:
        return "NA", "NA"
    metrics = qc_metrics(path)
    checked = int(metrics["rna_residual_checked_reads"])
    human = int(metrics["rna_residual_human_reads"])
    if human > checked:
        raise SystemExit(f"{path}: {human} residual-human reads exceed {checked} checked")
    return (f"{human / checked:.6f}" if checked else "NA"), checked

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
    p.add_argument("--rna-transcript-min-reads", type=int, default=3,
                   help="RNA reads a clb gene needs to count in clb_transcripts_detected")
    p.add_argument("--residual-human-qc", default=None,
                   help="rnaResidualHuman's QC fragment (rna_residual_human.py); --mode rna")
    p.add_argument("--depletion-qc", default=None,
                   help="mapReads' depletion QC fragment; the metagenome reads-per-million "
                        "denominator (reads_after_pangenome, else reads_after_t2t_phix)")
    p.add_argument("--fragment-variants", default=None,
                   help="compute_fragment_support_variants.py's TSV for this sample (V2): fills "
                        "pks_fragments / pks_reads_dedup / pks_fragments_dedup and tier_if_*")
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
    if a.rna_transcript_min_reads < 1:
        raise SystemExit("--rna-transcript-min-reads must be at least 1")
    if a.mode == "metagenome" and not a.depletion_qc:
        raise SystemExit("--mode metagenome needs --depletion-qc for pks_reads_per_million")
    depleted = depleted_reads(a.depletion_qc) if a.depletion_qc else None
    row = evidence_row(a.mode, a.sample, reads, genes, (b1, b2, b3), thresholds,
                       no_tiers=a.no_tiers, rna_library_type=a.rna_library_type,
                       depleted=depleted,
                       transcripts=(detected_genes(a.counts, a.rna_transcript_min_reads)
                                    if a.mode == "rna" else "NA"),
                       residual=(residual_human(a.residual_human_qc)
                                 if a.mode == "rna" else ("NA", "NA")))

    if a.mode == "wgs":
        row.update(robustness(row, a.fragment_variants))

    with Path(a.output).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, SCHEMAS[a.mode], delimiter="\t"); writer.writeheader()
        writer.writerow(row)
    if a.qc_fragment:
        with Path(a.qc_fragment).open("w", newline="") as handle:
            handle.write("Sample\tMetric\tValue\n")
            handle.write(f"{a.sample}\tevidence_outcome\t{evidence_outcome(row)}\n")
            if row.get("tier_robustness", "NA") != "NA":
                for field in ROBUSTNESS_FIELDS:
                    handle.write(f"{a.sample}\t{field}\t{row[field]}\n")


def read_variants(path):
    """{variant: (count, tier)} from compute_fragment_support_variants.py's report."""
    variants = {}
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 6 and fields[0] != "sample" and not fields[0].startswith("#"):
                variants[fields[1]] = (int(fields[2]), fields[5])
    return variants


def robustness(row, variants_path):
    """
    The V2 columns for a wgs row. The tier itself is never touched: read_evidence was
    already decided from raw reads. Untiered: NA. A tiered row with no reads has nothing to
    count another way (0, negative, stable); with reads but no variants file (a hand run
    without one), NA.
    """
    tier = row.get("read_evidence")
    if tier in (None, "untiered"):
        return {field: "NA" for field in ROBUSTNESS_FIELDS}
    if not variants_path:
        if int(row["pks_reads"]) == 0:
            out = {}
            for count_field, tier_field in ROBUSTNESS_VARIANTS.values():
                out.update({count_field: 0, tier_field: "negative"})
            return dict(out, tier_robustness="stable")
        return {field: "NA" for field in ROBUSTNESS_FIELDS}
    variants = read_variants(variants_path)
    raw = variants.get("raw")
    # One computation of the raw tier, not two that could drift: the variants script reads
    # the same counts, QC and island depth, so its raw row must be this row.
    if raw != (int(row["pks_reads"]), tier):
        raise SystemExit(f"{variants_path}: raw variant {raw} disagrees with read_evidence "
                         f"({row['pks_reads']}, {tier}); the two must read the same inputs")
    out = {}
    for variant, (count_field, tier_field) in ROBUSTNESS_VARIANTS.items():
        if variant not in variants:
            raise SystemExit(f"{variants_path}: no {variant} row")
        out[count_field], out[tier_field] = variants[variant]
    tiers = {out[tier_field] for _count, tier_field in ROBUSTNESS_VARIANTS.values()}
    out["tier_robustness"] = "stable" if tiers == {tier} else "unstable"
    return out


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
                 rna_library_type="unknown", depleted=None, transcripts="NA",
                 residual=("NA", "NA")):
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
        # transcripts: clb genes at >= --rna-transcript-min-reads; residual: (fraction,
        # checked reads) from residual_human(). NA when not supplied.
        return dict(sample=sample, mode=mode, clb_rna_reads=reads, rna_clb_genes_detected=genes,
                    rna_expression_breadth_1x=b1, rna_expression_breadth_2x=b2,
                    rna_expression_breadth_3x=b3, rna_library_type=rna_library_type,
                    clb_transcripts_detected=transcripts,
                    rna_residual_human_fraction=residual[0],
                    rna_residual_checked_reads=residual[1])
    raise ValueError(f"unknown mode {mode!r}; expected one of {', '.join(MODES)}")

if __name__ == "__main__": main()
