#!/usr/bin/env python3
"""Divergence from IHE3034 of the reads supporting the clb gene counts (A5).

Ludmil, 25 Sep A5: MAPQ >= 40 within IHE3034 says a read placed confidently *on this one
reference*; it says nothing about uniqueness across all microbes. The edit distance (NM)
of those reads can quantify how far they sit from IHE3034, and nothing more -- it is a
divergence metric relative to IHE3034, not evidence of organism identity or of microbial
specificity. Specificity is the competitive-mapping analysis (pks_specificity.nf), and is
stated relative to the panel it tested.

Reads measured: primary, mapped alignments at MAPQ >= --min-mapq that featureCounts'
--largestOverlap rule would assign to exactly one clb gene -- the reads behind
reads_clb_genes_align. A read overlapping two genes by the same number of bases is left
out, as featureCounts leaves it unassigned. Secondary and supplementary records are not
reads and are skipped. This only reads the alignment; no count, tier or call changes.

Per read:
    alignment columns  M/=/X + I + D CIGAR bases
    divergence         NM / alignment columns (0 = identical to IHE3034)
    read length        M/=/X + I + S + H CIGAR bases (the full read, clipping included)

Output 1, --qc-out: a Sample/Metric/Value fragment for build_qc_summary.py.

    divergence_ihe3034_reads          reads measured
    nm_ihe3034_median                 median NM
    divergence_ihe3034_median/_q1/_q3 median and quartiles of divergence
    divergence_ihe3034_frac_gt_2pct   fraction of reads with divergence > 2%
    divergence_ihe3034_frac_gt_5pct   fraction of reads with divergence > 5%
    clb_read_length_median/_q1/_q3/_min/_max

With no reads measured only divergence_ihe3034_reads (0) is written; the rest are NA in
the cohort table rather than an invented 0.

Output 2, --genes-out: one row per clb gene -- reads, median NM, median divergence.

Quartiles are linear interpolation between order statistics (numpy's default). The 2%
and 5% thresholds are compared exactly (NM * 100 > 2 * columns), not in floating point.

    samtools view SAMPLE.sorted.bam | ihe3034_divergence_qc.py --sample SAMPLE \\
        --annotation clb_genes.gff --qc-out SAMPLE.ihe3034_divergence.qc.tsv \\
        --genes-out SAMPLE.ihe3034_divergence_by_gene.tsv
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

CLB_GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
FLAG_SKIP = 0x4 | 0x100 | 0x800          # unmapped, secondary, supplementary
THRESHOLDS_PCT = (2, 5)

QC_METRICS = [
    "divergence_ihe3034_reads",
    "nm_ihe3034_median",
    "divergence_ihe3034_median",
    "divergence_ihe3034_q1",
    "divergence_ihe3034_q3",
    "divergence_ihe3034_frac_gt_2pct",
    "divergence_ihe3034_frac_gt_5pct",
    "clb_read_length_median",
    "clb_read_length_q1",
    "clb_read_length_q3",
    "clb_read_length_min",
    "clb_read_length_max",
]
GENE_FIELDS = ["Sample", "gene", "reads", "nm_ihe3034_median", "divergence_ihe3034_median"]


def load_genes(path):
    """clbA-clbS intervals from the GFF pksProfiler_align.nf builds: {contig: [(gene, s, e)]}."""
    genes = defaultdict(list)
    seen = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            name = attrs.get("Name", "")
            if name in CLB_GENES:
                genes[fields[0]].append((name, int(fields[3]), int(fields[4])))
                seen.add(name)
    missing = sorted(set(CLB_GENES) - seen)
    if missing:
        raise ValueError(f"{path}: annotation is missing {', '.join(missing)}")
    return genes


def parse_cigar(cigar):
    ops = CIGAR_RE.findall(cigar)
    if not ops or "".join(n + op for n, op in ops) != cigar:
        raise ValueError(f"malformed CIGAR {cigar!r}")
    return [(int(n), op) for n, op in ops]


def read_metrics(cigar):
    """(alignment columns, read length, reference span) of one CIGAR."""
    columns = length = span = 0
    for n, op in parse_cigar(cigar):
        if op in "M=X":
            columns += n; length += n; span += n
        elif op == "I":
            columns += n; length += n
        elif op == "D":
            columns += n; span += n
        elif op == "N":
            span += n
        elif op in "SH":
            length += n
    return columns, length, span


def assign_gene(genes, contig, start, end):
    """featureCounts --largestOverlap: the one gene with the most overlapping bases."""
    best, best_overlap, tied = None, 0, False
    for name, g_start, g_end in genes.get(contig, ()):
        overlap = min(end, g_end) - max(start, g_start) + 1
        if overlap <= 0:
            continue
        if overlap > best_overlap:
            best, best_overlap, tied = name, overlap, False
        elif overlap == best_overlap:
            tied = True
    return None if tied else best


def nm_tag(fields):
    for tag in fields[11:]:
        if tag.startswith("NM:i:"):
            return int(tag[5:])
    return None


def collect(lines, genes, min_mapq):
    """[(gene, nm, columns, read_length)] for the reads behind the clb counts."""
    reads = []
    for number, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            raise ValueError(f"SAM line {number}: {len(fields)} fields, expected at least 11")
        flag, pos, mapq, cigar = int(fields[1]), int(fields[3]), int(fields[4]), fields[5]
        if flag & FLAG_SKIP or mapq < min_mapq or cigar == "*":
            continue
        columns, length, span = read_metrics(cigar)
        gene = assign_gene(genes, fields[2], pos, pos + span - 1)
        if gene is None:
            continue
        nm = nm_tag(fields)
        if nm is None:
            raise ValueError(f"SAM line {number} ({fields[0]}): no NM tag; bowtie2 and "
                             "minimap2 always write one, so this is not their output")
        if columns == 0 or nm > columns:
            raise ValueError(f"SAM line {number} ({fields[0]}): NM {nm} with "
                             f"{columns} alignment columns")
        reads.append((gene, nm, columns, length))
    return reads


def quantile(sorted_values, q):
    """Linear interpolation between order statistics (numpy 'linear', R type 7)."""
    if not sorted_values:
        raise ValueError("quantile of no values")
    position = (len(sorted_values) - 1) * q
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    return sorted_values[lower] + (sorted_values[upper] - sorted_values[lower]) * (position - lower)


def fmt(value):
    value = float(value)
    return str(int(value)) if value.is_integer() else f"{value:.6f}"


def summarise(reads):
    """{metric: formatted value}; only the read count when there are no reads."""
    metrics = {"divergence_ihe3034_reads": str(len(reads))}
    if not reads:
        return metrics
    divergence = sorted(nm / columns for _g, nm, columns, _l in reads)
    nm_values = sorted(nm for _g, nm, _c, _l in reads)
    lengths = sorted(length for _g, _n, _c, length in reads)
    metrics["nm_ihe3034_median"] = fmt(quantile(nm_values, 0.5))
    metrics["divergence_ihe3034_median"] = fmt(quantile(divergence, 0.5))
    metrics["divergence_ihe3034_q1"] = fmt(quantile(divergence, 0.25))
    metrics["divergence_ihe3034_q3"] = fmt(quantile(divergence, 0.75))
    for pct in THRESHOLDS_PCT:
        above = sum(1 for _g, nm, columns, _l in reads if nm * 100 > pct * columns)
        metrics[f"divergence_ihe3034_frac_gt_{pct}pct"] = fmt(above / len(reads))
    metrics["clb_read_length_median"] = fmt(quantile(lengths, 0.5))
    metrics["clb_read_length_q1"] = fmt(quantile(lengths, 0.25))
    metrics["clb_read_length_q3"] = fmt(quantile(lengths, 0.75))
    metrics["clb_read_length_min"] = fmt(lengths[0])
    metrics["clb_read_length_max"] = fmt(lengths[-1])
    return metrics


def per_gene(reads):
    by_gene = defaultdict(list)
    for gene, nm, columns, _length in reads:
        by_gene[gene].append((nm, nm / columns))
    rows = []
    for gene in CLB_GENES:
        values = by_gene.get(gene, [])
        if values:
            rows.append((gene, len(values),
                         fmt(quantile(sorted(v[0] for v in values), 0.5)),
                         fmt(quantile(sorted(v[1] for v in values), 0.5))))
        else:
            rows.append((gene, 0, "NA", "NA"))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sample", required=True)
    parser.add_argument("--annotation", required=True, help="clbA-clbS gene GFF")
    parser.add_argument("--sam", default="-", help="SAM text (default: stdin)")
    parser.add_argument("--min-mapq", type=int, default=40)
    parser.add_argument("--qc-out", required=True)
    parser.add_argument("--genes-out", required=True)
    args = parser.parse_args(argv)

    genes = load_genes(args.annotation)
    handle = sys.stdin if args.sam == "-" else open(args.sam)
    try:
        reads = collect(handle, genes, args.min_mapq)
    finally:
        if handle is not sys.stdin:
            handle.close()

    metrics = summarise(reads)
    with open(args.qc_out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample", "Metric", "Value"])
        for metric in QC_METRICS:
            if metric in metrics:
                writer.writerow([args.sample, metric, metrics[metric]])
    with open(args.genes_out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(GENE_FIELDS)
        for row in per_gene(reads):
            writer.writerow([args.sample, *row])


if __name__ == "__main__":
    main()
