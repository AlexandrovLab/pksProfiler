#!/usr/bin/env python3
"""Residual-human check for RNA-seq, after every host-depletion pass.

Ludmil, 26 Sep mode memo section 4: record "total, unmapped, host-depleted, and residual-
human read counts". Total, unmapped and host-depleted are already QC columns; this is the
fourth. Definition (docs/running/wes_rna.md):

    Take a capped, systematic sample of the FINAL host-depleted reads (at most
    --rna_residual_human_sample_reads, default 200,000: every k-th read, k chosen so the
    sample fits the cap), realign it to the RNA splice index (--rna_splice_db) with
    minimap2 -ax <the human pass's splice preset> plus --rna_residual_human_args
    (default "-n 2 -m 20 -s 30": fewer minimizers, a lower chaining score and a lower
    alignment score than the preset's own floors), and count the reads with a primary
    alignment. residual-human fraction = mapped / sampled.

Why relaxed thresholds: re-running the depletion pass's own command on its own output
maps nothing by construction (a read it could place was removed). A human read survives
depletion by scoring under the preset's floors -- divergent, short-exon, junction-heavy --
and the relaxed check is exactly the search for those. It is a measurement, never a
filter: nothing is removed from the reads that go on to profiling. On simulated 100 bp
reads the relaxed floors recover about twice as many 12-16%-divergent reads as the preset
and placed none of 5,000 random reads on a 20 Mb random reference (the numbers are in
docs/running/wes_rna.md). Long reads
(--read_type ont|hifi) get the preset alone: those floors are on a different scale.

    rna_residual_human.py subsample --input reads.fastq.gz --cap N --output sample.fastq \
        --stats sample.tsv
    minimap2 -ax splice:sr ... index sample.fastq | \
        rna_residual_human.py count --sam - --sample S --stats sample.tsv --output qc.tsv

Written metrics:

    rna_residual_input_reads    final host-depleted reads the sample was drawn from
    rna_residual_checked_reads  reads in the capped sample
    rna_residual_human_reads    of those, with a primary alignment to the splice index
"""
import argparse
import gzip
import sys
from pathlib import Path

METRICS = ("rna_residual_input_reads", "rna_residual_checked_reads", "rna_residual_human_reads")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def count_fastq(path):
    with open_text(path) as fh:
        return sum(1 for _ in fh) // 4


def step_for(total, cap):
    """Every step-th read, so a sample of total reads holds at most cap of them."""
    if cap <= 0:
        raise ValueError("cap must be positive")
    return max(1, -(-total // cap))


def subsample(path, cap, output):
    """Systematic sample of at most `cap` records; returns (total, sampled)."""
    total = count_fastq(path)
    step = step_for(total, cap)
    sampled = 0
    with open_text(path) as fh, open(output, "w") as out:
        index = 0
        while sampled < cap:
            record = [fh.readline() for _ in range(4)]
            if not record[0]:
                break
            if index % step == 0:
                out.writelines(record)
                sampled += 1
            index += 1
    return total, sampled


def mapped_reads(sam_lines):
    """Distinct read names with a primary mapped record."""
    names = set()
    for line in sam_lines:
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t", 2)
        if len(fields) < 2:
            continue
        if int(fields[1]) & 2308 == 0:
            names.add(fields[0])
    return len(names)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("subsample")
    s.add_argument("--input", required=True)
    s.add_argument("--cap", required=True, type=int)
    s.add_argument("--output", required=True)
    s.add_argument("--stats", required=True)
    c = sub.add_parser("count")
    c.add_argument("--sam", default="-")
    c.add_argument("--sample", required=True)
    c.add_argument("--stats", required=True, help="the subsample step's --stats")
    c.add_argument("--output", required=True)
    a = ap.parse_args()

    if a.command == "subsample":
        total, sampled = subsample(a.input, a.cap, a.output)
        Path(a.stats).write_text(f"{total}\t{sampled}\n")
        return
    total, sampled = (int(x) for x in Path(a.stats).read_text().split())
    handle = sys.stdin if a.sam == "-" else open(a.sam)
    with handle:
        mapped = mapped_reads(handle)
    if mapped > sampled:
        raise SystemExit(f"{mapped} mapped reads exceed the {sampled} checked")
    with open(a.output, "w") as out:
        out.write("Sample\tMetric\tValue\n")
        out.write(f"{a.sample}\trna_residual_input_reads\t{total}\n")
        out.write(f"{a.sample}\trna_residual_checked_reads\t{sampled}\n")
        out.write(f"{a.sample}\trna_residual_human_reads\t{mapped}\n")


if __name__ == "__main__":
    main()
