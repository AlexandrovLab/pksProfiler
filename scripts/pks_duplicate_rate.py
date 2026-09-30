#!/usr/bin/env python3
"""Duplicate rate among pks-island reads, by position collapse (WES QC).

Ludmil, 26 Sep mode memo section 3: "WES duplicate behavior must be checked separately
from WGS." The input alignment's own 0x400 flags (wes_input_qc.py) describe the library;
they do not reach the pks reads, which are extracted to FASTQ, depleted and re-aligned
single-end to IHE3034, losing every flag on the way. So the pks-supporting reads are
collapsed here, on the profiling BAM, by the rule a duplicate marker applies to unpaired
reads:

    key = (reference, strand, unclipped 5' position)

A read is a duplicate when an earlier read had the same key. Mates are separate records
(the profiling alignment is single-end), so a duplicated pair counts twice, once per mate
-- the same as its two records would in the input alignment.

Input: `samtools view -F 2308 <profiling BAM> <contig>:<island start>-<island end>`, i.e.
primary mapped records overlapping the canonical island at the profiling MAPQ >= 40 (the
BAM is already MAPQ-filtered). Position collapse cannot tell a PCR duplicate from two
fragments that happen to start at the same base, so at high island depth it
overestimates duplication; off-target WES island depth is low, where the two coincide.

Written metrics (prefix --prefix, default wes_pks_island):

    <prefix>_reads       primary mapped island records
    <prefix>_duplicates  of those, sharing a key with an earlier record
"""
import argparse
import re
import sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")
DEFAULT_PREFIX = "wes_pks_island"
# The metrics written with the default prefix (tests/test_qc_denominator.py reads it).
METRICS = (f"{DEFAULT_PREFIX}_reads", f"{DEFAULT_PREFIX}_duplicates")


def five_prime(flag, pos, cigar):
    """Unclipped 5' reference position (1-based) of an aligned record."""
    ops = [(int(n), op) for n, op in CIGAR.findall(cigar)]
    if not ops:
        raise ValueError(f"unparseable CIGAR {cigar!r}")
    if flag & 16:
        ref_len = sum(n for n, op in ops if op in REF_CONSUMING)
        trailing = 0
        for n, op in reversed(ops):
            if op not in "SH":
                break
            trailing += n
        return pos + ref_len - 1 + trailing
    leading = 0
    for n, op in ops:
        if op not in "SH":
            break
        leading += n
    return pos - leading


def duplicate_counts(sam_lines):
    """(records, duplicates) over primary mapped records of SAM text lines."""
    seen = set()
    records = duplicates = 0
    for line in sam_lines:
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        flag = int(fields[1])
        if flag & 2308:   # unmapped, secondary, supplementary
            continue
        key = (fields[2], bool(flag & 16), five_prime(flag, int(fields[3]), fields[5]))
        records += 1
        if key in seen:
            duplicates += 1
        else:
            seen.add(key)
    return records, duplicates


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sam", default="-", help="SAM text (default: stdin)")
    ap.add_argument("--prefix", default=DEFAULT_PREFIX)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    handle = sys.stdin if a.sam == "-" else open(a.sam)
    with handle:
        records, duplicates = duplicate_counts(handle)
    with open(a.output, "w") as out:
        out.write("Sample\tMetric\tValue\n")
        out.write(f"{a.sample}\t{a.prefix}_reads\t{records}\n")
        out.write(f"{a.sample}\t{a.prefix}_duplicates\t{duplicates}\n")


if __name__ == "__main__":
    main()
