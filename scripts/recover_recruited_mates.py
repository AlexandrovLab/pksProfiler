#!/usr/bin/env python3
"""Recover both mates for read fragments with a credible local pks alignment."""
from __future__ import annotations
import argparse, csv, gzip, re

def base_and_mate(header):
    token = header.split()[0].lstrip("@")
    if token.endswith("/1") or token.endswith("/2"):
        return token[:-2], token[-1]
    fields = header.split()
    if len(fields) > 1 and fields[1][:1] in ("1", "2"):
        return token, fields[1][0]
    return token, "0"

def aligned_bases(cigar):
    return sum(int(n) for n, op in re.findall(r"(\d+)([M=X])", cigar))

def recruited_ids(sam, min_bases, min_identity):
    ids, alignments = set(), 0
    with open(sam) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip().split("\t")
            flag = int(fields[1])
            if flag & 4 or flag & 256 or flag & 2048:
                continue
            aligned = aligned_bases(fields[5])
            tags = {item.split(":", 2)[0]: item.split(":", 2)[2] for item in fields[11:] if item.count(":") >= 2}
            nm = int(tags.get("NM", aligned))
            identity = (aligned - nm) / aligned if aligned else 0
            if aligned >= min_bases and identity >= min_identity:
                ids.add(base_and_mate(fields[0])[0])
                alignments += 1
    return ids, alignments

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", required=True)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--single", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--min-aligned-bases", type=int, default=60)
    parser.add_argument("--min-identity", type=float, default=.90)
    args = parser.parse_args()
    wanted, alignments = recruited_ids(args.sam, args.min_aligned_bases, args.min_identity)
    found = {}
    with gzip.open(args.fastq, "rt") as handle:
        while True:
            record = [handle.readline() for _ in range(4)]
            if not record[0]:
                break
            base, mate = base_and_mate(record[0].rstrip())
            if base in wanted:
                found.setdefault(base, {})[mate] = "".join(record)
    pairs = singles = 0
    with gzip.open(args.r1, "wt") as r1, gzip.open(args.r2, "wt") as r2, gzip.open(args.single, "wt") as single:
        for mates in found.values():
            if "1" in mates and "2" in mates:
                r1.write(mates["1"])
                r2.write(mates["2"])
                pairs += 1
            else:
                for record in mates.values():
                    single.write(record)
                    singles += 1
    with open(args.stats, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows([["credible_alignments", alignments], ["recruited_fragment_ids", len(wanted)],
                          ["fragment_ids_found", len(found)], ["paired_fragments", pairs],
                          ["singleton_reads", singles], ["min_aligned_bases", args.min_aligned_bases],
                          ["min_identity", args.min_identity]])

if __name__ == "__main__":
    main()
