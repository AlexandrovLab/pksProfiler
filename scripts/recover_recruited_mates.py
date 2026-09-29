#!/usr/bin/env python3
"""Rebuild R1/R2 pairs (plus singletons) from a single interleaved or merged FASTQ.

Two uses, one pairing rule:

  recruitment (--sam given)   targeted assembly: only fragments with a credible local pks
                              alignment in the recruitment SAM are recovered, both mates.
  whole sample (--all-reads)  metagenome assembly (Ludmil, 26 Sep mode memo section 5):
                              every host-depleted read, reconstructed into R1/R2 plus
                              singletons, so MEGAHIT and read-to-contig alignment get pairs.
                              Depletion stays read-independent; pairs are rebuilt after it.

A mate is identified from the read name: a trailing /1 or /2 (what extractReads' samtools
fastq -N and filterReads write), else a Casava-style second field `1:...` / `2:...`, else
the read is unsuffixed (mate 0). A fragment with one /1 and one /2 is a pair; everything
else -- unsuffixed reads, orphans whose mate was depleted, and any further read claiming
a mate slot already taken -- is written as a singleton. No read is dropped or overwritten,
so R1 + R2 + singletons is exactly the input: assembly never narrows.

Whole-sample mode does not hold the sample in memory: records are keyed by fragment name
into one temporary file, sorted with an external `sort`, and paired in one streaming pass.

Mate-tag diagnostics (A7, Ludmil 25 Sep), in the stats table either way:
  mate_tag_1, mate_tag_2           reads carrying each suffix
  fragments_without_mate_tag       unsuffixed reads (normal for single-end input)
  orphan_reads                     a /1 or /2 whose mate is absent
  duplicate_mate_slot_overwrites   reads whose (name, mate) slot was already filled. These
                                   used to overwrite the earlier read; they are now kept as
                                   singletons, warned about when > 0, and fail the task
                                   above --max-duplicate-fraction.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter

CENSUS = ("mate_tag_1", "mate_tag_2", "fragments_without_mate_tag", "orphan_reads",
          "duplicate_mate_slot_overwrites")


def base_and_mate(header):
    token = header.split()[0].lstrip("@")
    if token.endswith("/1") or token.endswith("/2"):
        return token[:-2], token[-1]
    fields = header.split()
    if len(fields) > 1 and fields[1][:1] in ("1", "2") and fields[1][1:2] in ("", ":"):
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


def open_text(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t", **({"compresslevel": 1} if "w" in mode else {}))
    return open(path, mode)


def fastq_records(path):
    with open_text(path, "r") as handle:
        while True:
            record = [handle.readline() for _ in range(4)]
            if not record[0]:
                return
            if not record[3]:
                raise SystemExit(f"ERROR: {path}: truncated FASTQ record at {record[0].strip()!r}")
            yield record


class Pairer:
    """Writes fragments as pairs or singletons and keeps the mate-tag census."""

    def __init__(self, r1, r2, single):
        self.r1, self.r2, self.single = r1, r2, single
        self.census = Counter()

    def fragment(self, mates):
        """`mates` maps "1"/"2"/"0" to the records claiming that slot, in input order."""
        ones, twos, zeros = mates.get("1", []), mates.get("2", []), mates.get("0", [])
        self.census["mate_tag_1"] += len(ones)
        self.census["mate_tag_2"] += len(twos)
        self.census["fragments_without_mate_tag"] += len(zeros)
        self.census["duplicate_mate_slot_overwrites"] += (max(len(ones) - 1, 0) + max(len(twos) - 1, 0)
                                                          + max(len(zeros) - 1, 0))
        if ones and twos:
            self.r1.write(ones[0])
            self.r2.write(twos[0])
            self.census["paired_fragments"] += 1
            rest = ones[1:] + twos[1:]
        else:
            self.census["orphan_reads"] += len(ones[:1]) + len(twos[:1])
            rest = ones + twos
        for record in rest + zeros:
            self.single.write(record)
            self.census["singleton_reads"] += 1


def pair_recruited(fastq, wanted, pairer):
    found = {}
    for record in fastq_records(fastq):
        base, mate = base_and_mate(record[0].rstrip())
        if base in wanted:
            found.setdefault(base, {}).setdefault(mate, []).append("".join(record))
    for mates in found.values():
        pairer.fragment(mates)
    return len(found)


def pair_all(fastq, pairer, tmpdir=None, sort_memory="1G", threads=1):
    """Every read, paired by an external sort on the fragment name. Returns input reads."""
    total = 0
    with tempfile.TemporaryDirectory(dir=tmpdir) as work:
        keyed = os.path.join(work, "keyed.tsv")
        with open(keyed, "w") as out:
            for index, record in enumerate(fastq_records(fastq)):
                base, mate = base_and_mate(record[0].rstrip())
                fields = [line.rstrip("\n").replace("\t", " ") for line in record]
                # the index keeps input order within a slot, so "first" is well defined
                out.write(f"{base}\t{mate}\t{index:012d}\t" + "\t".join(fields) + "\n")
                total += 1
        sorted_path = os.path.join(work, "sorted.tsv")
        subprocess.run(["sort", "-t", "\t", "-k1,1", "-k2,2", "-k3,3", "-S", sort_memory,
                        f"--parallel={max(threads, 1)}", "-T", work, "-o", sorted_path, keyed],
                       check=True, env={**os.environ, "LC_ALL": "C"})
        os.unlink(keyed)
        current, mates = None, {}
        with open(sorted_path) as handle:
            for line in handle:
                base, mate, _index, header, seq, plus, qual = line.rstrip("\n").split("\t")
                if base != current:
                    if current is not None:
                        pairer.fragment(mates)
                    current, mates = base, {}
                mates.setdefault(mate, []).append(f"{header}\n{seq}\n{plus}\n{qual}\n")
            if current is not None:
                pairer.fragment(mates)
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sam", help="recruitment SAM: recover only these fragments")
    source.add_argument("--all-reads", action="store_true", help="pair every read (metagenome assembly)")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--single", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--min-aligned-bases", type=int, default=60)
    parser.add_argument("--min-identity", type=float, default=.90)
    parser.add_argument("--max-duplicate-fraction", type=float, default=0.001,
                        help="fail when duplicate_mate_slot_overwrites / reads exceeds this")
    parser.add_argument("--sort-memory", default="1G")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--tmpdir", default=None)
    args = parser.parse_args()

    with open_text(args.r1, "w") as r1, open_text(args.r2, "w") as r2, open_text(args.single, "w") as single:
        pairer = Pairer(r1, r2, single)
        if args.sam:
            wanted, alignments = recruited_ids(args.sam, args.min_aligned_bases, args.min_identity)
            found = pair_recruited(args.fastq, wanted, pairer)
            total = None
        else:
            total = pair_all(args.fastq, pairer, args.tmpdir, args.sort_memory, args.threads)

    census = pairer.census
    written = 2 * census["paired_fragments"] + census["singleton_reads"]
    if total is not None and written != total:
        raise SystemExit(f"ERROR: read {total} reads but wrote {written}: pairing must not "
                         "drop or duplicate a read")
    if args.sam:
        rows = [["credible_alignments", alignments], ["recruited_fragment_ids", len(wanted)],
                ["fragment_ids_found", found]]
    else:
        rows = [["input_reads", total]]
    rows += [["paired_fragments", census["paired_fragments"]], ["singleton_reads", census["singleton_reads"]]]
    if args.sam:
        rows += [["min_aligned_bases", args.min_aligned_bases], ["min_identity", args.min_identity]]
    rows += [[name, census[name]] for name in CENSUS]
    with open(args.stats, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)

    duplicates = census["duplicate_mate_slot_overwrites"]
    if duplicates:
        print(f"WARNING: {duplicates} read(s) claimed a mate slot already taken (duplicate read "
              "names); kept as singletons", file=sys.stderr)
        if written and duplicates / written > args.max_duplicate_fraction:
            raise SystemExit(f"ERROR: {duplicates} of {written} reads have duplicate names "
                             f"(> --max-duplicate-fraction {args.max_duplicate_fraction}); "
                             "the read names are malformed")
    if census["fragments_without_mate_tag"] and (census["mate_tag_1"] or census["mate_tag_2"]):
        print(f"WARNING: {census['fragments_without_mate_tag']} unsuffixed read(s) among "
              "mate-tagged ones; written as singletons", file=sys.stderr)


if __name__ == "__main__":
    main()
