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

Recruitment classes (A1/A2, Ludmil 25 Sep). The recruitment reference is the canonical
island plus a flank either side (indices/pks_island: 5 kb, header `island=... flank=...`);
--island-start/--island-end give the island's 1-based inclusive span in that reference.
Each primary alignment (Bowtie2 best-alignment reporting, no -k) is first held to
MAPQ >= --min-mapq, >= --min-aligned-bases aligned bases and >= --min-identity, then
placed by where it lies:

  boundary_spanning  crosses an island boundary (more than --boundary-tolerance-bp on each
                     side), or lies inside the island with its outer end within the
                     tolerance of a boundary and a soft clip longer than the tolerance
                     pointing out of the island. No query-coverage floor: a true junction
                     read is expected to be clipped where the host locus differs.
  island_internal    inside the island (within the tolerance) and query coverage
                     (query bases in the alignment / read length) >= --min-query-coverage.
                     A 150 bp read with a 60 bp local match (0.40) is rejected.
  flank_only         in a flank (overlapping the island by at most the tolerance) and query
                     coverage >= --min-query-coverage. A flank read clipped INTO the island
                     is flank-only: that is what a pks-negative relative sharing the flank
                     looks like, so it never counts as island support.
  rejected           everything else.

A fragment takes the best class among its mates (island_internal > boundary_spanning >
flank_only). island_supporting_fragments = island_internal + boundary_spanning drives the
"why did assembly fail?" diagnostic; flank-only fragments still go to the assemblers for
context unless --exclude-flank-only.

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
import contextlib
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


CLASSES = ("island_internal", "boundary_spanning", "flank_only")
RANK = {"island_internal": 3, "boundary_spanning": 2, "flank_only": 1, "rejected": 0}
REJECTIONS = ("low_mapq", "short_alignment", "low_identity", "low_query_coverage", "unmapped")
ISLAND_SUPPORTING = ("island_internal", "boundary_spanning")


def cigar_ops(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def aligned_bases(cigar):
    return sum(n for n, op in cigar_ops(cigar) if op in "M=X")


class RecruitmentRule:
    """The A1 recruitment rule. island_start/island_end: 1-based inclusive, in the
    coordinates of the island+flank recruitment reference."""

    def __init__(self, island_start, island_end, min_aligned_bases=60, min_identity=.90,
                 min_query_coverage=.80, min_mapq=20, boundary_tolerance_bp=10,
                 reference_length=None, include_flank_only=True):
        if not 1 <= island_start <= island_end:
            raise ValueError(f"bad island span {island_start}-{island_end}")
        if reference_length is not None and island_end > reference_length:
            raise ValueError(f"island end {island_end} beyond reference length {reference_length}")
        if not 0 < min_query_coverage <= 1:
            raise ValueError("min_query_coverage must be in (0,1]")
        if boundary_tolerance_bp < 0 or min_mapq < 0 or min_aligned_bases < 1:
            raise ValueError("tolerance and MAPQ must be >= 0, min aligned bases >= 1")
        self.island_start, self.island_end = island_start, island_end
        self.min_aligned_bases, self.min_identity = min_aligned_bases, min_identity
        self.min_query_coverage, self.min_mapq = min_query_coverage, min_mapq
        self.tolerance = boundary_tolerance_bp
        self.reference_length = reference_length
        self.include_flank_only = include_flank_only

    def assembled(self, cls):
        return cls in ISLAND_SUPPORTING or (cls == "flank_only" and self.include_flank_only)


def classify_alignment(pos, cigar, mapq, nm, rule):
    """(class, rejection reason or None) for one primary alignment. pos is SAM POS (1-based)."""
    ops = cigar_ops(cigar)
    if not ops:
        return "rejected", "unmapped"
    matched = sum(n for n, op in ops if op in "M=X")
    query_in_alignment = matched + sum(n for n, op in ops if op == "I")
    query_length = query_in_alignment + sum(n for n, op in ops if op in "SH")
    ref_length = sum(n for n, op in ops if op in "M=XDN")
    lead_clip = trail_clip = 0
    for n, op in ops:
        if op not in "SH":
            break
        lead_clip += n
    for n, op in reversed(ops):
        if op not in "SH":
            break
        trail_clip += n
    if mapq < rule.min_mapq:
        return "rejected", "low_mapq"
    if matched < rule.min_aligned_bases:
        return "rejected", "short_alignment"
    identity = (matched - nm) / matched if matched else 0
    if identity < rule.min_identity:
        return "rejected", "low_identity"
    coverage = query_in_alignment / query_length if query_length else 0

    ref_start, ref_end = pos, pos + ref_length - 1
    left, right, tol = rule.island_start, rule.island_end, rule.tolerance
    island_bp = max(0, min(ref_end, right) - max(ref_start, left) + 1)
    flank_bp = ref_length - island_bp
    if island_bp > tol and flank_bp > tol:
        return "boundary_spanning", None          # crosses a boundary, anchored both sides
    if island_bp > tol:                           # island side
        outward_left = abs(ref_start - left) <= tol and lead_clip > tol
        outward_right = abs(ref_end - right) <= tol and trail_clip > tol
        if outward_left or outward_right:
            return "boundary_spanning", None      # island end reached, host sequence differs
        if coverage >= rule.min_query_coverage:
            return "island_internal", None
        return "rejected", "low_query_coverage"
    if coverage >= rule.min_query_coverage:       # flank side, clipped-into-island included
        return "flank_only", None
    return "rejected", "low_query_coverage"


def check_reference(sq_lengths, rule):
    """The SAM must come from the reference the island coordinates describe."""
    if rule.reference_length is None:
        return
    if len(sq_lengths) != 1 or sq_lengths[0] != rule.reference_length:
        raise SystemExit(f"ERROR: recruitment SAM @SQ lengths {sq_lengths} do not match the "
                         f"expected island+flank reference length {rule.reference_length}; the "
                         "recruitment index and --island-start/--island-end disagree "
                         "(rebuild the index or fix pks_recruit_flank_bp)")


def classify_sam(sam, rule):
    """Per-fragment best class and the per-read census."""
    best, reads = {}, Counter()
    sq_lengths, checked = [], False
    # "-": the SAM on stdin (A8: targetedPksRecruit pipes bowtie2 straight in). It is read
    # exactly once, in order, so a stream is all this needs.
    with (contextlib.nullcontext(sys.stdin) if sam == "-" else open(sam)) as handle:
        for line in handle:
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    sq_lengths += [int(f[3:]) for f in line.rstrip("\n").split("\t") if f.startswith("LN:")]
                continue
            if not checked:
                check_reference(sq_lengths, rule)
                checked = True
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[1])
            if flag & 256 or flag & 2048:
                continue
            base = base_and_mate(fields[0])[0]
            if flag & 4:
                cls, reason = "rejected", "unmapped"
            else:
                tags = {item.split(":", 2)[0]: item.split(":", 2)[2] for item in fields[11:] if item.count(":") >= 2}
                matched = aligned_bases(fields[5])
                cls, reason = classify_alignment(int(fields[3]), fields[5], int(fields[4]),
                                                 int(tags.get("NM", matched)), rule)
            reads["aligned_reads" if reason != "unmapped" else "unmapped_reads"] += 1
            reads[f"reads_{cls}"] += 1
            if reason:
                reads[f"reads_rejected_{reason}"] += 1
            if base not in best or RANK[cls] > RANK[best[base]]:
                best[base] = cls
    if not checked:
        check_reference(sq_lengths, rule)
    return best, reads


def recruited_ids(sam, rule):
    """Fragment ids to assemble, and the number of accepted (non-rejected) alignments."""
    best, reads = classify_sam(sam, rule)
    ids = {base for base, cls in best.items() if rule.assembled(cls)}
    return ids, sum(reads[f"reads_{c}"] for c in CLASSES)


def recruitment_rows(best, reads, rule):
    fragments = Counter(best.values())
    rows = [["aligned_reads", reads["aligned_reads"]]]
    rows += [[f"reads_{c}", reads[f"reads_{c}"]] for c in CLASSES + ("rejected",)]
    rows += [[f"reads_rejected_{r}", reads[f"reads_rejected_{r}"]] for r in REJECTIONS if r != "unmapped"]
    rows += [[f"fragments_{c}", fragments[c]] for c in CLASSES + ("rejected",)]
    rows += [["island_supporting_fragments", sum(fragments[c] for c in ISLAND_SUPPORTING)],
             ["flank_only_fragments_assembled", fragments["flank_only"] if rule.include_flank_only else 0]]
    return rows


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
    source.add_argument("--sam", help="recruitment SAM, or - for stdin: recover only these fragments")
    source.add_argument("--all-reads", action="store_true", help="pair every read (metagenome assembly)")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--single", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--min-aligned-bases", type=int, default=60)
    parser.add_argument("--min-identity", type=float, default=.90)
    parser.add_argument("--min-query-coverage", type=float, default=.80,
                        help="island_internal / flank_only: aligned query bases / read length")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--boundary-tolerance-bp", type=int, default=10)
    parser.add_argument("--island-start", type=int,
                        help="island's first base, 1-based, in recruitment-reference coordinates")
    parser.add_argument("--island-end", type=int,
                        help="island's last base, 1-based inclusive, recruitment-reference coordinates")
    parser.add_argument("--reference-length", type=int,
                        help="expected @SQ LN of the recruitment reference (checked)")
    parser.add_argument("--exclude-flank-only", action="store_true",
                        help="do not pass flank-only fragments to the assemblers")
    parser.add_argument("--max-duplicate-fraction", type=float, default=0.001,
                        help="fail when duplicate_mate_slot_overwrites / reads exceeds this")
    parser.add_argument("--sort-memory", default="1G")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--tmpdir", default=None)
    args = parser.parse_args()
    rule = None
    if args.sam:
        if args.island_start is None or args.island_end is None:
            parser.error("--sam needs --island-start and --island-end")
        try:
            rule = RecruitmentRule(args.island_start, args.island_end, args.min_aligned_bases,
                                   args.min_identity, args.min_query_coverage, args.min_mapq,
                                   args.boundary_tolerance_bp, args.reference_length,
                                   not args.exclude_flank_only)
        except ValueError as error:
            parser.error(str(error))

    with open_text(args.r1, "w") as r1, open_text(args.r2, "w") as r2, open_text(args.single, "w") as single:
        pairer = Pairer(r1, r2, single)
        if args.sam:
            best, reads = classify_sam(args.sam, rule)
            wanted = {base for base, cls in best.items() if rule.assembled(cls)}
            alignments = sum(reads[f"reads_{c}"] for c in CLASSES)
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
        rows += recruitment_rows(best, reads, rule)
        rows += [["min_aligned_bases", args.min_aligned_bases], ["min_identity", args.min_identity],
                 ["min_query_coverage", args.min_query_coverage], ["min_mapq", args.min_mapq],
                 ["boundary_tolerance_bp", args.boundary_tolerance_bp],
                 ["include_flank_only", str(rule.include_flank_only).lower()],
                 ["island_start_in_reference", args.island_start],
                 ["island_end_in_reference", args.island_end],
                 ["reference_length", args.reference_length if args.reference_length else "NA"]]
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
