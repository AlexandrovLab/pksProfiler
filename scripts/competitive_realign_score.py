#!/usr/bin/env python3
"""Arm B of the read-level pks specificity check: dual-index competitive reassignment.

Ludmil's concern (September 2026): MAPQ >= 40 against IHE3034 alone shows confident
mapping to *that* reference, not uniqueness against the broader microbial universe.
clbB/clbC/clbH/clbI/clbJ/clbK/clbN/clbO carry NRPS/PKS megasynthase domains that are
conserved across unrelated bacterial secondary-metabolite clusters (the same domain-
homology failure mode M1 already caught once at the assembled-MAG-bin level, see
tests/test_pks_positive_definition.py -- clbB/clbK "hit every bin by domain homology").
This is that check at the read level, for tumour samples where evidence is too sparse
to ever assemble a bin.

Two independent index architectures must agree that a read belongs to the true
pks+ reference and not a decoy, because a single large concatenated decoy index can
distort MAPQ/AS purely from panel size and composition, independent of biology:

  B1  ONE concatenated index (true reference + full decoy panel), aligned with
      secondary alignments enabled so a read's best score against EACH side of the
      panel is visible even when its primary alignment landed on the other side.
      Per read: AS_best(true), AS_best(decoy), delta_AS = AS_best(true) - AS_best(decoy),
      and MAPQ as recomputed in this competitive context.

  B2  TWO separate indices -- true reference alone, decoy panel alone -- same
      aligner/parameters, each run independently. Per read: compare AS (higher wins),
      then NM (lower wins) as tiebreak. This is the AstraZeneca-NGS/disambiguate
      convention. An exact tie (equal AS and equal NM, or absent from a side of the
      panel matched by no alignment at all) is ambiguous, never force-assigned.

A read is "confirmed_true" only if B1 and B2 independently agree; "confirmed_decoy"
only if both agree the other way; anything else -- disagreement between arms, or an
arm being unable to call it -- is "ambiguous". Exact ties inside either arm are
ambiguous by construction (see classify_b1 / classify_b2), never forced toward
either side.

Input is SAM or BAM/CRAM (name-sorted or not; reads are grouped by QNAME in memory,
so this script is meant for the small per-sample supporting-read set already
narrowed by the pipeline's own MAPQ>=40-vs-IHE3034 filter, not a whole-genome BAM).
Secondary/supplementary records are read for B1 (so both sides of the panel are
seen per read) and ignored for B2 (each side already carries only its own contigs,
so the primary/best record per read is what the alignment produced).

    competitive_realign_score.py \\
        --b1-sam sampleX.b1.concatenated.sam \\
        --b2-true-sam sampleX.b2.true.sam \\
        --b2-decoy-sam sampleX.b2.decoy.sam \\
        --true-contigs true_contigs.txt \\
        --decoy-contigs decoy_contigs.txt \\
        --output sampleX.arm_b.tsv
"""
import argparse
import subprocess
import sys
from pathlib import Path

FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800

FIELDS = [
    "read_id", "b1_as_true", "b1_as_decoy", "b1_decoy_contig", "b1_delta_as", "b1_mapq",
    "b1_call", "b2_as_true", "b2_nm_true", "b2_as_decoy", "b2_nm_decoy", "b2_decoy_contig",
    "b2_call", "classification",
]


def validate_contig_lists(true_contigs, decoy_contigs):
    """A contig cannot be both -- that would make classify_b1/classify_b2 compare a
    reference against itself and always call it a tie."""
    overlap = true_contigs & decoy_contigs
    if overlap:
        raise ValueError(f"Contig(s) listed as both true and decoy: {sorted(overlap)}")


def read_contig_list(path):
    """One contig name per line; blank lines and '#' comments ignored."""
    names = set()
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                names.add(line.split()[0])
    if not names:
        raise ValueError(f"No contig names found in {path}")
    return names


def _open_text(path):
    """Yield SAM alignment (non-header) lines, converting BAM/CRAM via samtools if needed."""
    suffix = Path(path).suffix.lower()
    if suffix in (".bam", ".cram"):
        proc = subprocess.run(["samtools", "view", str(path)], check=True,
                               capture_output=True, text=True)
        return proc.stdout.splitlines()
    with Path(path).open() as handle:
        return [line.rstrip("\n") for line in handle if not line.startswith("@")]


def _tag(fields, tag, kind=int):
    prefix = f"{tag}:"
    for field in fields[11:]:
        if field.startswith(prefix):
            # e.g. AS:i:37 -> take the value after the second ':'
            value = field.split(":", 2)[2]
            return kind(value)
    return None


def parse_sam_records(path):
    """Yield (qname, flag, rname, mapq, AS, NM) for every alignment record."""
    for line in _open_text(path):
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        qname, flag, rname = fields[0], int(fields[1]), fields[2]
        mapq = int(fields[4])
        as_score = _tag(fields, "AS", int)
        nm = _tag(fields, "NM", int)
        yield qname, flag, rname, mapq, as_score, nm


def best_by_family(records, true_contigs, decoy_contigs):
    """B1: per read, the best AS seen against each family, plus the primary record's MAPQ
    and which specific decoy contig produced the best decoy-side score.

    records: iterable of (qname, flag, rname, mapq, as_score, nm) for ONE read.
    Returns dict(as_true, as_decoy, decoy_contig, mapq).
    """
    as_true = as_decoy = None
    decoy_contig = None
    mapq = None
    for _qname, flag, rname, rec_mapq, as_score, _nm in records:
        if flag & FLAG_UNMAPPED:
            continue
        if not (flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)):
            mapq = rec_mapq  # the primary record's MAPQ, whichever family it landed on
        if as_score is None:
            continue
        if rname in true_contigs:
            if as_true is None or as_score > as_true:
                as_true = as_score
        elif rname in decoy_contigs:
            if as_decoy is None or as_score > as_decoy:
                as_decoy = as_score
                decoy_contig = rname
    return {"as_true": as_true, "as_decoy": as_decoy, "decoy_contig": decoy_contig, "mapq": mapq}


def classify_b1(as_true, as_decoy, min_delta_as=0):
    """Which side the concatenated-index evidence favours. None means "no alignment
    to that family at all" -- distinct from an AS of 0, which cannot occur for a
    mapped read under minimap2/bowtie2 scoring, but is handled the same way regardless.
    """
    if as_true is None and as_decoy is None:
        return "ambiguous"
    if as_decoy is None:
        return "true"
    if as_true is None:
        return "decoy"
    delta = as_true - as_decoy
    if delta > min_delta_as:
        return "true"
    if delta < -min_delta_as:
        return "decoy"
    return "ambiguous"  # includes the exact tie (delta == 0) at min_delta_as == 0


def classify_b2(as_true, nm_true, as_decoy, nm_decoy):
    """Disambiguate convention: AS first (higher wins), NM as tiebreak (lower wins).
    An exact tie on both is ambiguous, never force-assigned to either side.
    """
    if as_true is None and as_decoy is None:
        return "ambiguous"
    if as_decoy is None:
        return "true"
    if as_true is None:
        return "decoy"
    if as_true > as_decoy:
        return "true"
    if as_true < as_decoy:
        return "decoy"
    # AS tie -> NM tiebreak
    if nm_true is None or nm_decoy is None:
        return "ambiguous"
    if nm_true < nm_decoy:
        return "true"
    if nm_true > nm_decoy:
        return "decoy"
    return "ambiguous"  # exact tie on both AS and NM


def final_classification(b1_call, b2_call):
    """A read is alignment-confirmed only if both independent arms agree."""
    if b1_call == "true" and b2_call == "true":
        return "confirmed_true"
    if b1_call == "decoy" and b2_call == "decoy":
        return "confirmed_decoy"
    return "ambiguous"


def group_by_read(records):
    grouped = {}
    for record in records:
        grouped.setdefault(record[0], []).append(record)
    return grouped


def best_single(records):
    """B2: the single alignment against a lone-family index. Prefer the primary
    record's AS/NM (a lone-family run should not normally emit secondaries, but if
    it does the primary is still the aligner's own best call)."""
    as_score = nm = None
    for _qname, flag, _rname, _mapq, rec_as, rec_nm in records:
        if flag & FLAG_UNMAPPED:
            continue
        if not (flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)):
            return rec_as, rec_nm
        if as_score is None or (rec_as is not None and rec_as > as_score):
            as_score, nm = rec_as, rec_nm
    return as_score, nm


def score_reads(b1_records, b2_true_records, b2_decoy_records, true_contigs, decoy_contigs,
                 min_delta_as=0):
    b1_by_read = group_by_read(b1_records)
    b2_true_by_read = group_by_read(b2_true_records)
    b2_decoy_by_read = group_by_read(b2_decoy_records)

    read_ids = sorted(set(b1_by_read) | set(b2_true_by_read) | set(b2_decoy_by_read))
    rows = []
    for read_id in read_ids:
        b1 = best_by_family(b1_by_read.get(read_id, []), true_contigs, decoy_contigs)
        b1_call = classify_b1(b1["as_true"], b1["as_decoy"], min_delta_as)
        delta_as = (None if b1["as_true"] is None or b1["as_decoy"] is None
                    else b1["as_true"] - b1["as_decoy"])

        as_true2, nm_true2 = best_single(b2_true_by_read.get(read_id, []))
        as_decoy2, nm_decoy2 = best_single(b2_decoy_by_read.get(read_id, []))
        decoy_contig2 = None
        for _qname, flag, rname, _mapq, _as, _nm in b2_decoy_by_read.get(read_id, []):
            if not (flag & FLAG_UNMAPPED):
                decoy_contig2 = rname
                break
        b2_call = classify_b2(as_true2, nm_true2, as_decoy2, nm_decoy2)

        rows.append({
            "read_id": read_id,
            "b1_as_true": b1["as_true"], "b1_as_decoy": b1["as_decoy"],
            "b1_decoy_contig": b1["decoy_contig"], "b1_delta_as": delta_as,
            "b1_mapq": b1["mapq"], "b1_call": b1_call,
            "b2_as_true": as_true2, "b2_nm_true": nm_true2,
            "b2_as_decoy": as_decoy2, "b2_nm_decoy": nm_decoy2,
            "b2_decoy_contig": decoy_contig2, "b2_call": b2_call,
            "classification": final_classification(b1_call, b2_call),
        })
    return rows


def write_rows(rows, output):
    with Path(output).open("w", newline="") as handle:
        handle.write("\t".join(FIELDS) + "\n")
        for row in rows:
            handle.write("\t".join("" if row[field] is None else str(row[field])
                                    for field in FIELDS) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--b1-sam", required=True,
                         help="SAM/BAM: reads aligned to the concatenated true+decoy index, "
                              "secondary alignments enabled (e.g. minimap2 -N 5)")
    parser.add_argument("--b2-true-sam", required=True,
                         help="SAM/BAM: the same reads aligned to the true-reference-only index")
    parser.add_argument("--b2-decoy-sam", required=True,
                         help="SAM/BAM: the same reads aligned to the decoy-panel-only index")
    parser.add_argument("--true-contigs", required=True,
                         help="text file, one true-reference contig name per line")
    parser.add_argument("--decoy-contigs", required=True,
                         help="text file, one decoy-panel contig name per line")
    parser.add_argument("--min-delta-as", type=float, default=0.0,
                         help="B1 margin (AS units) required to call true/decoy rather than "
                              "ambiguous; 0 calls any non-zero difference (default: 0)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    true_contigs = read_contig_list(args.true_contigs)
    decoy_contigs = read_contig_list(args.decoy_contigs)
    try:
        validate_contig_lists(true_contigs, decoy_contigs)
    except ValueError as exc:
        raise SystemExit(str(exc))

    rows = score_reads(
        list(parse_sam_records(args.b1_sam)),
        list(parse_sam_records(args.b2_true_sam)),
        list(parse_sam_records(args.b2_decoy_sam)),
        true_contigs, decoy_contigs, args.min_delta_as,
    )
    write_rows(rows, args.output)

    confirmed_true = sum(1 for r in rows if r["classification"] == "confirmed_true")
    confirmed_decoy = sum(1 for r in rows if r["classification"] == "confirmed_decoy")
    ambiguous = sum(1 for r in rows if r["classification"] == "ambiguous")
    print(f"reads={len(rows)} confirmed_true={confirmed_true} "
          f"confirmed_decoy={confirmed_decoy} ambiguous={ambiguous}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
