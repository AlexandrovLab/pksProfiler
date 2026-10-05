#!/usr/bin/env python3
"""V8b CARRIER ATTRIBUTION: which pks-positive organism most likely carries the island?

Ludmil, 25 Sep report, section 4:

  "V8b. Carrier attribution. Separately compare genuine pks reads against alternative
   pks-positive carriers/strains (pks+ E. coli, Klebsiella, Citrobacter, Enterobacter,
   etc.). This asks which organism is the likely carrier and should not be conflated
   with locus-detection specificity."

This script answers ONLY that question, and its vocabulary is deliberately disjoint from
V8a's. V8a (scripts/competitive_realign_score.py + scripts/specificity_decision.py)
decides whether a read belongs to the pks/clb locus at all, against pks-NEGATIVE decoys,
and says confirmed_true / confirmed_decoy / ambiguous and SURVIVES /
NOT_DISTINGUISHABLE_FROM_CROSS_MAPPING / UNDERPOWERED. Nothing here says any of those
words. Every member of the panel this script reads is pks-POSITIVE (see
ref/specificity_panels/v8b_carrier_panel.tsv and the island screen in
scripts/build_specificity_panels.py), so a read losing to another carrier is not evidence
against the locus -- it is evidence about the organism. The two outputs must never be
merged into one "specificity" number.

Method. The supporting reads are aligned once against the concatenated carrier panel with
secondary alignments enabled, so every carrier's best alignment score for a read is
visible even when the primary landed on another carrier. Per read:

    AS_best(carrier)         best alignment score against each carrier
    best                     the carrier with the highest AS_best
    delta_AS                 AS_best(best) - AS_best(runner-up, a DIFFERENT carrier)
    tied set                 every carrier within --tie-margin AS of the best

    attributed               exactly one carrier in the tied set and delta_AS >= --min-delta-as
    ambiguous_multi_carrier  otherwise (several carriers explain the read equally well)
    unmapped                 no alignment to any carrier
    off_map                  aligned, but only to contigs the carrier map does not name --
                             the index and the map are from different panel builds, and
                             the run is refused (see main) rather than reported as unmapped

Per sample, the attributed reads are counted per carrier and per genus. The sample-level
call is:

    CARRIER_ATTRIBUTED     >= --min-attributed-reads attributed reads AND the top carrier
                           holds >= --min-top-carrier-fraction of them
    CARRIER_GENUS_ONLY     the same fraction test passes at genus level but not at strain
                           level: the organism's genus is supported, the strain is not
    CARRIER_AMBIGUOUS      neither test passes although >= --min-attributed-reads reads
                           aligned to the panel: the read set does not separate the panel's
                           carriers (including the case where EVERY read ties across
                           genera, so nothing is attributed at all)
    CARRIER_UNDETERMINED   fewer than --min-attributed-reads reads aligned to the panel at
                           all: too little evidence to ask the question

5 Oct 2026 (V8 real-data pilot). CARRIER_UNDETERMINED used to mean "fewer than
--min-attributed-reads ATTRIBUTED reads", so a sample with hundreds of panel-aligned reads
that all tie across genera -- a property of the panel's near-identical islands, i.e. a
measured "cannot separate" -- reported the same verdict, and the same empty count
columns, as a sample with no reads and as an index/map mismatch. The summary now carries
n_mapped, n_off_map and tied_genus_sets so the three are told apart.

Strain-level and genus-level are reported separately and always both, because most short
reads over a mobile island shared at high identity across genera CANNOT separate strains,
and a strain-level answer printed without its genus-level context invites exactly the
over-reading this check exists to prevent. The reviewer's own note applies here: NM and
alignment-score margins are divergence metrics relative to the panel, never proof of
carrier identity.

    carrier_attribution.py --sample S --sam S.v8b_carrier.sam \\
        --carrier-map v8b_carrier_map.tsv \\
        --output-summary S.carrier_attribution.tsv \\
        --output-per-read S.carrier_per_read.tsv
"""
import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

FLAG_UNMAPPED = 0x4

ATTRIBUTED = "attributed"
AMBIGUOUS = "ambiguous_multi_carrier"
UNMAPPED = "unmapped"
OFF_MAP = "off_map"

CARRIER_ATTRIBUTED = "CARRIER_ATTRIBUTED"
CARRIER_GENUS_ONLY = "CARRIER_GENUS_ONLY"
CARRIER_AMBIGUOUS = "CARRIER_AMBIGUOUS"
CARRIER_UNDETERMINED = "CARRIER_UNDETERMINED"

PER_READ_FIELDS = [
    "read_id", "attribution", "best_carrier", "best_carrier_organism", "best_as",
    "runner_up_carrier", "runner_up_as", "delta_as", "carriers_within_margin",
    "tied_carriers", "best_genus",
]

SUMMARY_FIELDS = [
    "sample", "carrier_attribution", "panel_carriers", "n_reads", "n_attributed",
    "n_ambiguous_multi_carrier", "n_unmapped", "n_mapped", "n_off_map", "top_carrier", "top_carrier_organism",
    "top_carrier_reads", "top_carrier_fraction", "runner_up_carrier",
    "runner_up_carrier_reads", "top_genus", "top_genus_reads", "top_genus_fraction",
    "strain_level_resolved", "genus_level_resolved", "carrier_read_counts",
    "genus_read_counts", "tied_genus_sets", "min_delta_as", "tie_margin",
]


def read_carrier_map(path):
    """contig -> {carrier, organism, accession, genus}. Written by
    scripts/build_specificity_panels.py --panel v8b."""
    carriers = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            organism = row.get("carrier_organism", "")
            carriers[row["contig"]] = {
                "carrier": row["carrier_key"],
                "organism": organism,
                "accession": row.get("carrier_accession", ""),
                "genus": organism.split()[0] if organism else row["carrier_key"].split("_")[0],
            }
    if not carriers:
        raise ValueError(f"no contigs in {path}")
    return carriers


def parse_sam_records(path):
    """Yield (qname, flag, rname, AS) for every alignment record, SAM or BAM/CRAM."""
    suffix = Path(path).suffix.lower()
    if suffix in (".bam", ".cram"):
        proc = subprocess.run(["samtools", "view", str(path)], check=True,
                               capture_output=True, text=True)
        lines = proc.stdout.splitlines()
    else:
        with Path(path).open() as handle:
            lines = [line.rstrip("\n") for line in handle if not line.startswith("@")]
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        as_score = None
        for field in fields[11:]:
            if field.startswith("AS:"):
                as_score = int(field.split(":", 2)[2])
                break
        yield fields[0], int(fields[1]), fields[2], as_score


def best_per_carrier(records, carrier_map):
    """{carrier_key: best AS} for one read. Alignments to a contig that is not in the
    carrier map are ignored rather than guessed at: the map is the panel's definition."""
    best = {}
    for _qname, flag, rname, as_score in records:
        if flag & FLAG_UNMAPPED or as_score is None:
            continue
        entry = carrier_map.get(rname)
        if entry is None:
            continue
        key = entry["carrier"]
        if key not in best or as_score > best[key]:
            best[key] = as_score
    return best


def attribute_read(best, min_delta_as, tie_margin):
    """One read's carrier call. Returns a per-read row without the read id."""
    if not best:
        return {"attribution": UNMAPPED, "best_carrier": "", "best_as": "",
                "runner_up_carrier": "", "runner_up_as": "", "delta_as": "",
                "carriers_within_margin": 0, "tied_carriers": ""}
    ranked = sorted(best.items(), key=lambda kv: (-kv[1], kv[0]))
    top_carrier, top_as = ranked[0]
    runner_up, runner_up_as = (ranked[1] if len(ranked) > 1 else ("", None))
    delta = (top_as - runner_up_as) if runner_up_as is not None else None
    tied = [key for key, score in ranked if top_as - score <= tie_margin]
    resolved = len(tied) == 1 and (delta is None or delta >= min_delta_as)
    return {
        "attribution": ATTRIBUTED if resolved else AMBIGUOUS,
        "best_carrier": top_carrier,
        "best_as": top_as,
        "runner_up_carrier": runner_up,
        "runner_up_as": "" if runner_up_as is None else runner_up_as,
        "delta_as": "" if delta is None else delta,
        "carriers_within_margin": len(tied),
        "tied_carriers": ",".join(sorted(tied)),
    }


def attribute_reads(records, carrier_map, min_delta_as, tie_margin):
    grouped = defaultdict(list)
    for record in records:
        grouped[record[0]].append(record)
    rows = []
    for read_id in sorted(grouped):
        best = best_per_carrier(grouped[read_id], carrier_map)
        row = attribute_read(best, min_delta_as, tie_margin)
        if not best and any(not (flag & FLAG_UNMAPPED) and rname not in carrier_map
                            for _q, flag, rname, _as in grouped[read_id]):
            row["attribution"] = OFF_MAP
        row["read_id"] = read_id
        organism, genus = "", ""
        if row["best_carrier"]:
            for entry in carrier_map.values():
                if entry["carrier"] == row["best_carrier"]:
                    organism, genus = entry["organism"], entry["genus"]
                    break
        row["best_carrier_organism"] = organism
        row["best_genus"] = genus
        rows.append(row)
    return rows


def genus_of(carrier_map, carrier_key):
    for entry in carrier_map.values():
        if entry["carrier"] == carrier_key:
            return entry["genus"]
    return ""


def summarize(sample, rows, carrier_map, args):
    """Sample-level carrier attribution. Genus-level is computed from the tied sets, not
    from the strain call alone: reads that cannot separate two strains of one genus still
    support that genus, and dropping them would understate genus-level evidence."""
    n_reads = len(rows)
    attributed = [r for r in rows if r["attribution"] == ATTRIBUTED]
    ambiguous = [r for r in rows if r["attribution"] == AMBIGUOUS]
    unmapped = [r for r in rows if r["attribution"] == UNMAPPED]
    off_map = [r for r in rows if r["attribution"] == OFF_MAP]
    n_mapped = len(attributed) + len(ambiguous)

    carrier_counts = defaultdict(int)
    for row in attributed:
        carrier_counts[row["best_carrier"]] += 1

    genus_counts = defaultdict(int)
    tied_genus_sets = defaultdict(int)
    genus_attributed = 0
    for row in attributed + ambiguous:
        tied = [t for t in row["tied_carriers"].split(",") if t]
        if not tied:
            continue
        genera = {genus_of(carrier_map, key) for key in tied}
        if len(genera) == 1:
            genus_counts[genera.pop()] += 1
            genus_attributed += 1
        else:
            # Which genera a read could not separate: the reason a sample with many
            # aligned reads can still have empty genus counts.
            tied_genus_sets["+".join(sorted(genera))] += 1

    ranked_carriers = sorted(carrier_counts.items(), key=lambda kv: (-kv[1], kv[0]))
    ranked_genera = sorted(genus_counts.items(), key=lambda kv: (-kv[1], kv[0]))
    top_carrier, top_carrier_reads = (ranked_carriers[0] if ranked_carriers else ("", 0))
    runner_up, runner_up_reads = (ranked_carriers[1] if len(ranked_carriers) > 1 else ("", 0))
    top_genus, top_genus_reads = (ranked_genera[0] if ranked_genera else ("", 0))

    n_attributed = len(attributed)
    top_carrier_fraction = (top_carrier_reads / n_attributed) if n_attributed else 0.0
    top_genus_fraction = (top_genus_reads / genus_attributed) if genus_attributed else 0.0

    strain_ok = (n_attributed >= args.min_attributed_reads
                 and top_carrier_fraction >= args.min_top_carrier_fraction)
    genus_ok = (genus_attributed >= args.min_attributed_reads
                and top_genus_fraction >= args.min_top_carrier_fraction)

    if strain_ok:
        verdict = CARRIER_ATTRIBUTED
    elif genus_ok:
        verdict = CARRIER_GENUS_ONLY
    elif n_mapped < args.min_attributed_reads:
        verdict = CARRIER_UNDETERMINED
    else:
        verdict = CARRIER_AMBIGUOUS

    return {
        "sample": sample,
        "carrier_attribution": verdict,
        "panel_carriers": len({entry["carrier"] for entry in carrier_map.values()}),
        "n_reads": n_reads,
        "n_attributed": n_attributed,
        "n_ambiguous_multi_carrier": len(ambiguous),
        "n_unmapped": len(unmapped),
        "n_mapped": n_mapped,
        "n_off_map": len(off_map),
        "top_carrier": top_carrier,
        "top_carrier_organism": next(
            (e["organism"] for e in carrier_map.values() if e["carrier"] == top_carrier), ""),
        "top_carrier_reads": top_carrier_reads,
        "top_carrier_fraction": f"{top_carrier_fraction:.6f}",
        "runner_up_carrier": runner_up,
        "runner_up_carrier_reads": runner_up_reads,
        "top_genus": top_genus,
        "top_genus_reads": top_genus_reads,
        "top_genus_fraction": f"{top_genus_fraction:.6f}",
        "strain_level_resolved": strain_ok,
        "genus_level_resolved": genus_ok,
        "carrier_read_counts": ";".join(f"{k}:{v}" for k, v in ranked_carriers),
        "genus_read_counts": ";".join(f"{k}:{v}" for k, v in ranked_genera),
        "tied_genus_sets": ";".join(f"{k}:{v}" for k, v in
                                     sorted(tied_genus_sets.items(), key=lambda kv: (-kv[1], kv[0]))),
        "min_delta_as": args.min_delta_as,
        "tie_margin": args.tie_margin,
    }


def write_tsv(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--sam", required=True,
                   help="SAM/BAM: the supporting reads aligned to the concatenated V8b "
                        "carrier panel, secondary alignments enabled")
    p.add_argument("--carrier-map", required=True,
                   help="v8b_carrier_map.tsv: contig, carrier_key, carrier_organism, "
                        "carrier_accession")
    p.add_argument("--min-delta-as", type=float, default=1.0,
                   help="alignment-score margin over the runner-up carrier required to "
                        "attribute a read to one strain (default: 1)")
    p.add_argument("--tie-margin", type=float, default=0.0,
                   help="carriers within this many AS of the best are treated as tied "
                        "(default: 0, exact ties only)")
    p.add_argument("--min-attributed-reads", type=int, default=10,
                   help="below this many attributed reads the sample is "
                        "CARRIER_UNDETERMINED whatever the point estimate (default: 10)")
    p.add_argument("--min-top-carrier-fraction", type=float, default=0.8,
                   help="share of attributed reads the leading carrier (or genus) needs "
                        "for a resolved call (default: 0.8)")
    p.add_argument("--output-summary", required=True)
    p.add_argument("--output-per-read", required=True)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    carrier_map = read_carrier_map(args.carrier_map)
    rows = attribute_reads(list(parse_sam_records(args.sam)), carrier_map,
                            args.min_delta_as, args.tie_margin)
    summary = summarize(args.sample, rows, carrier_map, args)
    if summary["n_off_map"]:
        raise SystemExit(
            f"{args.sample}: {summary['n_off_map']} read(s) aligned only to contigs that "
            f"{args.carrier_map} does not name. The index and the carrier map are from "
            "different panel builds; rebuild both with scripts/build_specificity_panels.py "
            "--panel v8b rather than reporting those reads as unmapped.")
    write_tsv(args.output_per_read, PER_READ_FIELDS, rows)
    write_tsv(args.output_summary, SUMMARY_FIELDS, [summary])
    print(f"{args.sample}\t{summary['carrier_attribution']}"
          f"\ttop_carrier={summary['top_carrier'] or 'none'}"
          f" ({summary['top_carrier_reads']}/{summary['n_attributed']})"
          f"\ttop_genus={summary['top_genus'] or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
