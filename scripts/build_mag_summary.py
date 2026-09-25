#!/usr/bin/env python3
"""Build pks_mag_summary.tsv from per-sample MAG pipeline outputs."""

import argparse
import csv
import glob
import os
import sys
from pathlib import Path

from mag_utils import parse_hmmsearch_tblout, read_locus_evidence

ENTEROBACTERALES = "o__Enterobacterales"

# U1: a bin is a genome the pipeline actually reconstructed and can attribute pks
# evidence to; "unbinned" is the per-sample pool of contigs MetaBAT2 never placed in
# any bin. Confirmed missing before this: a contig can carry real pks-island alignment
# evidence and never surface anywhere once binning declines to place it -- common when
# abundance is low or the community is underpowered for composition-based binning --
# even though the read-level pksProfilerAlign/pksProfilerHMM lane may still show signal
# for the same sample. unit_type is how a reader (and build_master_summary.py) tells
# the two apart without guessing from bin_id alone: "positive signal in unbinned
# contigs, not attributable to any recovered genome" is a materially different claim
# from a bin-level positive call, and must never be conflated with one.
UNIT_TYPE_BIN = "bin"
UNIT_TYPE_UNBINNED = "unbinned"
UNBINNED_UNIT_ID = "unbinned"

FIELDNAMES = [
    "sample", "bin_id", "unit_type", "taxonomy", "completeness", "contamination",
    "genome_size", "contig_n50", "clb_genes_detected", "clb_genes", "best_evalue",
    "has_integrase", "has_transposase", "flanking_genes", "unexpected_taxon_flag",
    "locus_tier", "locus_genes_detected", "locus_breadth",
]


def is_unexpected_taxon(taxonomy):
    """True only when taxonomy is classified but not Enterobacterales.
    Unclassified bins are unknown, not anomalous."""
    return taxonomy not in ("unclassified", "") and ENTEROBACTERALES not in taxonomy


# CheckM2 writes these, capitalised exactly so -- see checkm2/predictQuality.py,
# which sets final_results['Genome_Size'] and final_results['Contig_N50'].
CHECKM2_REQUIRED = ("Name", "Completeness", "Contamination", "Genome_Size", "Contig_N50")


def parse_checkm2(tsv_path):
    """
    M4: this used `row.get("Genome_size", 0)` -- lowercase s, which CheckM2 has never
    written. The lookup never matched and the default turned a schema mismatch into a
    measurement: genome_size was 0 for all 8 bins of the v0.0.2 test, which also removed
    the cheapest sanity check on a bin (an E. coli MAG should be ~4.6 Mb, not 200 kb).

    A missing column now raises. A quality report we cannot read is a broken run, not a
    genome of size zero.
    """
    result = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        observed = list(reader.fieldnames or [])
        missing = [column for column in CHECKM2_REQUIRED if column not in observed]
        if missing:
            raise ValueError(
                f"CheckM2 report {tsv_path} is missing {', '.join(missing)}. "
                f"Columns present: {', '.join(observed) or 'none'}"
            )
        for row in reader:
            result[row["Name"]] = {
                "completeness": float(row["Completeness"]),
                "contamination": float(row["Contamination"]),
                "genome_size": int(row["Genome_Size"]),
                "contig_n50": int(row["Contig_N50"]),
            }
    return result


def parse_gtdbtk(summary_path):
    result = {}
    with open(summary_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            result[row["user_genome"]] = row.get("classification", "unclassified") or "unclassified"
    return result


def parse_context(context_path):
    has_int, has_tra, flank_set = False, False, set()
    with open(context_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            has_int |= row["has_integrase"].lower() == "true"
            has_tra |= row["has_transposase"].lower() == "true"
            if row.get("flanking_genes"):
                for g in row["flanking_genes"].split(";"):
                    if g.strip():
                        flank_set.add(g.strip())
    return has_int, has_tra, ";".join(sorted(flank_set))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkm2", required=True)
    parser.add_argument("--gtdbtk", required=True)
    parser.add_argument("--tblout_dir", required=True)
    parser.add_argument("--context_dir", required=True)
    parser.add_argument("--locus_dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    checkm2 = parse_checkm2(args.checkm2)
    gtdbtk = parse_gtdbtk(args.gtdbtk)
    locus_evidence = read_locus_evidence(args.locus_dir)

    rows = []
    for tblout_path in sorted(glob.glob(os.path.join(args.tblout_dir, "*.tblout"))):
        bin_id = Path(tblout_path).stem
        hits = parse_hmmsearch_tblout(tblout_path, args.evalue)
        if not hits:
            continue

        clb_genes_seen = {}
        for h in hits:
            gene = h["clb_gene"]
            clb_genes_seen[gene] = min(h["evalue"], clb_genes_seen.get(gene, float("inf")))
        clb_genes = sorted(clb_genes_seen.keys())
        best_evalue = min(clb_genes_seen.values())

        # A bin absent from the report is genuinely unmeasured; NA says so, where 0
        # would read as a measurement.
        qc = checkm2.get(bin_id, {"completeness": 0.0, "contamination": 0.0,
                                  "genome_size": "NA", "contig_n50": "NA"})
        taxonomy = gtdbtk.get(bin_id, "unclassified")
        unexpected = is_unexpected_taxon(taxonomy)
        # M1: the alignment-confirmed call for this bin -- see mag_utils.read_locus_evidence.
        # Absent means the bin was never aligned (should not happen once magBinLocusEvidence
        # runs for every bin), not that it failed the tier.
        locus = locus_evidence.get(bin_id, {"locus_tier": "NA", "locus_genes_detected": "NA",
                                             "locus_breadth": "NA"})

        context_path = os.path.join(args.context_dir, f"{bin_id}.context.tsv")
        if os.path.exists(context_path):
            has_int, has_tra, flanking = parse_context(context_path)
        else:
            has_int, has_tra, flanking = False, False, ""

        rows.append({
            "sample": args.sample,
            "bin_id": bin_id,
            "unit_type": UNIT_TYPE_BIN,
            "taxonomy": taxonomy,
            "completeness": qc["completeness"],
            "contamination": qc["contamination"],
            "genome_size": qc["genome_size"],
            "contig_n50": qc["contig_n50"],
            "clb_genes_detected": len(clb_genes),
            "clb_genes": ",".join(clb_genes),
            "best_evalue": best_evalue,
            "has_integrase": has_int,
            "has_transposase": has_tra,
            "flanking_genes": flanking,
            "unexpected_taxon_flag": unexpected,
            "locus_tier": locus["locus_tier"],
            "locus_genes_detected": locus["locus_genes_detected"],
            "locus_breadth": locus["locus_breadth"],
        })

    # U1: the unbinned pool is not driven by a tblout (no Prokka/hmmsearch runs on it --
    # see Modules/pks_mag.nf's pksMAG workflow), so it is never picked up by the tblout
    # loop above. Reported only when magBinLocusEvidence actually ran against it (i.e.
    # the sample had something left unbinned); a bin-shaped row with every bin-specific
    # field explicitly NA, not 0 or False, since none of those were measured for a pool
    # of contigs that is not a genome.
    unbinned = locus_evidence.get(UNBINNED_UNIT_ID)
    if unbinned is not None:
        rows.append({
            "sample": args.sample,
            "bin_id": UNBINNED_UNIT_ID,
            "unit_type": UNIT_TYPE_UNBINNED,
            "taxonomy": "NA",
            "completeness": "NA",
            "contamination": "NA",
            "genome_size": "NA",
            "contig_n50": "NA",
            "clb_genes_detected": "NA",
            "clb_genes": "",
            "best_evalue": "NA",
            "has_integrase": "NA",
            "has_transposase": "NA",
            "flanking_genes": "",
            "unexpected_taxon_flag": "NA",
            "locus_tier": unbinned["locus_tier"],
            "locus_genes_detected": unbinned["locus_genes_detected"],
            "locus_breadth": unbinned["locus_breadth"],
        })

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    bin_rows = sum(1 for r in rows if r["unit_type"] == UNIT_TYPE_BIN)
    if not rows:
        print(f"No pks+ bins found for {args.sample}", file=sys.stderr)
    else:
        note = " (plus an unbinned-pool row)" if unbinned is not None else ""
        print(f"Wrote {bin_rows} pks+ bins{note} to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
