#!/usr/bin/env python3
"""Summarize predicted prophages and pks-producer/lysogen co-occurrence."""

import argparse
import csv
import glob
import os
import re
from pathlib import Path

from mag_utils import parse_hmmsearch_tblout, SPECIFIC_CLB

PROPHAGE_FIELDS = [
    "sample", "bin_id", "host_taxonomy", "prophage_id", "host_contig",
    "start", "end", "length", "virus_score", "viral_taxonomy",
    "has_recA", "has_lexA", "clbS_like", "evidence_level",
]
INTERACTION_FIELDS = [
    "sample", "pks_producer_bin", "producer_taxonomy", "recipient_bin",
    "recipient_taxonomy", "recipient_prophage_count", "recipient_has_recA",
    "recipient_has_lexA", "recipient_clbS_like", "susceptibility_hypothesis",
    "interpretation",
]
MOBILITY_FIELDS = [
    "sample", "bin_id", "host_taxonomy", "distinct_clb_genes", "clb_genes",
    "provisional_pks_class", "biosynthetic_clb_detected", "specific_clb_detected",
    "clbP_detected", "clbS_detected", "clb_contig_count", "same_contig_clb_span",
    "nearby_integrase", "nearby_trna", "mobility_interpretation",
]
BIOSYNTHETIC_CLB = {"clbB", "clbC", "clbH", "clbI", "clbJ", "clbK", "clbN", "clbO"}


def parse_gtdbtk(path):
    with open(path) as fh:
        return {
            row["user_genome"]: row.get("classification", "") or "unclassified"
            for row in csv.DictReader(fh, delimiter="\t")
        }


def annotation_flags(gff_path):
    rec_a = lex_a = False
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = dict(re.findall(r"(\w+)=([^;]+)", fields[8]))
            gene = attrs.get("gene", "").lower()
            rec_a |= gene == "reca"
            lex_a |= gene == "lexa"
    return rec_a, lex_a


def clb_calls(tblout_path, evalue):
    genes = {hit["clb_gene"] for hit in parse_hmmsearch_tblout(tblout_path, evalue)}
    return genes, "clbS" in genes


def mobility_evidence(gff_path, tblout_path, evalue):
    hits = parse_hmmsearch_tblout(tblout_path, evalue)
    hit_tags = {hit["locus_tag"] for hit in hits}
    clb_cds, integrases, trnas = [], [], []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = dict(re.findall(r"(\w+)=([^;]+)", fields[8]))
            feature = {
                "contig": fields[0], "start": int(fields[3]), "end": int(fields[4]),
                "id": attrs.get("ID", attrs.get("locus_tag", "")),
                "text": (attrs.get("gene", "") + " " + attrs.get("product", "")).lower(),
            }
            if fields[2] == "CDS" and feature["id"] in hit_tags:
                clb_cds.append(feature)
            mobility_recombinase = any(term in feature["text"] for term in (
                "integrase", "site-specific recombinase", "phage recombinase",
            ))
            if fields[2] == "CDS" and mobility_recombinase:
                integrases.append(feature)
            if fields[2].lower() in ("trna", "rna") and "trna" in feature["text"]:
                trnas.append(feature)

    def close_to_clb(features, window):
        return any(
            feature["contig"] == clb["contig"]
            and feature["start"] <= clb["end"] + window
            and feature["end"] >= clb["start"] - window
            for feature in features for clb in clb_cds
        )

    contigs = {feature["contig"] for feature in clb_cds}
    span = ""
    if len(contigs) == 1 and clb_cds:
        span = max(x["end"] for x in clb_cds) - min(x["start"] for x in clb_cds) + 1
    return {
        "clb_contig_count": len(contigs), "same_contig_clb_span": span,
        "nearby_integrase": close_to_clb(integrases, 50000),
        "nearby_trna": close_to_clb(trnas, 10000),
    }

def read_prophages(path, min_length=3000, min_hallmarks=2):
    """Provirus intervals above the length and hallmark floors.

    T3: geNomad called 134 "viral contigs" on AA-3850, the top hits 369 bp with one
    gene and one hallmark. That count tracked how fragmented the assembly was, not
    biology. A provirus must now clear a length and a hallmark-count floor to be
    reported; a short interval with a single hallmark is assembly noise, and calling
    it a prophage put mobility claims on top of it.
    """
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("topology", "").lower() != "provirus":
                continue
            try:
                length_bp = int(float(row.get("length") or 0))
                hallmarks = int(float(row.get("n_hallmarks") or 0))
            except ValueError:
                continue
            if length_bp < min_length or hallmarks < min_hallmarks:
                continue
            match = re.fullmatch(r"(\d+)-(\d+)", row.get("coordinates", ""))
            if not match:
                continue
            start, end = sorted(map(int, match.groups()))
            seq_name = row.get("seq_name", "")
            rows.append({
                "prophage_id": seq_name,
                "host_contig": seq_name.split("|provirus_", 1)[0],
                "start": start,
                "end": end,
                "length": row.get("length", end - start + 1),
                "virus_score": row.get("virus_score", ""),
                "viral_taxonomy": row.get("taxonomy", ""),
            })
    return rows


def write_tsv(path, fields, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--gtdbtk", required=True)
    parser.add_argument("--gff-dir", required=True)
    parser.add_argument("--tblout-dir", required=True)
    parser.add_argument("--genomad-dir", required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--min-clb-genes", type=int, default=3)
    parser.add_argument("--min-provirus-length", type=int, default=3000)
    parser.add_argument("--min-provirus-hallmarks", type=int, default=2)
    parser.add_argument("--prophage-out", required=True)
    parser.add_argument("--interaction-out", required=True)
    parser.add_argument("--mobility-out", required=True)
    args = parser.parse_args()

    taxonomy = parse_gtdbtk(args.gtdbtk)
    bins = {}
    for tblout in glob.glob(os.path.join(args.tblout_dir, "*.tblout")):
        bin_id = Path(tblout).stem
        genes, clb_s = clb_calls(tblout, args.evalue)
        gff = os.path.join(args.gff_dir, f"{bin_id}.gff")
        rec_a, lex_a = annotation_flags(gff)
        summary = os.path.join(args.genomad_dir, f"{bin_id}.virus_summary.tsv")
        bins[bin_id] = {
            "taxonomy": taxonomy.get(bin_id, "unclassified"),
            "clb_genes": genes, "clbS_like": clb_s,
            "recA": rec_a, "lexA": lex_a, "prophages": read_prophages(summary, args.min_provirus_length,
                                        args.min_provirus_hallmarks),
            "mobility": mobility_evidence(gff, tblout, args.evalue),
        }

    inventory = []
    for bin_id, data in sorted(bins.items()):
        for phage in data["prophages"]:
            inventory.append({
                "sample": args.sample, "bin_id": bin_id,
                "host_taxonomy": data["taxonomy"], **phage,
                "has_recA": data["recA"], "has_lexA": data["lexA"],
                "clbS_like": data["clbS_like"],
                "evidence_level": "predicted_provirus",
            })

    interactions = []
    # M1: also require a specific, low-homology gene -- BIOSYNTHETIC_CLB alone is the
    # promiscuous megasynthase domains that a bare E-value cut lets any bacterium hit
    # (v0.0.2_functional_test_20260910, ERR525841: 7 of 8 bins otherwise qualified,
    # including two Bifidobacterium bins).
    producers = {
        key: val for key, val in bins.items()
        if len(val["clb_genes"]) >= args.min_clb_genes
        and val["clb_genes"] & BIOSYNTHETIC_CLB
        and val["clb_genes"] & SPECIFIC_CLB
    }
    recipients = {key: val for key, val in bins.items() if val["prophages"]}
    for producer_id, producer in sorted(producers.items()):
        for recipient_id, recipient in sorted(recipients.items()):
            protected = recipient["clbS_like"]
            sos_complete = recipient["recA"] and recipient["lexA"]
            hypothesis = "lower" if protected else ("plausible" if sos_complete else "uncertain")
            interactions.append({
                "sample": args.sample, "pks_producer_bin": producer_id,
                "producer_taxonomy": producer["taxonomy"],
                "recipient_bin": recipient_id,
                "recipient_taxonomy": recipient["taxonomy"],
                "recipient_prophage_count": len(recipient["prophages"]),
                "recipient_has_recA": recipient["recA"],
                "recipient_has_lexA": recipient["lexA"],
                "recipient_clbS_like": protected,
                "susceptibility_hypothesis": hypothesis,
                "interpretation": "co-occurrence_only_not_evidence_of_induction",
            })

    mobility_rows = []
    for bin_id, data in sorted(bins.items()):
        genes = data["clb_genes"]
        if not genes:
            continue
        count = len(genes)
        if genes == {"clbS"}:
            provisional = "clbS_only_possible_resistance"
        elif count < args.min_clb_genes:
            provisional = "isolated_or_low_evidence_clb_hits"
        elif not (genes & SPECIFIC_CLB):
            # M1: clears the gene-count floor on megasynthase homology alone
            # (clbB/clbC/clbH/clbI/clbJ/clbK/clbN/clbO) with none of the specific,
            # low-homology genes -- domain cross-reactivity, not a confirmed carrier.
            provisional = "megasynthase_only_low_specificity"
        elif count < 15:
            provisional = "partial_pks_candidate"
        elif count < 19:
            provisional = "near_complete_pks_candidate"
        else:
            provisional = "complete_pks_candidate"
        mobility = data["mobility"]
        hgt_markers = mobility["nearby_integrase"] or mobility["nearby_trna"]
        mobility_rows.append({
            "sample": args.sample, "bin_id": bin_id,
            "host_taxonomy": data["taxonomy"],
            "distinct_clb_genes": count,
            "clb_genes": ",".join(sorted(genes)),
            "provisional_pks_class": provisional,
            "biosynthetic_clb_detected": bool(genes & BIOSYNTHETIC_CLB),
            "specific_clb_detected": bool(genes & SPECIFIC_CLB),
            "clbP_detected": "clbP" in genes,
            "clbS_detected": "clbS" in genes,
            **mobility,
            "mobility_interpretation": (
                "HGT_supporting_features_not_direct_transfer_evidence"
                if hgt_markers else "no_local_HGT_marker_detected"
            ),
        })

    write_tsv(args.mobility_out, MOBILITY_FIELDS, mobility_rows)
    write_tsv(args.prophage_out, PROPHAGE_FIELDS, inventory)
    write_tsv(args.interaction_out, INTERACTION_FIELDS, interactions)


if __name__ == "__main__":
    main()
