#!/usr/bin/env python3
"""Extract genomic context around clb gene hits in a pks+ MAG bin or tumour contigs.

Reports, per clb hit: mobility markers in the flanking window (integrase,
transposase, tRNA) and whether the hit falls inside a geNomad-predicted provirus.
A tRNA in the flank matters because the canonical pks island integrates at a tRNA
locus; a clb gene inside a provirus interval is direct evidence the island is mobile."""

import argparse
import csv
import re
import sys

from mag_utils import parse_hmmsearch_tblout

_INTEGRASE_KW  = {"integrase", "recombinase", "resolvase"}
_TRANSPOSASE_KW = {"transposase", "insertion element", "is element"}


def parse_prokka_gff(gff_path):
    genes = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attr = dict(re.findall(r'(\w+)=([^;]+)', parts[8]))
            genes.append({
                "contig": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "locus_tag": attr.get("ID", attr.get("locus_tag", "")),
                "gene": attr.get("gene", ""),
                "product": attr.get("product", ""),
            })
    return genes


def parse_prokka_rna(gff_path):
    """tRNA features, which parse_prokka_gff skips because it keeps CDS only."""
    rnas = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2].lower() not in ("trna", "rna", "trna_gene"):
                continue
            if parts[2].lower() != "trna" and "trna" not in parts[8].lower():
                continue
            rnas.append({"contig": parts[0], "start": int(parts[3]), "end": int(parts[4])})
    return rnas


def parse_genomad_proviruses(path):
    """geNomad provirus intervals keyed by host contig.

    Follows build_community_prophage.py: keep rows whose topology is `provirus`, read the
    interval from `coordinates`, and recover the host contig by stripping the
    `|provirus_...` suffix from `seq_name`.
    """
    proviruses = {}
    if not path:
        return proviruses
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("topology") or "").lower() != "provirus":
                continue
            match = re.fullmatch(r"(\d+)-(\d+)", (row.get("coordinates") or "").strip())
            if not match:
                continue
            seq_name = row.get("seq_name", "")
            proviruses.setdefault(seq_name.split("|provirus_", 1)[0], []).append({
                "prophage_id": seq_name,
                "start": int(match.group(1)),
                "end": int(match.group(2)),
                "virus_score": row.get("virus_score", ""),
                "taxonomy": row.get("taxonomy", ""),
            })
    return proviruses


def overlapping_provirus(proviruses, contig, start, end):
    for provirus in proviruses.get(contig, []):
        if start <= provirus["end"] and end >= provirus["start"]:
            return provirus
    return None


def in_window(features, contig, pos, window):
    return any(
        f["contig"] == contig and f["start"] <= pos + window and f["end"] >= pos - window
        for f in features
    )


def get_flanking(genes, contig, pos, window=50000):
    return [
        g for g in genes
        if g["contig"] == contig
        and g["start"] <= pos + window
        and g["end"] >= pos - window
    ]


def _has_keyword(gene, keywords):
    text = (gene["gene"] + " " + gene["product"]).lower()
    return any(kw in text for kw in keywords)


def is_integrase(gene):
    return _has_keyword(gene, _INTEGRASE_KW)


def is_transposase(gene):
    return _has_keyword(gene, _TRANSPOSASE_KW)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--tblout", required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--window", type=int, default=50000)
    parser.add_argument("--genomad", default=None,
                        help="geNomad virus_summary.tsv; enables provirus overlap columns")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    genes = parse_prokka_gff(args.gff)
    rnas = parse_prokka_rna(args.gff)
    proviruses = parse_genomad_proviruses(args.genomad)
    hits = parse_hmmsearch_tblout(args.tblout, args.evalue)

    tag_to_gene = {g["locus_tag"]: g for g in genes}

    rows = []
    for hit in hits:
        anchor = tag_to_gene.get(hit["locus_tag"])
        mid = (anchor["start"] + anchor["end"]) // 2 if anchor else 0
        contig = anchor["contig"] if anchor else ""
        flanking = get_flanking(genes, contig, mid, args.window)
        has_int = any(is_integrase(g) for g in flanking)
        has_tra = any(is_transposase(g) for g in flanking)
        has_trna = in_window(rnas, contig, mid, args.window)
        provirus = overlapping_provirus(
            proviruses, contig,
            anchor["start"] if anchor else 0,
            anchor["end"] if anchor else 0,
        )
        flank_names = ";".join(
            g["gene"] or g["product"]
            for g in flanking
            if g["locus_tag"] != hit["locus_tag"]
        )
        rows.append({
            "locus_tag": hit["locus_tag"],
            "clb_gene": hit["clb_gene"],
            "evalue": hit["evalue"],
            "contig": contig,
            "has_integrase": has_int,
            "has_transposase": has_tra,
            "nearby_trna": has_trna,
            "in_prophage": provirus is not None,
            "prophage_id": provirus["prophage_id"] if provirus else "",
            "prophage_virus_score": provirus["virus_score"] if provirus else "",
            "prophage_taxonomy": provirus["taxonomy"] if provirus else "",
            "flanking_genes": flank_names,
        })

    fieldnames = ["locus_tag", "clb_gene", "evalue", "contig",
                  "has_integrase", "has_transposase", "nearby_trna",
                  "in_prophage", "prophage_id", "prophage_virus_score",
                  "prophage_taxonomy", "flanking_genes"]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} context rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
