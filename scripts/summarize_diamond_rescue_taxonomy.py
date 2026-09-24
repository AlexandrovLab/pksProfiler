#!/usr/bin/env python3
"""Join DIAMOND-rescued reads back to the organism the fast classifier called them.

select_diamond_reads.py decides which non-target reads to keep by DIAMOND identity
alone, then writes only the read ID onward -- it keeps neither `sseqid` (which clb gene
the read hit) nor the read's own krakenPrefilter taxid, so a rescued read's protein
homology and its community identity were computed and then both discarded. This
produces the table that was missing: per organism, per clb gene, how many rescued reads.
The organism is the krakenuniq per-read call from `krakenPrefilter` (Modules/pks_prefilter.nf),
resolved to its species ancestor with the same taxDB/NCBI-dump logic
build_clb_species_matrix.py already uses for the --pks_taxa lane -- reused here, not
re-derived, so "Unresolved_at_species" and "Unclassified" mean the same thing in both
tables.

Every row in this table is, by construction, a read the fast classifier's target-taxon
routing did *not* keep (krakenPrefilter only sends primary-taxon reads to profiling;
everything else lands in non_target, and this table is drawn from what DIAMOND rescued
out of that pile) -- so the table itself answers "was this caught by the fast classifier
or only the sensitive DIAMOND rescue?": if a read is in here, it was only DIAMOND.

    summarize_diamond_rescue_taxonomy.py --sample S1 --diamond diamond.tsv \\
        --kraken-output krakenuniq.output.txt --kraken-db DB --output out.tsv
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_clb_species_matrix import find_species_taxid, load_taxonomy, normalize_read_id

CLB_GENE_RE = re.compile(r"^clb[A-Za-z]$")

FIELDS = ["sample", "organism", "taxid", "clb_gene", "read_count"]


def clb_gene_from_sseqid(sseqid):
    """The reference protein headers are `id1|id2|clbX|score=...|hmmcov=...`."""
    for token in sseqid.split("|"):
        if CLB_GENE_RE.match(token):
            return token
    return "unknown"


def load_diamond_hits(path, max_evalue, min_aa, min_query_coverage):
    """Best passing hit per read: {read_id: clb_gene}.

    Same filters as select_diamond_reads.py's load_hits(), so a read counted here is a
    read select_diamond_reads.py actually rescued -- but sseqid is kept, and a read with
    more than one passing hit keeps the lowest-evalue one rather than the last one seen.
    """
    best = {}
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 8:
                raise ValueError(f"Malformed DIAMOND row {number}")
            query, sseqid, _pident, length, qlen, _slen, evalue, _bitscore = fields
            coverage = min(1.0, int(length) * 3 / int(qlen))
            if float(evalue) > max_evalue or int(length) < min_aa or coverage < min_query_coverage:
                continue
            read_id = normalize_read_id(query)
            evalue_f = float(evalue)
            if read_id not in best or evalue_f < best[read_id][1]:
                best[read_id] = (clb_gene_from_sseqid(sseqid), evalue_f)
    return {read_id: gene for read_id, (gene, _evalue) in best.items()}


def load_read_taxids(path):
    """{read_id: taxid} for classified reads, from krakenPrefilter's per-read output.

    Columns as in build_clb_species_matrix.read_kraken_classifications: status, read
    ID, assigned taxid.
    """
    taxids = {}
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed Kraken row {number}")
            if fields[0] != "C":
                continue
            taxid = fields[2].strip().split()[0]
            if taxid == "0":
                continue
            taxids[normalize_read_id(fields[1])] = taxid
    return taxids


def build_table(gene_by_read, taxid_by_read, parents, ranks, scientific_names):
    """{(organism, taxid, clb_gene): read_count}, one row per rescued read."""
    counts = defaultdict(int)
    species_cache = {}
    for read_id, gene in gene_by_read.items():
        taxid = taxid_by_read.get(read_id)
        if taxid is None:
            key = ("Unclassified", "0", gene)
        else:
            species_taxid = find_species_taxid(taxid, parents, ranks, species_cache)
            if species_taxid is None:
                key = ("Unresolved_at_species", "-1", gene)
            else:
                key = (scientific_names.get(species_taxid, f"taxid_{species_taxid}"),
                       species_taxid, gene)
        counts[key] += 1
    return counts


def write_table(sample, counts, output_path):
    with Path(output_path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t")
        writer.writeheader()
        for (organism, taxid, gene), count in sorted(
                counts.items(), key=lambda kv: (kv[0][0].lower(), kv[0][1], kv[0][2])):
            writer.writerow(dict(sample=sample, organism=organism, taxid=taxid,
                                 clb_gene=gene, read_count=count))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--diamond", required=True, type=Path)
    p.add_argument("--kraken-output", required=True, type=Path)
    p.add_argument("--kraken-db", required=True, type=Path)
    p.add_argument("--max-evalue", type=float, default=1e-5)
    p.add_argument("--min-aa", type=int, default=25)
    p.add_argument("--min-query-coverage", type=float, default=0.5)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()

    gene_by_read = load_diamond_hits(a.diamond, a.max_evalue, a.min_aa, a.min_query_coverage)
    if not gene_by_read:
        write_table(a.sample, {}, a.output)
        return

    taxid_by_read = load_read_taxids(a.kraken_output)
    parents, ranks, scientific_names = load_taxonomy(a.kraken_db)
    counts = build_table(gene_by_read, taxid_by_read, parents, ranks, scientific_names)
    write_table(a.sample, counts, a.output)


if __name__ == "__main__":
    main()
