#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


CLB_GENES = tuple(f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS")


def normalize_read_id(read_id):
    """Remove FASTA/FASTQ prefix and an optional mate suffix."""
    read_id = read_id.strip().split()[0]
    read_id = read_id.lstrip("@>")
    return re.sub(r"(?:/|\.)[12]$", "", read_id)


def parse_taxonomy(nodes_path, names_path):
    parents = {}
    ranks = {}
    scientific_names = {}

    with nodes_path.open() as handle:
        for line in handle:
            fields = [field.strip() for field in line.split("|")]

            if len(fields) < 3:
                continue

            taxid = fields[0]
            parents[taxid] = fields[1]
            ranks[taxid] = fields[2]

    with names_path.open() as handle:
        for line in handle:
            fields = [field.strip() for field in line.split("|")]

            if len(fields) < 4:
                continue

            taxid = fields[0]
            name = fields[1]
            name_class = fields[3]

            if name_class == "scientific name":
                scientific_names[taxid] = name

    return parents, ranks, scientific_names


def find_species_taxid(taxid, parents, ranks, cache):
    """Return the species ancestor for a taxonomy ID, if one exists."""
    if taxid in cache:
        return cache[taxid]

    original_taxid = taxid
    visited = set()

    while taxid and taxid not in visited:
        visited.add(taxid)

        if ranks.get(taxid) == "species":
            cache[original_taxid] = taxid
            return taxid

        parent = parents.get(taxid)

        if not parent or parent == taxid:
            break

        taxid = parent

    cache[original_taxid] = None
    return None


def read_kraken_classifications(
    kraken_path,
    parents,
    ranks,
):
    """
    Collect direct KrakenUniq classifications by read/template ID.

    If both mates are present, their classifications are collected under
    the same normalized read ID.
    """
    species_by_read = defaultdict(set)
    classified_reads = set()
    species_cache = {}

    with kraken_path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                continue

            status = fields[0]
            read_id = normalize_read_id(fields[1])
            taxid = fields[2].strip()

            if status != "C" or taxid == "0":
                continue

            classified_reads.add(read_id)

            species_taxid = find_species_taxid(
                taxid,
                parents,
                ranks,
                species_cache,
            )

            if species_taxid is not None:
                species_by_read[read_id].add(species_taxid)

    return species_by_read, classified_reads


def build_matrix(
    read_gene_path,
    species_by_read,
    classified_reads,
    scientific_names,
):
    matrix = defaultdict(lambda: defaultdict(int))

    with read_gene_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {"read_id", "Gene"}
        observed_columns = set(reader.fieldnames or [])

        if not required_columns.issubset(observed_columns):
            raise ValueError(
                f"{read_gene_path} must contain columns: "
                "read_id and Gene"
            )

        for row in reader:
            read_id = normalize_read_id(row["read_id"])
            gene = row["Gene"].strip()

            if gene not in CLB_GENES:
                continue

            species_taxids = species_by_read.get(read_id, set())

            if len(species_taxids) == 1:
                species_taxid = next(iter(species_taxids))
                species_name = scientific_names.get(
                    species_taxid,
                    f"taxid_{species_taxid}",
                )
                row_key = (species_name, species_taxid)

            elif len(species_taxids) > 1:
                # Conflicting species calls between records or mates
                row_key = ("Conflicting_species", "-1")

            elif read_id in classified_reads:
                # Classified by KrakenUniq, but not specifically enough
                # to resolve a species ancestor
                row_key = ("Unresolved_at_species", "-1")

            else:
                row_key = ("Unclassified", "0")

            matrix[row_key][gene] += 1

    return matrix


def write_matrix(matrix, output_path):
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Species", "TaxID", *CLB_GENES, "Total"])

        for species_name, taxid in sorted(
            matrix,
            key=lambda value: (value[0].lower(), value[1]),
        ):
            counts = matrix[(species_name, taxid)]
            gene_counts = [counts.get(gene, 0) for gene in CLB_GENES]

            writer.writerow(
                [
                    species_name,
                    taxid,
                    *gene_counts,
                    sum(gene_counts),
                ]
            )


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Join read-to-clb-gene assignments with direct KrakenUniq "
            "classifications and create a species-by-clb-gene matrix."
        )
    )

    parser.add_argument(
        "--read-gene",
        required=True,
        type=Path,
        help="TSV containing read_id and Gene columns.",
    )

    parser.add_argument(
        "--kraken-output",
        required=True,
        type=Path,
        help="Per-read KrakenUniq output file.",
    )

    parser.add_argument(
        "--taxonomy-dir",
        required=True,
        type=Path,
        help="Kraken database taxonomy directory.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output species-by-clb-gene TSV.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    nodes_path = args.taxonomy_dir / "nodes.dmp"
    names_path = args.taxonomy_dir / "names.dmp"

    if not args.read_gene.is_file():
        raise FileNotFoundError(
            f"Read-to-gene table not found: {args.read_gene}"
        )

    if not args.kraken_output.is_file():
        raise FileNotFoundError(
            f"KrakenUniq output not found: {args.kraken_output}"
        )

    if not nodes_path.is_file():
        raise FileNotFoundError(
            f"Taxonomy nodes file not found: {nodes_path}"
        )

    if not names_path.is_file():
        raise FileNotFoundError(
            f"Taxonomy names file not found: {names_path}"
        )

    parents, ranks, scientific_names = parse_taxonomy(
        nodes_path,
        names_path,
    )

    species_by_read, classified_reads = read_kraken_classifications(
        args.kraken_output,
        parents,
        ranks,
    )

    matrix = build_matrix(
        args.read_gene,
        species_by_read,
        classified_reads,
        scientific_names,
    )

    write_matrix(matrix, args.output)


if __name__ == "__main__":
    main()