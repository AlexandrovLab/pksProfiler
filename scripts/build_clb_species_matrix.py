#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


CLB_GENES = tuple(
    f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"
)


def normalize_read_id(read_id):
    """Normalize FASTA/FASTQ read identifiers for joining."""
    read_id = read_id.strip().split()[0]
    read_id = read_id.lstrip("@>")

    return re.sub(
        r"(?:/|\.)[12]$",
        "",
        read_id,
    )


def parse_ncbi_taxonomy(nodes_path, names_path):
    """Parse NCBI nodes.dmp and names.dmp files."""
    parents = {}
    ranks = {}
    scientific_names = {}

    with nodes_path.open() as handle:
        for line in handle:
            fields = [
                field.strip()
                for field in line.split("|")
            ]

            if len(fields) < 3:
                continue

            taxid = fields[0]
            parent_taxid = fields[1]
            rank = fields[2]

            parents[taxid] = parent_taxid
            ranks[taxid] = rank

    with names_path.open() as handle:
        for line in handle:
            fields = [
                field.strip()
                for field in line.split("|")
            ]

            if len(fields) < 4:
                continue

            taxid = fields[0]
            scientific_name = fields[1]
            name_class = fields[3]

            if name_class == "scientific name":
                scientific_names[taxid] = scientific_name

    return parents, ranks, scientific_names


def parse_krakenuniq_taxdb(taxdb_path):
    """
    Parse either supported KrakenUniq taxDB column order:

      taxid, parent taxid, scientific name, rank
      taxid, parent taxid, rank, scientific name
    """
    parents = {}
    ranks = {}
    scientific_names = {}

    known_ranks = {
        "no rank",
        "superkingdom",
        "kingdom",
        "subkingdom",
        "phylum",
        "subphylum",
        "class",
        "subclass",
        "infraclass",
        "cohort",
        "superorder",
        "order",
        "suborder",
        "infraorder",
        "parvorder",
        "superfamily",
        "family",
        "subfamily",
        "tribe",
        "subtribe",
        "genus",
        "subgenus",
        "species group",
        "species subgroup",
        "species",
        "subspecies",
        "varietas",
        "forma",
        "strain",
        "isolate",
        "clade",
    }

    with taxdb_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 4:
                continue

            taxid = fields[0].strip()
            parent_taxid = fields[1].strip()
            third_field = fields[2].strip()
            fourth_field = fields[3].strip()

            third_is_rank = third_field.lower() in known_ranks
            fourth_is_rank = fourth_field.lower() in known_ranks

            if fourth_is_rank and not third_is_rank:
                # TSCC layout: taxid, parent, scientific name, rank
                scientific_name = third_field
                rank = fourth_field
            elif third_is_rank and not fourth_is_rank:
                # Alternate layout: taxid, parent, rank, scientific name
                rank = third_field
                scientific_name = fourth_field
            elif fourth_is_rank:
                # Preserve the TSCC layout if both fields are ambiguous.
                scientific_name = third_field
                rank = fourth_field
            else:
                raise ValueError(
                    f"Cannot determine taxDB column order at "
                    f"{taxdb_path}:{line_number}: {line.rstrip()}"
                )

            parents[taxid] = parent_taxid
            ranks[taxid] = rank
            scientific_names[taxid] = scientific_name

    return parents, ranks, scientific_names


def load_taxonomy(kraken_db):
    """
    Load either an NCBI-style taxonomy directory or KrakenUniq taxDB.
    """
    taxonomy_dir = kraken_db / "taxonomy"
    nodes_path = taxonomy_dir / "nodes.dmp"
    names_path = taxonomy_dir / "names.dmp"
    taxdb_path = kraken_db / "taxDB"

    if nodes_path.is_file() and names_path.is_file():
        return parse_ncbi_taxonomy(
            nodes_path,
            names_path,
        )

    if taxdb_path.is_file():
        return parse_krakenuniq_taxdb(taxdb_path)

    raise FileNotFoundError(
        "No supported taxonomy files were found. Expected either "
        f"{nodes_path} and {names_path}, or {taxdb_path}."
    )


def find_species_taxid(
    taxid,
    parents,
    ranks,
    cache,
):
    """Find the species ancestor of a taxonomy ID."""
    if taxid in cache:
        return cache[taxid]

    original_taxid = taxid
    visited = set()

    while taxid and taxid not in visited:
        visited.add(taxid)

        if ranks.get(taxid, "").lower() == "species":
            cache[original_taxid] = taxid
            return taxid

        parent_taxid = parents.get(taxid)

        if not parent_taxid or parent_taxid == taxid:
            break

        taxid = parent_taxid

    cache[original_taxid] = None
    return None


def read_kraken_classifications(
    kraken_path,
    parents,
    ranks,
):
    """
    Collect direct KrakenUniq classifications by read/template ID.

    KrakenUniq output uses:
      column 1: classification status
      column 2: read identifier
      column 3: assigned taxonomy ID
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

            status = fields[0].strip()
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
                species_by_read[read_id].add(
                    species_taxid
                )

    return species_by_read, classified_reads


def build_matrix(
    read_gene_path,
    species_by_read,
    classified_reads,
    scientific_names,
):
    """
    Build a species-by-clb-gene count matrix.

    A read/template contributes once to the clb gene selected by
    the read-to-gene mapping.
    """
    matrix = defaultdict(
        lambda: defaultdict(int)
    )

    with read_gene_path.open() as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        expected_columns = {
            "read_id",
            "Gene",
        }

        observed_columns = set(
            reader.fieldnames or []
        )

        if not expected_columns.issubset(
            observed_columns
        ):
            raise ValueError(
                f"{read_gene_path} must contain columns "
                "named read_id and Gene."
            )

        for row in reader:
            read_id = normalize_read_id(
                row["read_id"]
            )
            gene = row["Gene"].strip()

            if gene not in CLB_GENES:
                continue

            species_taxids = species_by_read.get(
                read_id,
                set(),
            )

            if len(species_taxids) == 1:
                species_taxid = next(
                    iter(species_taxids)
                )

                species_name = scientific_names.get(
                    species_taxid,
                    f"taxid_{species_taxid}",
                )

                row_key = (
                    species_name,
                    species_taxid,
                )

            elif len(species_taxids) > 1:
                row_key = (
                    "Conflicting_species",
                    "-1",
                )

            elif read_id in classified_reads:
                row_key = (
                    "Unresolved_at_species",
                    "-1",
                )

            else:
                row_key = (
                    "Unclassified",
                    "0",
                )

            matrix[row_key][gene] += 1

    return matrix


def write_matrix(matrix, output_path):
    """Write the species-by-clb-gene matrix."""
    with output_path.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "Species",
                "TaxID",
                *CLB_GENES,
                "Total",
            ]
        )

        sorted_rows = sorted(
            matrix,
            key=lambda value: (
                value[0].lower(),
                value[1],
            ),
        )

        for species_name, taxid in sorted_rows:
            counts = matrix[
                (species_name, taxid)
            ]

            gene_counts = [
                counts.get(gene, 0)
                for gene in CLB_GENES
            ]

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
            "Join read-to-clb-gene assignments with direct "
            "KrakenUniq classifications and produce a "
            "species-by-clb-gene count matrix."
        )
    )

    parser.add_argument(
        "--read-gene",
        required=True,
        type=Path,
        help=(
            "TSV containing read_id and Gene columns."
        ),
    )

    parser.add_argument(
        "--kraken-output",
        required=True,
        type=Path,
        help="Per-read KrakenUniq output file.",
    )

    parser.add_argument(
        "--kraken-db",
        required=True,
        type=Path,
        help=(
            "KrakenUniq database containing either taxDB "
            "or taxonomy/nodes.dmp and taxonomy/names.dmp."
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help=(
            "Output species-by-clb-gene TSV file."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()

    if not args.read_gene.is_file():
        raise FileNotFoundError(
            "Read-to-gene table not found: "
            f"{args.read_gene}"
        )

    if not args.kraken_output.is_file():
        raise FileNotFoundError(
            "KrakenUniq output not found: "
            f"{args.kraken_output}"
        )

    if not args.kraken_db.is_dir():
        raise NotADirectoryError(
            "KrakenUniq database directory not found: "
            f"{args.kraken_db}"
        )

    parents, ranks, scientific_names = load_taxonomy(
        args.kraken_db
    )

    species_by_read, classified_reads = (
        read_kraken_classifications(
            args.kraken_output,
            parents,
            ranks,
        )
    )

    matrix = build_matrix(
        args.read_gene,
        species_by_read,
        classified_reads,
        scientific_names,
    )

    write_matrix(
        matrix,
        args.output,
    )


if __name__ == "__main__":
    main()