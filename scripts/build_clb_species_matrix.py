#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


CLB_GENES = tuple(
    f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"
)


def normalize_read_id(read_id):
    """Return the exact sequence identifier used for cross-file joins."""
    read_id = read_id.strip().split()[0]
    return read_id.lstrip("@>")


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


# Ranks seen in an NCBI-derived taxonomy. Only used to tell a rank column from a name
# column -- an unknown rank is fine, a column full of organism names is not.
KNOWN_RANKS = {
    "no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum",
    "subphylum", "superclass", "class", "subclass", "infraclass", "cohort",
    "superorder", "order", "suborder", "infraorder", "parvorder", "superfamily",
    "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "species group",
    "species subgroup", "species", "subspecies", "varietas", "forma", "strain",
    "clade", "serotype", "serogroup", "biotype", "genotype", "isolate", "morph",
    "section", "subsection", "series", "forma specialis", "pathogroup",
}


def assert_taxdb_column_order(rows, taxdb_path, sample_size=2000):
    """
    Fail if the rank is not where we read it from.

    F04 reported that this parser swapped rank and scientific name. It does not, for
    the KrakenUniq database in use: across its 2.5 million rows, field 4 holds a rank
    and field 3 never does. But the parser only ever *assumed* that, and a build that
    differed would have been mis-parsed in silence -- every organism named "species"
    or "no rank" -- which is the failure the finding describes even if the diagnosis
    was inverted. So the order is checked rather than assumed, on both candidates, and
    a file that matches neither is an error rather than a guess.
    """
    third = fourth = 0
    for fields in rows[:sample_size]:
        if len(fields) >= 4:
            third += fields[2].strip().lower() in KNOWN_RANKS
            fourth += fields[3].strip().lower() in KNOWN_RANKS
    if fourth >= max(third, 1):
        return
    if third > fourth:
        raise ValueError(
            f"{taxdb_path}: the rank looks like field 3, not field 4 "
            f"({third} of the first {min(len(rows), sample_size)} rows against "
            f"{fourth}). This parser reads field 3 as the scientific name and field 4 "
            "as the rank; that would mis-parse every organism here. Check the "
            "KrakenUniq build that produced this taxDB."
        )
    raise ValueError(
        f"{taxdb_path}: neither field 3 nor field 4 holds a recognisable taxonomic "
        f"rank in the first {min(len(rows), sample_size)} rows. This does not look "
        "like a KrakenUniq taxDB."
    )


def parse_krakenuniq_taxdb(taxdb_path):
    """
    Parse KrakenUniq taxDB.

    Columns, verified against krakenUniq_8_8_2023 and asserted at load:
      1. taxonomy ID
      2. parent taxonomy ID
      3. scientific name
      4. taxonomy rank
    """
    parents = {}
    ranks = {}
    scientific_names = {}

    with taxdb_path.open() as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle]

    assert_taxdb_column_order(rows, taxdb_path)

    for line_number, fields in enumerate(rows, start=1):
        if len(fields) < 4:
            continue

        taxid = fields[0].strip()
        parent_taxid = fields[1].strip()
        scientific_name = fields[2].strip()
        rank = fields[3].strip()

        if not taxid.isdigit():
            raise ValueError(
                f"Invalid taxid at {taxdb_path}:"
                f"{line_number}: {taxid!r}"
            )

        if not parent_taxid.isdigit():
            raise ValueError(
                f"Invalid parent taxid at {taxdb_path}:"
                f"{line_number}: {parent_taxid!r}"
            )

        if not scientific_name:
            raise ValueError(
                f"Missing scientific name at {taxdb_path}:"
                f"{line_number}"
            )

        if not rank:
            rank = "no rank"

        parents[taxid] = parent_taxid
        ranks[taxid] = rank
        scientific_names[taxid] = scientific_name

    return parents, ranks, scientific_names


def load_taxonomy(kraken_db):
    """
    Load the taxonomy that came with this database, or an NCBI dump if there is none.

    F05. The order used to be the other way round: an NCBI nodes.dmp/names.dmp pair
    beside the database won over its own taxDB. Kraken assigned these taxids by
    walking taxDB, so taxDB is the tree that explains them -- and it carries
    assembly- and sequence-level pseudo-taxids that the database build invents and
    NCBI has never heard of. Resolve one of those against an NCBI dump and there is
    no node, so there is no path up to a species: the read leaves no species row and
    nothing reports that it happened. Ludmil reproduced exactly that with a fixture.

    taxDB is also the snapshot the classification was made against. NCBI's taxonomy
    moves -- taxa merge, get renamed, get reclassified -- so a later dump disagrees
    with an earlier classification in places it will not point out.
    """
    taxonomy_dir = kraken_db / "taxonomy"
    nodes_path = taxonomy_dir / "nodes.dmp"
    names_path = taxonomy_dir / "names.dmp"
    taxdb_path = kraken_db / "taxDB"

    if taxdb_path.is_file():
        if nodes_path.is_file() and names_path.is_file():
            sys.stderr.write(
                f"[INFO] taxonomy: using {taxdb_path}, which came with this database. "
                f"An NCBI dump at {taxonomy_dir} is present and ignored -- it does not "
                "contain the database's own assembly and sequence taxids.\n")
        return parse_krakenuniq_taxdb(taxdb_path)

    if nodes_path.is_file() and names_path.is_file():
        sys.stderr.write(
            f"[WARN] taxonomy: no taxDB in {kraken_db}; falling back to the NCBI dump "
            f"at {taxonomy_dir}. Reads assigned to database-specific taxids will not "
            "resolve to a species.\n")
        return parse_ncbi_taxonomy(
            nodes_path,
            names_path,
        )

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
    Build a species-by-clb-gene direct-support matrix.

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
    """Write the species-by-clb-gene direct-support matrix."""
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
            "species-by-clb-gene direct-support matrix."
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
            "Output species-by-clb-gene support TSV file."
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
