#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path


CLB_GENES = tuple(
    f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"
)

GENE_COLUMNS = ("Sample", *CLB_GENES, "Total")
SPECIES_COLUMNS = ("Species", "TaxID", *CLB_GENES, "Total")
COMBINED_SPECIES_COLUMNS = (
    "Sample",
    *SPECIES_COLUMNS,
)


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(reader.fieldnames or ())
        rows = list(reader)

    return fieldnames, rows


def parse_nonnegative_integer(value, path, column):
    try:
        number = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"{path}: {column} must be an integer; observed {value!r}"
        ) from error

    if number < 0:
        raise ValueError(
            f"{path}: {column} must be nonnegative; observed {number}"
        )

    return number


def load_gene_support(paths):
    rows_by_sample = {}

    for path in sorted(paths):
        fieldnames, rows = read_tsv(path)

        if fieldnames != GENE_COLUMNS:
            raise ValueError(
                f"{path}: expected columns {GENE_COLUMNS}; "
                f"observed {fieldnames}"
            )

        if len(rows) != 1:
            raise ValueError(
                f"{path}: expected exactly one data row; "
                f"observed {len(rows)}"
            )

        row = rows[0]
        sample = row["Sample"].strip()

        if not sample:
            raise ValueError(f"{path}: Sample must not be empty")

        if sample in rows_by_sample:
            raise ValueError(f"Duplicate sample in gene support: {sample}")

        gene_total = 0

        for gene in CLB_GENES:
            gene_total += parse_nonnegative_integer(
                row[gene],
                path,
                gene,
            )

        reported_total = parse_nonnegative_integer(
            row["Total"],
            path,
            "Total",
        )

        if gene_total != reported_total:
            raise ValueError(
                f"{path}: gene sum {gene_total} does not equal "
                f"Total {reported_total}"
            )

        rows_by_sample[sample] = row

    return rows_by_sample


def sample_from_species_path(path):
    suffix = ".clb_species_support.tsv"

    if not path.name.endswith(suffix):
        raise ValueError(
            f"Cannot determine sample from species-support file: {path}"
        )

    return path.name[: -len(suffix)]


def load_species_support(paths):
    rows_by_sample = defaultdict(list)

    for path in sorted(paths):
        fieldnames, rows = read_tsv(path)

        if fieldnames != SPECIES_COLUMNS:
            raise ValueError(
                f"{path}: expected columns {SPECIES_COLUMNS}; "
                f"observed {fieldnames}"
            )

        sample = sample_from_species_path(path)

        if sample in rows_by_sample:
            raise ValueError(
                f"Duplicate species-support file for sample: {sample}"
            )

        # Ensure an empty, header-only matrix still registers the sample.
        rows_by_sample[sample] = []

        for row in rows:
            gene_total = 0

            for gene in CLB_GENES:
                gene_total += parse_nonnegative_integer(
                    row[gene],
                    path,
                    gene,
                )

            reported_total = parse_nonnegative_integer(
                row["Total"],
                path,
                "Total",
            )

            if gene_total != reported_total:
                raise ValueError(
                    f"{path}: row gene sum {gene_total} does not equal "
                    f"Total {reported_total}"
                )

            rows_by_sample[sample].append(row)

    return dict(rows_by_sample)


def validate_accounting(gene_rows, species_rows):
    gene_samples = set(gene_rows)
    species_samples = set(species_rows)

    if gene_samples != species_samples:
        raise ValueError(
            "Sample mismatch between overall and species support: "
            f"overall_only={sorted(gene_samples - species_samples)}, "
            f"species_only={sorted(species_samples - gene_samples)}"
        )

    for sample in sorted(gene_rows):
        overall = gene_rows[sample]
        rows = species_rows[sample]

        for gene in CLB_GENES:
            species_total = sum(int(row[gene]) for row in rows)
            overall_total = int(overall[gene])

            if species_total != overall_total:
                raise ValueError(
                    f"{sample}: {gene} species-support sum "
                    f"{species_total} does not equal overall support "
                    f"{overall_total}"
                )


def write_gene_output(rows_by_sample, output_path):
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=GENE_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()

        for sample in sorted(rows_by_sample):
            writer.writerow(rows_by_sample[sample])


def write_species_output(rows_by_sample, output_path):
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=COMBINED_SPECIES_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()

        for sample in sorted(rows_by_sample):
            for row in rows_by_sample[sample]:
                writer.writerow(
                    {
                        "Sample": sample,
                        **row,
                    }
                )


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Combine per-sample clb support tables and verify that "
            "species/support categories reconcile with overall counts."
        )
    )
    parser.add_argument(
        "--gene-files",
        nargs="+",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--species-files",
        nargs="+",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--gene-output",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--species-output",
        required=True,
        type=Path,
    )
    return parser.parse_args()


def main():
    args = parse_args()
    gene_rows = load_gene_support(args.gene_files)
    species_rows = load_species_support(args.species_files)
    validate_accounting(gene_rows, species_rows)
    write_gene_output(gene_rows, args.gene_output)
    write_species_output(species_rows, args.species_output)


if __name__ == "__main__":
    main()