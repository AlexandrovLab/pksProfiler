#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


CLB_GENES = tuple(
    f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"
)

SPECIES_COLUMNS = ("Species", "TaxID", *CLB_GENES, "Total")
OUTPUT_COLUMNS = ("Sample", *SPECIES_COLUMNS)


def read_species_support(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(reader.fieldnames or ())
        rows = list(reader)

    if fieldnames != SPECIES_COLUMNS:
        raise ValueError(
            f"{path}: expected columns {SPECIES_COLUMNS}; "
            f"observed {fieldnames}"
        )

    return rows


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


def sample_from_path(path):
    suffix = ".clb_species_support.tsv"

    if not path.name.endswith(suffix):
        raise ValueError(
            f"Cannot determine sample from species-support file: {path}"
        )

    sample = path.name[: -len(suffix)]

    if not sample:
        raise ValueError(f"Empty sample identifier in filename: {path}")

    return sample


def load_rows(paths):
    combined_rows = []
    observed_samples = set()

    for path in sorted(paths):
        sample = sample_from_path(path)

        if sample in observed_samples:
            raise ValueError(
                f"Duplicate species-support file for sample: {sample}"
            )

        observed_samples.add(sample)

        for row in read_species_support(path):
            gene_total = sum(
                parse_nonnegative_integer(row[gene], path, gene)
                for gene in CLB_GENES
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

            combined_rows.append(
                {
                    "Sample": sample,
                    **row,
                }
            )

    return combined_rows


def write_output(rows, output_path):
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_list_file(path):
    """One path per line. Blank lines ignored; every named file must be here."""
    paths = []
    for number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        name = line.strip()
        if not name:
            continue
        candidate = Path(name)
        if not candidate.exists():
            raise SystemExit(f"[ERROR] {path} line {number} names a file that is not here: {name}")
        paths.append(candidate)
    return paths


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Combine per-sample species-by-clb direct-support tables."
        )
    )
    # F02: one filename per input overflows the OS argument limit on a large cohort.
    # --species-files-from passes a single file listing them instead.
    parser.add_argument(
        "--species-files",
        nargs="+",
        default=[],
        type=Path,
    )
    parser.add_argument(
        "--species-files-from",
        default=None,
        type=Path,
        help="file of per-sample table paths, one per line (preferred)",
    )
    parser.add_argument(
        "--species-output",
        required=True,
        type=Path,
    )
    args = parser.parse_args()
    if bool(args.species_files) == bool(args.species_files_from):
        parser.error("give exactly one of --species-files or --species-files-from")
    if args.species_files_from:
        args.species_files = read_list_file(args.species_files_from)
    return args


def main():
    args = parse_args()
    rows = load_rows(args.species_files)
    write_output(rows, args.species_output)


if __name__ == "__main__":
    main()