#!/usr/bin/env python3
"""Merge per-sample clb count files into one Gene x Sample matrix.

Preferred form, and the one the pipeline uses:

    mergeGeneCounts.py --manifest counts.manifest.tsv --output merged.tsv

The manifest carries the sample identifier from the sample sheet next to the file it
produced, so no identifier is ever reconstructed from a filename. F03: names were
rebuilt with chained unanchored str.replace(), which collapsed the distinct sheet
identifiers `case` and `case.txt` onto one column, and the merge then exited 0 with a
duplicated column. Identifiers now arrive with the data and duplicates are refused.

The legacy positional form is kept for running the script by hand:

    mergeGeneCounts.py a.counts.txt b.counts.txt merged.tsv

It strips exactly one recognised terminal suffix -- not every occurrence anywhere in
the name -- and refuses duplicates in the same way.
"""
import math
import os
import sys

import pandas as pd

sys.path.insert(0, str(os.path.dirname(os.path.abspath(__file__))))
from clb_counts_schema import CountSchemaError, validate_clb_counts
from pandas.errors import EmptyDataError

EXPECTED_GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]

# Longest first: `.counts.txt` must win over `.txt`, or `case.counts.txt` loses only
# `.txt` and `case` and `case.txt` collide again by a different route.
KNOWN_SUFFIXES = (".counts.txt.gz", ".counts.tsv.gz", ".counts.txt", ".counts.tsv",
                  ".txt.gz", ".tsv.gz", ".txt", ".tsv")


def fail(message):
    raise ValueError(message)


def sample_from_filename(path):
    """Strip exactly one recognised terminal suffix. Anchored, and applied once."""
    base = os.path.basename(path)
    for suffix in KNOWN_SUFFIXES:
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def read_manifest(path):
    """Return [(sample_id, file)] from a two-column TSV, header optional."""
    try:
        lines = [line.rstrip("\n") for line in open(path) if line.strip()]
    except OSError as exc:
        fail(f"Cannot read manifest {path}: {exc}")
    if not lines:
        fail(f"Manifest is empty: {path}")

    if lines[0].split("\t")[:2] == ["sample_id", "file"]:
        lines = lines[1:]
    if not lines:
        fail(f"Manifest has a header and no rows: {path}")

    pairs = []
    for number, line in enumerate(lines, start=1):
        fields = line.split("\t")
        if len(fields) != 2:
            fail(f"Manifest line {number} is not two tab-separated fields: {line!r}")
        sample_id, counts_file = (field.strip() for field in fields)
        if not sample_id or not counts_file:
            fail(f"Manifest line {number} has an empty field: {line!r}")
        if not os.path.exists(counts_file):
            fail(f"Manifest line {number} names a file that is not here: {counts_file}")
        pairs.append((sample_id, counts_file))
    return pairs


def reject_duplicates(pairs):
    for label, index in (("sample identifier", 0), ("counts file", 1)):
        seen = set()
        duplicated = set()
        for pair in pairs:
            value = pair[index]
            if value in seen:
                duplicated.add(value)
            seen.add(value)
        if duplicated:
            fail(f"Duplicate {label}(s) in the merge input: {', '.join(sorted(duplicated))}. "
                 "Every column of the master table must belong to exactly one sample.")


def parse_args(argv):
    """Return ([(sample_id, file)], output_path)."""
    if "--manifest" in argv:
        import argparse
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--manifest", required=True)
        parser.add_argument("--output", required=True)
        args = parser.parse_args(argv)
        return read_manifest(args.manifest), args.output

    if len(argv) < 2:
        fail("usage: mergeGeneCounts.py --manifest M --output OUT"
             "   |   mergeGeneCounts.py COUNTS... OUT")
    return [(sample_from_filename(path), path) for path in argv[:-1]], argv[-1]


def read_counts(path, sample_name):
    """
    Return (counts_series, sample_name), failing on incomplete or malformed input.
    counts_series is a Series indexed by Geneid with the sample's counts. The name is
    given, never inferred from the file.
    """
    # Skip zero-byte
    try:
        if os.path.getsize(path) == 0:
            fail(f"Empty counts file: {path}")
    except OSError:
        fail(f"Cannot stat counts file: {path}")

    try:
        df = pd.read_csv(
            path, sep="\t", comment="#", low_memory=False,
            header=0, compression="infer"
        )
    except EmptyDataError:
        fail(f"Counts file contains no data: {path}")

    if df.empty:
        fail(f"Counts file contains no rows: {path}")

    # Expect at least the 6 annotation cols + 1 counts col
    if "Geneid" not in df.columns or df.shape[1] < 7:
        fail(f"Unexpected columns in {path}: {list(df.columns)}")

    # Determine counts column: last non-annotation column
    anno_cols = ["Geneid","Chr","Start","End","Strand","Length"]
    non_anno = [c for c in df.columns if c not in anno_cols]
    if not non_anno:
        fail(f"No counts column found in {path}")
    counts_col = non_anno[-1]

    # F06: one contract for both methods. This checked the gene set, duplicates and
    # non-negative finite values already, but accepted fractional counts -- and the
    # HMM merger beside it checked none of them.
    try:
        validated = validate_clb_counts(zip(df["Geneid"], df[counts_col]), path)
    except CountSchemaError as problem:
        fail(str(problem))

    counts = pd.Series(validated, name=sample_name, dtype="int64")
    counts = counts.reindex(EXPECTED_GENES)
    counts.name = sample_name

    return counts, sample_name

pairs, output_file = parse_args(sys.argv[1:])
reject_duplicates(pairs)

# Deterministic column order, independent of the order tasks happened to finish in.
pairs.sort(key=lambda pair: pair[0])

all_counts = []
used = []

for sample_name, path in pairs:
    counts, name = read_counts(path, sample_name)
    all_counts.append(counts)
    used.append(name)

if not all_counts:
    sys.stderr.write("[ERROR] No valid counts files found. Aborting.\n")
    sys.exit(1)

# Every input has already been proven to contain the same exact gene set.
merged_counts = pd.concat(all_counts, axis=1, join="inner").reindex(EXPECTED_GENES)

# Last line of defence: whatever route the names arrived by, the table may not ship
# with two columns of the same name.
if merged_counts.columns.duplicated().any():
    repeated = sorted(set(merged_counts.columns[merged_counts.columns.duplicated()]))
    fail(f"Refusing to write a table with duplicate sample column(s): {', '.join(repeated)}")

# Attach annotation from the first valid file
merged_counts.index.name = "Gene"
merged = merged_counts.reset_index()
merged.to_csv(output_file, sep="\t", index=False)


sys.stderr.write(f"[INFO] Merged {len(used)} samples -> {output_file}\n")
