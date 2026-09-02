#!/usr/bin/env python3
import math
import os
import sys

import pandas as pd
from pandas.errors import EmptyDataError

counts_files = sys.argv[1:-1]
output_file  = sys.argv[-1]

EXPECTED_GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]


def fail(message):
    raise ValueError(message)


def read_counts(path):
    """
    Return (counts_series, sample_name), failing on incomplete or malformed input.
    counts_series is a Series indexed by Geneid with the sample's counts.
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

    # Clean sample name
    base = os.path.basename(path)
    sample_name = (base
                   .replace(".counts.txt.gz","")
                   .replace(".counts.txt","")
                   .replace(".txt.gz","")
                   .replace(".txt",""))

    if df["Geneid"].duplicated().any():
        duplicates = sorted(df.loc[df["Geneid"].duplicated(), "Geneid"].unique())
        fail(f"Duplicate Geneid values in {path}: {', '.join(duplicates)}")

    observed = set(df["Geneid"])
    expected = set(EXPECTED_GENES)
    if observed != expected or len(df) != len(EXPECTED_GENES):
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        fail(f"Expected exactly clbA-clbS in {path}; missing={missing}, extra={extra}")

    numeric_counts = pd.to_numeric(df[counts_col], errors="raise")
    if numeric_counts.isna().any() or not numeric_counts.map(math.isfinite).all() or (numeric_counts < 0).any():
        fail(f"Counts must be finite non-negative numbers in {path}")

    counts = pd.Series(numeric_counts.values, index=df["Geneid"], name=sample_name)
    counts = counts.reindex(EXPECTED_GENES)
    counts.name = sample_name

    return counts, sample_name

# Read all usable files
all_counts = []
used = []

for p in counts_files:
    counts, name = read_counts(p)
    all_counts.append(counts)
    used.append(name)

if not all_counts:
    sys.stderr.write("[ERROR] No valid counts files found. Aborting.\n")
    sys.exit(1)

# Every input has already been proven to contain the same exact gene set.
merged_counts = pd.concat(all_counts, axis=1, join="inner").reindex(EXPECTED_GENES)

# Attach annotation from the first valid file
merged_counts.index.name = "Gene"
merged = merged_counts.reset_index()
merged.to_csv(output_file, sep="\t", index=False)


sys.stderr.write(f"[INFO] Merged {len(used)} samples -> {output_file}\n")
