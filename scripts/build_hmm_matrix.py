#!/usr/bin/env python3
"""Build a Gene x Sample count matrix from per-sample HMM count TSVs.

Preferred form, and the one the pipeline uses:

    build_hmm_matrix.py --manifest hmm.manifest.tsv --out matrix.tsv

The manifest carries the sample identifier from the sample sheet next to the file it
produced, so no identifier is reconstructed from a filename (F03). The legacy form

    build_hmm_matrix.py --inputs a.hmm_counts.tsv ... --out matrix.tsv --strip-suffix .hmm_counts.tsv

is kept for running the script by hand; it strips one anchored terminal suffix and
refuses duplicates in the same way.
"""
import argparse
import sys
import os
from pathlib import Path
from typing import Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from clb_counts_schema import CLB_GENES, CountSchemaError, validate_clb_counts


def fail(message):
    raise SystemExit(f"[ERROR] {message}")


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
    for label, index in (("sample identifier", 0), ("input file", 1)):
        seen = set()
        duplicated = set()
        for pair in pairs:
            value = pair[index]
            if value in seen:
                duplicated.add(value)
            seen.add(value)
        if duplicated:
            fail(f"Duplicate {label}(s): {', '.join(sorted(duplicated))}. "
                 "Every column of the matrix must belong to exactly one sample.")


def pairs_from_filenames(files, suffix_to_strip: Optional[str] = None):
    pairs = []
    for f in files:
        name = Path(f).name
        if suffix_to_strip and name.endswith(suffix_to_strip):
            sample = name[: -len(suffix_to_strip)]
        else:
            sample = Path(f).stem
        pairs.append((sample, str(f)))
    return pairs


def build_matrix(pairs) -> pd.DataFrame:
    if not pairs:
        fail("No input files provided.")
    reject_duplicates(pairs)

    # Deterministic column order, independent of task completion order.
    pairs = sorted(pairs, key=lambda pair: pair[0])

    columns = []
    for sample, f in pairs:
        # F06: the same contract the alignment merger holds its inputs to -- exactly
        # clbA-clbS, one row each, finite non-negative integers. This used to read the
        # counts as floats, group duplicate genes and sum them, and fill anything
        # missing with zero, so an unknown gene, a negative, a fraction or a
        # header-only sample became a plausible column.
        df = pd.read_csv(f, sep="\t", usecols=["Gene", "Count"], dtype=str)
        counts = validate_clb_counts(zip(df["Gene"], df["Count"]), f)
        s = pd.Series(counts, name=sample, dtype="int64")
        columns.append(s)

    matrix = pd.concat(columns, axis=1)

    # pd.concat is happy to produce two columns of the same name; the matrix is not.
    if matrix.columns.duplicated().any():
        repeated = sorted(set(matrix.columns[matrix.columns.duplicated()]))
        fail(f"Refusing to write a matrix with duplicate sample column(s): {', '.join(repeated)}")

    # Rows in the canonical gene order, columns sorted by sample.
    return matrix.reindex(CLB_GENES).reindex(sorted(matrix.columns), axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=None,
                    help="Two-column TSV: sample_id, file. Preferred; identifiers come from the sample sheet.")
    ap.add_argument("--inputs", nargs="+", default=None, help="Input TSV files with columns: Gene, Count")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--strip-suffix", default=None,
                    help="With --inputs: one anchored terminal suffix to strip to form the sample name")
    args = ap.parse_args()

    if bool(args.manifest) == bool(args.inputs):
        ap.error("give exactly one of --manifest or --inputs")

    pairs = (read_manifest(args.manifest) if args.manifest
             else pairs_from_filenames(args.inputs, suffix_to_strip=args.strip_suffix))

    try:
        matrix = build_matrix(pairs)
    except CountSchemaError as problem:
        raise SystemExit(f"[ERROR] {problem}")
    matrix.to_csv(args.out, sep="\t", index=True)


if __name__ == "__main__":
    main()
