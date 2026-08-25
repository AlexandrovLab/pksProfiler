#!/usr/bin/env python3
import pandas as pd
import sys, os
from pandas.errors import EmptyDataError

counts_files = sys.argv[1:-1]
output_file  = sys.argv[-1]

def read_counts(path):
    """
    Return (annot_df, counts_series, sample_name) or None if unusable.
    annot_df has ['Geneid','Chr','Start','End','Strand','Length'].
    counts_series is a Series indexed by Geneid with the sample's counts.
    """
    # Skip zero-byte
    try:
        if os.path.getsize(path) == 0:
            sys.stderr.write(f"[WARN] Empty file (0 bytes): {path}\n")
            return None
    except OSError:
        sys.stderr.write(f"[WARN] Cannot stat file, skipping: {path}\n")
        return None

    try:
        df = pd.read_csv(
            path, sep="\t", comment="#", low_memory=False,
            header=0, compression="infer"
        )
    except EmptyDataError:
        sys.stderr.write(f"[WARN] No data (only comments?): {path}\n")
        return None

    if df.empty:
        sys.stderr.write(f"[WARN] Empty dataframe after parsing: {path}\n")
        return None

    # Expect at least the 6 annotation cols + 1 counts col
    if "Geneid" not in df.columns or df.shape[1] < 7:
        # Some runs may lack headers; try force-naming if exactly 7+ columns
        if "Geneid" not in df.columns and df.shape[1] >= 7:
            df.columns = ["Geneid","Chr","Start","End","Strand","Length"] + \
                         [f"col{i}" for i in range(7, df.shape[1]+1)]
        else:
            sys.stderr.write(f"[WARN] Unexpected columns in {path}: {list(df.columns)}\n")
            return None

    # Determine counts column: last non-annotation column
    anno_cols = ["Geneid","Chr","Start","End","Strand","Length"]
    non_anno = [c for c in df.columns if c not in anno_cols]
    if not non_anno:
        sys.stderr.write(f"[WARN] No counts column found in {path}\n")
        return None
    counts_col = non_anno[-1]

    # Clean sample name
    base = os.path.basename(path)
    sample_name = (base
                   .replace(".counts.txt.gz","")
                   .replace(".counts.txt","")
                   .replace(".txt.gz","")
                   .replace(".txt",""))

    # Prepare outputs
    annot = df[anno_cols].copy()
    # if duplicates in Geneid exist, group/sum (featureCounts should not, but just in case)
    df = df[[ "Geneid", counts_col ]].copy()
    df = df.groupby("Geneid", as_index=False)[counts_col].sum()
    counts = df.set_index("Geneid")[counts_col]
    counts.name = sample_name

    return annot, counts, sample_name

# Read all usable files
ann_ref = None
all_counts = []
used = []

for p in counts_files:
    res = read_counts(p)
    if res is None:
        continue
    annot, counts, name = res
    if ann_ref is None:
        ann_ref = annot.drop_duplicates(subset=["Geneid"]).set_index("Geneid")
    all_counts.append(counts)
    used.append(name)

if not all_counts:
    sys.stderr.write("[ERROR] No valid counts files found. Aborting.\n")
    sys.exit(1)

# Outer-merge all counts on Geneid
merged_counts = pd.concat(all_counts, axis=1, join="outer")

# Attach annotation from the first valid file
merged_counts.index.name = "Gene"
merged = merged_counts.reset_index()
merged.to_csv(output_file, sep="\t", index=False)


sys.stderr.write(f"[INFO] Merged {len(used)} samples -> {output_file}\n")

