#!/usr/bin/env python3
import argparse
from pathlib import Path
import re

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# Bracken writes these. Nothing else is a Bracken report, and guessing which other
# column might be the abundance is how a taxonomy id gets plotted as a read count.
REQUIRED_COLUMNS = ("name", "new_est_reads")


def read_bracken_report(path: Path) -> pd.DataFrame:
    """
    Read a Bracken report, or refuse it.

    F11. This used to fall back to "a column that looks numeric" when new_est_reads
    was absent, and to the first column when name was. A file carrying only name and
    taxonomy_id was therefore accepted, and plotted Escherichia coli at 562 reads --
    562 being its taxonomy id. It also caught every parse error and returned an empty
    frame, which draws as "No taxa (empty or 0 reads)": a malformed file and a sample
    with nothing in it produced the same picture.

    An absent or empty file still reads as no taxa, because the workflow writes
    zero-byte reports deliberately for samples with no hits. A file with content that
    is not a Bracken report is an error.
    """
    if (not path.exists()) or path.stat().st_size == 0:
        return pd.DataFrame(columns=["taxon", "reads"])

    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception as problem:
        raise ValueError(f"{path}: not readable as a Bracken report: {problem}") from problem

    if df.empty:
        return pd.DataFrame(columns=["taxon", "reads"])

    columns = {str(c).lower(): c for c in df.columns}
    missing = [name for name in REQUIRED_COLUMNS if name not in columns]
    if missing:
        raise ValueError(
            f"{path}: not a Bracken report -- missing {', '.join(missing)}. "
            f"Columns present: {', '.join(str(c) for c in df.columns)}. "
            "Refusing to guess which column holds the abundance."
        )

    reads = pd.to_numeric(df[columns["new_est_reads"]], errors="coerce")
    if reads.notna().sum() == 0:
        raise ValueError(
            f"{path}: new_est_reads holds no numbers, so this is not an abundance "
            "column."
        )

    out = pd.DataFrame({"taxon": df[columns["name"]].astype(str), "reads": reads})
    out = out.dropna(subset=["taxon", "reads"])
    return out[out["reads"] > 0]


def barplot(ax, df: pd.DataFrame, title: str, top_n: int):
    ax.set_title(title)
    if df.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No taxa (empty or 0 reads)", ha="center", va="center")
        return

    df = df.sort_values("reads", ascending=False)
    if top_n > 0 and len(df) > top_n:
        df = df.iloc[:top_n].copy()

    ax.bar(df["taxon"], df["reads"])
    ax.set_ylabel("Estimated PKS-island reads")
    ax.set_xlabel("Taxon")
    ax.tick_params(axis="x", labelrotation=60)
    ax.margins(x=0.01)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, help="Sample ID for titles")
    ap.add_argument("--genus", required=True, help="Bracken genus report (*.bracken.G.report.txt)")
    ap.add_argument("--species", required=True, help="Bracken species report (*.bracken.S.report.txt)")
    ap.add_argument("--out", required=True, help="Output PDF")
    ap.add_argument("--top", type=int, default=20, help="Top N taxa to plot (0 = all)")
    args = ap.parse_args()

    sample = args.sample
    genus_path = Path(args.genus)
    species_path = Path(args.species)

    gdf = read_bracken_report(genus_path)
    sdf = read_bracken_report(species_path)

    out_pdf = Path(args.out)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(out_pdf) as pdf:
        # Genus page
        fig, ax = plt.subplots(figsize=(12, 6))
        barplot(ax, gdf, f"{sample} — Bracken Genus (top {args.top if args.top>0 else 'all'})", args.top)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Species page
        fig, ax = plt.subplots(figsize=(12, 6))
        barplot(ax, sdf, f"{sample} — Bracken Species (top {args.top if args.top>0 else 'all'})", args.top)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


if __name__ == "__main__":
    main()
