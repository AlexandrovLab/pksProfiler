#!/usr/bin/env python3
"""The pks island, the contigs that reconstruct it, and the read depth beneath them.

One figure, three stacked panels on a shared axis:

  island    the 19 clb genes as strand-aware arrows, filled by enzyme class
  contigs   one row per source, each aligned contig an arrow carrying its length
            and percent identity; the row label gives contig count and the share
            of the island that source recovers
  depth     read depth across the same window, with the median marked

The two sample types differ only in what the contig rows are:

    wgs          metaSPAdes and MEGAHIT, the two targeted assemblers. Agreement
                 between them is the evidence.
    metagenome   every recovered bin that carries island sequence, then the
                 unbinned pool. Island sequence inside a bin can be attributed to
                 a candidate organism; island sequence in the unbinned pool is
                 evidence of presence and says nothing about who carries it, so
                 the unbinned row is drawn in grey rather than a bin colour and
                 always sits last.

Drawing conventions follow how gene-cluster figures are normally published
(clinker, genoPlotR, EasyFig) and live in scripts/pks_figure_style.py, shared with
anything else that draws this island so the two cannot drift apart.

Alignments are filtered as the pipeline filters them -- at least 90% identity and
MAPQ at least 20 -- so the coverage printed here is what the pipeline reports
rather than a raw PAF span.

INPUTS are small summaries the pipeline already publishes: PAF alignments, a
bedGraph depth track. No contig or read sequence is opened.

    # targeted WGS
    python3 plot_island_reconstruction.py --mode wgs \
        --sample-dir results/by_sample/<sample> --gff ref/annotations/IHE3034.clbA-clbS.gff \
        --out <sample>.island_reconstruction.pdf

    # metagenome, every bin with island sequence plus the unbinned pool
    python3 plot_island_reconstruction.py --mode metagenome \
        --sample-dir results/by_sample/<sample> --gff ref/annotations/IHE3034.clbA-clbS.gff \
        --out <sample>.island_reconstruction.pdf
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import matplotlib.lines as mlines        # noqa: E402
import pks_figure_style as P             # noqa: E402

COVERAGE_COLOUR = "#51606e"
UNBINNED_COLOUR = "#9aa5b1"
# Bin colours, reused in order. The first is the same blue the WGS metaSPAdes row
# uses, so a reader moving between the two figures sees one palette.
BIN_COLOURS = ["#2f6f9f", "#1f9d8f", "#8c5aa8", "#d08c27",
               "#4a7f3e", "#7a6aa8", "#2f8f8a", "#a0705f"]

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.edgecolor": P.MUTED, "axes.labelcolor": P.INK, "text.color": P.INK,
    "xtick.color": P.MUTED, "ytick.color": P.MUTED,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.facecolor": "white", "savefig.bbox": "tight", "pdf.fonttype": 42,
})


# ---------------------------------------------------------------- inputs

def coverage(path, contig, island_start, island_end, pad):
    """Bin midpoints, depths, and the bp each bin covers inside the island.

    bedGraph bins are NOT uniform -- zero-depth stretches come out as a few wide
    intervals -- so breadth must be weighted by bp rather than counted as a
    fraction of bins, which badly overestimates it.
    """
    xs, ys, widths = [], [], []
    path = Path(path)
    if not path.exists():
        return xs, ys, widths
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] != contig:
                continue
            start, end = int(fields[1]) + 1, int(fields[2])
            if end < island_start - pad or start > island_end + pad:
                continue
            xs.append((start + end) / 2)
            ys.append(float(fields[3]))
            lo, hi = max(start, island_start), min(end, island_end)
            widths.append(hi - lo + 1 if hi >= lo else 0)
    return xs, ys, widths


def wgs_tracks(sample_dir):
    """metaSPAdes then MEGAHIT, whichever of their PAFs exist."""
    found = []
    for name, rel, colour in (
            ("metaspades", "contigs/metaspades/metaspades.vs_IHE3034.paf", "#2f6f9f"),
            ("megahit", "contigs/megahit/megahit.vs_IHE3034.paf", "#1f9d8f")):
        path = Path(sample_dir) / rel
        if path.exists():
            found.append((name, path, colour))
    return found


def metagenome_tracks(sample_dir):
    """Every bin carrying island sequence, then the unbinned pool last.

    Bins are ordered by how much of the island they recover, so the strongest
    candidate carrier is the first row a reader meets. The unbinned pool is always
    last and always grey: it is not a genome and must not read as one.
    """
    directory = Path(sample_dir) / "genomes" / "locus_alignment"
    if not directory.is_dir():
        return []
    bins, unbinned = [], []
    for path in sorted(directory.glob("*.vs_IHE3034.paf")):
        unit = path.name[: -len(".vs_IHE3034.paf")]
        (unbinned if re.search(r"unbinn?ed", unit, re.I) else bins).append((unit, path))
    return bins, unbinned


def explicit_tracks(specs, island_start, island_end):
    """Tracks from `--paf label=path`, in the order given.

    A Nextflow task stages its inputs flat, so a process cannot rely on the directory
    layout --sample-dir discovery walks; it passes what it has instead. The ordering
    rule is the same either way: anything whose label names the unbinned pool is drawn
    grey and pushed last, because it is not a genome and the figure must not let it
    read as one.
    """
    named, pooled = [], []
    for spec in specs:
        if "=" not in spec:
            sys.exit(f"--paf expects label=path, got {spec!r}")
        label, path = spec.split("=", 1)
        hits = P.paf_hits(path, island_start, island_end)
        if not hits:
            continue
        (pooled if re.search(r"unbinn?ed", label, re.I) else named).append((label, hits))
    tracks = [(label, hits, BIN_COLOURS[i % len(BIN_COLOURS)])
              for i, (label, hits) in enumerate(named)]
    tracks += [(label, hits, UNBINNED_COLOUR) for label, hits in pooled]
    return tracks


# ---------------------------------------------------------------- figure

def build(args, genes):
    sample_dir = Path(args.sample_dir)
    island_start, island_end = args.island_start, args.island_end
    length = island_end - island_start + 1
    pad = int(length * 0.035)
    xlim = (island_start - pad, island_end + pad)

    if args.paf:
        tracks = explicit_tracks(args.paf, island_start, island_end)
        if not tracks:
            sys.exit("no --paf alignment passes the identity and MAPQ filters")
    elif args.mode == "wgs":
        sources = wgs_tracks(sample_dir)
        if not sources:
            sys.exit(f"{sample_dir} has no targeted-assembly PAF under contigs/; was "
                     "this sample run with --tumor_targeted_assembly?")
        tracks = [(name, P.paf_hits(path, island_start, island_end), colour)
                  for name, path, colour in sources]
    else:
        result = metagenome_tracks(sample_dir)
        if not result:
            sys.exit(f"{sample_dir} has no genomes/locus_alignment/*.vs_IHE3034.paf; was "
                     "this run --enable_mags, and did it reach "
                     "alignMagBinToCanonicalReference?")
        bins, unbinned = result
        scored = []
        for unit, path in bins:
            hits = P.paf_hits(path, island_start, island_end)
            if hits:
                scored.append((unit, hits, P.union_bp(hits, island_start, island_end)))
        scored.sort(key=lambda t: -t[2])
        tracks = [(unit, hits, BIN_COLOURS[i % len(BIN_COLOURS)])
                  for i, (unit, hits, _bp) in enumerate(scored)]
        for unit, path in unbinned:
            hits = P.paf_hits(path, island_start, island_end)
            if hits:
                # metabat2 names the pool "bin.unbinned"; "bin.unbinned contigs" reads as
                # though it were a bin, which is the one thing this row must not suggest.
                tracks.append(("unbinned contigs", hits, UNBINNED_COLOUR))
        if not tracks:
            sys.exit("no bin or unbinned pool of this sample has an island alignment "
                     "passing the identity and MAPQ filters")

    depth_path = (Path(args.coverage) if Path(args.coverage).is_absolute()
                  else sample_dir / args.coverage)
    xs, ys, _widths = coverage(depth_path, args.contig,
                               island_start, island_end, pad)

    rows = sum(max(1, len(hits)) for _, hits, _ in tracks)
    figure = plt.figure(figsize=(7.4, 1.78 + 0.15 * rows))
    band = max(0.38, 0.19 * rows)
    # top=0.90 rather than 0.775: with the title and the statistics line gone there is
    # one short line of header to clear, not three.
    grid = figure.add_gridspec(3, 1, height_ratios=[0.26, band, 0.68], hspace=0.035,
                               left=0.185, right=0.735, top=0.90, bottom=0.165)

    def side(ax, text, y=0.5, x=-0.012):
        ax.text(x, y, text, transform=ax.transAxes, ha="right", va="center",
                fontsize=6.9, color=P.INK, fontweight="medium", linespacing=1.25)

    # --- island -------------------------------------------------------------
    ax_island = figure.add_subplot(grid[0])
    P.gene_track(ax_island, genes, y=0.0, h=0.74, label_min_bp=850, fontsize=5.8)
    ax_island.set_xlim(*xlim)
    ax_island.set_ylim(-0.42, 1.18)
    ax_island.axis("off")
    side(ax_island, f"pks island\n{len(genes)} clb genes", y=0.33)

    # --- contigs ------------------------------------------------------------
    ax_contigs = figure.add_subplot(grid[1], sharex=ax_island)
    ax_contigs.set_facecolor("#f4f2ef")
    for spine in ax_contigs.spines.values():
        spine.set_visible(False)
    y, ticks, labels = 0.0, [], []
    for name, hits, colour in tracks:
        ax_contigs.hlines(y, island_start, island_end, color="#dcd7d0", lw=1.1, zorder=1)
        for index, hit in enumerate(hits):
            yy = y - index * 0.44
            left = max(hit["ts"], xlim[0])
            right = min(hit["te"], xlim[1])
            P.arrow(ax_contigs, left, right, yy, 0.32, colour, hit["strand"],
                    head_bp=min(850, (right - left) * 0.3), ec="white", lw=0.4)
            text = f"{hit['aligned']/1000:.1f} kb · {hit['identity']*100:.1f}%"
            if (right - left) > length * 0.24:
                # Wide enough to carry the label inside, which is where it reads best.
                ax_contigs.text((left + right) / 2, yy, text, va="center", ha="center",
                                fontsize=5.7, color="white", fontweight="bold", zorder=5)
            elif right + length * (0.012 + 0.009 * len(text)) < island_end:
                ax_contigs.text(right + length * 0.012, yy, text, va="center", ha="left",
                                fontsize=5.7, color=P.INK)
            else:
                # A contig ending near the right edge would push its label off the axis
                # and into the legend, so it goes on the other side instead. The width
                # estimate above is deliberately generous: a label that collides with
                # the legend is far worse than one placed left a little early.
                ax_contigs.text(left - length * 0.012, yy, text, va="center", ha="right",
                                fontsize=5.7, color=P.INK)
        span = (len(hits) - 1) * 0.50
        recovered = P.union_bp(hits, island_start, island_end)
        ticks.append(y - span / 2)
        labels.append(f"{name}\n{len(hits)} contig{'s' if len(hits) != 1 else ''} "
                      f"· {recovered/length*100:.1f}%")
        y -= span + 0.70
    ax_contigs.set_yticks(ticks)
    ax_contigs.set_yticklabels(labels, fontsize=6.0)
    ax_contigs.set_xlim(*xlim)
    ax_contigs.set_ylim(y + 0.46, 0.40)
    ax_contigs.tick_params(axis="y", length=0)
    ax_contigs.tick_params(axis="x", labelbottom=False, length=0)

    # --- depth --------------------------------------------------------------
    ax_depth = figure.add_subplot(grid[2], sharex=ax_island)
    if xs:
        ax_depth.fill_between(xs, ys, step="mid", color=COVERAGE_COLOUR, alpha=0.26, lw=0)
        ax_depth.plot(xs, ys, drawstyle="steps-mid", color=COVERAGE_COLOUR, lw=0.5)
        median = sorted(ys)[len(ys) // 2]
        ax_depth.axhline(median, color="#b5462f", lw=0.7, ls=(0, (4, 2)), zorder=4)
        ax_depth.set_ylim(0, max(ys) * 1.08)
        ax_depth.text(1.006, median / (max(ys) * 1.12), f"median\n{median:.0f}×",
                      transform=ax_depth.transAxes, va="center", ha="left", fontsize=5.9,
                      color="#b5462f", fontweight="bold", linespacing=1.2)
    else:
        ax_depth.text(0.5, 0.5, "no depth track found", transform=ax_depth.transAxes,
                      ha="center", va="center", fontsize=6.2, color=P.MUTED, style="italic")
        ax_depth.set_ylim(0, 1)
    for boundary in (island_start, island_end):
        ax_depth.axvline(boundary, color=P.MUTED, lw=0.5, ls=":")
    ax_depth.set_xlim(*xlim)
    ax_depth.grid(axis="y", color=P.GRID, lw=0.45)
    ax_depth.set_axisbelow(True)
    ax_depth.tick_params(axis="y", labelsize=6.2)
    steps = max(1, int(round(length / 10000)))
    xticks = [island_start + i * 10000 for i in range(steps)] + [island_end]
    ax_depth.set_xticks(xticks)
    ax_depth.set_xticklabels([f"{i*10}" for i in range(steps)] + [f"{length/1000:.1f}"],
                             fontsize=6.5)
    ax_depth.set_xlabel("position within the pks island (kb)", fontsize=6.6, labelpad=1)
    side(ax_depth, "read depth\n(×, 50 bp bins)", y=0.5, x=-0.085)

    handles = [mlines.Line2D([], [], marker="s", ls="", ms=4.6,
                             mfc=P.CLASS_COL[key], mec="none", label=P.CLASS_LABEL[key])
               for key in P.CLASS_ORDER]
    figure.legend(handles=handles, loc="center left", bbox_to_anchor=(0.755, 0.52),
                  ncol=1, fontsize=5.8, frameon=False, handletextpad=0.4,
                  labelspacing=0.42, borderpad=0.0)

    # One line above the figure: the sample, nothing else. No title and no summary
    # statistics -- every number they would repeat is already in the panel row labels
    # (contig count, share of the island) or on the depth axis (median), and a figure
    # going into a deck or a manuscript carries its caption there rather than baked in.
    # The numbers still go to stdout, where they are useful for checking a run.
    if args.label:
        figure.text(0.185, 0.962, args.label, ha="left", va="top",
                    fontsize=7.6, color=P.MUTED)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(out, dpi=args.dpi)
    for extra in args.also or []:
        figure.savefig(out.with_suffix(f".{extra}"), dpi=args.dpi)
    plt.close(figure)

    for name, hits, _ in tracks:
        recovered = P.union_bp(hits, island_start, island_end)
        lowest = min((hit["identity"] for hit in hits), default=0)
        print(f"  {name:24s} {len(hits):2d} contigs  {recovered:6,d} bp  "
              f"{recovered/length*100:6.2f}%  min identity {lowest*100:.1f}%")
    if xs:
        print(f"  depth median {sorted(ys)[len(ys)//2]:.0f}x  max {max(ys):.0f}x")
    print(f"\n  wrote {out}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--mode", choices=["wgs", "metagenome"], default=None,
                        help="which layout to discover tracks from under --sample-dir. Not needed when --paf is given.")
    parser.add_argument("--paf", action="append", default=None, metavar="LABEL=PATH",
                        help="a contig track, given explicitly and repeatable, in the order to draw. Overrides discovery; this is what a pipeline task uses, since it stages its inputs flat.")
    parser.add_argument("--gff", required=True)
    parser.add_argument("--out", required=True, help="output figure (.pdf or .png)")
    parser.add_argument("--also", nargs="*", default=None,
                        help="extra formats to write beside --out, e.g. --also png")
    parser.add_argument("--label", default=None,
                        help="the one line above the figure. Defaults to the sample "
                             "directory's name; pass an empty string for no line at all, "
                             "which is what a figure for a cohort whose identifiers are "
                             "controlled should use.")
    parser.add_argument("--coverage", default="alignment/coverage_50bp_bins.bedGraph",
                        help="depth track. Relative to --sample-dir, or an absolute path when a task stages it flat.")
    parser.add_argument("--contig", default="NC_017628.1")
    parser.add_argument("--island-start", type=int, default=2193827)
    parser.add_argument("--island-end", type=int, default=2244594)
    parser.add_argument("--dpi", type=int, default=400)
    args = parser.parse_args(argv)

    if not args.paf and not args.mode:
        sys.exit("give either --paf label=path (repeatable) or --mode for discovery")
    if not Path(args.sample_dir).is_dir():
        sys.exit(f"not a directory: {args.sample_dir}")
    if args.label is None:
        args.label = Path(args.sample_dir).name
    genes = P.read_gff_genes(args.gff)
    if not genes:
        sys.exit(f"no genes read from {args.gff}")
    return build(args, genes)


if __name__ == "__main__":
    sys.exit(main())
