#!/usr/bin/env python3
"""Shared drawing conventions for pks island figures.

Follows how gene-cluster comparison figures are normally drawn in the literature
(clinker / clustermap.js, Gilchrist & Chooi 2021; genoPlotR; EasyFig):

  - genes as strand-aware arrows, not plain rectangles
  - fill colour encodes functional class, with a legend
  - aligned segments carry their percent identity
  - a scale bar instead of a full numeric axis on the track panels

Colour assignments for the 19 clb genes follow the standard description of the island:
three PKS (clbC, clbI, clbO), three NRPS (clbH, clbJ, clbN), two hybrid NRPS/PKS
(clbB, clbK), a regulator (clbR), a self-resistance peptidase (clbS), and the remaining
accessory / tailoring / editing proteins.
"""
import re
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D

CLASS = {
    "clbC": "PKS", "clbI": "PKS", "clbO": "PKS",
    "clbH": "NRPS", "clbJ": "NRPS", "clbN": "NRPS",
    "clbB": "hybrid", "clbK": "hybrid",
    "clbR": "regulator", "clbS": "resistance",
}
CLASS_COL = {
    "PKS": "#2f6f9f", "NRPS": "#8c5aa8", "hybrid": "#1f9d8f",
    "regulator": "#d08c27", "resistance": "#b5462f", "accessory": "#aab4bd",
}
CLASS_ORDER = ["PKS", "NRPS", "hybrid", "accessory", "regulator", "resistance"]
CLASS_LABEL = {
    "PKS": "polyketide synthase", "NRPS": "non-ribosomal peptide synthase",
    "hybrid": "hybrid NRPS/PKS", "accessory": "accessory / tailoring / editing",
    "regulator": "regulator (clbR)", "resistance": "self-resistance (clbS)",
}

INK, MUTED, GRID = "#1f2b38", "#5b6b7a", "#dde3e8"


def gene_class(name):
    return CLASS.get(name, "accessory")


def read_gff_genes(path):
    """(name, start, end, strand) for every gene, in genomic order."""
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"gene=([^;]+)", f[8])
        if m:
            out.append((m.group(1), int(f[3]), int(f[4]), f[6]))
    return sorted(out, key=lambda g: g[1])


def arrow(ax, x0, x1, y, h, color, strand="+", lw=0.5, ec="white", zorder=3,
          min_head_frac=0.30, head_bp=None):
    """A gene/feature arrow. Head length adapts so short genes stay recognisable."""
    span = x1 - x0
    if span <= 0:
        return
    head = head_bp if head_bp is not None else span * 0.35
    head = min(head, span * min_head_frac if span * min_head_frac > 0 else head)
    if strand == "-":
        pts = [(x1, y + h / 2), (x0 + head, y + h / 2), (x0, y),
               (x0 + head, y - h / 2), (x1, y - h / 2)]
    else:
        pts = [(x0, y - h / 2), (x1 - head, y - h / 2), (x1, y),
               (x1 - head, y + h / 2), (x0, y + h / 2)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor=ec,
                         linewidth=lw, zorder=zorder))


def gene_track(ax, genes, y=0.0, h=0.55, label_min_bp=900, fontsize=6.0,
               label_rotation=0, head_bp=1100):
    """Draw the clb genes as coloured, strand-aware arrows."""
    for name, s, e, strand in genes:
        arrow(ax, s, e, y, h, CLASS_COL[gene_class(name)], strand, head_bp=head_bp)
        if (e - s) >= label_min_bp:
            ax.text((s + e) / 2, y + h * 0.78, name.replace("clb", ""),
                    ha="center", va="bottom", fontsize=fontsize, color=INK,
                    rotation=label_rotation)


def gene_legend(ax, loc="upper left", fontsize=5.9, ncol=3, bbox=None, present=None):
    keys = [k for k in CLASS_ORDER if present is None or k in present]
    h = [Line2D([], [], marker="s", ls="", ms=5, mfc=CLASS_COL[k], mec="none",
                label=CLASS_LABEL[k]) for k in keys]
    return ax.legend(handles=h, loc=loc, fontsize=fontsize, frameon=False, ncol=ncol,
                     handletextpad=0.35, columnspacing=1.1, borderpad=0.1,
                     labelspacing=0.3, bbox_to_anchor=bbox)


def scale_bar(ax, x0, y, length_bp, label, color=INK, fontsize=6.2, lw=1.1):
    ax.plot([x0, x0 + length_bp], [y, y], color=color, lw=lw,
            solid_capstyle="butt", clip_on=False, zorder=6)
    for xx in (x0, x0 + length_bp):
        ax.plot([xx, xx], [y - 0.045, y + 0.045], color=color, lw=lw,
                clip_on=False, zorder=6)
    ax.text(x0 + length_bp / 2, y + 0.085, label, ha="center", va="bottom",
            fontsize=fontsize, color=color, clip_on=False)


def paf_hits(path, island_start, island_end, min_identity=0.90, min_mapq=20,
             min_aligned=0):
    """PAF records overlapping the island, filtered the way the pipeline filters them.

    The pipeline's summarize_tumor_pks_contigs.py keeps hits at >=90% identity and
    MAPQ >=20; applying the same thresholds here keeps the figure consistent with the
    numbers the pipeline reports.
    """
    out = []
    try:
        txt = open(path).read()
    except OSError:
        return out
    for line in txt.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        try:
            qs, qe = int(f[2]) + 1, int(f[3])
            ts, te = int(f[7]) + 1, int(f[8])
            nmatch, blocklen, mapq = int(f[9]), int(f[10]), int(f[11])
        except ValueError:
            continue
        if blocklen <= 0 or te < island_start or ts > island_end:
            continue
        ident = nmatch / blocklen
        if ident < min_identity or mapq < min_mapq:
            continue
        aligned = min(te, island_end) - max(ts, island_start) + 1
        if aligned < min_aligned:
            continue
        out.append({"contig": f[0], "qlen": int(f[1]), "qs": qs, "qe": qe,
                    "ts": ts, "te": te, "strand": f[4], "identity": ident,
                    "mapq": mapq, "aligned": aligned, "blocklen": blocklen})
    return sorted(out, key=lambda r: r["ts"])


def union_bp(hits, start, end):
    cov = set()
    for h in hits:
        a, b = max(h["ts"], start), min(h["te"], end)
        if b >= a:
            cov.update(range(a, b + 1))
    return len(cov)
