#!/usr/bin/env python3
"""Draw the pipeline schematic as a publication figure.

    python3 workflow_logo/make_figure.py

Two things drove the layout. One box per process gave thirty near-empty rectangles and a
figure twice as tall as it was wide, so processes are grouped into the four phases each
lane actually has and the box is sized to hold a list rather than one line of text. And
the trunk runs across the page instead of down it, because four sequential steps with no
branching do not need four rows.

The result is close to square and the type is large enough to read at the width a README
displays it, which was the other complaint.

Conventions, the ones journal figures tend to share:

  * one type family, four sizes, and no size used for two different jobs
  * hairline rules; no shadow, no gradient, no bevel
  * lane identity carried by position AND colour, so the figure survives greyscale and
    colour-vision deficiency
  * dashed outline means optional or conditional, everywhere, with no exceptions

Palette checked with the dataviz validator (categorical, light surface): worst adjacent
CVD separation dE 16.6 deutan, 23.1 normal vision, all four checks pass.
"""
from pathlib import Path
from xml.sax.saxutils import escape

# ---------------------------------------------------------------- design tokens
INK     = "#161c21"
BODY    = "#39454e"
MUTED   = "#6f7a83"
RULE    = "#d2dade"
SURFACE = "#ffffff"
PAGE    = "#fcfcfb"
SHARED  = "#44525c"
METAGEN = "#2166ac"
TUMOUR  = "#b2182b"
GATE    = "#b8892b"

FAM  = "Helvetica Neue, Helvetica, Arial, sans-serif"
MONO = "ui-monospace, SFMono-Regular, Menlo, Consolas, monospace"

SZ_LANE  = 13.5   # lane headers
SZ_PHASE = 13     # phase headers
SZ_STEP  = 12.5   # process names
SZ_NOTE  = 10.5   # annotations
SZ_LEG   = 10.5   # legend

W_BOLD, W_REG = 0.552, 0.497

# ---------------------------------------------------------------- page geometry
W       = 930
MARGIN  = 30
GUTTER  = 46
COL     = (W - 2 * MARGIN - GUTTER) // 2      # 412
LEFT_X  = MARGIN
RIGHT_X = MARGIN + COL + GUTTER

PAD      = 14
HEAD_H   = 27     # phase header strip
STEP_H   = 21     # one process line
NOTE_H   = 16
PHASE_GAP = 22
RADIUS   = 3

_out = []
def emit(s): _out.append(s)


def txt(x, y, s, size, fill, family=FAM, weight="400", anchor="start", tracking=None):
    tr = f' letter-spacing="{tracking}"' if tracking else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{family}" font-size="{size}" '
            f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{tr}>'
            f'{escape(s)}</text>')


def check(s, size, ratio, width, where):
    if len(s) * size * ratio > width:
        print(f"  ! may overflow in {where}: {s[:58]}")


def phase(x, y, w, title, note, steps, colour, dashed=False):
    """A phase: a coloured header strip over a list of the processes it contains.

    steps: (name, optional trailing tag, is_gate, is_optional)
    """
    h = HEAD_H + (NOTE_H if note else 0) + len(steps) * STEP_H + 12
    dash = ' stroke-dasharray="5 3"' if dashed else ""
    emit(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{RADIUS}" fill="{SURFACE}" '
         f'stroke="{colour}" stroke-width="1"{dash}/>')
    # header strip
    emit(f'<path d="M{x} {y + RADIUS} a{RADIUS} {RADIUS} 0 0 1 {RADIUS} -{RADIUS} '
         f'h{w - 2 * RADIUS} a{RADIUS} {RADIUS} 0 0 1 {RADIUS} {RADIUS} '
         f'v{HEAD_H - RADIUS} h-{w} z" fill="{colour}"/>')
    emit(txt(x + PAD, y + 18.5, title, SZ_PHASE, "#ffffff", weight="700"))
    check(title, SZ_PHASE, W_BOLD, w - 2 * PAD, title)

    yy = y + HEAD_H
    if note:
        yy += NOTE_H
        emit(txt(x + PAD, yy - 4, note, SZ_NOTE, MUTED))
        check(note, SZ_NOTE, W_REG, w - 2 * PAD, title)

    for name, tag, is_gate, is_opt in steps:
        yy += STEP_H
        base = yy - 5
        mark_x = x + PAD + 3
        if is_gate:
            emit(f'<path d="M{mark_x} {base - 7} l4.5 4.5 l-4.5 4.5 l-4.5 -4.5 z" '
                 f'fill="{GATE}"/>')
        elif is_opt:
            emit(f'<circle cx="{mark_x}" cy="{base - 3.5}" r="3.4" fill="none" '
                 f'stroke="{colour}" stroke-width="1.2" stroke-dasharray="2 1.6"/>')
        else:
            emit(f'<circle cx="{mark_x}" cy="{base - 3.5}" r="3.4" fill="{colour}"/>')
        fill = "#8a6519" if is_gate else BODY
        emit(txt(x + PAD + 15, base, name, SZ_STEP, fill,
                 family=MONO if not is_gate else FAM,
                 weight="600" if is_gate else "500"))
        if tag:
            emit(txt(x + w - PAD, base, tag, SZ_NOTE, MUTED, anchor="end"))
        check(name + "   " + (tag or ""), SZ_STEP, W_REG, w - 2 * PAD - 15, name)
    return y + h


def arrow(x, y1, y2, colour):
    emit(f'<path d="M{x} {y1} V{y2 - 6}" stroke="{colour}" stroke-width="1.4" '
         f'marker-end="url(#tip)"/>')


# ---------------------------------------------------------------- content
TRUNK = [("extractReads", "unmapped reads"),
         ("filterReads", "fastp"),
         ("mapReads", "GRCh38 → T2T+PhiX")]

METAGENOME = [
    ("ASSEMBLE AND BIN", "de novo from the complete stream, never a filtered subset", [
        ("megahitAssemble", None, False, False),
        ("contig-count gate", "0 contigs = a true negative", True, False),
        ("alignToContigs → jgiContigDepths", None, False, False),
        ("metabat2Bin", "fans out per bin", False, False)]),
    ("CHARACTERISE EACH BIN", "annotation runs first, so locus_tags match downstream", [
        ("prokkaAnnotate", None, False, False),
        ("checkm2Predict · gtdbtkClassify", "completeness, taxonomy", False, False),
        ("hmmsearchClb · genomadProphages", "clb HMM, proviruses", False, False),
        ("pks-positive bin filter", "distinct clb genes", True, False)]),
    ("PLACE IT IN CONTEXT", "what sits beside the island, and in which organism", [
        ("extractGenomicContext", "±50 kb, integrase, tRNA", False, False),
        ("magSummaryTable", None, False, False)]),
]

TUMOUR_LANE = [
    ("PROFILE AND SCORE", "counts the island; does not yet try to rebuild it", [
        ("prefilter", "optional, narrows counting only", False, True),
        ("pksProfilerAlign", "Bowtie2 → featureCounts", False, False),
        ("classifyTumorPksEvidence", "reads, genes, breadth", False, False),
        ("evidence-tier gate", "broad ≥ 30/8/7.5%", True, False)]),
    ("REBUILD THE ISLAND", "returns to the complete stream first — assembly never narrows", [
        ("re-join MAPPED_READS", "not the profiled subset", False, False),
        ("targetedPksRecruit", "≥60 bp, ≥90% identity", False, False),
        ("targetedPksMegahit ∥ …Spades", "two assemblers", False, False),
        ("alignTargetedContigs", "vs the IHE3034 island", False, False),
        ("summarizeTargetedPksEvidence", "assembler agreement", False, False)]),
    ("OPTIONAL CONTEXT", "off unless asked for; needs genomad_db", [
        ("tumour contig context", "prokka → clb HMM → geNomad", False, True)]),
]


def render():
    emit(f'<rect width="{W}" height="{{H}}" fill="{PAGE}"/>')
    full = W - 2 * MARGIN

    # ------------------------------------------------ trunk, across the page
    y = 30
    emit(txt(MARGIN, y, "SHARED TRUNK — every sample, every mode", SZ_LANE, SHARED,
             weight="700", tracking="1.2"))
    y += 12
    emit(f'<path d="M{MARGIN} {y} H{W - MARGIN}" stroke="{SHARED}" stroke-width="2"/>')
    y += 16

    n = len(TRUNK) + 1
    step_gap = 26
    bw = (full - (n - 1) * step_gap) / n
    bh = 46
    cells = [("Sequencing input", "BAM · CRAM · FASTQ")] + TRUNK
    for i, (name, note) in enumerate(cells):
        bx = MARGIN + i * (bw + step_gap)
        emit(f'<rect x="{bx}" y="{y}" width="{bw}" height="{bh}" rx="{RADIUS}" '
             f'fill="{SURFACE}" stroke="{RULE}" stroke-width="1"/>')
        emit(f'<path d="M{bx + 0.5} {y + RADIUS} v{bh - 2 * RADIUS}" stroke="{SHARED}" '
             f'stroke-width="3" stroke-linecap="round"/>')
        emit(txt(bx + 12, y + 20, name, SZ_STEP, INK,
                 family=FAM if i == 0 else MONO, weight="600"))
        emit(txt(bx + 12, y + 36, note, SZ_NOTE, MUTED))
        check(name, SZ_STEP, W_BOLD, bw - 20, name)
        if i < n - 1:
            ax = bx + bw
            emit(f'<path d="M{ax + 5} {y + bh / 2} H{ax + step_gap - 6}" stroke="{RULE}" '
                 f'stroke-width="1.4" marker-end="url(#tip)"/>')
    y += bh + 16

    # the branch point
    emit(f'<rect x="{MARGIN}" y="{y}" width="{full}" height="44" rx="{RADIUS}" '
         f'fill="#eef2f4" stroke="{SHARED}" stroke-width="1.6"/>')
    emit(txt(MARGIN + PAD, y + 19, "MAPPED_READS", SZ_STEP, INK, family=MONO, weight="700"))
    emit(txt(MARGIN + PAD + 128, y + 19,
             "the complete host-depleted stream — both lanes branch from here, and the "
             "tumour lane returns to it", SZ_NOTE, MUTED))
    emit(txt(MARGIN + PAD, y + 35,
             "profiling narrows what is counted · assembly always reads the whole stream",
             SZ_NOTE, SHARED))
    branch_y = y + 44

    # ------------------------------------------------ two lanes
    lane_top = branch_y + 74
    for cx, colour, cond in ((LEFT_X + COL / 2, METAGEN, "sample_type = metagenome"),
                             (RIGHT_X + COL / 2, TUMOUR, "sample_type = tumor_wgs")):
        emit(f'<path d="M{W / 2} {branch_y} V{branch_y + 20} H{cx} V{lane_top - 34}" '
             f'fill="none" stroke="{colour}" stroke-width="1.4" marker-end="url(#tip)"/>')
        pill_w = len(cond) * SZ_NOTE * W_REG + 26
        emit(f'<rect x="{cx - pill_w / 2}" y="{branch_y + 24}" width="{pill_w}" height="19" '
             f'rx="9.5" fill="{PAGE}" stroke="{colour}" stroke-width="1"/>')
        emit(txt(cx, branch_y + 37.5, cond, SZ_NOTE, colour, family=MONO, anchor="middle"))

    for x, label, colour in ((LEFT_X, "METAGENOME LANE", METAGEN),
                             (RIGHT_X, "TUMOUR LANE", TUMOUR)):
        emit(txt(x, lane_top - 12, label, SZ_LANE, colour, weight="700", tracking="1.2"))
        emit(f'<path d="M{x} {lane_top - 5} H{x + COL}" stroke="{colour}" stroke-width="2"/>')

    ends = []
    for x, phases, colour in ((LEFT_X, METAGENOME, METAGEN),
                              (RIGHT_X, TUMOUR_LANE, TUMOUR)):
        yy = lane_top + 12
        for i, (title, note, steps) in enumerate(phases):
            optional = title == "OPTIONAL CONTEXT"
            bottom = phase(x, yy, COL, title, note, steps, colour, dashed=optional)
            if i < len(phases) - 1:
                arrow(x + COL / 2, bottom, bottom + PHASE_GAP, colour)
            yy = bottom + PHASE_GAP
        ends.append(yy - PHASE_GAP)

    # ------------------------------------------------ shared tail
    low = max(ends)
    tail_y = low + 50
    for x, colour in ((LEFT_X + COL / 2, METAGEN), (RIGHT_X + COL / 2, TUMOUR)):
        emit(f'<path d="M{x} {ends[0] if x < W / 2 else ends[1]} V{low + 24} H{W / 2} '
             f'V{tail_y - 6}" fill="none" stroke="{colour}" stroke-width="1.4" '
             f'marker-end="url(#tip)"/>')
    emit(f'<rect x="{MARGIN}" y="{tail_y}" width="{full}" height="46" rx="{RADIUS}" '
         f'fill="{SURFACE}" stroke="{SHARED}" stroke-width="1.6"/>')
    emit(f'<path d="M{MARGIN + 0.5} {tail_y + RADIUS} v{46 - 2 * RADIUS}" stroke="{SHARED}" '
         f'stroke-width="3" stroke-linecap="round"/>')
    emit(txt(MARGIN + PAD, tail_y + 20, "strainTyping", SZ_STEP, INK, family=MONO,
             weight="700"))
    emit(txt(MARGIN + PAD + 106, tail_y + 20,
             "MLST → sequence type · clonal complex · phylogroup", SZ_NOTE, MUTED))
    emit(txt(MARGIN + PAD, tail_y + 36,
             "assembled units from either lane — a partial allele profile is reported as "
             "insufficient_loci, never guessed", SZ_NOTE, MUTED))

    # ------------------------------------------------ legend
    ly = tail_y + 46 + 32
    emit(f'<path d="M{MARGIN} {ly - 18} H{W - MARGIN}" stroke="{RULE}" stroke-width="1"/>')
    x = MARGIN
    for kind, label in (("run", "runs by default"), ("opt", "optional"),
                        ("gate", "decision or filter")):
        if kind == "gate":
            emit(f'<path d="M{x + 4} {ly - 8} l4.5 4.5 l-4.5 4.5 l-4.5 -4.5 z" fill="{GATE}"/>')
        elif kind == "opt":
            emit(f'<circle cx="{x + 4}" cy="{ly - 3.5}" r="3.4" fill="none" stroke="{MUTED}" '
                 f'stroke-width="1.2" stroke-dasharray="2 1.6"/>')
        else:
            emit(f'<circle cx="{x + 4}" cy="{ly - 3.5}" r="3.4" fill="{MUTED}"/>')
        emit(txt(x + 15, ly, label, SZ_LEG, MUTED))
        x += 15 + len(label) * SZ_LEG * W_REG + 30
    emit(txt(W - MARGIN, ly, "pksProfiler v0.0.2", SZ_LEG, MUTED, anchor="end",
             weight="600"))
    return ly + 22


def main():
    height = render()
    defs = ('<defs><marker id="tip" viewBox="0 0 8 8" refX="6.6" refY="4" markerWidth="5.5" '
            'markerHeight="5.5" orient="auto-start-reverse">'
            '<path d="M0 1.3 L7 4 L0 6.7 z" fill="context-stroke" stroke="none"/></marker></defs>')
    body = "\n".join(_out).replace("{H}", str(height))
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{height}" '
           f'viewBox="0 0 {W} {height}" role="img" '
           f'aria-label="pksProfiler workflow: a shared trunk of read extraction, quality '
           f'filtering and human-read removal produces MAPPED_READS, the complete '
           f'host-depleted stream, which branches into a metagenome lane that assembles and '
           f'bins the community and a tumour lane that scores read evidence and then returns '
           f'to the complete stream to reassemble the island; assembled output from either '
           f'lane is strain typed">\n{defs}\n{body}\n</svg>\n')
    out = Path(__file__).resolve().parent / "v2.svg"
    out.write_text(svg)
    print(f"wrote {out}  {W}x{height}  {len(svg) / 1024:.1f} KB")


if __name__ == "__main__":
    main()
