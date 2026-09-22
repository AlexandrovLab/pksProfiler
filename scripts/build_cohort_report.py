#!/usr/bin/env python3
"""One page for the whole cohort: clb genes, breadth, evidence tier, contig support.

The pipeline writes a directory per sample. At two thousand samples that is two
thousand directories and no view of the cohort, which is the thing anyone reads first.

This reads what the pipeline already produced and arranges it:

    cohort/gene_counts/pks.gene.counts.align.txt        clb reads per gene, per sample
    by_sample/<sample>/read_evidence.tsv                tier, reads, genes, breadth
    by_sample/<sample>/contigs/final_evidence/final_pks_evidence.tsv  contig support

It derives nothing. The tier comes from the pipeline's own classification and the
breadth from the depth it already computed; the hand-built snapshot this replaces
(pksProfiler_analysis/mutographs_evidence_snapshots/build_mutographs_evidence.py)
recomputed both, with the v0.0.1 thresholds, and so could disagree with the run it
was describing.

The layout is that snapshot's, kept deliberately: same blue buckets for the gene
matrix, same log-scaled read bar, same three breadth bars, same evidence badge. Its
`contig pending` placeholder is now filled in, in the row where it always belonged.

    build_cohort_report.py --results results --output cohort/pks_cohort_report.html
"""
import argparse
import csv
import html
import math
from pathlib import Path

GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]

TIERS = ["negative", "localized_indeterminate", "multi_gene", "broad_island",
         "extensive_island"]
TIER_RANK = {name: index for index, name in enumerate(TIERS)}

# The snapshot's palette. Blue buckets for the gene matrix (a magnitude, so one hue
# light to dark), greys for breadth, and an ordered ramp for the evidence badge. The
# snapshot had four classes; v0.0.2 has five, so multi_gene takes a step between
# broad and indeterminate rather than borrowing another hue.
BUCKETS = ((10, "#2166ac"), (3, "#67a9cf"), (1, "#d1e5f0"))
EMPTY = "#f2f2f2"
TIER_FILL = {"extensive_island": "#9e2a2b", "broad_island": "#e09f3e",
             "multi_gene": "#f0c987", "localized_indeterminate": "#8c8c8c",
             "negative": "#d9d9d9", "not_classified": "#d9d9d9"}
TIER_INK = {"multi_gene": "#3b2a12", "negative": "#4d4d4d",
            "not_classified": "#4d4d4d"}
TIER_LABEL = {"extensive_island": "Extensive", "broad_island": "Broad",
              "multi_gene": "Multi-gene", "localized_indeterminate": "Indeterminate",
              "negative": "Negative", "not_classified": "Not classified"}
AGREEMENT_FILL = {"concordant": "#4d7c5f", "discordant": "#9e2a2b",
                  "single_assembler_only": "#e09f3e", "no_contig_support": "#bdbdbd"}
# The badge is 88px at 9px Arial; `single_assembler_only` is 21 characters and would
# spill out of it. The full value stays in the TSV and in the hover text.
AGREEMENT_LABEL = {"concordant": "Concordant", "discordant": "Discordant",
                   "single_assembler_only": "One assembler",
                   "no_contig_support": "No support"}

# geometry, following the snapshot and extended rightwards for the contig columns
LEFT, GENE_X, CELL, ROW_H, TOP = 15, 205, 15, 13, 142
READS_X, B1_X, CLASS_X = 510, 630, 725
CONTIG_X, COV_X, AGREE_X, FIG_X = 825, 925, 1045, 1145
WIDTH = 1270


def bucket(count):
    for floor, colour in BUCKETS:
        if count >= floor:
            return colour
    return EMPTY


def read_tsv(path):
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gene_counts(matrix_path):
    rows = read_tsv(matrix_path)
    if not rows:
        return {}
    samples = [c for c in rows[0].keys() if c != "Gene"]
    counts = {s: {g: 0 for g in GENES} for s in samples}
    for row in rows:
        if row.get("Gene") not in GENES:
            continue
        for sample in samples:
            try:
                counts[sample][row["Gene"]] = int(float(row[sample] or 0))
            except (TypeError, ValueError):
                counts[sample][row["Gene"]] = 0
    return counts


def number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def collect(results):
    counts = gene_counts(results / "cohort/gene_counts/pks.gene.counts.align.txt")
    by_sample = results / "by_sample"
    samples = sorted({p.name for p in by_sample.iterdir() if p.is_dir()} | set(counts)) \
        if by_sample.is_dir() else sorted(counts)

    records = []
    for sample in samples:
        evidence = read_tsv(by_sample / sample / "read_evidence.tsv")
        evidence = evidence[0] if evidence else {}
        # pks_targeted.nf publishes this to contigs/final_evidence/, not contigs/.
        # Reading the wrong path did not error -- read_tsv returned nothing, `assembled`
        # went False, and the report simply showed no contig columns for samples that had
        # assembled perfectly well. Silent, and only visible by noticing an absence.
        contigs = read_tsv(
            by_sample / sample / "contigs/final_evidence/final_pks_evidence.tsv")
        contigs = contigs[0] if contigs else {}
        records.append({
            "sample": sample,
            "genes": counts.get(sample, {g: 0 for g in GENES}),
            "tier": evidence.get("read_evidence", "not_classified"),
            "pks_reads": int(number(evidence.get("pks_reads"))),
            "genes_detected": int(number(evidence.get("clb_genes_detected"))),
            "breadth_1x": number(evidence.get("island_breadth_1x")),
            "breadth_2x": number(evidence.get("island_breadth_2x")),
            "breadth_3x": number(evidence.get("island_breadth_3x")),
            "assembled": bool(contigs),
            "structural": contigs.get("final_structural_evidence", ""),
            "agreement": contigs.get("assembler_agreement", ""),
            "megahit_contigs": int(number(contigs.get("megahit_supporting_contigs"))),
            "metaspades_contigs": int(number(contigs.get("metaspades_supporting_contigs"))),
            "megahit_cov": number(contigs.get("megahit_reference_coverage")),
            "metaspades_cov": number(contigs.get("metaspades_reference_coverage")),
            "figures": [p.name for p in
                        sorted((by_sample / sample / "figures").glob("*"))
                        if p.is_file() and not p.name.startswith(".")],
        })

    records.sort(key=lambda r: (TIER_RANK.get(r["tier"], -1), r["genes_detected"],
                                r["pks_reads"], r["breadth_1x"]), reverse=True)
    return records


def esc(value):
    return html.escape(str(value), quote=True)


def figure_anchor(sample, filename):
    """Stable id shared by the row link and the panel it reveals."""
    safe = "".join(c if c.isalnum() else "-" for c in f"{sample}-{filename}")
    return f"fig-{safe}"


def svg_matrix(records):
    height = TOP + ROW_H * len(records) + 40
    maxlog = max([math.log10(1 + r["pks_reads"]) for r in records] or [1]) or 1

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" '
           f'viewBox="0 0 {WIDTH} {height}" role="img" '
           f'aria-label="clb gene detection, evidence tier and contig support per sample">',
           '<rect width="100%" height="100%" fill="white"/>',
           '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111}'
           '.h{font-size:13px;font-weight:bold}.n{font-size:9px}'
           '.note{font-size:10px;fill:#444}'
           '.row:hover .zebra{fill:#eef4fb}'
           '.cell{stroke:white}.cell:hover{stroke:#111;stroke-width:1.2}</style>',
           '<text x="15" y="22" class="h" font-size="16">pks cohort report</text>',
           f'<text x="15" y="42" class="note">{len(records)} samples, ordered by evidence '
           'tier, then genes detected, then reads</text>']

    legend = [("#2166ac", "≥10 reads/gene", 15), ("#67a9cf", "3–9", 126),
              ("#d1e5f0", "1–2", 192), (EMPTY, "0", 258)]
    for colour, label, x in legend:
        stroke = ' stroke="#ccc"' if colour == EMPTY else ""
        out.append(f'<rect x="{x}" y="55" width="11" height="11" fill="{colour}"{stroke}/>'
                   f'<text x="{x + 16}" y="64" class="n">{label}</text>')

    out.append('<text x="15" y="100" class="h">Sample</text>')
    for index, gene in enumerate(GENES):
        x = GENE_X + index * CELL + 10
        out.append(f'<text x="{x}" y="103" class="n" transform="rotate(-55 {x} 103)">{gene}</text>')
    for label, x in (("PKS reads", READS_X), ("Breadth ≥1×", B1_X),
                     ("Read evidence", CLASS_X), ("Contigs", CONTIG_X),
                     ("Island recovered", COV_X), ("Assemblers", AGREE_X),
                     ("Figures", FIG_X)):
        out.append(f'<text x="{x}" y="100" class="h">{label}</text>')
    out.append(f'<line x1="15" y1="113" x2="{WIDTH - 15}" y2="113" stroke="black"/>')

    for index, record in enumerate(records):
        y = TOP + index * ROW_H
        tip = (f'{record["sample"]} · {TIER_LABEL.get(record["tier"], record["tier"])} '
               f'· {record["pks_reads"]} reads across {record["genes_detected"]} genes '
               f'· breadth {100 * record["breadth_1x"]:.1f}% at ≥1×, '
               f'{100 * record["breadth_2x"]:.1f}% at ≥2×, '
               f'{100 * record["breadth_3x"]:.1f}% at ≥3×')
        # <title> as well as data-tip. JupyterLab and other HTML viewers render files
        # in a sandboxed iframe with JavaScript disabled, and a hover that only works
        # in a full browser is a hover that does not work.
        out.append(f'<g class="row" data-tip="{esc(tip)}"><title>{esc(tip)}</title>')
        zebra = "#fafafa" if index % 2 else "white"
        out.append(f'<rect class="zebra" x="1" y="{y - ROW_H + 2}" width="{WIDTH - 2}" '
                   f'height="{ROW_H}" fill="{zebra}"/>')
        out.append(f'<text x="{LEFT}" y="{y}" class="n">{esc(record["sample"])}</text>')

        for gene_index, gene in enumerate(GENES):
            value = record["genes"].get(gene, 0)
            cell_tip = f'{record["sample"]} · {gene} · {value} reads'
            out.append(f'<rect class="cell" x="{GENE_X + gene_index * CELL}" y="{y - 9}" '
                       f'width="13" height="11" fill="{bucket(value)}" '
                       f'data-tip="{esc(cell_tip)}"><title>{esc(cell_tip)}</title></rect>')

        bar = 78 * math.log10(1 + record["pks_reads"]) / maxlog
        out.append(f'<rect x="{READS_X}" y="{y - 7}" width="{bar:.1f}" height="7" fill="#4d4d4d"/>'
                   f'<text x="{READS_X + 105}" y="{y}" class="n" text-anchor="end">'
                   f'{record["pks_reads"]}</text>')

        # Only >=1x is drawn. The three breadths are nested by definition and only
        # >=1x feeds the tier; >=2x and >=3x are depth detail, so they live in the
        # hover text and the TSV rather than as two shorter copies of the same bar.
        pct = 100 * record["breadth_1x"]
        out.append(f'<rect x="{B1_X}" y="{y - 7}" width="{pct * 0.55:.1f}" height="7" '
                   f'fill="#636363"/><text x="{B1_X + 78}" y="{y}" class="n" '
                   f'text-anchor="end">{pct:.1f}%</text>')

        tier = record["tier"]
        out.append(f'<rect x="{CLASS_X}" y="{y - 9}" width="88" height="11" rx="2" '
                   f'fill="{TIER_FILL.get(tier, "#d9d9d9")}"/>'
                   f'<text x="{CLASS_X + 44}" y="{y}" class="n" text-anchor="middle" '
                   f'style="fill:{TIER_INK.get(tier, "white")}">'
                   f'{esc(TIER_LABEL.get(tier, tier))}</text>')

        if record["assembled"]:
            best = max(record["megahit_cov"], record["metaspades_cov"])
            contig_tip = (f'{record["sample"]} · MEGAHIT {record["megahit_contigs"]} contigs, '
                          f'{100 * record["megahit_cov"]:.1f}% · metaSPAdes '
                          f'{record["metaspades_contigs"]} contigs, '
                          f'{100 * record["metaspades_cov"]:.1f}% · '
                          f'{record["structural"] or "no call"} \u00b7 '
                          f'{record["agreement"] or "no agreement call"}')
            out.append(f'<g data-tip="{esc(contig_tip)}">'
                       f'<title>{esc(contig_tip)}</title>'
                       f'<text x="{CONTIG_X}" y="{y}" class="n">'
                       f'MH {record["megahit_contigs"]} / SP {record["metaspades_contigs"]}</text>'
                       f'<rect x="{COV_X}" y="{y - 7}" width="{100 * best * 0.9:.1f}" '
                       f'height="7" fill="#2166ac"/>'
                       f'<text x="{COV_X + 110}" y="{y}" class="n" text-anchor="end">'
                       f'{100 * best:.1f}%</text>'
                       f'<rect x="{AGREE_X}" y="{y - 9}" width="88" height="11" rx="2" '
                       f'fill="{AGREEMENT_FILL.get(record["agreement"], "#bdbdbd")}"/>'
                       f'<text x="{AGREE_X + 44}" y="{y}" class="n" text-anchor="middle" '
                       f'style="fill:white">'
                       f'{esc(AGREEMENT_LABEL.get(record["agreement"], record["agreement"] or "n/a"))}'
                       f'</text></g>')
        else:
            out.append(f'<text x="{CONTIG_X}" y="{y}" class="n" style="fill:#999">'
                       'assembly not run</text>')

        # Figures open from the row itself. The href is relative to cohort/, and the
        # full path is in the link's title: a sandboxed viewer may refuse to navigate,
        # and then at least the hover says where the file is.
        offset = 0
        for name in record["figures"]:
            label = ("coverage" if "circos" in name
                     else "contigs" if "contig" in name
                     else name.rsplit(".", 1)[0])
            href = f'../by_sample/{record["sample"]}/figures/{name}'
            anchor = figure_anchor(record["sample"], name)
            out.append(f'<a href="#{esc(anchor)}">'
                       f'<title>{esc(href)}</title>'
                       f'<text x="{FIG_X + offset}" y="{y}" class="n" '
                       f'style="fill:#2166ac;text-decoration:underline">{esc(label)}</text>'
                       f'</a>')
            offset += 8 + 5 * len(label)
        if not record["figures"]:
            out.append(f'<text x="{FIG_X}" y="{y}" class="n" style="fill:#bbb">none</text>')
        out.append("</g>")

    note = ("Breadth is the pipeline's own ≥1×/≥2×/≥3× over the "
            "island at MAPQ ≥40. Tiers and contig calls are read from the run's own "
            "outputs, not recomputed here. Tiers rank how strong the evidence is; they are "
            "not a validated positive/negative test.")
    out.append(f'<text x="15" y="{height - 16}" class="note">{note}</text>')
    out.append("</svg>")
    return "\n".join(out)


def figure_panels(records):
    """Hidden panels revealed by :target when a row's figure link is clicked.

    They sit directly under the table and only one is ever open, so this is an inline
    viewer rather than a list of everything. No script: JupyterLab and similar viewers
    sandbox the page and run none, and a report whose figures only open in a full
    browser does not open its figures.
    """
    panels = []
    for record in records:
        for name in record["figures"]:
            anchor = figure_anchor(record["sample"], name)
            href = f'../by_sample/{record["sample"]}/figures/{name}'
            if name.endswith(".svg"):
                view = (f'<object type="image/svg+xml" data="{esc(href)}" width="100%">'
                        f'<a href="{esc(href)}">{esc(name)}</a></object>')
            else:
                # PDFs do not embed reliably, least of all sandboxed.
                view = (f'<p><a href="{esc(href)}">Open {esc(name)}</a>'
                        f' &mdash; PDFs cannot be shown inline here.</p>')
            panels.append(
                f'<div class="figpanel" id="{esc(anchor)}">'
                f'<div class="fighead">{esc(record["sample"])} &middot; {esc(name)}'
                f'<a class="figclose" href="#top">close</a></div>{view}</div>')
    return "".join(panels)


def page(records, results):
    tiers = {tier: sum(1 for r in records if r["tier"] == tier) for tier in TIERS}
    unclassified = sum(1 for r in records if r["tier"] not in TIER_RANK)
    assembled = sum(1 for r in records if r["assembled"])

    tiles = "".join(
        f'<div class="tile"><div class="n">{count}</div>'
        f'<div class="l">{esc(TIER_LABEL[tier])}</div></div>'
        for tier, count in sorted(tiers.items(), key=lambda kv: -TIER_RANK[kv[0]]))
    if unclassified:
        tiles += (f'<div class="tile"><div class="n">{unclassified}</div>'
                  f'<div class="l">Not classified</div></div>')
    tiles += (f'<div class="tile"><div class="n">{assembled}</div>'
              f'<div class="l">with contig analysis</div></div>')

    no_assembly = ("<p class='muted'>Targeted assembly did not run for any sample "
                   "(<code>--tumor_targeted_assembly</code>), so the contig columns are "
                   "empty.</p>" if assembled == 0 else "")

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>pks cohort report</title>
<style>
  body {{ margin:0; padding:28px 16px; background:#f7f7f5; color:#111;
          font:14px/1.5 Arial, Helvetica, sans-serif; }}
  main {{ max-width:1460px; margin:0 auto; }}
  h1 {{ font-size:20px; margin:0 0 4px; }}
  .muted {{ color:#666; font-size:12px; }}
  .tiles {{ display:flex; flex-wrap:wrap; gap:8px; margin:16px 0; }}
  .tile {{ background:#fff; border:1px solid #ddd; border-radius:6px;
           padding:8px 14px; min-width:118px; }}
  .tile .n {{ font-size:19px; font-weight:bold; }}
  .tile .l {{ font-size:11px; color:#666; }}
  .scroll {{ overflow-x:auto; background:#fff; border:1px solid #ddd;
             border-radius:6px; padding:6px; }}
  #tip {{ position:fixed; pointer-events:none; opacity:0; transition:opacity .08s;
          background:#111; color:#fff; font-size:12px; padding:6px 9px;
          border-radius:4px; max-width:420px; z-index:10; }}
  code {{ font-size:12px; }}
  /* :target, not a click handler -- sandboxed viewers run no JavaScript. */
  .figpanel {{ display:none; }}
  .figpanel:target {{ display:block; background:#fff; border:1px solid #ddd;
                      border-radius:6px; margin:8px 0; padding:10px; }}
  .fighead {{ font-size:12px; color:#666; margin-bottom:8px; }}
  .figclose {{ float:right; }}
  .figpanel object {{ border:1px solid #eee; border-radius:4px; min-height:340px; }}
</style></head>
<body><main id="top">
  <h1>pks cohort report</h1>
  <p class="muted">{len(records)} samples &middot; {esc(results)}</p>
  <div class="tiles">{tiles}</div>
  {no_assembly}
  <div class="scroll">{svg_matrix(records)}</div>
  {figure_panels(records)}
  <h2 style="font-size:15px;margin:26px 0 6px">Where these numbers come from</h2>
  <ul class="muted">
    <li><code>cohort/gene_counts/pks.gene.counts.align.txt</code> &mdash; the gene matrix</li>
    <li><code>by_sample/&lt;sample&gt;/read_evidence.tsv</code> &mdash; tier, reads, genes, breadth</li>
    <li><code>by_sample/&lt;sample&gt;/contigs/final_evidence/final_pks_evidence.tsv</code> &mdash; contig columns</li>
    <li><code>RUN_REPORT.txt</code> &mdash; the code, parameters and references behind them</li>
    <li><code>by_sample/&lt;sample&gt;/figures/</code> &mdash; what the Figures links open;
        they are relative, so this page has to stay in <code>cohort/</code></li>
  </ul>
  <p class="muted">This page reads and arranges; it derives nothing.</p>
</main>
<div id="tip"></div>
<script>
  // Hover anywhere with a data-tip: a gene cell, a whole row, the contig group.
  var tip = document.getElementById('tip');
  document.addEventListener('mouseover', function (event) {{
    var node = event.target.closest('[data-tip]');
    if (!node) {{ tip.style.opacity = 0; return; }}
    tip.textContent = node.getAttribute('data-tip');
    tip.style.opacity = 1;
  }});
  document.addEventListener('mousemove', function (event) {{
    var pad = 14;
    var x = Math.min(event.clientX + pad, window.innerWidth - tip.offsetWidth - 8);
    var y = Math.min(event.clientY + pad, window.innerHeight - tip.offsetHeight - 8);
    tip.style.left = x + 'px';
    tip.style.top = y + 'px';
  }});
  document.addEventListener('mouseleave', function () {{ tip.style.opacity = 0; }});
</script>
</body></html>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--table", default=None, type=Path)
    args = parser.parse_args()

    records = collect(args.results)
    if not records:
        raise SystemExit(f"[ERROR] no samples found under {args.results}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(page(records, args.results))

    if args.table:
        fields = ["sample", *GENES, "pks_reads", "genes_detected", "island_breadth_1x",
                  "island_breadth_2x", "island_breadth_3x", "read_evidence",
                  "megahit_supporting_contigs", "metaspades_supporting_contigs",
                  "megahit_reference_coverage", "metaspades_reference_coverage",
                  "assembler_agreement", "final_structural_evidence"]
        with args.table.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader()
            for record in records:
                writer.writerow({
                    "sample": record["sample"], **record["genes"],
                    "pks_reads": record["pks_reads"],
                    "genes_detected": record["genes_detected"],
                    "island_breadth_1x": f"{record['breadth_1x']:.6f}",
                    "island_breadth_2x": f"{record['breadth_2x']:.6f}",
                    "island_breadth_3x": f"{record['breadth_3x']:.6f}",
                    "read_evidence": record["tier"],
                    "megahit_supporting_contigs": record["megahit_contigs"],
                    "metaspades_supporting_contigs": record["metaspades_contigs"],
                    "megahit_reference_coverage": f"{record['megahit_cov']:.6f}",
                    "metaspades_reference_coverage": f"{record['metaspades_cov']:.6f}",
                    "assembler_agreement": record["agreement"],
                    "final_structural_evidence": record["structural"],
                })

    print(f"[cohort report] {len(records)} samples -> {args.output}")


if __name__ == "__main__":
    main()
