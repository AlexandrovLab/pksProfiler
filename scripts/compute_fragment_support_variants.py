#!/usr/bin/env python3
"""Four countings of "how much evidence supports a pks-positive call", and
whether the frozen tier cascade is stable across all of them.

Why this exists. The PI's concern: the pipeline's read-evidence tier
(multi_gene / broad_island / extensive_island) is currently driven by a single
number -- raw reads assigned to a clb gene -- and that number can move for
reasons that have nothing to do with how much of the island is really there.
Two mates of one physical DNA fragment are each their own row in this BAM (see
below), so a fragment can contribute two reads instead of one; PCR/optical
duplicates from library prep can contribute an arbitrary number of copies of
the same original fragment. Before trusting a tier call for a sample already
flagged pks-positive, this recomputes the same evidence four ways and checks
the tier does not move:

  1. raw           what pksProfiler_align.nf already publishes -- unread,
                   unchanged, not a second implementation that could drift
                   from the pipeline's own number.
  2. dup_marked    the same reads, with samtools-markdup-flagged PCR/optical
                   duplicates (SAM flag 0x400) excluded before gene
                   assignment and before breadth.
  3. fragment      reads grouped into distinct DNA fragments (read ID with a
                   trailing /1 or /2 stripped) before counting -- two mates of
                   one fragment landing in the same gene collapse to one unit
                   of evidence, not two.
  4. dup_collapsed fragment grouping AND duplicate exclusion together.

Why fragments are not "pairs", and why featureCounts -p / --countReadPairs is
not the answer. pksProfiler_align.nf aligns with `bowtie2 ... -U`, single-end
mode: filterReads.nf has already merged R1/R2 into one interleaved FASTQ
stream, tagging every read ID that does not already end in /1 or /2 with a
mate suffix (see Modules/filter_reads.nf's awk block). Every mate is aligned
and reported as its own independent single-end alignment record; the BAM
carries no PAIRED or PROPER_PAIR flag information at all. "Unique fragment"
here can only mean "distinct read ID with a trailing /1 or /2 stripped" -- a
string operation on the read name, not a SAM pairing flag, and not something
featureCounts' pairing mode can compute, because there is no real pairing
information left for it to use.

samtools markdup on single-end-style input. samtools markdup ordinarily sits
downstream of collate + fixmate, which use paired-end mate information (ms/MC
tags) to match up a read with its mate before deciding duplication. That
machinery has nothing to attach to here: FLAG 0x1 (paired) is never set on
these alignments, so fixmate has no mate to reconcile and collate has nothing
to regroup by name for. Tested directly against a small coordinate-sorted BAM
with no fixmate/collate step at all (see tests/test_fragment_support_variants.py):
samtools markdup marks single (unpaired) duplicates correctly, purely from
each record's own reference position, orientation and CIGAR -- exactly the
"same start, same strand, same length" rule that matters for two independent
reads or fragment mates that happen to have landed on an identical footprint.
No fixmate step is run, and none is needed: fixmate exists to make a read's
mate's coordinates visible to markdup, and there is no real mate here to make
visible. Running fixmate on unpaired-flagged records would not add real
pairing information; it would just be a no-op step that risks masking the
fact that this data has none.

One rule, reused, not reimplemented. Which clb gene a read belongs to is
scripts/assign_reads_to_genes.py's assign() -- the same canonical
largest-overlap-with-alphabetical-tiebreak rule the taxonomy lane already
uses, imported here directly. Tier classification is
scripts/classify_tumor_pks_evidence.classify() -- the same frozen threshold
cascade the pipeline already uses, imported here directly. Breadth uses the
same samtools depth invocation and canonical island coordinates
(params.pks_contig/pks_start_1based/pks_end_1based) as pksProfiler_align.nf
and pks_targeted.nf.

A fragment whose two mates disagree on gene assignment. assign() decides gene
membership per read (a mate IS a read, in this single-end BAM), not per
fragment -- there is no "fragment vs gene" overlap to measure, only "read vs
gene". So when a fragment's two mates land in different genes, each mate's own
assignment stands, and the fragment counts as one unit of *fragment-level*
evidence for EACH gene at least one of its reads was assigned to. A fragment
split across two genes is not force-assigned to a single "winning" gene (there
is no principled way to pick one over the other -- unlike overlap arbitration
within one read, which assign() already resolves) and is not dropped from
either gene's evidence (that would throw away real support for both). The one
consequence to note when reading the report: per-gene fragment counts can sum
to MORE than the total distinct-fragment count (variant 1's per-gene reads
never do this, because assign() gives each read at most one gene). The total
fragment figure fed to the tier cascade is always the number of DISTINCT
fragment IDs assigned to any gene at all, counted once each, never the sum of
the per-gene fragment counts.

This is a targeted, on-demand check for samples ALREADY flagged pks-positive,
not a new pipeline stage: main.nf does not run this for every sample, only a
reviewer decides to run it for one. See the CLI example at the bottom of this
docstring.

    python3 compute_fragment_support_variants.py \\
      --sample SAMPLE01 \\
      --bam    SAMPLE01.sorted.bam \\
      --counts SAMPLE01.counts.txt \\
      --qc     SAMPLE01.alignment.qc.tsv \\
      --coverage SAMPLE01.coverage.txt \\
      --annotation clb_genes.gff \\
      --output SAMPLE01.fragment_support_variants.tsv

--bam/--counts/--qc/--coverage are exactly pksProfiler_align.nf's own
published outputs for that sample (${sample_dir}/${sampleID}/alignment/ for
the BAM+qc.tsv, ${sample_dir}/${sampleID}/ for counts.txt+coverage.txt) --
nothing here is regenerated for variant 1. --annotation is the same
--pks_genome_annotation GFF the pipeline already uses. samtools and bedtools
default to whatever is on PATH; pass --samtools/--bedtools to pin the exact
binaries from conda_envs/pks_align_env.yml (samtools 1.21, bedtools 2.31.1)
if PATH does not already resolve to them.
"""
import argparse
import importlib.util
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


classify_tumor = _load("classify_tumor_pks_evidence", "classify_tumor_pks_evidence.py")
assign_reads = _load("assign_reads_to_genes", "assign_reads_to_genes.py")

VARIANTS = ("raw", "dup_marked", "fragment", "dup_collapsed")

FROZEN_THRESHOLDS = {
    "multi_gene":       (5, 3, 0.01),
    "broad_island":     (30, 8, 0.075),
    "extensive_island": (100, 10, 0.15),
}


def fragment_id(read_id):
    """Strip a trailing /1 or /2 mate suffix. filter_reads.nf guarantees every
    read ID has one (added if the raw FASTQ header lacked it), so this is the
    one and only definition of "which fragment does this read belong to"."""
    if read_id.endswith("/1") or read_id.endswith("/2"):
        return read_id[:-2]
    return read_id


def run(cmd, **kwargs):
    kwargs.setdefault("check", True)
    return subprocess.run(cmd, **kwargs)


def build_clb_bed(annotation, out_bed, min_overlap_script_dir=ROOT):
    run([sys.executable, str(min_overlap_script_dir / "clb_gene_bed.py"),
         "--annotation", str(annotation), "--output", str(out_bed)],
        capture_output=True, text=True)


def bam_to_gene_assignment(bam, clb_bed, bedtools, min_overlap):
    """bedtools bamtobed | bedtools intersect -wo | assign_reads_to_genes.assign()
    -- exactly extractPksIslandReads's own pipeline in Modules/pks_taxa.nf,
    reused here rather than reimplemented."""
    bamtobed = run([bedtools, "bamtobed", "-i", str(bam)], capture_output=True, text=True)
    intersect = run([bedtools, "intersect", "-a", "-", "-b", str(clb_bed), "-wo"],
                     input=bamtobed.stdout, capture_output=True, text=True)
    rows = list(assign_reads.parse(intersect.stdout.splitlines(keepends=True),
                                    read_name_field=3, gene_field=9, overlap_field=-1))
    return assign_reads.assign(rows, min_overlap=min_overlap)


def duplicate_read_ids(bam, samtools):
    """Read IDs samtools markdup flagged 0x400 (PCR/optical duplicate)."""
    result = run([samtools, "view", "-f", "0x400", str(bam)], capture_output=True, text=True)
    return {line.split("\t", 1)[0] for line in result.stdout.splitlines() if line}


def markdup(bam, out_bam, samtools):
    """samtools markdup directly on the coordinate-sorted BAM.

    No collate/fixmate step: this BAM has no real pairing (see module
    docstring) so there is no mate information for fixmate to compute, and
    every record here is FLAG-unpaired, which is exactly what markdup needs to
    treat each alignment as an independent single read rather than requiring
    the ms/MC tags fixmate would otherwise add for paired data."""
    run([samtools, "markdup", "-s", str(bam), str(out_bam)], capture_output=True, text=True)
    run([samtools, "index", str(out_bam)], capture_output=True, text=True)


def depth_breadth(bam, contig, start, end, island_length, samtools, extra_flags=()):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".depth.tsv", delete=False) as handle:
        depth_path = handle.name
    cmd = [samtools, "depth", "-aa", "-s", "-Q", "40",
           "-r", f"{contig}:{start}-{end}", *extra_flags, str(bam)]
    result = run(cmd, capture_output=True, text=True)
    Path(depth_path).write_text(result.stdout)
    try:
        return classify_tumor.breadths(depth_path, island_length)[0]
    finally:
        Path(depth_path).unlink(missing_ok=True)


def gene_counts_from_assignment(assignment, exclude_ids=frozenset(), by_fragment=False):
    """{gene: count} from an assign() result, optionally excluding read IDs
    (duplicate exclusion) and optionally counting distinct fragments touching
    a gene rather than reads assigned to it (fragment grouping). Also returns
    the total count fed to the tier cascade: reads (len of the filtered
    assignment) for read-level variants, or the number of DISTINCT fragment
    IDs assigned to ANY gene for fragment-level variants -- never the sum of
    the per-gene fragment counts, which can double-count a fragment split
    across two genes (see module docstring)."""
    filtered = {read_id: gene_overlap for read_id, gene_overlap in assignment.items()
                if read_id not in exclude_ids}

    per_gene = defaultdict(set) if by_fragment else defaultdict(int)
    for read_id, (gene, _overlap) in filtered.items():
        if by_fragment:
            per_gene[gene].add(fragment_id(read_id))
        else:
            per_gene[gene] += 1

    if by_fragment:
        gene_counts = {gene: len(frags) for gene, frags in per_gene.items()}
        total = len({fragment_id(read_id) for read_id in filtered})
    else:
        gene_counts = dict(per_gene)
        total = len(filtered)

    return gene_counts, total


def genes_detected(gene_counts):
    return sum(1 for count in gene_counts.values() if count >= 1)


def compute_variants(args):
    thresholds = {
        "multi_gene":       (args.multi_gene_reads, args.multi_gene_genes, args.multi_gene_breadth),
        "broad_island":     (args.broad_island_reads, args.broad_island_genes, args.broad_island_breadth),
        "extensive_island": (args.extensive_island_reads, args.extensive_island_genes, args.extensive_island_breadth),
    }
    for lower, higher in (("multi_gene", "broad_island"), ("broad_island", "extensive_island")):
        if any(h < l for l, h in zip(thresholds[lower], thresholds[higher])):
            raise SystemExit(f"tier thresholds must be non-decreasing: {lower} exceeds {higher}")

    island_length = args.pks_end - args.pks_start + 1

    # --- Variant 1: raw. Read the pipeline's own published files, unchanged. ---
    qc_metrics = classify_tumor.qc_metrics(args.qc)
    raw_reads = int(float(qc_metrics.get("reads_clb_genes_align", 0)))
    raw_genes = classify_tumor.detected_genes(args.counts)
    raw_breadth = classify_tumor.breadths(args.coverage, island_length)[0]
    raw_tier = classify_tumor.classify(raw_reads, raw_genes, raw_breadth, thresholds)

    results = {
        "raw": {"count": raw_reads, "genes_detected": raw_genes,
                "breadth_1x": raw_breadth, "tier": raw_tier},
    }

    # --- Shared read-to-gene assignment (variants 2/3/4): the one canonical rule. ---
    with tempfile.TemporaryDirectory() as tmp:
        clb_bed = Path(tmp) / "clb_genes.bed"
        build_clb_bed(args.annotation, clb_bed)
        assignment = bam_to_gene_assignment(args.bam, clb_bed, args.bedtools, args.min_overlap)

        markdup_bam = Path(tmp) / "markdup.bam"
        markdup(args.bam, markdup_bam, args.samtools)
        dup_ids = duplicate_read_ids(markdup_bam, args.samtools)

        dup_breadth = depth_breadth(markdup_bam, args.pks_contig, args.pks_start, args.pks_end,
                                     island_length, args.samtools, extra_flags=("-G", "1024"))

    # --- Variant 2: dup_marked. Same reads, duplicates (0x400) excluded. ---
    dm_genes, dm_reads = gene_counts_from_assignment(assignment, exclude_ids=dup_ids, by_fragment=False)
    dm_genes_detected = genes_detected(dm_genes)
    dm_tier = classify_tumor.classify(dm_reads, dm_genes_detected, dup_breadth, thresholds)
    results["dup_marked"] = {"count": dm_reads, "genes_detected": dm_genes_detected,
                              "breadth_1x": dup_breadth, "tier": dm_tier}

    # --- Variant 3: fragment. Same assignment, no dedup, grouped by fragment.
    # Breadth is unchanged from variant 1: fragment grouping is a counting
    # choice over already-aligned reads, not a change to which bases are
    # covered. ---
    fr_genes, fr_total = gene_counts_from_assignment(assignment, exclude_ids=frozenset(), by_fragment=True)
    fr_genes_detected = genes_detected(fr_genes)
    fr_tier = classify_tumor.classify(fr_total, fr_genes_detected, raw_breadth, thresholds)
    results["fragment"] = {"count": fr_total, "genes_detected": fr_genes_detected,
                            "breadth_1x": raw_breadth, "tier": fr_tier}

    # --- Variant 4: dup_collapsed. Fragment grouping AND duplicate exclusion.
    # Breadth is the same -G 1024 depth as variant 2: fragment grouping does
    # not change which bases are covered, only dedup does. ---
    dc_genes, dc_total = gene_counts_from_assignment(assignment, exclude_ids=dup_ids, by_fragment=True)
    dc_genes_detected = genes_detected(dc_genes)
    dc_tier = classify_tumor.classify(dc_total, dc_genes_detected, dup_breadth, thresholds)
    results["dup_collapsed"] = {"count": dc_total, "genes_detected": dc_genes_detected,
                                 "breadth_1x": dup_breadth, "tier": dc_tier}

    return results


def verdict(results):
    baseline = results["raw"]["tier"]
    diverged = [name for name in VARIANTS[1:] if results[name]["tier"] != baseline]
    if not diverged:
        return "PASS", []
    notes = [f"{name}: {baseline} -> {results[name]['tier']}" for name in diverged]
    return "FLAG", notes


def write_report(path, sample, results):
    verdict_label, notes = verdict(results)
    with open(path, "w") as handle:
        handle.write("sample\tvariant\treads_or_fragments\tgenes_detected\tisland_breadth_1x\ttier\n")
        for name in VARIANTS:
            row = results[name]
            handle.write(f"{sample}\t{name}\t{row['count']}\t{row['genes_detected']}\t"
                         f"{row['breadth_1x']:.6f}\t{row['tier']}\n")
        handle.write("\n")
        handle.write(f"# verdict\t{verdict_label}\n")
        if notes:
            for note in notes:
                handle.write(f"# diverged\t{note}\n")
        else:
            handle.write("# diverged\tnone -- tier is stable across all four variants\n")
    return verdict_label, notes


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--bam", required=True, help="post-alignment sorted BAM, pre-dedup (pksProfiler_align.nf's own ${sampleID}.sorted.bam)")
    p.add_argument("--counts", required=True, help="published ${sampleID}.counts.txt (featureCounts --largestOverlap)")
    p.add_argument("--qc", required=True, help="published ${sampleID}.alignment.qc.tsv")
    p.add_argument("--coverage", required=True, help="published ${sampleID}.coverage.txt (samtools depth -aa -s -Q40 over the island)")
    p.add_argument("--annotation", required=True, help="--pks_genome_annotation GFF (clbA-clbS source)")
    p.add_argument("--output", required=True)
    p.add_argument("--min-overlap", type=int, default=1, help="matches params.min_gene_overlap_bp")
    p.add_argument("--pks-contig", default="NC_017628.1")
    p.add_argument("--pks-start", type=int, default=2193827)
    p.add_argument("--pks-end", type=int, default=2244594)
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--bedtools", default="bedtools")
    p.add_argument("--multi-gene-reads", type=int, default=5)
    p.add_argument("--multi-gene-genes", type=int, default=3)
    p.add_argument("--multi-gene-breadth", type=float, default=.01)
    p.add_argument("--broad-island-reads", type=int, default=30)
    p.add_argument("--broad-island-genes", type=int, default=8)
    p.add_argument("--broad-island-breadth", type=float, default=.075)
    p.add_argument("--extensive-island-reads", type=int, default=100)
    p.add_argument("--extensive-island-genes", type=int, default=10)
    p.add_argument("--extensive-island-breadth", type=float, default=.15)
    args = p.parse_args()

    results = compute_variants(args)
    verdict_label, notes = write_report(args.output, args.sample, results)
    print(f"{args.sample}: {verdict_label}" + ("; " + "; ".join(notes) if notes else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
