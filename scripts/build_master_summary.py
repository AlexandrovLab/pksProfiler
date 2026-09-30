#!/usr/bin/env python3
"""Join every per-sample result into one table.

pksProfiler writes each stage's results separately, which is right for the pipeline and
wrong for reading. This walks a finished (or partly finished) results directory and emits
one row per sample carrying whatever each stage produced, with `NA` where a stage did not
run. Stages that were never enabled simply contribute no columns.

Run it after a pipeline run:

    python3 scripts/build_master_summary.py --results results \\
        --output results/cohort/pks.master_summary.tsv

It only reads published output, so it is safe to re-run at any time and does not require
re-running the pipeline.

Inside the pipeline it does not read the published tree at all (Ludmil, 2026-09-25 report
P0-5b): Nextflow publishes asynchronously, so a published copy could still be missing when
the masterSummary task ran, and a file left at the same --outdir by an earlier run could be
read in its place. The task stages every file it reads as a declared input and passes

    --inputs-manifest pks.master_summary.inputs.tsv    sample<TAB>kind<TAB>path

instead of --results. The two are mutually exclusive; --results stays for hand runs and
simply builds the same manifest by looking where each kind is published (PUBLISHED below).
"""
import argparse
import csv
import glob
import os
import sys
from pathlib import Path
from collections import defaultdict

from mag_utils import POSITIVE_LOCUS_TIERS
from build_qc_summary import WES_QC_COLUMNS

# Column groups, in the order they appear. Each group is skipped entirely when the stage
# that produces it left no files behind.
QC_COLUMNS = ["input_alignment_records", "total_primary_reads", "extracted_unmapped_reads",
              "extracted_low_score_mapped_reads", "extracted_reads",
              "filter_input_reads", "reads_after_fastp", "reads_after_hg38", "reads_after_t2t_phix",
              "reads_mapped_ihe3034_mapq40", "num_clb_genes_align", "reads_clb_genes_align",
              "num_clb_genes_hmm", "reads_clb_genes_hmm"]
# pks_reads_per_million / host_depleted_reads: metagenome only (26 Sep memo section 5), NA for wgs.
EVIDENCE_COLUMNS = ["read_evidence", "pks_reads_per_million", "pks_reads", "host_depleted_reads",
                    "clb_genes_detected", "island_breadth_1x"]
# Ludmil, 26 Sep mode memo sections 3, 4 and 6: WES and RNA evidence keep their own
# column names (classify_tumor_pks_evidence.py SCHEMAS) and are never folded into the WGS
# tier columns above; `mode` and `evidence_label` say which kind of evidence a row is.
MODE_COLUMNS = ["mode", "evidence_label"]
# GAMEPLAN rc1: how the sample ended (build_qc_summary.py terminal_state), right after the
# mode so the row's outcome is read before any of its numbers.
TERMINAL_COLUMNS = ["terminal_state"]
WES_EVIDENCE_COLUMNS = ["clb_wes_reads", "wes_clb_genes_detected", "wes_island_breadth_1x"]
RNA_EVIDENCE_COLUMNS = ["clb_rna_reads", "rna_clb_genes_detected", "rna_expression_breadth_1x",
                        "rna_library_type",
                        # memo section 4 RNA QC: clb genes at >= N reads, residual human
                        "clb_transcripts_detected", "rna_residual_human_fraction"]
# The mode each evidence kind implies, for a file written before `mode` was a column.
KIND_MODE = {"wes_evidence": "wes", "rna_evidence": "rna"}
CONTIG_COLUMNS = ["final_structural_evidence", "assembler_agreement",
                  "megahit_reference_covered_bp", "megahit_reference_coverage",
                  "megahit_supporting_contigs", "metaspades_reference_covered_bp",
                  "metaspades_reference_coverage", "metaspades_supporting_contigs",
                  "recruited_fragment_ids", "paired_fragments",
                  # Ludmil, 25 Sep A2: recruitment classes; the diagnostic reads
                  # island-supporting (island_internal + boundary_spanning) only.
                  "island_supporting_fragments", "boundary_spanning_fragments",
                  "flank_only_fragments", "assembly_input_diagnostic"]
MAG_COLUMNS = ["mag_bins_total", "mag_bins_pks_positive", "pks_mag_taxonomy",
               "pks_mag_completeness", "pks_mag_clb_genes",
               # U1: the per-sample pool of contigs MetaBAT2 never placed in any bin,
               # kept as its own pair of columns rather than folded into mag_bins_*
               # above -- a positive here means "this sample carries clb sequence
               # outside any recovered genome," not "one more recovered genome is
               # pks-positive," and conflating the two would misreport which claim
               # the evidence actually supports.
               "mag_unbinned_locus_tier", "mag_unbinned_pks_positive"]
COMMUNITY_COLUMNS = ["community_prophages_total", "community_pks_producers",
                     "community_neighbours_assessed", "community_neighbours_with_prophage",
                     "island_nearby_integrase", "island_nearby_trna", "island_in_prophage"]
TYPING_COLUMNS = ["typing_status", "ST", "clonal_complex", "phylogroup"]
# Ludmil, 25 Sep report section 4. V8a asks whether the sequence belongs to the pks/clb
# locus (challenged against pks-NEGATIVE decoys); V8b asks which pks-POSITIVE organism
# carries it. The report says they "should not be conflated", so they are two column
# groups with disjoint names and disjoint vocabularies, and a run may have either, both
# or neither. `locus_specificity_arms` records whether Arm A ran, because SURVIVES from
# one arm and SURVIVES from two are not the same claim.
V8A_COLUMNS = ["locus_specificity", "locus_specificity_arms", "locus_specificity_retention",
               "locus_specificity_recomputed_tier", "locus_specificity_flagged_absorbers"]
V8B_COLUMNS = ["carrier_attribution", "carrier_top", "carrier_top_organism",
               "carrier_top_fraction", "carrier_top_genus", "carrier_runner_up"]

TRUE = {"true", "yes", "1", "t"}

# P0-5b: every kind this table reads, and where it is published under by_sample/<sample>/
# (per-sample) or the results root (cohort-level). main.nf tags each staged file with one
# of these kinds; --results globs these paths for a hand run.
PUBLISHED = {
    "read_evidence": "read_evidence.tsv",
    "wes_evidence": "wes_evidence.tsv",
    "rna_evidence": "rna_evidence.tsv",
    "final_evidence": "contigs/final_evidence/final_pks_evidence.tsv",
    "mag_summary": "genomes/pks_mag_summary.tsv",
    "prophage_inventory": "community/community_prophage_inventory.tsv",
    "community_interactions": "community/pks_community_interactions.tsv",
    "island_mobility": "community/pks_island_mobility.tsv",
    # No current process publishes this (the tumour contig-context lane was removed);
    # kept so a hand run over an older tree still reads it.
    "contig_mobility": "community/contig_pks_mobility.tsv",
    "strain_types": "strain/strain_types.tsv",
    # Ludmil V8a and V8b: two different questions, two kinds, two column groups. Never
    # merged into one "specificity" column -- see V8A_COLUMNS/V8B_COLUMNS below.
    "locus_specificity": "locus_specificity/locus_specificity_summary.tsv",
    "carrier_attribution": "carrier_attribution/carrier_attribution.tsv",
}
COHORT_PUBLISHED = {"qc_summary": "cohort/qc/pks.qc.summary.tsv"}


def inputs_from_manifest(path):
    """{kind: [(sample, path), ...]} from a staged-input manifest; paths relative to it.

    Kinds this table does not read (the cohort report's, staged by the same channels)
    are skipped. Cohort-level kinds carry an empty sample.
    """
    base = Path(path).parent
    inputs = defaultdict(list)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            kind, rel = row.get("kind") or "", row.get("path") or ""
            if rel and (kind in PUBLISHED or kind in COHORT_PUBLISHED):
                inputs[kind].append((row.get("sample") or "", str(base / rel)))
    return inputs


def inputs_from_results(results):
    """The same shape, by looking where each kind is published. Hand runs only."""
    inputs = defaultdict(list)
    for kind, rel in COHORT_PUBLISHED.items():
        path = os.path.join(results, rel)
        if os.path.exists(path):
            inputs[kind].append(("", path))
    for kind, rel in PUBLISHED.items():
        for path in sorted(glob.glob(f"{results}/by_sample/*/{rel}")):
            inputs[kind].append((sample_of(path), path))
    return inputs


def rows(path, delim="\t"):
    try:
        with open(path) as fh:
            yield from csv.DictReader(fh, delimiter=delim)
    except (OSError, csv.Error) as exc:
        print(f"  warning: could not read {path}: {exc}", file=sys.stderr)


def first(path):
    for row in rows(path):
        return row
    return {}


def sample_of(path, row=None):
    """The sample a published file belongs to.

    Files under by_sample/<sample>/ no longer carry a <sample>. filename prefix, because
    the directory already says which sample they are. Every published table also has a
    `sample` column; prefer it when present, since it survives a file being copied out.
    """
    if row and row.get("sample"):
        return row["sample"]
    parts = os.path.normpath(path).split(os.sep)
    if "by_sample" in parts:
        i = parts.index("by_sample")
        if i + 1 < len(parts):
            return parts[i + 1]
    return None


def read_qc(inputs, out):
    if not inputs["qc_summary"]:
        return False
    for _, path in inputs["qc_summary"]:
        for row in rows(path):
            sample = row.get("Sample")
            if sample:
                for c in QC_COLUMNS + TERMINAL_COLUMNS:
                    out[sample][c] = row.get(c, "NA")
    return True


def read_wes_qc(inputs, out):
    """Memo section 3: the WES QC columns, present only when the QC table carries them."""
    found = False
    for _, path in inputs["qc_summary"]:
        for row in rows(path):
            sample = row.get("Sample")
            if sample and WES_QC_COLUMNS[0] in row:
                found = True
                for c in WES_QC_COLUMNS:
                    out[sample][c] = row.get(c, "NA")
    return found


def has_terminal_state(inputs):
    return any("terminal_state" in first(path) for _, path in inputs["qc_summary"])


def evidence_label(mode, read_evidence=None):
    """What a row's evidence may be read as. Mirrored in build_cohort_report.py."""
    if mode == "wgs":
        return ("WGS evidence (not tiered: long-read)" if read_evidence == "untiered"
                else "WGS genomic evidence tier")
    return {"metagenome": "Metagenome read evidence (not tiered)",
            "wes": "WES off-target DNA evidence (not tiered)",
            "rna": "RNA expression evidence (not tiered)"}.get(mode, "NA")


def read_mode(inputs, out):
    """`mode` and `evidence_label` from whichever evidence file a sample has."""
    found = False
    for kind in ("read_evidence", "wes_evidence", "rna_evidence"):
        for sample, path in inputs[kind]:
            row = first(path)
            sample = sample or row.get("sample")
            if row and sample:
                found = True
                mode = row.get("mode") or KIND_MODE.get(kind, "NA")
                out[sample]["mode"] = mode
                out[sample]["evidence_label"] = evidence_label(mode, row.get("read_evidence"))
    return found


def read_mode_evidence(kind, columns):
    def reader(inputs, out):
        found = False
        for sample, path in inputs[kind]:
            row = first(path)
            sample = sample or row.get("sample")
            if row and sample:
                found = True
                for c in columns:
                    out[sample][c] = row.get(c, "NA")
        return found
    return reader


def read_locus_specificity(inputs, out):
    """V8a only. `verdict` is the locus-detection call; nothing here names an organism."""
    found = False
    for sample, path in inputs["locus_specificity"]:
        row = first(path)
        sample = sample or (row.get("sample") if row else None)
        if not (row and sample):
            continue
        found = True
        out[sample]["locus_specificity"] = row.get("verdict", "NA")
        out[sample]["locus_specificity_arms"] = row.get("arms_used", "NA")
        out[sample]["locus_specificity_retention"] = row.get("retention", "NA")
        out[sample]["locus_specificity_recomputed_tier"] = row.get("recomputed_tier", "NA")
        out[sample]["locus_specificity_flagged_absorbers"] = row.get("flagged_absorbers", "")
    return found


def read_carrier_attribution(inputs, out):
    """V8b only. Strain and genus are both carried: a strain-level carrier printed
    without its genus-level context invites exactly the over-reading V8b guards against
    (Ludmil, same report: a divergence metric is not proof of carrier identity)."""
    found = False
    for sample, path in inputs["carrier_attribution"]:
        row = first(path)
        sample = sample or (row.get("sample") if row else None)
        if not (row and sample):
            continue
        found = True
        out[sample]["carrier_attribution"] = row.get("carrier_attribution", "NA")
        out[sample]["carrier_top"] = row.get("top_carrier", "NA")
        out[sample]["carrier_top_organism"] = row.get("top_carrier_organism", "NA")
        out[sample]["carrier_top_fraction"] = row.get("top_carrier_fraction", "NA")
        out[sample]["carrier_top_genus"] = row.get("top_genus", "NA")
        out[sample]["carrier_runner_up"] = row.get("runner_up_carrier", "NA")
    return found


read_evidence = read_mode_evidence("read_evidence", EVIDENCE_COLUMNS)
read_wes_evidence = read_mode_evidence("wes_evidence", WES_EVIDENCE_COLUMNS)
read_rna_evidence = read_mode_evidence("rna_evidence", RNA_EVIDENCE_COLUMNS)


def read_contigs(inputs, out):
    found = False
    for sample, path in inputs["final_evidence"]:
        row = first(path)
        sample = sample or row.get("sample")
        if row and sample:
            found = True
            for c in CONTIG_COLUMNS:
                out[sample][c] = row.get(c, "NA")
    return found


# The columns read_mags() needs from pks_mag_summary.tsv (build_mag_summary.py's
# FIELDNAMES). tests/test_master_summary.py builds that file with the real producer,
# so a rename on either side fails a test instead of silently zeroing the table.
MAG_SUMMARY_REQUIRED = ("unit_type", "taxonomy", "completeness",
                        "locus_tier", "locus_genes_detected", "locus_breadth")
LOCUS_TIER_ORDER = ("multi_gene", "broad_island", "extensive_island")


def int_or_zero(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def float_or_zero(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def locus_rank(b):
    """Strongest locus evidence first: tier, then genes detected, then breadth."""
    tier = b.get("locus_tier")
    return (LOCUS_TIER_ORDER.index(tier) if tier in LOCUS_TIER_ORDER else -1,
            int_or_zero(b.get("locus_genes_detected")),
            float_or_zero(b.get("locus_breadth")))


def read_mags(inputs, out):
    found = False
    for sample, path in inputs["mag_summary"]:
        all_rows = list(rows(path))
        if not all_rows:
            continue
        sample = sample or sample_of(path, all_rows[0])
        if not sample:
            continue
        found = True

        # U1: a row absent "unit_type" (older runs, or fixtures predating that column)
        # is a bin -- only the unbinned pseudo-unit row ever sets it to "unbinned".
        bins = [b for b in all_rows if b.get("unit_type", "bin") != "unbinned"]
        unbinned = next((b for b in all_rows if b.get("unit_type") == "unbinned"), None)

        # Ludmil, 2026-09-25 report P0-2: this read distinct_clb_genes/clb_genes and
        # classification/host_taxonomy -- columns build_mag_summary.py no longer writes
        # (clb_genes is now a comma-joined name list, so int() failed and returned 0),
        # and called a bin positive on HMM gene count. Every current-schema bin therefore
        # read as negative: mag_bins_pks_positive 0, best None, taxonomy NA. Positivity is
        # the alignment-confirmed canonical-locus tier, the same rule the MAG lane and
        # the cohort report use; a file without that column is refused, not guessed at.
        missing = [c for c in MAG_SUMMARY_REQUIRED if bins and c not in bins[0]]
        if missing:
            sys.exit(f"ERROR: {path} lacks {', '.join(missing)}: it predates the "
                     "locus-tier MAG schema (scripts/build_mag_summary.py). Rebuild it "
                     "rather than summarising it with the wrong positivity rule.")

        positive = [b for b in bins if b.get("locus_tier") in POSITIVE_LOCUS_TIERS]
        best = max(positive, key=locus_rank) if positive else None
        out[sample]["mag_bins_total"] = str(len(bins))
        out[sample]["mag_bins_pks_positive"] = str(len(positive))
        out[sample]["pks_mag_taxonomy"] = (best or {}).get("taxonomy") or "NA"
        out[sample]["pks_mag_completeness"] = (best or {}).get("completeness") or "NA"
        out[sample]["pks_mag_clb_genes"] = str(int_or_zero(best.get("locus_genes_detected"))) if best else "0"

        # U1: kept apart from the bin columns above by construction -- see MAG_COLUMNS.
        if unbinned is not None:
            tier = unbinned.get("locus_tier", "NA")
            out[sample]["mag_unbinned_locus_tier"] = tier
            out[sample]["mag_unbinned_pks_positive"] = "yes" if tier in POSITIVE_LOCUS_TIERS else "no"
    return found


def read_community(inputs, out):
    found = False
    for sample, path in inputs["prophage_inventory"]:
        sample = sample or sample_of(path)
        if not sample:
            continue
        found = True
        out[sample]["community_prophages_total"] = str(len(list(rows(path))))
    for sample, path in inputs["community_interactions"]:
        pairs = list(rows(path))
        sample = sample or sample_of(path, pairs[0] if pairs else None)
        if not sample:
            continue
        found = True
        out[sample]["community_pks_producers"] = str(len({p.get("pks_producer_bin") for p in pairs if p.get("pks_producer_bin")}))
        out[sample]["community_neighbours_assessed"] = str(len(pairs))
        def has_prophage(p):
            try:
                return int(p.get("recipient_prophage_count") or 0) > 0
            except ValueError:
                return False
        out[sample]["community_neighbours_with_prophage"] = str(sum(1 for p in pairs if has_prophage(p)))
    for sample, path in inputs["island_mobility"]:
        mob = list(rows(path))
        sample = sample or sample_of(path, mob[0] if mob else None)
        if not sample:
            continue
        found = True
        anyof = lambda col: "yes" if any((m.get(col) or "").strip().lower() in TRUE for m in mob) else "no"
        out[sample]["island_nearby_integrase"] = anyof("nearby_integrase")
        out[sample]["island_nearby_trna"] = anyof("nearby_trna")
    # tumour contigs report the same three flags per clb hit
    for sample, path in inputs["contig_mobility"]:
        hits = list(rows(path))
        if not hits:
            continue
        sample = sample or sample_of(path, hits[0])
        if not sample:
            continue
        found = True
        anyof = lambda col: "yes" if any((h.get(col) or "").strip().lower() in TRUE for h in hits) else "no"
        out[sample]["island_nearby_integrase"] = anyof("has_integrase")
        out[sample]["island_nearby_trna"] = anyof("nearby_trna")
        out[sample]["island_in_prophage"] = anyof("in_prophage")
    return found


def read_typing(inputs, out):
    found = False
    for sample, path in inputs["strain_types"]:
        units = list(rows(path))
        if not units:
            continue
        sample = sample or sample_of(path, units[0])
        if not sample:
            continue
        found = True
        typed = [u for u in units if u.get("status") == "typed"]
        pick = typed[0] if typed else units[0]
        out[sample]["typing_status"] = pick.get("status", "NA")
        for c in ("ST", "clonal_complex", "phylogroup"):
            out[sample][c] = pick.get(c, "NA")
    return found


def main():
    ap = argparse.ArgumentParser()
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--inputs-manifest",
                        help="sample<TAB>kind<TAB>path, one row per staged input; what the "
                             "pipeline passes (Ludmil, 2026-09-25 report P0-5b)")
    source.add_argument("--results", help="a pksProfiler results directory (hand runs)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--sample-type", default=None,
                    help="the run's --sample_type; the mode of every sample that has no "
                         "evidence file of its own (Ludmil, 26 Sep mode memo section 6)")
    ap.add_argument("--expected-samples", default=None,
                    help="one sample ID per line; samples outside this list are "
                         "dropped even if their directory is still on disk")
    a = ap.parse_args()

    inputs = (inputs_from_manifest(a.inputs_manifest) if a.inputs_manifest
              else inputs_from_results(a.results))
    out = defaultdict(dict)
    stages = [
        ("read profiling", QC_COLUMNS, read_qc),
        ("WES QC", WES_QC_COLUMNS, read_wes_qc),
        ("evidence tier", EVIDENCE_COLUMNS, read_evidence),
        ("WES evidence", WES_EVIDENCE_COLUMNS, read_wes_evidence),
        ("RNA evidence", RNA_EVIDENCE_COLUMNS, read_rna_evidence),
        ("contig reassembly", CONTIG_COLUMNS, read_contigs),
        ("draft genomes", MAG_COLUMNS, read_mags),
        ("community context", COMMUNITY_COLUMNS, read_community),
        ("strain typing", TYPING_COLUMNS, read_typing),
        ("V8a locus specificity", V8A_COLUMNS, read_locus_specificity),
        ("V8b carrier attribution", V8B_COLUMNS, read_carrier_attribution),
    ]
    columns = ["sample"]
    # Mode first, so a table mixing modalities can never be read without it.
    if read_mode(inputs, out) or a.sample_type:
        columns += MODE_COLUMNS
    if has_terminal_state(inputs):
        columns += TERMINAL_COLUMNS
    for name, cols, fn in stages:
        present = fn(inputs, out)
        print(f"  {name:20} {'found' if present else 'not run — columns omitted'}", file=sys.stderr)
        if present:
            columns += cols

    # With --results this walks by_sample/ by path, so a directory left over from an
    # older run is indistinguishable from one this run produced (the same class of bug
    # as build_cohort_report.py's, finding 5). The pipeline's manifest cannot name one,
    # but the filter stays: it is what the task was given as the run's sample list.
    if a.expected_samples:
        with open(a.expected_samples) as fh:
            expected = {line.strip() for line in fh if line.strip()}
        out = {sample: cols for sample, cols in out.items() if sample in expected}

    if not out:
        sys.exit(f"No per-sample results found in {a.inputs_manifest or a.results}")
    if a.sample_type:
        for cols in out.values():
            cols.setdefault("mode", a.sample_type)
            cols.setdefault("evidence_label", evidence_label(cols["mode"], cols.get("read_evidence")))

    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    with open(a.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", restval="NA", extrasaction="ignore")
        w.writeheader()
        for sample in sorted(out):
            w.writerow({"sample": sample, **out[sample]})
    print(f"{len(out)} samples x {len(columns)} columns -> {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
