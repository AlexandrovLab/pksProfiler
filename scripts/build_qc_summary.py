#!/usr/bin/env python3
"""Combine per-stage pksProfiler QC fragments into one cohort table.

The table is built from whatever fragments exist. A sample that failed mid-run, or
never ran at all, gets a row recording that rather than taking the whole cohort
table down with it -- losing 96 samples' QC because one hit a time cap is the
failure this guards against. Corruption is still fatal: conflicting values for the
same metric, counts that increase down the depletion chain, and out-of-range gene
counts all mean a bug, not a missing sample.
"""

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


# F09: one column, one quantity, and nothing stands in for anything else.
#
# `input_reads` used to mean "primary records in the input alignment, or for FASTQ the
# reads entering fastp", and when the alignment count was missing -- which it always
# was, because nothing emitted it -- it silently became the post-extraction count. A
# column that is sometimes the library and sometimes what survived extraction cannot
# show how much of the library was lost, which is why F01 was invisible here.
#
# Each of these is now written by exactly one stage, under its own name, and is absent
# rather than substituted when that stage did not run.
# A5: written by ihe3034_divergence_qc.py. The read count is an integer like every other
# count; the rest are decimals, and NA when no read was measured.
DIVERGENCE_COUNT = "divergence_ihe3034_reads"
# metric -> upper bound (None: any non-negative value)
DIVERGENCE_DECIMALS = {
    "nm_ihe3034_median": None,
    "divergence_ihe3034_median": 1.0,
    "divergence_ihe3034_q1": 1.0,
    "divergence_ihe3034_q3": 1.0,
    "divergence_ihe3034_frac_gt_2pct": 1.0,
    "divergence_ihe3034_frac_gt_5pct": 1.0,
    "clb_read_length_median": None,
    "clb_read_length_q1": None,
    "clb_read_length_q3": None,
    "clb_read_length_min": None,
    "clb_read_length_max": None,
}
DIVERGENCE_COLUMNS = [DIVERGENCE_COUNT, *DIVERGENCE_DECIMALS]

OUTPUT_COLUMNS = [
    "Sample",
    "input_alignment_records",   # from the index: records, secondary/supplementary included
    # E2: the same index counts by input-reference category (summarize_idxstats.py).
    # Mapped records per category, then unmapped records with and without a position.
    "input_mapped_canonical",
    "input_mapped_mito",
    "input_mapped_alt",
    "input_mapped_unplaced",
    "input_mapped_decoy",
    "input_mapped_viral",
    "input_mapped_other",
    "input_unmapped_placed",
    "input_unmapped_unplaced",
    "total_primary_reads",       # exact primary count; only with --exact_input_counts
    # E1 (Ludmil, 25 Sep): extraction takes FLAG 4 records, plus -- when
    # extract_low_score_mode is on -- primary mapped records with AS below
    # --extract_min_as. extracted_reads is the sum, what extraction wrote, and what enters
    # fastp (filter_input_reads); it is the extraction count in the attrition chain.
    "extracted_unmapped_reads",  # FLAG 4 records extracted (the only class before E1)
    "extracted_low_score_mapped_reads",  # primary mapped records with AS < extract_min_as
    "extracted_reads",           # everything extraction wrote: the two above, summed
    "extract_low_score_no_as_records",  # mapped primary records with no AS tag, never added
    "extract_low_score_mode",    # on_auto_bwa_mem / off_auto_not_bwa_mem / on_requested / off_requested
    "filter_input_reads",        # what entered fastp
    "reads_after_fastp",
    "reads_after_hg38",
    "reads_after_t2t_phix",
    "reads_after_pangenome",
    # P1 (25 Sep): the metagenome prefilter lane. Kraken routing of the depleted reads,
    # DIAMOND rescue of the non-target remainder, and the merged candidate set.
    "prefilter_input_reads",
    "prefilter_primary_reads",
    "prefilter_non_target_reads",
    "diamond_input_reads",
    "diamond_rescue_fragments",
    "diamond_rescued_reads",
    "candidate_primary_reads",
    "candidate_rescued_reads",
    "candidate_total_reads",
    # reads_mapped_ihe3034_mapq40 is deliberately NOT a column here. It is reads mapping
    # anywhere on the 5.1 Mb IHE3034 genome, of which the island is 1%, so beside the
    # clb columns it invites the reading "lots of reads, no island" as though the two
    # were the same measurement. It is still written per sample in
    # by_sample/<sample>/alignment/*.alignment.qc.tsv -- exported, as F09 asked, and
    # still used below to bound the clb count -- just not shown in the cohort table.
    "num_clb_genes_align",
    "reads_clb_genes_align",
    # Ludmil, 25 Sep A5: divergence FROM IHE3034 of the reads behind reads_clb_genes_align
    # (ihe3034_divergence_qc.py). A QC distribution, not a call and not a specificity
    # measure: a read near 0 is close to this one reference, which says nothing about
    # what else it would also match. Specificity is the competitive-mapping lane.
    *DIVERGENCE_COLUMNS,
    "num_clb_genes_hmm",
    "reads_clb_genes_hmm",
    # F07: reads whose best hit tied across two clb models and were therefore counted
    # for neither. Reported rather than left as a gap in the counts.
    "hmm_ambiguous_reads",
    # F10: why a sample is absent from the species table. no_pks_reads,
    # below_rank_threshold, no_species_identified, species_identified -- or NA when
    # taxonomy did not run at all.
    "taxonomy_status",
    "clb_species_reported",
    "status",
    # GAMEPLAN rc1 (release mechanics): how the sample ended -- exactly one of
    # completed_positive / completed_negative / completed_untiered / failed_at_<stage> /
    # excluded_<reason>. See terminal_state() below.
    "terminal_state",
]

# Ludmil, 26 Sep mode memo sections 3 and 4: mode-specific QC columns. They are added to
# the table (just before `status`) only when the run expects the mode's QC stage, so a WGS
# or metagenome table is unchanged. Definitions: docs/running/wes_rna.md.
#
# WES. The input alignment's own duplicate marking and capture overlap (wes_input_qc.py),
# and duplicates among the pks-island reads by position collapse (pks_duplicate_rate.py).
# The two rates are different quantities and are never combined.
WES_QC_COLUMNS = [
    "wes_input_primary_mapped",           # primary mapped records in the input alignment
    "wes_input_duplicate_fraction",       # of those, flagged 0x400 by the submitter
    "wes_capture_on_target_fraction",     # of those, overlapping --wes_capture_bed
    "wes_capture_off_target_fraction",    # 1 - on-target; capture_bed_not_provided without a BED
    "wes_pks_island_reads",               # primary island reads in the profiling BAM (MAPQ >= 40)
    "wes_pks_island_duplicate_fraction",  # of those, same (contig, strand, unclipped 5') as another
]
# RNA. A capped sample of the final host-depleted reads realigned to the splice index
# under relaxed floors (rna_residual_human.py); the fraction that still maps.
RNA_QC_COLUMNS = [
    "rna_residual_checked_reads",
    "rna_residual_human_reads",
    "rna_residual_human_fraction",
]
CAPTURE_NOT_PROVIDED = "capture_bed_not_provided"

# L4 (Ludmil, 1 Oct): of the reads E1's low-score rule added at extraction
# (extracted_low_score_mapped_reads), how many are still there after each stage: counted by
# name in mapReads (fastp, hg38, T2T+phiX) and lowScoreRescueProfile (in the profiling BAM,
# then featureCounts-assigned to a clb gene). Only samples that rescued a read have them;
# the columns appear, after extract_low_score_mode, only when some sample does, so a run
# with the rule off (or FASTQ input) writes the same table as before.
LOW_SCORE_RESCUE_COLUMNS = [
    "low_score_rescued_after_fastp",
    "low_score_rescued_after_hg38",
    "low_score_rescued_after_t2t_phix",
    "low_score_rescued_profiled",
    "low_score_rescued_clb_assigned",
]
# V2 robustness (decision A, 1 Oct): the tier the same evidence would get if counted as
# DNA fragments, with duplicates removed, or both -- scripts/compute_fragment_support_variants.py
# run inside classifyPksReadEvidence. read_evidence itself is still the raw-read tier. Only
# tiered (short-read wgs) samples have them; the columns appear, before `status`, only then.
TIER_ROBUSTNESS_COLUMNS = [
    "pks_fragments",
    "tier_if_fragments",
    "pks_reads_dedup",
    "tier_if_dedup",
    "pks_fragments_dedup",
    "tier_if_fragments_dedup",
    "tier_robustness",
]
OPTIONAL_COLUMNS = LOW_SCORE_RESCUE_COLUMNS + TIER_ROBUSTNESS_COLUMNS
TIER_VALUES = {"negative", "localized_indeterminate", "multi_gene", "broad_island",
               "extensive_island"}

# Pipeline stages in execution order, and the QC metrics whose presence shows the stage
# finished for a sample. A failed sample is failed_at_<the first expected stage it is
# missing>: stages downstream of a failure never run, so the first gap is the failure.
# alignment_profiling and hmm_profiling run side by side; alignment is listed first.
STAGE_METRICS = {
    "extraction": ("extracted_unmapped_reads",),
    "wes_input_qc": ("wes_input_primary_mapped", "wes_capture_status"),
    "read_qc": ("filter_input_reads", "reads_after_fastp"),
    "host_depletion": ("reads_after_hg38", "reads_after_t2t_phix"),
    "rna_residual_human": ("rna_residual_checked_reads", "rna_residual_human_reads"),
    "prefilter": ("prefilter_input_reads",),
    "diamond_rescue": ("diamond_input_reads",),
    "candidate_merge": ("candidate_total_reads",),
    "alignment_profiling": ("reads_clb_genes_align", "num_clb_genes_align"),
    "hmm_profiling": ("reads_clb_genes_hmm", "num_clb_genes_hmm"),
    "wes_pks_duplicates": ("wes_pks_island_reads", "wes_pks_island_duplicates"),
    "evidence": ("evidence_outcome",),
    "taxonomy": ("taxonomy_status",),
}
STAGE_ORDER = tuple(STAGE_METRICS)
# evidence_outcome (classify_tumor_pks_evidence.py) -> completed state. Only tiered
# short-read wgs is ever positive or negative; wes, rna, metagenome and long-read wgs
# evidence is untiered.
COMPLETED_STATE = {"positive": "completed_positive", "negative": "completed_negative",
                   "untiered": "completed_untiered"}
EXCLUSION_REASON = re.compile(r"^[a-z0-9][a-z0-9_]*$")

# complete       every required metric present
# incomplete     some required metrics missing -- the sample ran partially
# no_qc_produced expected from the sample sheet, but emitted no QC at all
# Required means "every run produces this". The alignment-only metrics are not here:
# a FASTQ run has no input alignment, and demanding one would mark every FASTQ sample
# incomplete.
REQUIRED = [
    "filter_input_reads",
    "reads_after_fastp",
    "reads_after_hg38",
    "reads_after_t2t_phix",
]


# Metrics whose value is an outcome rather than a count, with the outcomes they may
# take. F10: the taxonomy lane reports why a sample has no species, and "no species"
# is not a number.
TEXT_METRICS = {
    "taxonomy_status": {"no_pks_reads", "below_rank_threshold",
                        "no_species_identified", "species_identified"},
    "evidence_outcome": set(COMPLETED_STATE),
    # E1: scripts/header_provenance.py LOW_SCORE_STATES.
    "extract_low_score_mode": {"on_auto_bwa_mem", "off_auto_not_bwa_mem",
                               "on_requested", "off_requested"},
    "wes_capture_status": {"capture_bed_provided", CAPTURE_NOT_PROVIDED},
    "tier_if_fragments": TIER_VALUES,
    "tier_if_dedup": TIER_VALUES,
    "tier_if_fragments_dedup": TIER_VALUES,
    "tier_robustness": {"stable", "unstable"},
}


def fraction(part, whole):
    return "NA" if part is None or not whole else f"{part / whole:.6f}"


def mode_qc_values(sample, values):
    """The WES and RNA QC columns for one sample, with their consistency checks."""
    def bounded(part_name, whole_name):
        part, whole = values.get(part_name), values.get(whole_name)
        if part is not None and whole is not None and part > whole:
            raise ValueError(f"Impossible QC counts for {sample}: {part_name} ({part}) "
                             f"exceeds {whole_name} ({whole})")
        return part, whole

    out = {}
    duplicates, mapped = bounded("wes_input_primary_duplicates", "wes_input_primary_mapped")
    on_target, _ = bounded("wes_capture_on_target", "wes_input_primary_mapped")
    out["wes_input_primary_mapped"] = mapped
    out["wes_input_duplicate_fraction"] = fraction(duplicates, mapped)
    status = values.get("wes_capture_status")
    if status == CAPTURE_NOT_PROVIDED:
        out["wes_capture_on_target_fraction"] = CAPTURE_NOT_PROVIDED
        out["wes_capture_off_target_fraction"] = CAPTURE_NOT_PROVIDED
    else:
        on = fraction(on_target, mapped)
        out["wes_capture_on_target_fraction"] = on
        out["wes_capture_off_target_fraction"] = "NA" if on == "NA" else f"{1 - float(on):.6f}"
    pks_duplicates, pks_reads = bounded("wes_pks_island_duplicates", "wes_pks_island_reads")
    out["wes_pks_island_reads"] = pks_reads
    out["wes_pks_island_duplicate_fraction"] = fraction(pks_duplicates, pks_reads)

    human, checked = bounded("rna_residual_human_reads", "rna_residual_checked_reads")
    bounded("rna_residual_checked_reads", "rna_residual_input_reads")
    # The check samples the final host-depleted reads, so its population is exactly
    # that count: after the pangenome passes when they ran, else after T2T + phiX.
    final = values.get("reads_after_pangenome", values.get("reads_after_t2t_phix"))
    drawn = values.get("rna_residual_input_reads")
    if drawn is not None and final is not None and drawn != final:
        raise ValueError(f"Impossible QC counts for {sample}: rna_residual_input_reads ({drawn}) "
                         f"is not the final host-depleted count ({final})")
    out["rna_residual_checked_reads"] = checked
    out["rna_residual_human_reads"] = human
    out["rna_residual_human_fraction"] = fraction(human, checked)
    return out


def mode_qc_columns(stages):
    """The mode QC columns this run's table carries, from its expected stages."""
    columns = []
    if any(stage in stages for stage in ("wes_input_qc", "wes_pks_duplicates")):
        columns += WES_QC_COLUMNS
    if "rna_residual_human" in stages:
        columns += RNA_QC_COLUMNS
    return columns


def table_columns(stages, metrics=None):
    """
    OUTPUT_COLUMNS with the run's mode QC columns just before `status`, and the optional
    L4 / V2 groups when any sample carries them (L4 after extract_low_score_mode, V2
    before the mode QC columns).
    """
    present = {metric for values in (metrics or {}).values() for metric in values}
    rescue = LOW_SCORE_RESCUE_COLUMNS if present & set(LOW_SCORE_RESCUE_COLUMNS) else []
    robustness = TIER_ROBUSTNESS_COLUMNS if present & set(TIER_ROBUSTNESS_COLUMNS) else []
    after = OUTPUT_COLUMNS.index("extract_low_score_mode") + 1
    at = OUTPUT_COLUMNS.index("status")
    return (OUTPUT_COLUMNS[:after] + rescue + OUTPUT_COLUMNS[after:at] + robustness
            + mode_qc_columns(stages) + OUTPUT_COLUMNS[at:])


def read_list_file(path):
    """One path per line. Blank lines ignored; every named file must be here."""
    paths = []
    for number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        name = line.strip()
        if not name:
            continue
        candidate = Path(name)
        if not candidate.exists():
            raise SystemExit(f"[ERROR] {path} line {number} names a file that is not here: {name}")
        paths.append(candidate)
    return paths


def parse_args():
    parser = argparse.ArgumentParser()
    # F02: a cohort of tens of thousands overflows the OS argument limit when every
    # input is its own argument. --inputs-file passes one filename instead.
    parser.add_argument("--inputs", nargs="*", default=[], type=Path)
    parser.add_argument("--inputs-file", default=None, type=Path,
                        help="file of QC fragment paths, one per line (preferred)")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--expected-samples",
        type=Path,
        help="file of sample IDs, one per line, that entered the run; any that "
             "produced no QC fragment are still given a row",
    )
    parser.add_argument(
        "--expected-stages", default=None,
        help="comma-separated stages this run should have produced QC for, from: "
             + ", ".join(STAGE_ORDER) + ". main.nf passes them; without it they are "
             "inferred as read_qc, host_depletion plus every stage any sample reported")
    parser.add_argument(
        "--excluded-samples", type=Path, default=None,
        help="sample<TAB>reason per line (reason: lower-case, digits, underscores); those "
             "samples end as excluded_<reason> whatever QC they have")
    args = parser.parse_args()
    if args.inputs and args.inputs_file:
        parser.error("give --inputs or --inputs-file, not both")
    if args.inputs_file:
        args.inputs = read_list_file(args.inputs_file)
    return args


def decimal_metric(metric, text, path):
    """A5 decimals: finite, non-negative, and within the metric's bound."""
    try:
        value = float(text)
    except ValueError:
        raise ValueError(f"Malformed {metric} in {path}: {text!r}") from None
    bound = DIVERGENCE_DECIMALS[metric]
    if not math.isfinite(value) or value < 0 or (bound is not None and value > bound):
        raise ValueError(f"Out-of-range {metric} in {path}: {text!r}")
    return value


def format_decimal(value):
    return str(int(value)) if float(value).is_integer() else f"{value:.6f}"


def check_divergence(sample, mapped):
    """A5: the divergence fragment must describe a subset of the reads it sits beside,
    and its own quantiles must be ordered."""
    measured = mapped.get(DIVERGENCE_COUNT)
    ceiling = mapped.get("reads_mapped_ihe3034_mapq40")
    if measured is not None and ceiling is not None and measured > ceiling:
        raise ValueError(f"Impossible QC counts for {sample}: {DIVERGENCE_COUNT} ({measured}) "
                         f"exceeds reads_mapped_ihe3034_mapq40 ({ceiling})")
    decimals = [m for m in DIVERGENCE_DECIMALS if mapped.get(m) is not None]
    if decimals and not measured:
        raise ValueError(f"Impossible QC values for {sample}: {', '.join(decimals)} "
                         f"reported with {DIVERGENCE_COUNT}={measured}")
    for chain in (("divergence_ihe3034_q1", "divergence_ihe3034_median", "divergence_ihe3034_q3"),
                  ("divergence_ihe3034_frac_gt_5pct", "divergence_ihe3034_frac_gt_2pct"),
                  ("clb_read_length_min", "clb_read_length_q1", "clb_read_length_median",
                   "clb_read_length_q3", "clb_read_length_max")):
        present = [(m, mapped[m]) for m in chain if mapped.get(m) is not None]
        for (low_name, low), (high_name, high) in zip(present, present[1:]):
            if low > high:
                raise ValueError(f"Impossible QC values for {sample}: {low_name} ({low}) "
                                 f"exceeds {high_name} ({high})")


def load_fragments(paths):
    metrics = defaultdict(dict)

    for path in paths:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["Sample", "Metric", "Value"]:
                raise ValueError(
                    f"Unexpected QC header in {path}: {reader.fieldnames}"
                )

            for row in reader:
                sample = row["Sample"].strip()
                metric = row["Metric"].strip()
                value_text = row["Value"].strip()

                if not sample or not metric:
                    raise ValueError(f"Malformed QC row in {path}: {row}")

                # Every metric is a count except the few that report an outcome. Those
                # are checked against their allowed values rather than waved through,
                # so a typo in a status is still an error.
                if metric in TEXT_METRICS:
                    if value_text not in TEXT_METRICS[metric]:
                        raise ValueError(
                            f"Unknown {metric} in {path}: {value_text!r}. "
                            f"Expected one of {', '.join(sorted(TEXT_METRICS[metric]))}")
                elif metric in DIVERGENCE_DECIMALS:
                    value = decimal_metric(metric, value_text, path)
                elif not value_text.isdigit():
                    raise ValueError(f"Malformed QC row in {path}: {row}")

                if metric in TEXT_METRICS:
                    value = value_text
                elif metric not in DIVERGENCE_DECIMALS:
                    value = int(value_text)
                previous = metrics[sample].get(metric)
                if previous is not None and previous != value:
                    raise ValueError(
                        f"Conflicting {metric} values for {sample}: "
                        f"{previous} and {value}"
                    )
                metrics[sample][metric] = value

    return metrics


def output_row(sample, values):
    mapped = {
        "Sample": sample,
        "input_alignment_records": values.get("input_alignment_records"),
        "total_primary_reads": values.get("total_primary_reads"),
        "extracted_unmapped_reads": values.get("extracted_unmapped_reads"),
        "extracted_low_score_mapped_reads": values.get("extracted_low_score_mapped_reads"),
        "extracted_reads": values.get("extracted_reads"),
        "extract_low_score_no_as_records": values.get("extract_low_score_no_as_records"),
        "extract_low_score_mode": values.get("extract_low_score_mode"),
        "filter_input_reads": values.get("filter_input_reads"),
        # A6: named for its MAPQ >= 40 filter (was reads_mapped_ihe3034).
        "reads_mapped_ihe3034_mapq40": values.get("reads_mapped_ihe3034_mapq40"),
        "reads_after_fastp": values.get("reads_after_fastp"),
        "reads_after_hg38": values.get("reads_after_hg38"),
        "reads_after_t2t_phix": values.get("reads_after_t2t_phix"),
        "reads_after_pangenome": values.get("reads_after_pangenome"),
        "num_clb_genes_align": values.get("num_clb_genes_align"),
        "reads_clb_genes_align": values.get("reads_clb_genes_align"),
        "hmm_ambiguous_reads": values.get("hmm_ambiguous_reads"),
        "taxonomy_status": values.get("taxonomy_status"),
        "clb_species_reported": values.get("clb_species_reported"),
        "num_clb_genes_hmm": values.get("num_clb_genes_hmm"),
        "reads_clb_genes_hmm": values.get("reads_clb_genes_hmm"),
    }
    for column in OUTPUT_COLUMNS:
        if column.startswith(("input_mapped_", "input_unmapped_", "prefilter_",
                              "diamond_", "candidate_")) or column in DIVERGENCE_COLUMNS:
            mapped[column] = values.get(column)

    missing = [column for column in REQUIRED if mapped[column] is None]
    mapped["status"] = "incomplete" if missing else "complete"

    # Reads can only be lost down this chain. Records come first because they count
    # secondary and supplementary alignments as well, so they are >= primary reads.
    #
    # E1: the extraction step in the chain is extracted_reads, everything extraction wrote
    # (FLAG 4 + low-score mapped), since that is what enters fastp. A fragment written
    # before E1 has only extracted_unmapped_reads, and then that was everything extracted.
    extraction_count = ("extracted_reads" if mapped["extracted_reads"] is not None
                        else "extracted_unmapped_reads")
    count_order = [
        "input_alignment_records",
        "total_primary_reads",
        extraction_count,
        "filter_input_reads",
        "reads_after_fastp",
        "reads_after_hg38",
        "reads_after_t2t_phix",
        "reads_after_pangenome",
        "reads_mapped_ihe3034_mapq40",
        "reads_clb_genes_align",
    ]
    observed = [
        (column, mapped[column])
        for column in count_order
        if mapped[column] is not None
    ]
    for (upstream_name, upstream), (downstream_name, downstream) in zip(
        observed, observed[1:]
    ):
        if downstream > upstream:
            raise ValueError(
                f"Impossible QC counts for {sample}: {downstream_name} "
                f"({downstream}) exceeds {upstream_name} ({upstream})"
            )

    hmm_reads = mapped["reads_clb_genes_hmm"]
    if mapped["reads_after_pangenome"] is not None:
        depleted_metric = "reads_after_pangenome"
    else:
        depleted_metric = "reads_after_t2t_phix"

    depleted_reads = mapped[depleted_metric]
    # Both operands have to exist. The ordering loop above filters None out of its own
    # comparisons and this check did not, so a sample missing its depletion metric raised
    # TypeError instead of being reported. That is reachable in production: with
    # sample_failure_strategy=ignore a failed filterReads leaves exactly this gap, and the
    # crash then takes down masterQCSummary and with it the whole cohort's QC table --
    # the precise outcome 'ignore' exists to prevent. A sample missing the metric is
    # already flagged by the `status` column; there is simply nothing to compare.
    if hmm_reads is not None and depleted_reads is not None and hmm_reads > depleted_reads:
        raise ValueError(
            f"Impossible QC counts for {sample}: reads_clb_genes_hmm "
            f"({hmm_reads}) exceeds {depleted_metric} ({depleted_reads})"
        )

    # E2 / P1: parts must add up to their totals. Each group is checked only when every
    # member is present, so a partially-run sample is reported, not rejected.
    for total_name, parts in (
        ("input_alignment_records",
         [c for c in OUTPUT_COLUMNS if c.startswith(("input_mapped_", "input_unmapped_"))]),
        ("prefilter_input_reads", ["prefilter_primary_reads", "prefilter_non_target_reads"]),
        ("candidate_total_reads", ["candidate_primary_reads", "candidate_rescued_reads"]),
        # E1: extraction's two classes add up to what it wrote.
        ("extracted_reads", ["extracted_unmapped_reads", "extracted_low_score_mapped_reads"]),
    ):
        values_present = [mapped.get(p) for p in parts]
        total = mapped.get(total_name)
        if total is None or any(v is None for v in values_present):
            continue
        if sum(values_present) != total:
            raise ValueError(
                f"Impossible QC counts for {sample}: {' + '.join(parts)} "
                f"({sum(values_present)}) != {total_name} ({total})")
    for rescued, pool in (("diamond_rescued_reads", "diamond_input_reads"),
                          ("diamond_input_reads", "prefilter_non_target_reads")):
        if mapped.get(rescued) is not None and mapped.get(pool) is not None \
                and mapped[rescued] > mapped[pool]:
            raise ValueError(f"Impossible QC counts for {sample}: {rescued} "
                             f"({mapped[rescued]}) exceeds {pool} ({mapped[pool]})")

    check_divergence(sample, mapped)

    # L4: a rescued read can only be lost down the chain, which starts with what extraction
    # added. Two chains: reads through depletion, then alignment records in the profiling
    # BAM (a long read's supplementary alignment is a second record of one read).
    reads_chain = ["extracted_low_score_mapped_reads"] + LOW_SCORE_RESCUE_COLUMNS[:3]
    for chain in (reads_chain, LOW_SCORE_RESCUE_COLUMNS[3:]):
        present = [(name, values[name]) for name in chain if values.get(name) is not None]
        for (upstream_name, upstream), (downstream_name, downstream) in zip(present, present[1:]):
            if downstream > upstream:
                raise ValueError(f"Impossible QC counts for {sample}: {downstream_name} "
                                 f"({downstream}) exceeds {upstream_name} ({upstream})")
    for column in OPTIONAL_COLUMNS:
        mapped[column] = values.get(column)

    for field in ("num_clb_genes_align", "num_clb_genes_hmm"):
        genes_detected = mapped[field]
        if genes_detected is not None and not 0 <= genes_detected <= 19:
            raise ValueError(f"Invalid {field} for {sample}: {genes_detected}")

    for column in DIVERGENCE_DECIMALS:
        if mapped[column] is not None:
            mapped[column] = format_decimal(mapped[column])

    mapped.update(mode_qc_values(sample, values))
    return {
        column: ("NA" if mapped.get(column) is None else mapped[column])
        for column in OUTPUT_COLUMNS + WES_QC_COLUMNS + RNA_QC_COLUMNS + OPTIONAL_COLUMNS
    }


def parse_stages(text):
    stages = [s.strip() for s in (text or "").split(",") if s.strip()]
    unknown = [s for s in stages if s not in STAGE_METRICS]
    if unknown:
        raise ValueError(f"Unknown stage(s) {', '.join(unknown)}; expected: {', '.join(STAGE_ORDER)}")
    return [s for s in STAGE_ORDER if s in stages]


def infer_stages(metrics):
    """For hand runs without --expected-stages: the always-run stages plus any stage
    at least one sample reported. A stage every sample failed cannot be inferred."""
    seen = {metric for values in metrics.values() for metric in values}
    return [stage for stage in STAGE_ORDER
            if stage in ("read_qc", "host_depletion")
            or any(m in seen for m in STAGE_METRICS[stage])]


def terminal_state(values, stages, excluded_reason=None):
    """Exactly one of completed_positive / completed_negative / completed_untiered /
    failed_at_<stage> / excluded_<reason>.

    `values` is the sample's raw metrics (empty for a sample that emitted no QC). A run
    whose stages include no evidence step (HMM-only profiling) has nothing to tier, so a
    sample that finished every stage is completed_untiered.
    """
    if excluded_reason:
        return f"excluded_{excluded_reason}"
    for stage in stages:
        if any(values.get(metric) is None for metric in STAGE_METRICS[stage]):
            return f"failed_at_{stage}"
    outcome = values.get("evidence_outcome")
    return COMPLETED_STATE.get(outcome, "completed_untiered")


def read_excluded(path):
    if path is None:
        return {}
    excluded = {}
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 2 or not fields[0].strip() or not EXCLUSION_REASON.match(fields[1].strip()):
            raise ValueError(f"{path} line {number}: expected sample<TAB>reason with a reason of "
                             f"lower-case letters, digits and underscores, got {line!r}")
        excluded[fields[0].strip()] = fields[1].strip()
    return excluded


def write_summary(metrics, output_path, expected=(), stages=None, excluded=None):
    excluded = excluded or {}
    if stages is None:
        stages = infer_stages(metrics)
    samples = sorted(set(metrics) | set(expected) | set(excluded))
    columns = table_columns(stages, metrics)
    with output_path.open("w", newline="") as handle:
        # extrasaction: the row carries metrics the table does not show, such as
        # reads_mapped_ihe3034_mapq40, which is still used for the monotonic check above.
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for sample in samples:
            if sample in metrics:
                row = output_row(sample, metrics[sample])
            else:
                row = {column: "NA" for column in columns}
                row["Sample"] = sample
                row["status"] = "no_qc_produced"
            row["terminal_state"] = terminal_state(metrics.get(sample, {}), stages,
                                                   excluded.get(sample))
            writer.writerow(row)


def read_expected(path):
    if path is None:
        return []
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def main():
    args = parse_args()
    metrics = load_fragments(args.inputs)
    expected = read_expected(args.expected_samples)
    excluded = read_excluded(args.excluded_samples)

    if not metrics and not expected and not excluded:
        raise ValueError("No QC records were supplied")

    stages = parse_stages(args.expected_stages) if args.expected_stages is not None else None
    write_summary(metrics, args.output, expected, stages, excluded)


if __name__ == "__main__":
    main()
