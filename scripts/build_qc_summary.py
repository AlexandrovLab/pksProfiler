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
OUTPUT_COLUMNS = [
    "Sample",
    "input_alignment_records",   # from the index: records, secondary/supplementary included
    "total_primary_reads",       # exact primary count; only with --exact_input_counts
    "extracted_unmapped_reads",  # what extraction wrote
    "filter_input_reads",        # what entered fastp
    "reads_after_fastp",
    "reads_after_hg38",
    "reads_after_t2t_phix",
    "reads_after_pangenome",
    # reads_mapped_ihe3034 is deliberately NOT a column here. It is reads mapping
    # anywhere on the 5.1 Mb IHE3034 genome, of which the island is 1%, so beside the
    # clb columns it invites the reading "lots of reads, no island" as though the two
    # were the same measurement. It is still written per sample in
    # by_sample/<sample>/alignment/*.alignment.qc.tsv -- exported, as F09 asked, and
    # still used below to bound the clb count -- just not shown in the cohort table.
    "num_clb_genes_align",
    "reads_clb_genes_align",
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
]

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
}


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
    args = parser.parse_args()
    if args.inputs and args.inputs_file:
        parser.error("give --inputs or --inputs-file, not both")
    if args.inputs_file:
        args.inputs = read_list_file(args.inputs_file)
    return args


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
                elif not value_text.isdigit():
                    raise ValueError(f"Malformed QC row in {path}: {row}")

                value = value_text if metric in TEXT_METRICS else int(value_text)
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
        "filter_input_reads": values.get("filter_input_reads"),
        "reads_mapped_ihe3034": values.get("reads_mapped_ihe3034"),
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

    missing = [column for column in REQUIRED if mapped[column] is None]
    mapped["status"] = "incomplete" if missing else "complete"

    # Reads can only be lost down this chain. Records come first because they count
    # secondary and supplementary alignments as well, so they are >= primary reads.
    count_order = [
        "input_alignment_records",
        "total_primary_reads",
        "extracted_unmapped_reads",
        "filter_input_reads",
        "reads_after_fastp",
        "reads_after_hg38",
        "reads_after_t2t_phix",
        "reads_after_pangenome",
        "reads_mapped_ihe3034",
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

    for field in ("num_clb_genes_align", "num_clb_genes_hmm"):
        genes_detected = mapped[field]
        if genes_detected is not None and not 0 <= genes_detected <= 19:
            raise ValueError(f"Invalid {field} for {sample}: {genes_detected}")

    return {
        column: ("NA" if mapped.get(column) is None else mapped[column])
        for column in OUTPUT_COLUMNS
    }


def write_summary(metrics, output_path, expected=()):
    samples = sorted(set(metrics) | set(expected))
    with output_path.open("w", newline="") as handle:
        # extrasaction: the row carries metrics the table does not show, such as
        # reads_mapped_ihe3034, which is still used for the monotonic check above.
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for sample in samples:
            if sample in metrics:
                writer.writerow(output_row(sample, metrics[sample]))
            else:
                row = {column: "NA" for column in OUTPUT_COLUMNS}
                row["Sample"] = sample
                row["status"] = "no_qc_produced"
                writer.writerow(row)


def read_expected(path):
    if path is None:
        return []
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def main():
    args = parse_args()
    metrics = load_fragments(args.inputs)
    expected = read_expected(args.expected_samples)

    if not metrics and not expected:
        raise ValueError("No QC records were supplied")

    write_summary(metrics, args.output, expected)


if __name__ == "__main__":
    main()
