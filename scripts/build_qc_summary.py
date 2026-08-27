#!/usr/bin/env python3
"""Combine per-stage pksProfiler QC fragments into one cohort table."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


OUTPUT_COLUMNS = [
    "Sample",
    "input_reads",
    "unmapped_reads",
    "reads_after_fastp",
    "reads_after_hg38",
    "reads_after_t2t_phix",
    "reads_mapping_ihe3034",
    "reads_mapping_clb",
    "clb_genes_detected",
    "reads_passing_hmm_threshold",
    "clb_genes_detected_hmm",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


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

                if not sample or not metric or not value_text.isdigit():
                    raise ValueError(f"Malformed QC row in {path}: {row}")

                value = int(value_text)
                previous = metrics[sample].get(metric)
                if previous is not None and previous != value:
                    raise ValueError(
                        f"Conflicting {metric} values for {sample}: "
                        f"{previous} and {value}"
                    )
                metrics[sample][metric] = value

    if not metrics:
        raise ValueError("No QC records were supplied")

    return metrics


def output_row(sample, values):
    filter_input = values.get("filter_input_reads")
    input_reads = values.get("bam_input_primary_records", filter_input)
    unmapped_reads = values.get("extracted_unmapped_reads", filter_input)

    mapped = {
        "Sample": sample,
        "input_reads": input_reads,
        "unmapped_reads": unmapped_reads,
        "reads_after_fastp": values.get("reads_after_fastp"),
        "reads_after_hg38": values.get("reads_after_hg38"),
        "reads_after_t2t_phix": values.get("reads_after_t2t_phix"),
        "reads_mapping_ihe3034": values.get("reads_mapping_ihe3034"),
        "reads_mapping_clb": values.get("reads_mapping_clb"),
        "clb_genes_detected": values.get("clb_genes_detected"),
        "reads_passing_hmm_threshold": values.get(
            "reads_passing_hmm_threshold"
        ),
        "clb_genes_detected_hmm": values.get("clb_genes_detected_hmm"),
    }

    required = [
        "input_reads",
        "unmapped_reads",
        "reads_after_fastp",
        "reads_after_hg38",
        "reads_after_t2t_phix",
    ]
    missing = [column for column in required if mapped[column] is None]
    if missing:
        raise ValueError(
            f"Missing required QC metric(s) for {sample}: {', '.join(missing)}"
        )

    count_chains = [
        [
            "input_reads",
            "unmapped_reads",
            "reads_after_fastp",
            "reads_after_hg38",
            "reads_after_t2t_phix",
            "reads_mapping_ihe3034",
            "reads_mapping_clb",
        ],
        ["reads_after_t2t_phix", "reads_passing_hmm_threshold"],
    ]
    for count_order in count_chains:
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

    for field in ("clb_genes_detected", "clb_genes_detected_hmm"):
        genes_detected = mapped[field]
        if genes_detected is not None and not 0 <= genes_detected <= 19:
            raise ValueError(f"Invalid {field} for {sample}: {genes_detected}")

    return {
        column: ("NA" if mapped[column] is None else mapped[column])
        for column in OUTPUT_COLUMNS
    }


def write_summary(metrics, output_path):
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for sample in sorted(metrics):
            writer.writerow(output_row(sample, metrics[sample]))


def main():
    args = parse_args()
    metrics = load_fragments(args.inputs)
    write_summary(metrics, args.output)


if __name__ == "__main__":
    main()
