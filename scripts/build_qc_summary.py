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


OUTPUT_COLUMNS = [
    "Sample",
    "input_reads",
    "unmapped_reads",
    "reads_after_fastp",
    "reads_after_hg38",
    "reads_after_t2t_phix",
    "reads_after_pangenome",
    "num_clb_genes_align",
    "reads_clb_genes_align",
    "num_clb_genes_hmm",
    "reads_clb_genes_hmm",
    "status",
]

# complete       every required metric present
# incomplete     some required metrics missing -- the sample ran partially
# no_qc_produced expected from the sample sheet, but emitted no QC at all
REQUIRED = [
    "input_reads",
    "unmapped_reads",
    "reads_after_fastp",
    "reads_after_hg38",
    "reads_after_t2t_phix",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="*", default=[], type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--expected-samples",
        type=Path,
        help="file of sample IDs, one per line, that entered the run; any that "
             "produced no QC fragment are still given a row",
    )
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
        "reads_after_pangenome": values.get("reads_after_pangenome"),
        "num_clb_genes_align": values.get("num_clb_genes_align"),
        "reads_clb_genes_align": values.get("reads_clb_genes_align"),
        "num_clb_genes_hmm": values.get("num_clb_genes_hmm"),
        "reads_clb_genes_hmm": values.get("reads_clb_genes_hmm"),
    }

    missing = [column for column in REQUIRED if mapped[column] is None]
    mapped["status"] = "incomplete" if missing else "complete"

    count_order = [
        "input_reads",
        "unmapped_reads",
        "reads_after_fastp",
        "reads_after_hg38",
        "reads_after_t2t_phix",
        "reads_after_pangenome",
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
    if hmm_reads is not None and hmm_reads > depleted_reads:
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
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
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
