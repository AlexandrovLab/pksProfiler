#!/usr/bin/env python3
"""Assert on tests/check_e2e_execution.sh's E1 run (bwa mem -T 0 style BAM/CRAM input).

Ludmil, 25 Sep E1. Three samples built from the same records
(tests/fixtures/generate_low_score_bam_fixtures.py):
    --bwa-bam   @PG bwa mem -T 0, BAM   -> auto on: FLAG 4 + AS < 30 mapped reads
    --bwa-cram  the same, as CRAM       -> auto on, decoded with --cram_reference
    --star-bam  @PG STAR, BAM           -> auto off: FLAG 4 only
Checked: the extracted FASTQ read for read (names, /1 /2, sequence in sequencing
orientation), the QC counts and their add-up, the mode metric, the per-sample provenance
columns, the cohort summary and RUN_REPORT line, the run-log line, and that the recovered
pks reads reach the clb counts. Exits 1 on the first failed group.
"""
import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

FAILURES = []


def check(label, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {label}" + ("" if ok else f" ({detail})"))
    if not ok:
        FAILURES.append(label)


def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fastq_records(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    return sorted([lines[i][1:].split()[0], lines[i + 1], lines[i + 3]] for i in range(0, len(lines), 4))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--expected", required=True, type=Path)
    ap.add_argument("--log", required=True, type=Path)
    ap.add_argument("--bwa-bam", required=True)
    ap.add_argument("--bwa-cram", required=True)
    ap.add_argument("--star-bam", required=True)
    a = ap.parse_args()
    expected = json.loads(a.expected.read_text())
    samples = {a.bwa_bam: ("on", "on_auto_bwa_mem", "yes", "0"),
               a.bwa_cram: ("on", "on_auto_bwa_mem", "yes", "0"),
               a.star_bam: ("off", "off_auto_not_bwa_mem", "no", "NA")}
    qc = {r["Sample"]: r for r in read_tsv(a.results / "cohort" / "qc" / "pks.qc.summary.tsv")}
    log = a.log.read_text()

    for sample, (rule, state, bwa, score) in samples.items():
        exp = expected[rule]
        fq = a.results / "by_sample" / sample / "intermediates" / f"{sample}.UNMAPPED.fastq.gz"
        check(f"{sample}: extracted FASTQ exists (save_intermediates)", fq.exists())
        if fq.exists():
            got = fastq_records(fq)
            check(f"{sample}: extracted FASTQ is exactly the {rule} rule's {len(exp['fastq'])} reads, "
                  "same names, /1 /2 and orientation as samtools fastq -N",
                  got == exp["fastq"], f"got {[r[0] for r in got]}")
        row = qc.get(sample, {})
        for column in ("extracted_unmapped_reads", "extracted_low_score_mapped_reads", "extracted_reads"):
            check(f"{sample}: {column} == {exp[column]}", row.get(column) == str(exp[column]),
                  f"got {row.get(column)!r}")
        no_as = "NA" if exp["no_as"] is None else str(exp["no_as"])
        check(f"{sample}: extract_low_score_no_as_records == {no_as} (counted, never added)",
              row.get("extract_low_score_no_as_records") == no_as,
              f"got {row.get('extract_low_score_no_as_records')!r}")
        check(f"{sample}: extract_low_score_mode == {state}", row.get("extract_low_score_mode") == state,
              f"got {row.get('extract_low_score_mode')!r}")
        check(f"{sample}: filter_input_reads == extracted_reads (what fastp received)",
              row.get("filter_input_reads") == str(exp["extracted_reads"]), f"got {row.get('filter_input_reads')!r}")
        check(f"{sample}: reads_clb_genes_align == {exp['pks_reads']} (recovered pks reads reach the counts)",
              row.get("reads_clb_genes_align") == str(exp["pks_reads"]), f"got {row.get('reads_clb_genes_align')!r}")
        check(f"{sample}: num_clb_genes_align == {len(exp['pks_genes'])} ({', '.join(exp['pks_genes'])})",
              row.get("num_clb_genes_align") == str(len(exp["pks_genes"])),
              f"got {row.get('num_clb_genes_align')!r}")
        prov = a.results / "by_sample" / sample / "provenance" / "input_header.tsv"
        prow = read_tsv(prov)[0] if prov.exists() else {}
        check(f"{sample}: input_header.tsv aligner_bwa_mem={bwa} bwa_min_score={score} "
              f"extract_low_score_mapped={state}",
              (prow.get("aligner_bwa_mem"), prow.get("bwa_min_score"), prow.get("extract_low_score_mapped"))
              == (bwa, score, state), f"got {prow}")
        check(f"{sample}: input_header.tsv carries no command line",
              prov.exists() and "synthetic_hg38.fa" not in prov.read_text() and "CL:" not in prov.read_text())
        check(f"{sample}: the run log says which rule it got",
              f"extractReads {sample}: low-score mapped extraction {state}" in log)

    summary = a.results / "cohort" / "provenance" / "input_header_summary.tsv"
    states = {r["key"]: r["n_files"] for r in read_tsv(summary) if r["section"] == "extract_low_score"} \
        if summary.exists() else {}
    check("input_header_summary.tsv counts 2 on_auto_bwa_mem and 1 off_auto_not_bwa_mem",
          states == {"on_auto_bwa_mem": "2", "off_auto_not_bwa_mem": "1"}, f"got {states}")
    report = a.results / "RUN_REPORT.txt"
    text = report.read_text() if report.exists() else ""
    check("RUN_REPORT.txt shows the low-score extraction states and parameters",
          "low-score mapped extraction" in text and "extract_low_score_mapped" in text
          and "extract_min_as" in text)

    if FAILURES:
        print(f"{len(FAILURES)} E1 check(s) failed", file=sys.stderr)
        sys.exit(1)
    print("E1 low-score extraction checks passed")


if __name__ == "__main__":
    main()
