#!/usr/bin/env python3
"""Assert on the real output of tests/check_e2e_execution.sh's Nextflow run.

Every value asserted here was independently worked out by hand against the exact
fixtures tests/fixtures/generate_e2e_fixtures.py builds (see that script's LOCI
and negative-read generation) using the real conda-cached bowtie2/minimap2/
samtools binaries, before this was wired into a Nextflow run -- so a failure here
means the pipeline disagrees with a real, independently reproduced tool run, not
that the expected numbers were guessed.

Exits 1 with a message naming the first failed check; exits 0 and prints one
"ok" line per check on success.
"""
import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

# From tests/fixtures/generate_e2e_fixtures.py: 4 loci, one read pair (2 reads) each.
POSITIVE_EXPECTED_GENES = {"clbS", "clbN", "clbD", "clbA"}
POSITIVE_EXPECTED_READS = 2 * len(POSITIVE_EXPECTED_GENES)  # 8
POSITIVE_MIN_BREADTH = 0.01   # multi_gene floor
POSITIVE_MAX_BREADTH = 0.075  # must stay below broad_island, or the fixture drifted
NEGATIVE_PAIRS = 4
NEGATIVE_READS = 2 * NEGATIVE_PAIRS  # 8

FAILURES = []


def check(label, condition, detail=""):
    if condition:
        print(f"  ok    {label}")
    else:
        message = f"  FAIL  {label}" + (f" ({detail})" if detail else "")
        print(message)
        FAILURES.append(label)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mate_suffix_census(fastq_gz_path):
    """Return {base_read_id: sorted(['1','2', ...])} for a bgzip/gzip FASTQ."""
    bases = Counter()
    suffix_seen = {}
    with gzip.open(fastq_gz_path, "rt") as handle:
        for i, line in enumerate(handle):
            if i % 4 != 0:
                continue
            read_id = line.rstrip("\n")[1:].split(" ")[0]
            if "/" not in read_id:
                raise ValueError(f"{fastq_gz_path}: read {read_id!r} has no mate suffix")
            base, _, mate = read_id.rpartition("/")
            bases[base] += 1
            suffix_seen.setdefault(base, set()).add(mate)
    return bases, suffix_seen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--positive-sample", default="e2e_positive")
    parser.add_argument("--negative-sample", default="e2e_negative")
    parser.add_argument(
        "--expect-extracted-unmapped-reads", type=int, default=None,
        help="for a BAM/CRAM input run: both samples' extract.qc.tsv "
             "extracted_unmapped_reads must equal this (extractReads' own "
             "samtools-fastq-derived count -- the F01 read-routing metric).")
    args = parser.parse_args()

    results = args.results
    pos, neg = args.positive_sample, args.negative_sample

    # ---------- read_evidence.tsv (tumour tier classification) ----------
    pos_evidence_path = results / "by_sample" / pos / "read_evidence.tsv"
    neg_evidence_path = results / "by_sample" / neg / "read_evidence.tsv"
    check("positive sample's read_evidence.tsv exists", pos_evidence_path.exists())
    check("negative sample's read_evidence.tsv exists", neg_evidence_path.exists())
    if pos_evidence_path.exists():
        row = read_tsv(pos_evidence_path)[0]
        check("positive tier is multi_gene (real bowtie2 alignment + featureCounts)",
              row["read_evidence"] == "multi_gene", f"got {row['read_evidence']!r}")
        check("positive pks_reads == 8",
              row["pks_reads"] == str(POSITIVE_EXPECTED_READS), f"got {row['pks_reads']!r}")
        check("positive clb_genes_detected == 4",
              row["clb_genes_detected"] == str(len(POSITIVE_EXPECTED_GENES)),
              f"got {row['clb_genes_detected']!r}")
        breadth = float(row["island_breadth_1x"])
        check("positive island_breadth_1x is between the multi_gene and broad_island floors",
              POSITIVE_MIN_BREADTH <= breadth < POSITIVE_MAX_BREADTH, f"got {breadth}")
    if neg_evidence_path.exists():
        row = read_tsv(neg_evidence_path)[0]
        check("negative tier is negative (0 clb reads, real bowtie2 alignment)",
              row["read_evidence"] == "negative", f"got {row['read_evidence']!r}")
        check("negative pks_reads == 0", row["pks_reads"] == "0", f"got {row['pks_reads']!r}")

    # ---------- counts.txt (featureCounts per-gene matrix) ----------
    pos_counts_path = results / "by_sample" / pos / "counts.txt"
    if pos_counts_path.exists():
        genes_with_reads = set()
        total = 0
        with pos_counts_path.open() as handle:
            for line in handle:
                if line.startswith("#") or line.startswith("Geneid"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if not fields[0].startswith("clb"):
                    continue
                count = int(float(fields[-1] or 0))
                total += count
                if count > 0:
                    genes_with_reads.add(fields[0])
        check("counts.txt assigns reads to exactly the 4 simulated genes",
              genes_with_reads == POSITIVE_EXPECTED_GENES, f"got {sorted(genes_with_reads)}")
        check("counts.txt clb column total matches read_evidence's pks_reads",
              total == POSITIVE_EXPECTED_READS, f"got {total}")
    else:
        check("positive sample's counts.txt exists", False)

    # ---------- cohort QC summary ----------
    qc_path = results / "cohort" / "qc" / "pks.qc.summary.tsv"
    check("cohort QC summary exists", qc_path.exists())
    if qc_path.exists():
        rows = {row["Sample"]: row for row in read_tsv(qc_path)}
        for sample, expected_reads in ((pos, POSITIVE_EXPECTED_READS), (neg, NEGATIVE_READS)):
            row = rows.get(sample)
            check(f"{sample} has a QC row", row is not None)
            if row is None:
                continue
            check(f"{sample} status is complete", row["status"] == "complete", f"got {row['status']!r}")
            check(f"{sample} filter_input_reads == {expected_reads}",
                  row["filter_input_reads"] == str(expected_reads),
                  f"got {row['filter_input_reads']!r}")
            after_fastp = int(row["reads_after_fastp"])
            check(f"{sample} reads_after_fastp > 0 and <= filter_input_reads",
                  0 < after_fastp <= expected_reads, f"got {after_fastp}")
            # u-1 regression guard, cohort-level view: neither synthetic host reference
            # is related to either read set, so host depletion must remove nothing.
            check(f"{sample} reads_after_hg38 == reads_after_fastp (no spurious host match)",
                  row["reads_after_hg38"] == str(after_fastp), f"got {row['reads_after_hg38']!r}")
            check(f"{sample} reads_after_t2t_phix == reads_after_fastp (no spurious host match)",
                  row["reads_after_t2t_phix"] == str(after_fastp), f"got {row['reads_after_t2t_phix']!r}")
        check("positive reads_clb_genes_align == 8",
              rows.get(pos, {}).get("reads_clb_genes_align") == str(POSITIVE_EXPECTED_READS))
        check("negative reads_clb_genes_align == 0",
              rows.get(neg, {}).get("reads_clb_genes_align") == "0")

        # F01, from Ludmil's 2026-09-19 report: `samtools fastq` silently dropped a
        # read category depending on which -o/-0/-1/-2/-s routing flags were passed.
        # tests/test_extract_read_routing.py checks the invocation statically; this is
        # the same claim checked empirically, against extractReads' real output, for a
        # real BAM input run.
        if args.expect_extracted_unmapped_reads is not None:
            for sample in (pos, neg):
                row = rows.get(sample, {})
                check(f"{sample} extracted_unmapped_reads == {args.expect_extracted_unmapped_reads} "
                      "(F01: samtools fastq routed every category to the one output stream)",
                      row.get("extracted_unmapped_reads") == str(args.expect_extracted_unmapped_reads),
                      f"got {row.get('extracted_unmapped_reads')!r}")

    # ---------- u-1: both mates of every fragment survive host depletion ----------
    # This is the empirical version of the mate_pair_check_20260923 investigation:
    # that investigation found no live bug, but had no real pipeline run to check
    # against. This does. save_intermediates=true (set by check_e2e_execution.sh)
    # publishes exactly the file checked here.
    for sample, expected_pairs in ((pos, len(POSITIVE_EXPECTED_GENES)), (neg, NEGATIVE_PAIRS)):
        depleted_path = results / "by_sample" / sample / "intermediates" / f"{sample}.host_depleted.fastq.gz"
        check(f"{sample} host_depleted.fastq.gz (save_intermediates) exists", depleted_path.exists())
        if not depleted_path.exists():
            continue
        bases, suffixes = mate_suffix_census(depleted_path)
        check(f"{sample}: every original fragment ({expected_pairs}) appears exactly twice after host depletion",
              len(bases) == expected_pairs and set(bases.values()) == {2},
              f"got {len(bases)} distinct fragments, counts {sorted(set(bases.values()))}")
        check(f"{sample}: every fragment kept BOTH /1 and /2 (u-1 mate-suffix survival)",
              all(mates == {"1", "2"} for mates in suffixes.values()),
              f"got {[(k, sorted(v)) for k, v in suffixes.items() if v != {'1', '2'}]}")

    # ---------- cohort-wide roll-ups exist and mention both samples ----------
    for name in ("pks_cohort_report.tsv", "pks.master_summary.tsv"):
        path = results / "cohort" / name
        check(f"cohort/{name} exists", path.exists())
        if path.exists():
            rows = read_tsv(path)
            sample_column = "sample" if "sample" in (rows[0].keys() if rows else []) else "Sample"
            samples_seen = {row.get(sample_column) for row in rows}
            check(f"cohort/{name} has a row for both samples",
                  {pos, neg} <= samples_seen, f"got {samples_seen}")

    print()
    if FAILURES:
        print(f"e2e execution assertions: {len(FAILURES)} FAILED: {', '.join(FAILURES)}", file=sys.stderr)
        return 1
    print("e2e execution assertions: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
