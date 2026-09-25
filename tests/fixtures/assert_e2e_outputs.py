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

# u-3 follow-up: cohort-scaling fixtures (tests/fixtures/generate_e2e_fixtures.py's
# BROAD_WINDOWS/EXTENSIVE_WINDOWS/BORDERLINE_WINDOWS). Every value below was worked out
# the same way the positive/negative ones above were -- a real bowtie2 --very-sensitive
# --no-unal | samtools view -q40 | sort, real featureCounts --largestOverlap, real
# samtools depth over the 50,768 bp island, and the real
# scripts/classify_tumor_pks_evidence.py -- run by hand against these exact fixtures
# before this was wired into a cohort Nextflow run. See that script for the tier floors
# (multi_gene/broad_island/extensive_island/localized_indeterminate).
TIER_FIXTURE_EXPECTATIONS = {
    # tier name -> (reads_clb_genes_align, clb_genes_detected, (min_breadth, max_breadth))
    "broad_island":              (34, 8, (0.075, 0.15)),
    "extensive_island":          (106, 12, (0.15, 1.0)),
    "localized_indeterminate":   (2, 1, (0.0, 0.01)),
}

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


def parse_cohort_samples(spec):
    """'"name:tier,name:tier"' -> [(name, tier), ...]. Empty string -> []."""
    if not spec:
        return []
    pairs = []
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        name, _, tier = item.partition(":")
        if not name or not tier:
            raise SystemExit(f"--cohort-samples entry {item!r} is not 'name:tier'")
        pairs.append((name, tier))
    return pairs


def check_hmm_evidence(results, sample, expected_reads, expected_genes, expected_gene_names=None):
    """pksProfiler_hmm.nf's own outputs: hmm_counts.tsv and hmm.qc.tsv.

    A separate code path from the bowtie2 lane's read_evidence.tsv/counts.txt --
    pksProfilerHMM/hmm_best_hit.py had no real execution coverage anywhere in this
    repo before u-3's follow-up (this run).
    """
    counts_path = results / "by_sample" / sample / "hmm" / "hmm_counts.tsv"
    qc_path = results / "by_sample" / sample / "hmm" / "hmm.qc.tsv"
    check(f"{sample} hmm_counts.tsv exists", counts_path.exists())
    check(f"{sample} hmm.qc.tsv exists", qc_path.exists())
    if counts_path.exists():
        genes_with_hits = set()
        total = 0
        for row in read_tsv(counts_path):
            count = int(row["Count"])
            total += count
            if count > 0:
                genes_with_hits.add(row["Gene"])
        check(f"{sample} hmm_counts.tsv total assigned reads == {expected_reads}",
              total == expected_reads, f"got {total}")
        check(f"{sample} hmm_counts.tsv detects exactly {expected_genes} gene(s)",
              len(genes_with_hits) == expected_genes, f"got {sorted(genes_with_hits)}")
        if expected_gene_names is not None:
            check(f"{sample} hmm_counts.tsv hits exactly {sorted(expected_gene_names)}",
                  genes_with_hits == expected_gene_names, f"got {sorted(genes_with_hits)}")
    if qc_path.exists():
        metrics = {row["Metric"]: row["Value"] for row in read_tsv(qc_path)}
        check(f"{sample} reads_clb_genes_hmm == {expected_reads}",
              metrics.get("reads_clb_genes_hmm") == str(expected_reads),
              f"got {metrics.get('reads_clb_genes_hmm')!r}")
        check(f"{sample} num_clb_genes_hmm == {expected_genes}",
              metrics.get("num_clb_genes_hmm") == str(expected_genes),
              f"got {metrics.get('num_clb_genes_hmm')!r}")
        check(f"{sample} hmm_ambiguous_reads == 0 (no cross-gene tie in a synthetic fixture)",
              metrics.get("hmm_ambiguous_reads") == "0",
              f"got {metrics.get('hmm_ambiguous_reads')!r}")


def check_read_evidence_tier(results, sample, expected_tier, expected_reads=None,
                              expected_genes=None, expected_breadth_range=None):
    """The bowtie2 lane's per-sample tier call (classify_tumor_pks_evidence.py's real
    output for this exact fixture), used both for the original positive/negative pair
    and for u-3's broad_island/extensive_island/localized_indeterminate additions."""
    path = results / "by_sample" / sample / "read_evidence.tsv"
    check(f"{sample}'s read_evidence.tsv exists", path.exists())
    if not path.exists():
        return
    row = read_tsv(path)[0]
    check(f"{sample} tier is {expected_tier}",
          row["read_evidence"] == expected_tier, f"got {row['read_evidence']!r}")
    if expected_reads is not None:
        check(f"{sample} pks_reads == {expected_reads}",
              row["pks_reads"] == str(expected_reads), f"got {row['pks_reads']!r}")
    if expected_genes is not None:
        check(f"{sample} clb_genes_detected == {expected_genes}",
              row["clb_genes_detected"] == str(expected_genes), f"got {row['clb_genes_detected']!r}")
    if expected_breadth_range is not None:
        breadth = float(row["island_breadth_1x"])
        lo, hi = expected_breadth_range
        check(f"{sample} island_breadth_1x is in [{lo}, {hi})",
              lo <= breadth < hi, f"got {breadth}")


def check_cohort_membership(results, expected_samples, strict=False):
    """Every cohort-wide roll-up carries exactly the expected sample set.

    With strict=True (a full, named cohort -- u-3's 5-sample run), membership is
    exact set equality AND one row per sample: a 2-sample check can't tell "both
    samples are in there" apart from "both samples are in there twice, and so is a
    stray third row" -- the same dedup/sorting boundary a 2-sample cohort can't
    exercise at all. With strict=False (kept for the original 2-sample runs above),
    it is the original, weaker containment check.
    """
    for name in ("pks_cohort_report.tsv", "pks.master_summary.tsv"):
        path = results / "cohort" / name
        check(f"cohort/{name} exists", path.exists())
        if not path.exists():
            continue
        rows = read_tsv(path)
        sample_column = "sample" if "sample" in (rows[0].keys() if rows else []) else "Sample"
        samples_seen = [row.get(sample_column) for row in rows]
        if strict:
            check(f"cohort/{name} has exactly one row per expected sample, no extras",
                  sorted(samples_seen) == sorted(expected_samples),
                  f"got {sorted(samples_seen)}, expected {sorted(expected_samples)}")
        else:
            check(f"cohort/{name} has a row for both samples",
                  set(expected_samples) <= set(samples_seen), f"got {set(samples_seen)}")

    if strict:
        # masterTableAlign (Gene x Sample matrix): membership is column headers, not rows.
        gene_matrix_path = results / "cohort" / "gene_counts" / "pks.gene.counts.align.txt"
        check("cohort/gene_counts/pks.gene.counts.align.txt exists", gene_matrix_path.exists())
        if gene_matrix_path.exists():
            with gene_matrix_path.open(newline="") as handle:
                header = next(csv.reader(handle, delimiter="\t"))
            columns = set(header[1:])
            check("pks.gene.counts.align.txt has exactly one column per expected sample",
                  columns == set(expected_samples), f"got {sorted(columns)}")

        qc_path = results / "cohort" / "qc" / "pks.qc.summary.tsv"
        if qc_path.exists():
            rows = {row["Sample"] for row in read_tsv(qc_path)}
            check("cohort QC summary has exactly one row per expected sample",
                  rows == set(expected_samples), f"got {sorted(rows)}")


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


def read_names(fastq_gz_path):
    """Return the list of read ids (QNAME, whitespace-truncated) in a gzip/bgzip FASTQ."""
    names = []
    with gzip.open(fastq_gz_path, "rt") as handle:
        for i, line in enumerate(handle):
            if i % 4 != 0:
                continue
            names.append(line.rstrip("\n")[1:].split(" ")[0])
    return names


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
    parser.add_argument(
        "--profiling-method", choices=["bowtie2", "hmm"], default="bowtie2",
        help="which profiling lane produced these results: the bowtie2/featureCounts "
             "read_evidence.tsv + counts.txt lane (default), or pksProfilerHMM's "
             "hmm_counts.tsv + hmm.qc.tsv lane (u-3 follow-up).")
    parser.add_argument(
        "--cohort-samples", default="",
        help="'name:tier,name:tier,...' for samples beyond --positive-sample/"
             "--negative-sample (u-3 follow-up's broad_island/extensive_island/"
             "localized_indeterminate cohort-scaling fixtures). Their tier is checked "
             "exactly, and cohort-wide roll-ups are checked for exact membership "
             "(one row per sample, none missing, none duplicated) across the full "
             "named set instead of the weaker 'contains both' check used without it.")
    parser.add_argument(
        "--single-end", action="store_true",
        help="the positive/negative samples are single-end FASTQ (fastq1 only, no "
             "fastq2). filterReads never mate-suffixes a single-file input (see "
             "Modules/filter_reads.nf's input_list.size()==1 branch), so the u-1 "
             "mate-pair survival check does not apply; a simpler single-read survival "
             "check runs instead. Use --positive-expected-reads/--positive-expected-"
             "genes/--positive-min-breadth/--positive-max-breadth/--positive-gene-"
             "names/--negative-expected-reads to describe a single-end fixture, which "
             "carries half the bases (and so a different tier and gene set) of its "
             "paired-end namesake built from the same window positions.")
    parser.add_argument("--positive-expected-reads", type=int, default=None)
    parser.add_argument("--positive-expected-genes", type=int, default=None)
    parser.add_argument("--positive-min-breadth", type=float, default=None)
    parser.add_argument("--positive-max-breadth", type=float, default=None)
    parser.add_argument(
        "--positive-gene-names", default=None,
        help="comma-separated clb gene names the positive sample must detect "
             "(defaults to the paired-end fixture's clbS,clbN,clbD,clbA).")
    parser.add_argument("--negative-expected-reads", type=int, default=None)
    parser.add_argument("--positive-expected-tier", default="multi_gene")
    parser.add_argument(
        "--positive-source-fastq", default=None,
        help="--single-end only: the exact FASTQ staged as the positive sample's "
             "fastq1, so post-host-depletion read identity can be checked against it.")
    parser.add_argument(
        "--negative-source-fastq", default=None,
        help="--single-end only: same as --positive-source-fastq, for the negative sample.")
    args = parser.parse_args()

    results = args.results
    pos, neg = args.positive_sample, args.negative_sample
    extra_samples = parse_cohort_samples(args.cohort_samples)
    all_samples = [pos, neg] + [name for name, _tier in extra_samples]

    # Overridable positive/negative expectations -- default to the original paired-end
    # fixture's numbers so every existing call site (BAM/CRAM/HMM/cohort) is unaffected.
    pos_reads = args.positive_expected_reads if args.positive_expected_reads is not None else POSITIVE_EXPECTED_READS
    pos_genes = args.positive_expected_genes if args.positive_expected_genes is not None else len(POSITIVE_EXPECTED_GENES)
    pos_min_breadth = args.positive_min_breadth if args.positive_min_breadth is not None else POSITIVE_MIN_BREADTH
    pos_max_breadth = args.positive_max_breadth if args.positive_max_breadth is not None else POSITIVE_MAX_BREADTH
    pos_gene_names = (set(args.positive_gene_names.split(",")) if args.positive_gene_names
                      else POSITIVE_EXPECTED_GENES)
    neg_reads = args.negative_expected_reads if args.negative_expected_reads is not None else NEGATIVE_READS

    if args.profiling_method == "bowtie2":
        # ---------- read_evidence.tsv (tumour tier classification) ----------
        check_read_evidence_tier(results, pos, args.positive_expected_tier,
                                  expected_reads=pos_reads,
                                  expected_genes=pos_genes,
                                  expected_breadth_range=(pos_min_breadth, pos_max_breadth))
        check_read_evidence_tier(results, neg, "negative", expected_reads=0)

        for sample, tier in extra_samples:
            expected_reads, expected_genes, breadth_range = TIER_FIXTURE_EXPECTATIONS[tier]
            check_read_evidence_tier(results, sample, tier,
                                      expected_reads=expected_reads,
                                      expected_genes=expected_genes,
                                      expected_breadth_range=breadth_range)

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
            check(f"counts.txt assigns reads to exactly the {len(pos_gene_names)} simulated genes",
                  genes_with_reads == pos_gene_names, f"got {sorted(genes_with_reads)}")
            check("counts.txt clb column total matches read_evidence's pks_reads",
                  total == pos_reads, f"got {total}")
        else:
            check("positive sample's counts.txt exists", False)
    else:
        # ---------- pksProfilerHMM: hmm_counts.tsv + hmm.qc.tsv ----------
        # Same fixture reads as the bowtie2 lane above (see check_e2e_execution.sh),
        # run for real through nhmmscan/hmm_best_hit.py instead -- a code path this
        # repo shipped since v0.0.1 with zero prior execution coverage.
        check_hmm_evidence(results, pos, expected_reads=POSITIVE_EXPECTED_READS,
                            expected_genes=len(POSITIVE_EXPECTED_GENES),
                            expected_gene_names=POSITIVE_EXPECTED_GENES)
        check_hmm_evidence(results, neg, expected_reads=0, expected_genes=0,
                            expected_gene_names=set())

    # ---------- cohort QC summary ----------
    qc_path = results / "cohort" / "qc" / "pks.qc.summary.tsv"
    check("cohort QC summary exists", qc_path.exists())
    if qc_path.exists():
        rows = {row["Sample"]: row for row in read_tsv(qc_path)}
        expected_reads_by_sample = {pos: pos_reads, neg: neg_reads}
        for sample, tier in extra_samples:
            expected_reads_by_sample[sample] = TIER_FIXTURE_EXPECTATIONS[tier][0]

        for sample, expected_reads in expected_reads_by_sample.items():
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

        if args.profiling_method == "bowtie2":
            check(f"positive reads_clb_genes_align == {pos_reads}",
                  rows.get(pos, {}).get("reads_clb_genes_align") == str(pos_reads))
            check("negative reads_clb_genes_align == 0",
                  rows.get(neg, {}).get("reads_clb_genes_align") == "0")
            for sample, tier in extra_samples:
                expected_reads = TIER_FIXTURE_EXPECTATIONS[tier][0]
                check(f"{sample} reads_clb_genes_align == {expected_reads}",
                      rows.get(sample, {}).get("reads_clb_genes_align") == str(expected_reads),
                      f"got {rows.get(sample, {}).get('reads_clb_genes_align')!r}")

        # F01, from Ludmil's 2026-09-19 report: `samtools fastq` silently dropped a
        # read category depending on which -o/-0/-1/-2/-s routing flags were passed.
        # tests/test_extract_read_routing.py checks the invocation statically; this is
        # the same claim checked empirically, against extractReads' real output, for a
        # real BAM/CRAM input run.
        if args.expect_extracted_unmapped_reads is not None:
            for sample in (pos, neg):
                row = rows.get(sample, {})
                check(f"{sample} extracted_unmapped_reads == {args.expect_extracted_unmapped_reads} "
                      "(F01: samtools fastq routed every category to the one output stream)",
                      row.get("extracted_unmapped_reads") == str(args.expect_extracted_unmapped_reads),
                      f"got {row.get('extracted_unmapped_reads')!r}")

    if args.single_end:
        # ---------- single-end: every read survives host depletion, none duplicated ----------
        # filterReads never mate-suffixes a single-file input (Modules/filter_reads.nf's
        # input_list.size()==1 branch skips the awk tagging block entirely and calls
        # fastp -i directly), so there is no mate to lose and mate_suffix_census's "every
        # read id has a /1 or /2" assumption does not hold here -- it would raise on the
        # first untagged read id. This is the single-end analogue: read identity and
        # count survive host depletion intact, checked directly against the exact FASTQ
        # bytes staged for the run rather than a fixture-side read count, so a real
        # duplication or drop is caught even if it happened to preserve the total.
        expected_names_by_sample = {pos: args.positive_source_fastq, neg: args.negative_source_fastq}
        for sample, source_fastq in expected_names_by_sample.items():
            depleted_path = results / "by_sample" / sample / "intermediates" / f"{sample}.host_depleted.fastq.gz"
            check(f"{sample} host_depleted.fastq.gz (save_intermediates) exists", depleted_path.exists())
            if not depleted_path.exists() or not source_fastq:
                continue
            source_names = read_names(Path(source_fastq))
            depleted_names = read_names(depleted_path)
            check(f"{sample}: every single-end read ({len(source_names)}) survives host depletion, "
                  "none dropped, none duplicated",
                  sorted(depleted_names) == sorted(source_names),
                  f"source had {len(source_names)} reads, depleted has {len(depleted_names)}; "
                  f"missing={sorted(set(source_names) - set(depleted_names))}, "
                  f"extra={sorted(set(depleted_names) - set(source_names))}")
    else:
        # ---------- u-1: both mates of every fragment survive host depletion ----------
        # This is the empirical version of the mate_pair_check_20260923 investigation:
        # that investigation found no live bug, but had no real pipeline run to check
        # against. This does. save_intermediates=true (set by check_e2e_execution.sh)
        # publishes exactly the file checked here. Independent of --profiling-method: host
        # depletion runs upstream of the align/hmm branch either way.
        expected_pairs_by_sample = {pos: pos_reads // 2, neg: neg_reads // 2}
        for sample, tier in extra_samples:
            expected_pairs_by_sample[sample] = TIER_FIXTURE_EXPECTATIONS[tier][0] // 2

        for sample, expected_pairs in expected_pairs_by_sample.items():
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

    # ---------- cohort-wide roll-ups exist and carry every sample ----------
    check_cohort_membership(results, all_samples, strict=bool(extra_samples))

    print()
    if FAILURES:
        print(f"e2e execution assertions: {len(FAILURES)} FAILED: {', '.join(FAILURES)}", file=sys.stderr)
        return 1
    print("e2e execution assertions: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
