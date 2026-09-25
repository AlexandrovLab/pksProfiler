#!/usr/bin/env python3
"""Assert on the real output of tests/check_malformed_input_execution.sh's run.

Three samples, one cohort: a healthy BAM (real reads, real synthetic host references,
identical construction to check_e2e_execution.sh's "positive" fixture) plus two
deliberately malformed BAMs (see tests/fixtures/generate_malformed_input_fixtures.py).
`params.sample_failure_strategy` defaults to 'ignore' (main.nf), so a bad sample must
not silently succeed, must not be reported as though it completed, and must not take
the rest of the cohort down with it -- all three are checked here, against the real
files this run produced, not a static read of the source.

This is also the live half of Ludmil's finding-8 acceptance test
(tests/test_sample_stage_failure_policy.py's own docstring: "actually forcing a task to
fail needs a live Nextflow run, out of scope for a unit test").
"""
import argparse
import csv
import glob
import sys
from pathlib import Path

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path,
                         help="the -work-dir this run used, to resolve RUN_REPORT.txt's "
                              "'stderr: <work-dir>/<hash>*/.command.err' pointers")
    parser.add_argument("--healthy-sample", required=True)
    parser.add_argument("--truncated-sample", required=True)
    parser.add_argument("--bad-magic-sample", required=True)
    args = parser.parse_args()

    results = args.results
    healthy = args.healthy_sample
    malformed = (args.truncated_sample, args.bad_magic_sample)

    # ---------- the healthy sample: unaffected by its malformed cohort-mates ----------
    qc_path = results / "cohort" / "qc" / "pks.qc.summary.tsv"
    check("cohort QC summary exists", qc_path.exists())
    rows = {}
    if qc_path.exists():
        rows = {row["Sample"]: row for row in read_tsv(qc_path)}
        check(f"{healthy} has a QC row", healthy in rows)
        if healthy in rows:
            check(f"{healthy} status is complete (a bad cohort-mate must not affect it)",
                  rows[healthy]["status"] == "complete", f"got {rows[healthy]['status']!r}")
            # healthy.bam is built from generate_e2e_fixtures.py's "positive" fixture
            # (4 loci, one read pair each) aligned against the synthetic hg38
            # reference -- the same construction and the same expected count as
            # check_e2e_execution.sh's own BAM-input positive run.
            check(f"{healthy} extracted_unmapped_reads == 8",
                  rows[healthy].get("extracted_unmapped_reads") == "8",
                  f"got {rows[healthy].get('extracted_unmapped_reads')!r}")

        # ---------- the malformed samples: recorded, not silently dropped or faked ----------
        for sample in malformed:
            check(f"{sample} has a QC row (F09: it must be visible, not silently absent)",
                  sample in rows)
            if sample in rows:
                check(f"{sample} status is NOT complete (its extraction genuinely failed)",
                      rows[sample]["status"] != "complete", f"got {rows[sample]['status']!r}")

    for sample in malformed:
        read_evidence = results / "by_sample" / sample / "read_evidence.tsv"
        check(f"{sample} has no read_evidence.tsv (no fabricated tier for a sample "
              "whose extraction never produced reads)", not read_evidence.exists())

    # ---------- RUN_REPORT.txt: the failure is recorded, with a pointer to its own log ----------
    report_path = results / "RUN_REPORT.txt"
    check("RUN_REPORT.txt exists", report_path.exists())
    if not report_path.exists():
        print()
        if FAILURES:
            print(f"malformed-input assertions: {len(FAILURES)} FAILED: {', '.join(FAILURES)}", file=sys.stderr)
            return 1
        return 0

    report = report_path.read_text()
    check("RUN_REPORT.txt has a WHAT FAILED section", "WHAT FAILED" in report)

    stderr_globs = []
    if "WHAT FAILED" in report:
        failed_section = report[report.index("WHAT FAILED"):]
        check("WHAT FAILED names extractReads",
              "extractReads" in failed_section, "no extractReads failure listed")
        # Two independent malformed samples -> two independent failed tasks, not one
        # error suppressing or masking the other.
        check("WHAT FAILED lists two separate failed tasks (one per malformed sample)",
              failed_section.count("extractReads") >= 2,
              f"found {failed_section.count('extractReads')} extractReads mentions")
        for line in failed_section.splitlines():
            line = line.strip()
            if line.startswith("stderr:"):
                stderr_globs.append(line.split("stderr:", 1)[1].strip())

    check("RUN_REPORT.txt points to at least 2 .command.err logs for the failed tasks",
          len(stderr_globs) >= 2, f"found {len(stderr_globs)}: {stderr_globs}")

    # ---------- the diagnostics themselves are the real, specific tool messages ----------
    expected_substrings = {
        "truncated": "EOF block",
        "bad_magic": "neither BAM nor CRAM",
    }
    found_text = ""
    for pattern in stderr_globs:
        for match in glob.glob(pattern):
            found_text += Path(match).read_text(errors="replace") + "\n"

    for label, substring in expected_substrings.items():
        check(f"a failed task's .command.err contains the {label} diagnostic ({substring!r})",
              substring in found_text, f"not found in {len(stderr_globs)} referenced log(s)")

    print()
    if FAILURES:
        print(f"malformed-input assertions: {len(FAILURES)} FAILED: {', '.join(FAILURES)}", file=sys.stderr)
        return 1
    print("malformed-input assertions: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
