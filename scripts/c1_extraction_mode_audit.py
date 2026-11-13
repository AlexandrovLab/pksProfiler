#!/usr/bin/env python3
"""C1: which extraction mode did each cohort's samples actually resolve to?

Ludmil's v3 code review, first retrospective check:

    Count `extract_low_score_mode` across every cohort QC summary. If every relevant WGS
    sample is `on_auto_bwa_mem` or `on_requested`, the lossy auto-fallback did not affect
    those samples. If any sample is `off_auto_not_bwa_mem`, review its upstream aligner
    provenance before interpreting it.

WHY IT MATTERS. scripts/header_provenance.py:154 reads

    return "on_auto_bwa_mem" if bwa_mem == "yes" else "off_auto_not_bwa_mem"

and bwa_mem has THREE values -- yes, no, unconfirmed. A BAM carrying a real bwa @PG whose
command line was stripped or rewritten resolves to `unconfirmed` and is then treated
exactly like a non-BWA aligner: it takes the FLAG4-only path and never gets low-AS mapped
reads rescued. Fewer reads reach the classifier, so reads, genes, breadth and tier can all
come out lower. That is a sensitivity loss, not a false positive.

Two things this answers:

  1. Did it fire at all?  Any `off_auto_not_bwa_mem` sample took the narrower path.
  2. Was a cohort MIXED?  Mixed resolution within one cohort is the dangerous case: some
     samples were more sensitive than others, and a prevalence difference between cohorts
     (or between arms) can then be an artefact of extraction rather than biology.

`off_requested` is NOT a finding: that is someone asking for off explicitly. Only the
`auto` path can fail closed without anyone deciding to.

PRIVACY: reads cohort QC summaries, which are derived from protected data, and writes ONLY
counts per cohort per mode -- no sample identifier, no row. Run it yourself.

    python3 scripts/c1_extraction_mode_audit.py \
        --qc NAME=/path/to/cohort/qc/pks.qc.summary.tsv [--qc NAME=... ...] \
        [--out c1_extraction_modes.tsv]

Point --qc at as many cohorts as you have; the label before `=` is what appears in the
output, so use the cohort name, never a path that leaks one.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# scripts/header_provenance.py LOW_SCORE_STATES
RESCUED = ("on_auto_bwa_mem", "on_requested")
NOT_RESCUED = ("off_auto_not_bwa_mem", "off_requested")
ORDER = ("on_auto_bwa_mem", "on_requested", "off_auto_not_bwa_mem", "off_requested")
COLUMN = "extract_low_score_mode"


def read_modes(path):
    """[mode, ...] for every sample row that carries the column."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    if not rows:
        return None, "empty QC summary"
    if COLUMN not in rows[0]:
        # A run from before the column existed, or a mode that does not extract.
        return None, f"no {COLUMN} column (pre-2.1 run, or a FASTQ-only cohort)"
    return [(r.get(COLUMN) or "").strip() for r in rows if (r.get(COLUMN) or "").strip()], None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc", action="append", required=True, metavar="NAME=PATH",
                    help="cohort label and its cohort/qc/pks.qc.summary.tsv; repeatable")
    ap.add_argument("--out", default="c1_extraction_modes.tsv")
    a = ap.parse_args()

    rows_out, notes, any_mixed, any_lossy = [], [], False, False

    for spec in a.qc:
        if "=" not in spec:
            sys.exit(f"--qc needs NAME=PATH (got {spec!r})")
        name, _, path = spec.partition("=")
        name = name.strip()
        if not Path(path).exists():
            notes.append(f"{name}: QC summary not found at that path; skipped")
            continue
        modes, why = read_modes(path)
        if modes is None:
            notes.append(f"{name}: {why}")
            continue

        counts = {m: modes.count(m) for m in ORDER if modes.count(m)}
        other = sorted(set(modes) - set(ORDER))
        for m in other:
            counts[m] = modes.count(m)
        n = len(modes)
        distinct = len(counts)
        mixed = distinct > 1
        lossy = counts.get("off_auto_not_bwa_mem", 0)
        any_mixed = any_mixed or mixed
        any_lossy = any_lossy or bool(lossy)

        for mode in list(ORDER) + other:
            if mode in counts:
                rows_out.append([name, mode, counts[mode], f"{100*counts[mode]/n:.1f}", n,
                                 "MIXED" if mixed else "uniform"])
        if other:
            notes.append(f"{name}: unexpected mode value(s) {other}")

    if not rows_out:
        sys.exit("no cohort produced any rows; check the --qc paths")

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["# C1: resolved extraction mode per cohort. Counts only, no identifiers."])
        w.writerow([f"# off_auto_not_bwa_mem is the fail-closed path (header_provenance.py:154);"])
        w.writerow([f"# off_requested is a deliberate choice and is not a finding."])
        w.writerow(["cohort", "extract_low_score_mode", "samples", "percent", "cohort_n", "uniformity"])
        w.writerows(rows_out)

    width = max(len(r[0]) for r in rows_out)
    for cohort, mode, n, pct, total, uni in rows_out:
        flag = "  <-- fail-closed path" if mode == "off_auto_not_bwa_mem" else ""
        print(f"  {cohort:<{width}}  {mode:<22} {n:>5} / {total:<5} {pct:>5}%  {uni}{flag}")
    for n in notes:
        print(f"  NOTE {n}")

    print()
    if any_mixed:
        print("  MIXED RESOLUTION FOUND. Samples within one cohort took different extraction")
        print("  paths, so their sensitivity is not comparable. Do not read a prevalence")
        print("  difference involving that cohort as biology until this is resolved.")
    elif any_lossy:
        print("  Uniform, but on the fail-closed path. No within-cohort confounding, though")
        print("  low-AS mapped reads were never rescued for those samples; check their")
        print("  upstream aligner provenance before comparing against a rescued cohort.")
    else:
        print("  Every sample resolved to a rescuing mode. The auto-fallback did not affect")
        print("  these results.")
    print(f"\n  wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
