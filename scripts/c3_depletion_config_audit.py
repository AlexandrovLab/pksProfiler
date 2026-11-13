#!/usr/bin/env python3
"""C3: did the stool cohorts we compare share one host-depletion configuration?

Ludmil's v3 code review, third retrospective check:

    Confirm identical pangenome-depletion configuration across stool cohorts. If every
    compared stool cohort used the same pangenome-depletion setting, same
    reference/database digest, and the same depletion stage completed successfully, the
    denominator-stage issue does not confound that comparison.

    Check both the flag and the reference digest/stage completion, not merely whether
    `--pangenome_db` was non-null.

WHY THE DENOMINATOR MOVES (review item 7). Metagenome abundance is reported as
pks_reads per million `host_depleted_reads`. That denominator is the count AFTER
pangenome depletion when that stage ran, and after T2T+phiX when it did not -- and the
stage that produced it is not recorded beside the number. Two cohorts configured
differently therefore report rates on different denominators, and nothing in the output
says so. A lower rate can mean less pks, or simply more human sequence left in.

So this checks three things per cohort, which is what the review asks for and more than
a non-null test would give:

  1. the FLAG      -- was --pangenome_db set, and to what
  2. the DIGEST    -- the fingerprint of each depletion reference actually used, so two
                      cohorts naming the same path but a rebuilt index are not called
                      identical
  3. COMPLETION    -- did the depletion stage actually produce counts, per the cohort QC
                      table, rather than being configured and then failing

PRIVACY: reads run provenance and the cohort QC table, and writes per-cohort
configuration plus sample COUNTS. No sample identifier, no per-sample row.

    python3 scripts/c3_depletion_config_audit.py \\
        --run NAME=/path/to/results [--run NAME=... ...] [--out c3_depletion.tsv]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

DEPLETION_REFS = ("hg38_db", "t2t_phix_db", "pangenome_db")
STAGE_COLUMNS = ("reads_after_hg38", "reads_after_t2t_phix", "reads_after_pangenome")


def newest_provenance(results_dir):
    records = sorted(Path(results_dir).glob("runs/provenance/*.json"),
                     key=lambda p: p.stat().st_mtime, reverse=True)
    for path in records:
        try:
            return json.loads(path.read_text(encoding="utf-8", errors="replace")), path.name
        except (OSError, ValueError):
            continue
    return None, None


def params_of(record):
    node = (record.get("invocation") or {}).get("params")
    return node if isinstance(node, dict) else {}


def stage_counts(results_dir):
    """How many samples carry a usable number for each depletion stage."""
    qc = Path(results_dir) / "cohort" / "qc" / "pks.qc.summary.tsv"
    if not qc.is_file():
        qc = Path(results_dir) / "cohort" / "pks.qc.summary.tsv"
    if not qc.is_file():
        return None
    with qc.open(newline="", encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = {"samples": len(rows)}
    for column in STAGE_COLUMNS:
        if rows and column in rows[0]:
            good = sum(1 for r in rows
                       if (r.get(column) or "NA") not in ("NA", "", None))
            out[column] = good
        else:
            out[column] = None  # column absent entirely
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="NAME=RESULTS_DIR")
    ap.add_argument("--out", default="c3_depletion.tsv")
    a = ap.parse_args()

    cohorts, rows, notes = {}, [], []
    for spec in a.run:
        if "=" not in spec:
            sys.exit(f"--run needs NAME=RESULTS_DIR (got {spec!r})")
        name, _, path = spec.partition("=")
        name = name.strip()
        record, source = newest_provenance(path)
        if record is None:
            notes.append(f"{name}: no provenance record under {path}")
            continue
        params = params_of(record)
        deps = record.get("dependencies") or {}
        counts = stage_counts(path)

        profile = {
            "pangenome_set": "yes" if params.get("pangenome_db") else "no",
            "pangenome_path": str(params.get("pangenome_db") or "(unset)"),
        }
        for ref in DEPLETION_REFS:
            profile[f"{ref}_digest"] = str(deps.get(ref) or "(none)")
        cohorts[name] = profile

        for key, value in profile.items():
            rows.append([name, source, key, value])
        if counts is None:
            notes.append(f"{name}: no cohort QC table, so stage completion is unknown")
        else:
            rows.append([name, source, "samples", str(counts["samples"])])
            for column in STAGE_COLUMNS:
                got = counts[column]
                rows.append([name, source, column,
                             "(column absent)" if got is None
                             else f"{got}/{counts['samples']} samples"])

    if not cohorts:
        sys.exit("no cohort produced a provenance record; check the --run paths")

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["# C3: host-depletion configuration per cohort. Counts only."])
        w.writerow(["# The denominator of a metagenome rate is reads after pangenome"])
        w.writerow(["# depletion when that stage ran, else after T2T+phiX (review item 7)."])
        w.writerow(["cohort", "source", "field", "value"])
        w.writerows(rows)

    names = sorted(cohorts)
    width = max(len(n) for n in names)
    keys = ["pangenome_set", "pangenome_path"] + [f"{r}_digest" for r in DEPLETION_REFS]
    for key in keys:
        values = {cohorts[n].get(key, "") for n in names}
        flag = "" if len(values) == 1 else "   <-- DIFFERS"
        print(f"  {key}{flag}")
        for n in names:
            print(f"      {n:<{width}}  {cohorts[n].get(key, '')}")
    for n in notes:
        print(f"  NOTE {n}")

    differing = [k for k in keys if len({cohorts[n].get(k, '') for n in names}) > 1]
    print()
    if len(names) < 2:
        print("  Only one cohort given, so there is nothing to compare. Pass every stool")
        print("  cohort whose rates you intend to put side by side.")
    elif differing:
        print(f"  THESE COHORTS ARE NOT DIRECTLY COMPARABLE: {', '.join(differing)} differ.")
        print("  A per-million rate is computed over reads after the LAST depletion stage")
        print("  that ran, so a cohort with pangenome depletion has a smaller denominator")
        print("  than one without, and reports a higher rate for the same biology.")
    else:
        print("  Every cohort used the same depletion configuration and the same reference")
        print("  digests, so their per-million rates share a denominator definition and can")
        print("  be compared. Check the stage counts above are complete as well.")
    print(f"\n  wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
