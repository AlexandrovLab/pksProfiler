#!/usr/bin/env python3
"""C2: what parameter values did the runs that produced our results actually use?

Ludmil's v3 code review, second retrospective check:

    Confirm `min_gene_overlap_bp` and tier thresholds in provenance for result-producing
    runs. If `min_gene_overlap_bp` remained at 1 bp and tier thresholds remained at
    shipped defaults, the overlap-plumbing defect did not alter already-produced primary
    results.

WHY. Review item 12: Modules/pksProfiler_align.nf passes `--largestOverlap` to
featureCounts but never passes `--minOverlap`, so params.min_gene_overlap_bp is declared
and then not plumbed through. At the shipped default of 1 bp that is harmless, because
featureCounts' own default is also 1 -- the parameter and the behaviour agree by
coincidence. Any OTHER value silently does nothing, and a run would report an overlap
rule it did not apply.

The tier thresholds matter for the same reason in reverse: a result is only comparable to
another result if both were called with the same ladder. These are read from what the run
recorded, not from the current source, which is the whole point -- the question is what
the run DID, not what today's defaults are.

This reads run provenance (parameters, references, versions). It does not open any
sequencing or clinical file, and prints one row per parameter per run.

    python3 scripts/c2_result_parameter_audit.py \\
        --run NAME=/path/to/results [--run NAME=... ...] [--out c2_parameters.tsv]

Point --run at a results directory; it finds runs/provenance/*.json, newest first, and
falls back to RUN_REPORT.txt if no JSON is there.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

# name -> the value shipped as default. A run matching all of these was called with the
# ladder v1.0.0 documents; anything else has to be stated wherever its numbers are used.
EXPECTED = {
    "min_gene_overlap_bp": "1",
    "tumor_multi_gene_min_pks_reads": "5",
    "tumor_multi_gene_min_clb_genes": "3",
    "tumor_multi_gene_min_breadth": "0.01",
    "tumor_broad_island_min_pks_reads": "30",
    "tumor_broad_island_min_clb_genes": "8",
    "tumor_broad_island_min_breadth": "0.075",
    "tumor_extensive_island_min_pks_reads": "100",
    "tumor_extensive_island_min_clb_genes": "10",
    "tumor_extensive_island_min_breadth": "0.15",
}
# min_gene_overlap_bp is the one where a non-default value is not merely a difference but
# a silent no-op, because the value never reaches featureCounts (review item 12).
NOT_PLUMBED = {"min_gene_overlap_bp"}


def same_number(a, b):
    try:
        return abs(float(a) - float(b)) < 1e-12
    except (TypeError, ValueError):
        return str(a).strip() == str(b).strip()


def from_json(path):
    try:
        record = json.loads(Path(path).read_text(encoding="utf-8", errors="replace"))
    except (OSError, ValueError):
        return None
    # finalize_provenance writes them under invocation.params; the older spellings are
    # kept as fallbacks so a record from an earlier release still reads.
    params = None
    for path in (("invocation", "params"), ("invocation", "parameters"),
                 ("params",), ("parameters",)):
        node = record
        for key in path:
            node = node.get(key) if isinstance(node, dict) else None
            if node is None:
                break
        if isinstance(node, dict) and node:
            params = node
            break
    if not isinstance(params, dict):
        return None
    return {k: ("" if v is None else str(v)) for k, v in params.items()}


def from_run_report(path):
    """Fallback: the PARAMETERS section prints `  name  value` pairs."""
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    found = {}
    for name in EXPECTED:
        m = re.search(rf"^\s*{re.escape(name)}\s+(\S+)\s*$", text, re.M)
        if m:
            found[name] = m.group(1)
    return found or None


def read_params(results_dir):
    """(params, source) for the newest provenance record under a results directory."""
    root = Path(results_dir)
    jsons = sorted(root.glob("runs/provenance/*.json"),
                   key=lambda p: p.stat().st_mtime if p.exists() else 0, reverse=True)
    for j in jsons:
        params = from_json(j)
        if params:
            return params, f"provenance/{j.name}"
    report = root / "RUN_REPORT.txt"
    if report.is_file():
        params = from_run_report(report)
        if params:
            return params, "RUN_REPORT.txt"
    return None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="NAME=RESULTS_DIR")
    ap.add_argument("--out", default="c2_parameters.tsv")
    a = ap.parse_args()

    rows, notes = [], []
    any_drift, any_noop, any_missing = False, False, False

    for spec in a.run:
        if "=" not in spec:
            sys.exit(f"--run needs NAME=RESULTS_DIR (got {spec!r})")
        name, _, path = spec.partition("=")
        name = name.strip()
        params, source = read_params(path)
        if params is None:
            notes.append(f"{name}: no provenance record found under {path}")
            continue

        for key, expected in EXPECTED.items():
            if key not in params:
                rows.append([name, source, key, "(not recorded)", expected, "NOT RECORDED"])
                any_missing = True
                continue
            actual = params[key]
            if same_number(actual, expected):
                rows.append([name, source, key, actual, expected, "default"])
            else:
                verdict = "NON-DEFAULT (NOT PLUMBED)" if key in NOT_PLUMBED else "NON-DEFAULT"
                rows.append([name, source, key, actual, expected, verdict])
                any_drift = True
                if key in NOT_PLUMBED:
                    any_noop = True

    if not rows:
        sys.exit("no run produced any parameters; check the --run paths")

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["# C2: parameter values the result-producing runs actually used."])
        w.writerow(["# Read from run provenance, not from the current source."])
        w.writerow(["run", "source", "parameter", "recorded_value", "shipped_default", "verdict"])
        w.writerows(rows)

    width = max(len(r[0]) for r in rows)
    for run, _src, key, actual, expected, verdict in rows:
        if verdict == "default":
            continue
        print(f"  {run:<{width}}  {key:<38} {actual:<10} (default {expected})  {verdict}")
    for n in notes:
        print(f"  NOTE {n}")

    print()
    runs = len({r[0] for r in rows})
    if not (any_drift or any_missing):
        print(f"  All {len(EXPECTED)} parameters sat at their shipped defaults in every one of")
        print(f"  the {runs} run(s) checked. The overlap-plumbing defect (item 12) did not alter")
        print("  these results, and every run used the same tier ladder.")
    else:
        if any_noop:
            print("  A run set min_gene_overlap_bp away from 1 but featureCounts never received")
            print("  it, so the overlap rule it reports is not the rule it applied. Those counts")
            print("  need re-deriving once the parameter is plumbed through.")
        if any_drift:
            print("  A tier threshold differs from the shipped ladder. Results from that run are")
            print("  not directly comparable with results called under the defaults.")
        if any_missing:
            print("  A parameter was not recorded at all, so what that run used cannot be")
            print("  established from its own provenance.")
    print(f"\n  wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
