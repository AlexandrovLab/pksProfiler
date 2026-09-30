#!/usr/bin/env python3
"""F16 tool versions, asserted on a real run: every tool-running task recorded what it ran,
from its own pinned conda env, and RUN_REPORT.txt compared it and found no mismatch.

Shared by assert_mode_outputs.py, assert_e2e_outputs.py and assert_normal_parity.py, and
runnable on its own:

    assert_tool_versions.py --results DIR --work-dir DIR
"""
import argparse
import csv
import json
import sys
from pathlib import Path

GOOD = ("ok", "unpinned", "version_unread")


def check_tool_versions(results, work_dir, check):
    """`check(label, condition, detail)` is the caller's reporter."""
    results, work_dir = Path(results), Path(work_dir)
    records = sorted(work_dir.glob("*/*/.tool_versions.tsv"))
    check("tasks wrote .tool_versions.tsv records", len(records) > 0, f"{len(records)} found")
    bad, outside, tools = [], [], {}
    for record in records:
        with record.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                tools.setdefault(row["tool"], set()).add(row["version"])
                if row["status"] not in GOOD:
                    bad.append(f"{row['process']}:{row['tool']}={row['status']}")
                prefix = row.get("conda_prefix") or ""
                if not prefix or not row["path"].startswith(prefix.rstrip("/") + "/"):
                    outside.append(f"{row['process']}:{row['tool']} {row['path']}")
    check("every recorded tool is its pinned build (status ok/unpinned)", not bad, ", ".join(bad[:5]))
    check("every recorded tool ran from its task's own conda env", not outside, ", ".join(outside[:5]))
    if "minimap2" in tools:
        check("minimap2 that ran is the pinned 2.31", tools["minimap2"] == {"2.31"},
              f"got {sorted(tools['minimap2'])}")

    report = results / "RUN_REPORT.txt"
    text = report.read_text() if report.exists() else ""
    check("RUN_REPORT.txt has TOOL VERSIONS (as run)", "TOOL VERSIONS (as run)" in text)
    check("RUN_REPORT.txt has no TOOL VERSION MISMATCH", "TOOL VERSION MISMATCH" not in text)
    sidecars = sorted((results / "runs" / "provenance").glob("*.json"),
                      key=lambda p: p.stat().st_mtime)
    versions = json.loads(sidecars[-1].read_text()).get("tool_versions", {}) if sidecars else {}
    check("provenance record lists no tool-version problems", versions.get("problems") == [],
          f"got {versions.get('problems')!r}")
    check("provenance record read every task's record",
          versions.get("tasks_recorded", 0) >= len(records) > 0,
          f"{versions.get('tasks_recorded')} of {len(records)}")
    return tools


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    a = ap.parse_args()
    failures = []

    def check(label, condition, detail=""):
        print(f"  {'ok  ' if condition else 'FAIL'}  {label}" + ("" if condition or not detail else f" ({detail})"))
        if not condition:
            failures.append(label)

    check_tool_versions(a.results, a.work_dir, check)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
