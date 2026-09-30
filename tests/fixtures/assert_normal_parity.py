#!/usr/bin/env python3
"""A matched normal is processed exactly like its tumour: asserted on one real run.

Ludmil, 26 Sep mode memo: "Process matched normals identically. Matched normal/blood
provides a built-in specificity control and should use the same extraction, depletion,
profiling, and QC." tests/test_matched_normal_parity.py proves it by construction (the
pipeline carries no tumour/normal role at all); this is the execution half.

tests/check_mode_execution.sh (mode normal_parity) runs ONE synthetic BAM under two
sample names, one labelled as a tumour and one as its matched normal, in the same run.
Then, for every per-sample process:

  * both samples ran it, and no process ran for one and not the other;
  * the rendered task script (.command.sh) is byte-identical once the sample name is
    replaced -- same tools, flags, presets, indexes and thresholds;
  * the tool-version record (.tool_versions.tsv) is identical -- same binaries;

and every per-sample output (extraction/filter/depletion/alignment QC, gene counts,
coverage, read evidence, the host-depleted reads published by --save_intermediates) is
identical modulo the sample name, as are their rows of the cohort QC table and master
summary. Nothing here reads anything but this run's synthetic fixtures.

    assert_normal_parity.py --results DIR --work-dir DIR --tumour T --normal N
"""
import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assert_tool_versions import check_tool_versions  # noqa: E402

FAILURES = []
PLACEHOLDER = "<SAMPLE>"
# Text outputs compared modulo the sample name. Everything else (BAM/BAI, PDF, gzip) is
# compared by digest, except files that embed a timestamp or an absolute path.
TEXT_SUFFIXES = (".tsv", ".txt", ".bedgraph", ".csv")
NOT_COMPARED = (".html", ".json", ".pdf", ".svg", ".png")


def check(label, condition, detail=""):
    print(f"  {'ok  ' if condition else 'FAIL'}  {label}" + ("" if condition or not detail else f" ({detail})"))
    if not condition:
        FAILURES.append(label)


def tasks(work_dir):
    """{(process, tag): task dir} from each task's .command.run header."""
    found = {}
    for run in work_dir.glob("*/*/.command.run"):
        for line in run.read_text(errors="replace").splitlines()[:6]:
            m = re.match(r"^### name: '(.+?)(?: \((.+)\))?'$", line)
            if m:
                found[(m.group(1), m.group(2))] = run.parent
                break
    return found


def neutral(text, sample):
    return text.replace(sample, PLACEHOLDER)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def outputs(sample_dir, sample):
    """{relative path with the sample name neutralised: file}."""
    return {neutral(str(p.relative_to(sample_dir)), sample): p
            for p in sample_dir.rglob("*") if p.is_file()}


def rows(path, key):
    with path.open(newline="") as fh:
        return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--tumour", required=True)
    ap.add_argument("--normal", required=True)
    a = ap.parse_args()
    tumour, normal = a.tumour, a.normal
    print(f"[normal_parity] {tumour} vs {normal}")

    # ---- the same processes, with the same rendered scripts ----
    found = tasks(a.work_dir)
    by_sample = {s: {proc: d for (proc, tag), d in found.items() if tag == s} for s in (tumour, normal)}
    t_procs, n_procs = set(by_sample[tumour]), set(by_sample[normal])
    check("processes ran for the tumour", len(t_procs) >= 5, f"{sorted(t_procs)}")
    check("the normal ran exactly the tumour's processes", t_procs == n_procs,
          f"tumour only {sorted(t_procs - n_procs)}, normal only {sorted(n_procs - t_procs)}")
    for stage in ("extractReads", "filterReads", "mapReads", "pksProfilerAlign", "classifyPksReadEvidence"):
        check(f"{stage} ran for both", stage in t_procs and stage in n_procs)
    for proc in sorted(t_procs & n_procs):
        t_dir, n_dir = by_sample[tumour][proc], by_sample[normal][proc]
        t_sh = neutral((t_dir / ".command.sh").read_text(errors="replace"), tumour)
        n_sh = neutral((n_dir / ".command.sh").read_text(errors="replace"), normal)
        check(f"{proc}: task script identical modulo the sample name", t_sh == n_sh)
        t_rec, n_rec = t_dir / ".tool_versions.tsv", n_dir / ".tool_versions.tsv"
        if t_rec.exists() or n_rec.exists():
            check(f"{proc}: same tool binaries and versions",
                  t_rec.exists() and n_rec.exists() and t_rec.read_text() == n_rec.read_text())

    # ---- the same outputs ----
    t_out = outputs(a.results / "by_sample" / tumour, tumour)
    n_out = outputs(a.results / "by_sample" / normal, normal)
    check("the normal published exactly the tumour's per-sample files", set(t_out) == set(n_out),
          f"tumour only {sorted(set(t_out) - set(n_out))}, normal only {sorted(set(n_out) - set(t_out))}")
    compared = 0
    for rel in sorted(set(t_out) & set(n_out)):
        t_file, n_file = t_out[rel], n_out[rel]
        if rel.endswith(NOT_COMPARED) or rel.endswith(".bam") or rel.endswith(".bai"):
            # HTML/JSON/figures carry timestamps or paths; the BAM header carries the
            # sample-named file in its @PG command line. Their content is covered by the
            # counts, coverage and evidence tables compared below.
            continue
        if rel.endswith(TEXT_SUFFIXES):
            same = neutral(t_file.read_text(errors="replace"), tumour) == \
                   neutral(n_file.read_text(errors="replace"), normal)
        else:
            same = digest(t_file) == digest(n_file)
        check(f"{rel} identical", same)
        compared += 1
    check("compared the depletion, profiling and evidence outputs", compared >= 5, f"{compared} files")
    for needle in ("depletion.qc.tsv", "counts.txt", "read_evidence.tsv", "host_depleted.fastq.gz"):
        check(f"a {needle} was compared", any(needle in rel for rel in t_out))

    # ---- the same QC ----
    qc = rows(a.results / "cohort/qc/pks.qc.summary.tsv", "Sample")
    t_qc = {k: v for k, v in qc.get(tumour, {}).items() if k != "Sample"}
    n_qc = {k: v for k, v in qc.get(normal, {}).items() if k != "Sample"}
    check("cohort QC rows identical apart from the sample name", t_qc == n_qc and bool(t_qc),
          f"differs in {sorted(k for k in t_qc if t_qc.get(k) != n_qc.get(k))}")
    check("the pair has real pks evidence (not a vacuous comparison)",
          int(t_qc.get("reads_clb_genes_align") or 0) > 0, f"got {t_qc.get('reads_clb_genes_align')!r}")
    master = rows(a.results / "cohort/pks.master_summary.tsv", "sample")
    t_m = {k: v for k, v in master.get(tumour, {}).items() if k != "sample"}
    n_m = {k: v for k, v in master.get(normal, {}).items() if k != "sample"}
    check("master summary rows identical apart from the sample name", t_m == n_m and bool(t_m),
          f"differs in {sorted(k for k in t_m if t_m.get(k) != n_m.get(k))}")

    # ---- F16 ----
    check_tool_versions(a.results, a.work_dir, check)

    print()
    if FAILURES:
        print(f"[normal_parity] {len(FAILURES)} FAILED: {', '.join(FAILURES)}", file=sys.stderr)
        return 1
    print("[normal_parity] all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
