#!/usr/bin/env python3
"""Assert on one real run of tests/check_mode_execution.sh (26 Sep mode memo section 6).

For each mode: which branch ran (the evidence file of that mode and no other), which
outputs appeared, what was intentionally disabled (no tier outside short-read wgs, no
targeted assembly outside it, no MAGs), the depletion preset the host pass really used,
and the terminal state the QC summary gives the sample.

    assert_mode_outputs.py --mode wgs|wes|rna|metagenome|wgs_ont|rna_ont \
        --results DIR --work-dir DIR --sample S --fixtures DIR/mode_fixtures.json [--log F]
"""
import argparse
import csv
import json
import sys
from pathlib import Path

FAILURES = []


def check(label, condition, detail=""):
    print(f"  {'ok  ' if condition else 'FAIL'}  {label}" + ("" if condition or not detail else f" ({detail})"))
    if not condition:
        FAILURES.append(label)


def one_row(path):
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else {}


def rows_by(path, key):
    with path.open(newline="") as fh:
        return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}


def host_commands(work_dir):
    """The rendered mapReads scripts of this run (pipeline code, not data)."""
    return [p.read_text(errors="replace") for p in work_dir.glob("*/*/.command.sh")
            if "reads_after_hg38" in p.read_text(errors="replace")]


EXPECT = {
    #  mode       evidence file         tier/label field       expected terminal state
    "wgs":        ("read_evidence.tsv", "multi_gene",          "completed_positive"),
    "wes":        ("wes_evidence.tsv",  None,                  "completed_untiered"),
    "rna":        ("rna_evidence.tsv",  None,                  "completed_untiered"),
    "metagenome": ("read_evidence.tsv", "untiered",            "completed_untiered"),
    "wgs_ont":    ("read_evidence.tsv", "untiered",            "completed_untiered"),
    "rna_ont":    ("rna_evidence.tsv",  None,                  "completed_untiered"),
}
EVIDENCE_FILES = ("read_evidence.tsv", "wes_evidence.tsv", "rna_evidence.tsv")
READS_FIELD = {"read_evidence.tsv": "pks_reads", "wes_evidence.tsv": "clb_wes_reads",
               "rna_evidence.tsv": "clb_rna_reads"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=sorted(EXPECT))
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fixtures", required=True, type=Path)
    ap.add_argument("--log", type=Path, default=None)
    a = ap.parse_args()
    fx = json.loads(a.fixtures.read_text())
    mode, sample = a.mode, a.sample
    sdir = a.results / "by_sample" / sample
    evidence_file, tier, terminal = EXPECT[mode]
    print(f"[{mode}] {sample}")

    # ---- the branch: this mode's evidence file, and no other ----
    ev_path = sdir / evidence_file
    check(f"{evidence_file} exists", ev_path.exists())
    for other in EVIDENCE_FILES:
        if other != evidence_file:
            check(f"no {other} (other modes' schema not written)", not (sdir / other).exists())
    ev = one_row(ev_path) if ev_path.exists() else {}
    expected_mode = mode.split("_")[0]
    check(f"evidence mode column is {expected_mode}", ev.get("mode") == expected_mode, f"got {ev.get('mode')!r}")
    if tier is not None:
        check(f"read_evidence is {tier}", ev.get("read_evidence") == tier, f"got {ev.get('read_evidence')!r}")

    reads = {"wgs": fx["short_positive_reads"], "wes": fx["short_positive_reads"],
             "metagenome": fx["short_positive_reads"], "rna": fx["rna_short_pks_reads"],
             "wgs_ont": fx["long_pks_reads"], "rna_ont": fx["long_pks_reads"]}[mode]
    got = ev.get(READS_FIELD[evidence_file])
    check(f"{READS_FIELD[evidence_file]} == {reads}", got == str(reads), f"got {got!r}")

    # ---- disabled on purpose ----
    # wgs: the positive fixture is multi_gene, below --tumor_contig_tiers, so it is not
    # assembled either; every other mode never assembles.
    check("no targeted assembly (contigs/ absent)", not (sdir / "contigs").exists())
    check("no MAG lane (genomes/ absent)", not (sdir / "genomes").exists())
    if mode == "metagenome":
        depleted = ev.get("host_depleted_reads")
        check("host_depleted_reads == 8 (nothing in the fixture is host)", depleted == "8", f"got {depleted!r}")
        check("pks_reads_per_million == 1000000.0000",
              ev.get("pks_reads_per_million") == "1000000.0000", f"got {ev.get('pks_reads_per_million')!r}")
    if mode == "rna":
        check("rna_library_type recorded (polyA)", ev.get("rna_library_type") == "polyA",
              f"got {ev.get('rna_library_type')!r}")

    # ---- depletion: what was removed, and with which preset ----
    qc = rows_by(a.results / "cohort/qc/pks.qc.summary.tsv", "Sample").get(sample, {})
    check("QC status complete", qc.get("status") == "complete", f"got {qc.get('status')!r}")
    check(f"terminal_state {terminal}", qc.get("terminal_state") == terminal, f"got {qc.get('terminal_state')!r}")
    host_in, host_out = {
        "wgs": (fx["short_positive_reads"], fx["short_positive_reads"]),
        "wes": (fx["short_positive_reads"], fx["short_positive_reads"]),
        "metagenome": (fx["short_positive_reads"], fx["short_positive_reads"]),
        "rna": (fx["rna_short_input"], fx["rna_short_pks_reads"]),
        "wgs_ont": (fx["wgs_ont_input"], fx["long_pks_reads"]),
        "rna_ont": (fx["rna_ont_input"], fx["long_pks_reads"]),
    }[mode]
    check(f"filter_input_reads == {host_in}", qc.get("filter_input_reads") == str(host_in),
          f"got {qc.get('filter_input_reads')!r}")
    check(f"reads_after_hg38 == {host_out} (every human read removed, every pks read kept)",
          qc.get("reads_after_hg38") == str(host_out), f"got {qc.get('reads_after_hg38')!r}")
    commands = host_commands(a.work_dir)
    check("found the rendered host-depletion command", len(commands) >= 1, f"{len(commands)} found")
    human_preset = {"wgs": "-ax sr ", "wes": "-ax sr ", "metagenome": "-ax sr ",
                    "rna": "-ax splice:sr ", "wgs_ont": "-ax map-ont ", "rna_ont": "-ax splice "}[mode]
    check(f"human pass used minimap2 {human_preset.strip()}", any(human_preset in c for c in commands))
    if mode in ("rna", "rna_ont"):
        check("human pass got the junction BED (-j)", any("-j " in c for c in commands))
    if mode.endswith("_ont"):
        check("no -ax sr call for long reads", not any("-ax sr " in c for c in commands))
        check("no --frag for long reads", not any("--frag" in c for c in commands))

    # ---- master summary carries the mode ----
    master = rows_by(a.results / "cohort/pks.master_summary.tsv", "sample").get(sample, {})
    check(f"master summary mode is {expected_mode}", master.get("mode") == expected_mode,
          f"got {master.get('mode')!r}")
    check("master summary terminal_state matches", master.get("terminal_state") == terminal,
          f"got {master.get('terminal_state')!r}")
    if mode == "wgs_ont" and a.log:
        check("the run logged that targeted assembly was switched off",
              "targeted recruitment and assembly switched off" in a.log.read_text(errors="replace"))

    print()
    if FAILURES:
        print(f"[{mode}] {len(FAILURES)} FAILED: {', '.join(FAILURES)}", file=sys.stderr)
        return 1
    print(f"[{mode}] all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
