#!/usr/bin/env python3
"""Assert on one real run of tests/check_mode_execution.sh (26 Sep mode memo section 6).

For each mode: which branch ran (the evidence file of that mode and no other), which
outputs appeared, what was intentionally disabled (no tier outside short-read wgs, no
targeted assembly outside it, no MAGs), the depletion preset the host pass really used,
and the terminal state the QC summary gives the sample.

    assert_mode_outputs.py --mode wgs|wgs_assembly|wes|rna|metagenome|wgs_ont|rna_ont \
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
    "wgs_assembly": ("read_evidence.tsv", "extensive_island",  "completed_positive"),
    "wes":        ("wes_evidence.tsv",  None,                  "completed_untiered"),
    "rna":        ("rna_evidence.tsv",  None,                  "completed_untiered"),
    "metagenome": ("read_evidence.tsv", "untiered",            "completed_untiered"),
    "wgs_ont":    ("read_evidence.tsv", "untiered",            "completed_untiered"),
    "rna_ont":    ("rna_evidence.tsv",  None,                  "completed_untiered"),
}
EVIDENCE_FILES = ("read_evidence.tsv", "wes_evidence.tsv", "rna_evidence.tsv")
READS_FIELD = {"read_evidence.tsv": "pks_reads", "wes_evidence.tsv": "clb_wes_reads",
               "rna_evidence.tsv": "clb_rna_reads"}


def metrics(path):
    with path.open(newline="") as fh:
        return {r["metric"]: r["value"] for r in csv.DictReader(fh, delimiter="\t")}


def check_targeted_recruitment(sdir, work_dir, fx):
    """Ludmil, 25 Sep A1/A2: recruitment ran with the new rule and reports its classes.

    The fixture's extensive-island pairs are all island-internal; its flank pairs are
    flank-only; each boundary pair's best mate crosses a boundary; the chimeric pairs
    (60 bp island + 90 bp random per mate) are rejected on query coverage."""
    r = fx["recruitment"]
    internal, boundary, flank = (r["island_internal_fragments"], r["boundary_spanning_fragments"],
                                 r["flank_only_fragments"])
    supporting = internal + boundary
    stats_path = sdir / "contigs/recruitment/pks_recruitment.tsv"
    check("pks_recruitment.tsv exists", stats_path.exists())
    stats = metrics(stats_path) if stats_path.exists() else {}
    for name, want in (("fragments_island_internal", internal), ("fragments_boundary_spanning", boundary),
                       ("fragments_flank_only", flank), ("fragments_rejected", r["rejected_fragments"]),
                       ("island_supporting_fragments", supporting),
                       ("recruited_fragment_ids", supporting + flank),
                       ("paired_fragments", supporting + flank),
                       ("reads_rejected", r["rejected_reads"]),
                       ("reads_rejected_low_query_coverage", r["rejected_reads"]),
                       ("min_mapq", "20"), ("min_query_coverage", "0.8"),
                       ("boundary_tolerance_bp", "10"), ("island_start_in_reference", "5001"),
                       ("island_end_in_reference", "55768"), ("reference_length", "60768")):
        check(f"recruitment {name} == {want}", stats.get(name) == str(want), f"got {stats.get(name)!r}")
    recruit_cmds = [p.read_text(errors="replace") for p in work_dir.glob("*/*/.command.sh")
                    if "recover_recruited_mates.py" in p.read_text(errors="replace")
                    and "--sam" in p.read_text(errors="replace")]
    check("found the rendered recruitment command", len(recruit_cmds) == 1, f"{len(recruit_cmds)} found")
    check("recruitment Bowtie2 has no -k 1", bool(recruit_cmds) and all(" -k 1 " not in c for c in recruit_cmds))
    final_path = sdir / "contigs/final_evidence/final_pks_evidence.tsv"
    check("final_pks_evidence.tsv exists", final_path.exists())
    final = one_row(final_path) if final_path.exists() else {}
    for name, want in (("island_supporting_fragments", supporting),
                       ("boundary_spanning_fragments", boundary), ("flank_only_fragments", flank)):
        check(f"final evidence {name} == {want}", final.get(name) == str(want), f"got {final.get(name)!r}")
    check("assembly_input_diagnostic reads island support",
          final.get("assembly_input_diagnostic") in ("island_fragments_no_contig", "island_contigs_recovered"),
          f"got {final.get('assembly_input_diagnostic')!r}")


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

    reads = {"wgs": fx["short_positive_reads"], "wgs_assembly": None, "wes": fx["short_positive_reads"],
             "metagenome": fx["short_positive_reads"], "rna": fx["rna_short_pks_reads"],
             "wgs_ont": fx["long_pks_reads"], "rna_ont": fx["long_pks_reads"]}[mode]
    got = ev.get(READS_FIELD[evidence_file])
    if reads is not None:   # wgs_assembly: the tier is the point, not the exact count
        check(f"{READS_FIELD[evidence_file]} == {reads}", got == str(reads), f"got {got!r}")

    # ---- disabled on purpose ----
    # wgs: the positive fixture is multi_gene, below --tumor_contig_tiers, so it is not
    # assembled either; every other mode never assembles. wgs_assembly widens the tiers.
    if mode == "wgs_assembly":
        check_targeted_recruitment(sdir, a.work_dir, fx)
    else:
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
        "wgs_assembly": (fx["recruitment"]["input_reads"], fx["recruitment"]["input_reads"]),
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
    human_preset = {"wgs": "-ax sr ", "wgs_assembly": "-ax sr ", "wes": "-ax sr ", "metagenome": "-ax sr ",
                    "rna": "-ax splice:sr ", "wgs_ont": "-ax map-ont ", "rna_ont": "-ax splice "}[mode]
    check(f"human pass used minimap2 {human_preset.strip()}", any(human_preset in c for c in commands))
    if mode in ("rna", "rna_ont"):
        check("human pass got the junction BED (-j)", any("-j " in c for c in commands))
    if mode.endswith("_ont"):
        check("no -ax sr call for long reads", not any("-ax sr " in c for c in commands))
        check("no --frag for long reads", not any("--frag" in c for c in commands))

    # ---- master summary carries the mode ----
    master = rows_by(a.results / "cohort/pks.master_summary.tsv", "sample").get(sample, {})
    if mode == "wgs_assembly":
        r = fx["recruitment"]
        fragments = str(r["island_internal_fragments"] + r["boundary_spanning_fragments"])
        check(f"master summary island_supporting_fragments == {fragments}",
              master.get("island_supporting_fragments") == fragments,
              f"got {master.get('island_supporting_fragments')!r}")
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
