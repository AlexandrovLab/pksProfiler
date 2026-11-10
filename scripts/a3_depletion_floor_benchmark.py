#!/usr/bin/env python3
"""Tier A 3: what does the host-depletion removal floor actually discard?

From Ludmil's v3 code review:

    Host depletion inherits an unbenchmarked minimap2 removal floor. The pipeline keeps
    only minimap2-unmapped reads (`samtools fastq -f 4`). Therefore any read that
    minimap2 reports as mapped under the selected short-read preset is removed, without
    a pksProfiler-specific benchmark of MAPQ, identity, aligned fraction, or alignment
    score. This is conceptually analogous to the upstream weak-human-hit problem that
    motivated low-AS rescue: a weak human hit can hide a bacterial read. It applies to
    every cohort and remains unquantified.

    Fix: Benchmark true human reads, true pks reads, related bacterial sequence, and
    low-complexity reads against MAPQ, alignment score, identity, and aligned fraction.
    Set a defensible removal rule only after that benchmark.

TWO THINGS THIS MEASURES
------------------------
1. THE ASYMMETRY. For aligned input the pipeline rescues reads the UPSTREAM aligner
   mapped weakly (AS below the low-score threshold), on the reasoning that a weak human
   hit can hide a bacterial read. Its own depletion passes apply no such rule: `-f 4`
   discards anything minimap2 aligns at all, however poorly.

2. THE COMPOUNDING, which the review understates. The floor is not applied once. In
   Modules/map_reads.nf it is applied by the human pass, then again by the T2T+phiX
   pass, and then once per pangenome index in a loop. Each pass independently re-exposes
   every surviving read to the same unbenchmarked rule, so a small per-pass false-removal
   rate compounds across passes. This script therefore models the stages IN SERIES:
   stage k sees only the reads that survived stages 1..k-1, exactly as the pipeline does.

For a read whose origin is known NOT to be human, any alignment at all is a read the
pipeline throws away: a false removal, and a sensitivity loss that is invisible
downstream because the read never reaches profiling.

PRIVACY
-------
Uses only public reference sequence and reads simulated from it. It opens no sample
file, and reports only distributions and counts.

USAGE
-----
    python3 scripts/a3_depletion_floor_benchmark.py \
        --index $DB_ROOT/human-GRC-db.mmi \
        --index $DB_ROOT/human-GCA-phix-db.mmi \
        --pks-fasta indices/pks_island/pks_island.fna \
        --bacteria-fasta indices/GCF_000025745.1/..._genomic.fna \
        --human-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        --out validation/a3_depletion_floor.tsv

Pass --index once per depletion stage, in pipeline order. Without --human-fasta the
true-human arm is skipped and the run measures false removal only, which is the
sensitivity-critical half but not a complete removal rule.
"""
from __future__ import annotations

import argparse
import csv
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_PRESET = "sr"      # HostDepletion.DNA_PRESET['short']
READ_LEN = 150
SUBS = {"A": "C", "C": "A", "G": "T", "T": "G"}
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

# Every read is named "<class>|<n>" so that all classes go through ONE minimap2
# invocation per stage and are split apart afterwards. The human index is ~9 GB, so
# loading it once per stage instead of once per class is the whole cost of the run.
NAME_SEP = "|"

METRICS = ("mapq", "AS", "identity", "aligned_fraction")


# ---------------------------------------------------------------- read simulation

def load_sequence(path, limit_bp=None):
    """Concatenated sequence of a public reference FASTA, upper-cased, N dropped."""
    chunks, total = [], 0
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            piece = line.strip().upper().replace("N", "")
            if not piece:
                continue
            chunks.append(piece)
            total += len(piece)
            if limit_bp and total >= limit_bp:
                break
    return "".join(chunks)


def simulate(sequence, count, divergence, rng, label):
    """`count` reads of READ_LEN, each base substituted at rate `divergence`."""
    reads = []
    if len(sequence) <= READ_LEN:
        return reads
    for index in range(count):
        start = rng.randrange(0, len(sequence) - READ_LEN)
        bases = list(sequence[start:start + READ_LEN])
        for position, base in enumerate(bases):
            if base in SUBS and rng.random() < divergence:
                bases[position] = SUBS[base]
        reads.append((f"{label}{NAME_SEP}{index}", "".join(bases)))
    return reads


def low_complexity(count, rng):
    """Homopolymer and short-repeat reads: the classic spurious-alignment class."""
    reads = []
    for index in range(count):
        unit = "".join(rng.choice("ACGT") for _ in range((index % 3) + 1))
        reads.append((f"low_complexity{NAME_SEP}{index}", (unit * READ_LEN)[:READ_LEN]))
    return reads


def build_classes(args, rng):
    classes = []
    pks = load_sequence(args.pks_fasta)
    if not pks:
        sys.exit(f"no sequence read from {args.pks_fasta}")
    # A carrier in a real cohort is not identical to the reference island, so the
    # divergent arms are the ones that matter for sensitivity.
    for divergence, name in ((0.00, "pks_identical"),
                             (0.02, "pks_2pct_divergent"),
                             (0.05, "pks_5pct_divergent")):
        classes.append((name, simulate(pks, args.reads, divergence, rng, name)))
    classes.append(("low_complexity", low_complexity(args.reads, rng)))
    if args.bacteria_fasta:
        other = load_sequence(args.bacteria_fasta)
        classes.append(("related_bacteria",
                        simulate(other, args.reads, 0.0, rng, "related_bacteria")))
    if args.human_fasta:
        # 0 means load the whole assembly, so reads are drawn genome-wide rather than
        # from the first few megabases of chr1.
        human = load_sequence(args.human_fasta, limit_bp=args.human_limit_bp or None)
        classes.append(("true_human",
                        simulate(human, args.reads, 0.0, rng, "true_human")))
    return [(label, reads) for label, reads in classes if reads]


# ---------------------------------------------------------------- alignment

def parse_tags(fields):
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def as_int(value):
    if value is None:
        return None
    return int(value) if value.lstrip("+-").isdigit() else None


def run_stage(reads, index, preset, threads, minimap2, workdir):
    """One depletion pass. Returns {read_name: record} for every primary record.

    `record['mapped']` is True for a read this pass would DISCARD, because the pipeline
    keeps only the unmapped (`samtools fastq -f 4 -F 2304`).
    """
    fastq = Path(workdir) / "stage.fq"
    with fastq.open("w") as handle:
        for name, sequence in reads:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")

    sam = Path(workdir) / "stage.sam"
    with sam.open("w") as out:
        done = subprocess.run(
            [minimap2, "-2", "-ax", preset, "--secondary=no",
             "-t", str(threads), str(index), str(fastq)],
            stdout=out, stderr=subprocess.DEVNULL, check=False)
    if done.returncode != 0:
        sys.exit(f"minimap2 exited {done.returncode} aligning against {index}")

    records = {}
    with sam.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 0x100 or flag & 0x800:      # -F 2304, exactly as the pipeline does
                continue
            tags = parse_tags(fields[11:])
            mismatches = as_int(tags.get("NM"))
            aligned_bp = (sum(int(n) for n, op in CIGAR.findall(fields[5]) if op in "M=X")
                          if fields[5] != "*" else 0)
            records[fields[0]] = {
                "mapped": not flag & 4,           # -f 4 keeps ONLY the unmapped
                "mapq": int(fields[4]),
                "AS": as_int(tags.get("AS")),
                "identity": (1 - mismatches / aligned_bp)
                            if (mismatches is not None and aligned_bp) else None,
                "aligned_fraction": aligned_bp / READ_LEN if aligned_bp else 0.0,
            }
    return records


def run_serial(classes, indexes, preset, threads, minimap2, workdir):
    """Apply every depletion stage IN SERIES, as Modules/map_reads.nf does.

    Stage k sees only the reads that survived stages 1..k-1. Returns a list of
    per-stage dicts {class: {"n_in", "removed": [records]}}.
    """
    surviving = [(label, list(reads)) for label, reads in classes]
    stages = []
    for index in indexes:
        live = [(label, reads) for label, reads in surviving if reads]
        if not live:
            break
        flat = [read for _, reads in live for read in reads]
        records = run_stage(flat, index, preset, threads, minimap2, workdir)

        stage = {}
        carried = []
        for label, reads in surviving:
            removed, kept = [], []
            for name, sequence in reads:
                record = records.get(name)
                if record and record["mapped"]:
                    removed.append(record)
                else:
                    kept.append((name, sequence))
            stage[label] = {"n_in": len(reads), "removed": removed}
            carried.append((label, kept))
        stages.append((Path(index).name, stage))
        surviving = carried
    return stages, dict(surviving)


# ---------------------------------------------------------------- reporting

def quantiles(records, metric):
    values = sorted(r[metric] for r in records if r[metric] is not None)
    if not values:
        return "", "", ""
    return (f"{values[0]:.4g}",
            f"{values[len(values) // 2]:.4g}",
            f"{values[-1]:.4g}")


def build_rows(classes, stages, surviving):
    rows = []
    starting = {label: len(reads) for label, reads in classes}
    for stage_name, stage in stages:
        for label, _ in classes:
            entry = stage.get(label)
            if entry is None:
                continue
            removed = entry["removed"]
            row = {
                "stage": stage_name,
                "class": label,
                "reads_in": entry["n_in"],
                "removed_here": len(removed),
                "percent_of_stage_input": (f"{100 * len(removed) / entry['n_in']:.2f}"
                                           if entry["n_in"] else ""),
            }
            for metric in METRICS:
                low, mid, high = quantiles(removed, metric)
                row[f"{metric}_min"] = low
                row[f"{metric}_median"] = mid
                row[f"{metric}_max"] = high
            rows.append(row)

    for label, _ in classes:
        start = starting[label]
        left = len(surviving.get(label, []))
        row = {
            "stage": "ALL_STAGES",
            "class": label,
            "reads_in": start,
            "removed_here": start - left,
            "percent_of_stage_input": f"{100 * (start - left) / start:.2f}" if start else "",
        }
        for metric in METRICS:
            row[f"{metric}_min"] = row[f"{metric}_median"] = row[f"{metric}_max"] = ""
        rows.append(row)
    return rows


def write_table(path, rows, indexes):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("# Tier A 3: what the -f 4 host-depletion floor removes, by class and stage.\n")
        handle.write("# Simulated reads through the pipeline's own minimap2 depletion command.\n")
        handle.write("# Stages are applied IN SERIES: each sees only the survivors of the last.\n")
        handle.write("# For any non-human class, removed_here is a FALSE removal.\n")
        handle.write(f"# stages: {', '.join(Path(i).name for i in indexes)}\n")
        writer = csv.DictWriter(handle, list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def report(rows, classes, human_arm_present, total_pipeline_stages):
    header = (f"  {'stage':<26}{'class':<22}{'in':>7}{'rm':>7}{'%':>8}"
              f"{'MAPQ med':>10}{'AS med':>9}{'ident med':>11}")
    print(header)
    print("  " + "-" * (len(header) - 2))
    last_stage = None
    for row in rows:
        if row["stage"] != last_stage:
            if last_stage is not None:
                print()
            last_stage = row["stage"]
        print(f"  {row['stage'][:25]:<26}{row['class']:<22}{row['reads_in']:>7}"
              f"{row['removed_here']:>7}{row['percent_of_stage_input']:>8}"
              f"{row['mapq_median']:>10}{row['AS_median']:>9}{row['identity_median']:>11}")
    print()

    totals = {r["class"]: r for r in rows if r["stage"] == "ALL_STAGES"}
    pks_rows = [r for name, r in totals.items() if name.startswith("pks")]
    worst = max(pks_rows, key=lambda r: int(r["removed_here"]), default=None)

    if worst and int(worst["removed_here"]):
        pct = float(worst["percent_of_stage_input"])
        print("  FALSE REMOVAL of genuine pks reads, across the benchmarked stages:")
        print(f"    {worst['removed_here']} of {worst['reads_in']} "
              f"({pct:.2f}%) in class {worst['class']}.")
        print("    Those reads never reach profiling, so the loss is invisible downstream.")
        if total_pipeline_stages and total_pipeline_stages > len(
                {r['stage'] for r in rows if r['stage'] != 'ALL_STAGES'}):
            benchmarked = len({r["stage"] for r in rows if r["stage"] != "ALL_STAGES"})
            per_stage = 1 - (1 - pct / 100) ** (1 / benchmarked) if benchmarked else 0
            projected = 1 - (1 - per_stage) ** total_pipeline_stages
            print()
            print(f"    Benchmarked {benchmarked} stage(s); a full run with the pangenome")
            print(f"    database enabled applies the SAME rule {total_pipeline_stages} times"
                  " in series.")
            print(f"    At the measured per-stage rate ({100 * per_stage:.3f}%), the projected")
            print(f"    cumulative false-removal over {total_pipeline_stages} stages is"
                  f" {100 * projected:.2f}%.")
            print("    That projection assumes stages are independent, which is an upper")
            print("    bound: reads removed by a human index are correlated across human")
            print("    references. It is an extrapolation, not a measurement.")
    elif worst:
        print("  No simulated pks read was removed by the benchmarked stage(s),")
        print("  so the current floor costs no measurable pks sensitivity here.")

    human = totals.get("true_human")
    if human_arm_present and human:
        print()
        print(f"  True human reads removed: {human['removed_here']} of "
              f"{human['reads_in']} ({human['percent_of_stage_input']}%). This is the")
        print("  floor working as intended; a removal rule is defensible where the human")
        print("  and pks metric distributions separate.")
    elif not human_arm_present:
        print()
        print("  NOTE: no --human-fasta was supplied, so there is no true-human row to")
        print("        compare against. A threshold change needs both sides.")


# ---------------------------------------------------------------- entry point

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--index", action="append", dest="indexes", default=[],
                        help="a minimap2 .mmi depletion stage; repeat in pipeline order")
    parser.add_argument("--human-index", action="append", dest="indexes",
                        help=argparse.SUPPRESS)   # back-compatible alias
    parser.add_argument("--pks-fasta", required=True,
                        help="pks island or pks-carrying genome: the true-pks source")
    parser.add_argument("--human-fasta", default=None,
                        help="human reference to simulate TRUE human reads from")
    parser.add_argument("--bacteria-fasta", default=None,
                        help="a related bacterial genome, e.g. the E. coli backbone")
    parser.add_argument("--preset", default=DEFAULT_PRESET,
                        help="minimap2 preset; must match the pipeline's (default: sr)")
    parser.add_argument("--reads", type=int, default=5000, help="reads per class")
    parser.add_argument("--human-limit-bp", type=int, default=0,
                        help="cap on human sequence loaded; 0 loads the whole assembly")
    parser.add_argument("--total-pipeline-stages", type=int, default=0,
                        help="how many times a full run applies the floor (human + T2T + "
                             "one per pangenome index); used only to project the "
                             "compounded rate beyond the stages benchmarked here")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--minimap2", default=shutil.which("minimap2") or "minimap2")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--tmpdir", default=None,
                        help="where to stage the FASTQ and SAM (default: system temp)")
    parser.add_argument("--out", default="a3_depletion_floor.tsv")
    args = parser.parse_args(argv)

    if not args.indexes:
        parser.error("at least one --index is required")

    rng = random.Random(args.seed)
    classes = build_classes(args, rng)
    if not classes:
        sys.exit("no class produced reads")

    with tempfile.TemporaryDirectory(dir=args.tmpdir) as workdir:
        stages, surviving = run_serial(classes, args.indexes, args.preset,
                                       args.threads, args.minimap2, workdir)

    rows = build_rows(classes, stages, surviving)
    write_table(args.out, rows, args.indexes)
    report(rows, classes, bool(args.human_fasta), args.total_pipeline_stages)
    print(f"\n  wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
