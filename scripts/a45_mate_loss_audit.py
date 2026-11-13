#!/usr/bin/env python3
"""Tier A 4+5: how many mates does extraction actually lose, and how human are they?

From Ludmil's v3 code review (items 4 and 5, one defect in two branches): extraction
keeps a read only if IT qualifies -- unmapped (`-f 4`), or mapped below the AS rescue
threshold when rescue is on. Its mate is never considered, so a fragment with one
qualifying read and one better-scoring human-mapped mate is truncated at the door, and
`recover_recruited_mates.py` cannot restore sequence that never entered the extracted
stream. The prescribed fix (QNAME recovery routed through host depletion) costs an
extra decode pass over every alignment, so before building it: measure.

TWO QUESTIONS, ONE SCRIPT

  1. The census. Of the fragments extraction keeps anything from, how many lose their
     other mate? (both_mates_extracted vs orphan_fragments)

  2. The character of what is lost. For each orphan, the mate that stayed behind is by
     definition a primary mapped record. Its alignment score says what we gave up:
     a mate at AS >= 100 is a confident human match that depletion would re-remove
     anyway -- recovering it buys nothing. A mate at AS 30-49 sits one step above the
     rescue line, in exactly the gray zone that motivated low-AS rescue. A mate with
     MAPQ 0 mapped somewhere ambiguous. The histogram of lost-mate AS and the MAPQ-0
     fraction are the numbers that decide whether items 4-5 need a feature or a
     documented, quantified bound.

HOW IT WORKS

Everything comes from the original BAM/CRAM -- the extracted FASTQ is not needed (it
is only retained under --save_intermediates). Pass 1 streams the alignment with the
extraction's own record filter (-F 2304: secondaries and supplementaries out, QCFAIL
and duplicates through) and replays the extraction predicate on each primary record,
writing the names of qualifying reads. Pass 2 fetches every primary record bearing
those names (`samtools view -N`), so each candidate fragment's two sides can be
compared; grouping is done with a disk-based sort, so memory stays flat no matter the
cohort. Choose --low-score-mode to match how the run was actually extracted (C1 showed
one mode per cohort; the provenance report says which).

PRIVACY

This reads patient-derived alignments, so it is written for Ammal to run, not for
Claude. The default output is aggregate only: counts, rates and histogram bins, no
read name, no sequence, no per-read row. --per-fragment exists for local debugging and
its output stays on the cluster: it carries read names and belongs in no report, deck,
commit or email.

    python3 scripts/a45_mate_loss_audit.py \
        --alignment /path/sample.bam \
        --low-score-mode on --min-as 30 \
        --out validation/a45_mate_loss.SAMPLE.tsv

    # CRAM needs the same reference the run used:
    #   --cram-reference /path/Homo_sapiens_assembly38.alt.masked.fasta
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# The extraction's record filter, exactly as Modules/extract_reads.nf applies it:
# secondaries and supplementaries excluded, QCFAIL and duplicates passed (item 10
# documents that choice; this audit must mirror it, not improve on it).
EXCLUDE_FLAGS = "2304"

# Lost-mate AS bins. The 30-49 band is the gray zone: one step above the rescue
# threshold bwa itself uses as its "too weak to report" default.
AS_BINS = [
    ("AS_under_30", lambda score: score < 30),
    ("AS_30_49", lambda score: 30 <= score < 50),
    ("AS_50_99", lambda score: 50 <= score < 100),
    ("AS_100_plus", lambda score: score >= 100),
]

FIELDS = [
    "sample", "low_score_mode", "min_as",
    "primary_records_seen", "qualifying_reads",
    "fragments_touched", "both_mates_extracted", "orphan_fragments",
    "orphan_fraction_of_fragments", "single_end_or_mateless",
    "lost_mate_no_AS", "lost_mate_AS_under_30", "lost_mate_AS_30_49",
    "lost_mate_AS_50_99", "lost_mate_AS_100_plus",
    "lost_mate_mapq0", "lost_mate_mapq0_fraction",
    "lost_mate_gray_zone_fraction",
]


def run(cmd, **kwargs):
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, **kwargs)


def tag_value(fields, tag):
    prefix = tag + ":"
    for field in fields[11:]:
        if field.startswith(prefix):
            value = field.split(":", 2)[2]
            try:
                return int(float(value))
            except ValueError:
                return None
    return None


def qualifies(flag, score, low_score_on, min_as):
    """The extraction predicate from Modules/extract_reads.nf.

    Unmapped always qualifies. With rescue on, a mapped record with AS below the
    threshold qualifies too; a mapped record with NO AS tag is dropped by the rescue
    awk (counted, not extracted), so it does not qualify either way.
    """
    if flag & 4:
        return True
    if low_score_on and score is not None and score < min_as:
        return True
    return False


def segment_of(flag):
    """1 for a first-segment record, 2 for last, 0 for anything else."""
    first, last = bool(flag & 64), bool(flag & 128)
    if first and not last:
        return 1
    if last and not first:
        return 2
    return 0


def view_cmd(samtools, alignment, cram_reference, names_file=None):
    cmd = [samtools, "view", "-F", EXCLUDE_FLAGS]
    if cram_reference:
        cmd += ["-T", cram_reference]
    if names_file:
        cmd += ["-N", names_file]
    cmd.append(alignment)
    return cmd


def pass_one(args, names_path):
    """Stream every primary record; write the names of qualifying reads."""
    seen = 0
    qualifying = 0
    proc = run(view_cmd(args.samtools, args.alignment, args.cram_reference))
    with open(names_path, "w", encoding="utf-8") as names:
        last_name = None
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            seen += 1
            flag = int(fields[1])
            if qualifies(flag, tag_value(fields, "AS"), args.low_score_on, args.min_as):
                qualifying += 1
                # Both mates often qualify back to back in name order; this cheap
                # dedup keeps the -N file smaller. samtools -N dedups anyway.
                if fields[0] != last_name:
                    names.write(fields[0] + "\n")
                    last_name = fields[0]
    if proc.wait() != 0:
        sys.exit(f"samtools view failed on {args.alignment} (pass 1)")
    return seen, qualifying


def pass_two(args, names_path, workdir):
    """Fetch every primary record of the candidate names; group per fragment on disk."""
    rows_path = Path(workdir) / "candidates.tsv"
    proc = run(view_cmd(args.samtools, args.alignment, args.cram_reference,
                        names_file=str(names_path)))
    with rows_path.open("w", encoding="utf-8") as rows:
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            score = tag_value(fields, "AS")
            keep = qualifies(flag, score, args.low_score_on, args.min_as)
            rows.write("\t".join((
                fields[0], str(segment_of(flag)), "1" if keep else "0",
                "NA" if score is None else str(score), fields[4],
            )) + "\n")
    if proc.wait() != 0:
        sys.exit(f"samtools view -N failed on {args.alignment} (pass 2)")

    sorted_path = Path(workdir) / "candidates.sorted.tsv"
    env = dict(os.environ, LC_ALL="C")
    with sorted_path.open("w", encoding="utf-8") as out:
        done = subprocess.run(["sort", "-k1,1", "-T", workdir, str(rows_path)],
                              stdout=out, env=env)
    if done.returncode != 0:
        sys.exit("sort failed grouping candidate fragments")
    rows_path.unlink()
    return sorted_path


def classify(sorted_path, per_fragment_path=None):
    """Walk the name-grouped records; count fragment fates and bin lost mates."""
    counts = {
        "fragments_touched": 0, "both_mates_extracted": 0, "orphan_fragments": 0,
        "single_end_or_mateless": 0, "lost_mate_no_AS": 0, "lost_mate_mapq0": 0,
    }
    for name, _ in AS_BINS:
        counts[f"lost_mate_{name}"] = 0

    per_fragment = open(per_fragment_path, "w", encoding="utf-8") \
        if per_fragment_path else None
    if per_fragment:
        per_fragment.write("# PER-FRAGMENT, carries read names. Keep on the cluster.\n")
        per_fragment.write("fragment\tfate\tlost_mate_AS\tlost_mate_MAPQ\n")

    def settle(name, records):
        counts["fragments_touched"] += 1
        kept = [r for r in records if r["keep"]]
        lost = [r for r in records if not r["keep"]]
        if not lost and len(records) >= 2:
            counts["both_mates_extracted"] += 1
            fate, lost_as, lost_mapq = "both_extracted", "", ""
        elif len(records) == 1:
            counts["single_end_or_mateless"] += 1
            fate, lost_as, lost_mapq = "single_end", "", ""
        elif kept and lost:
            counts["orphan_fragments"] += 1
            mate = lost[0]
            fate, lost_as, lost_mapq = "orphan", mate["as"], mate["mapq"]
            if mate["as"] == "NA":
                counts["lost_mate_no_AS"] += 1
            else:
                score = int(mate["as"])
                for bin_name, in_bin in AS_BINS:
                    if in_bin(score):
                        counts[f"lost_mate_{bin_name}"] += 1
                        break
            if mate["mapq"] == "0":
                counts["lost_mate_mapq0"] += 1
        else:
            # Every record qualifies but none is a clean mate pair (e.g. two
            # same-segment records). Count it as touched and move on.
            fate, lost_as, lost_mapq = "irregular", "", ""
        if per_fragment:
            per_fragment.write(f"{name}\t{fate}\t{lost_as}\t{lost_mapq}\n")

    current_name, records = None, []
    with sorted_path.open(encoding="utf-8") as handle:
        for line in handle:
            name, segment, keep, score, mapq = line.rstrip("\n").split("\t")
            if name != current_name:
                if current_name is not None:
                    settle(current_name, records)
                current_name, records = name, []
            records.append({"segment": segment, "keep": keep == "1",
                            "as": score, "mapq": mapq})
    if current_name is not None:
        settle(current_name, records)
    if per_fragment:
        per_fragment.close()
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", required=True, help="the original BAM/CRAM")
    parser.add_argument("--sample", default=None,
                        help="label for the output row (default: file stem)")
    parser.add_argument("--low-score-mode", choices=["on", "off"], required=True,
                        help="how the run extracted: 'on' = unmapped plus AS<threshold "
                             "rescue, 'off' = unmapped only. Must match the run "
                             "(by_sample/<s>/provenance/input_header.tsv says which).")
    parser.add_argument("--min-as", type=int, default=30,
                        help="the run's --extract_min_as (default 30)")
    parser.add_argument("--cram-reference", default=None,
                        help="reference FASTA for CRAM decode; must be the one the "
                             "run used")
    parser.add_argument("--samtools", default=shutil.which("samtools") or "samtools")
    parser.add_argument("--tmpdir", default=None,
                        help="where the name list and sort spill live; needs room for "
                             "one line per qualifying read")
    parser.add_argument("--per-fragment", default=None,
                        help="also write one row per fragment HERE. Carries read "
                             "names: keep it on the cluster.")
    parser.add_argument("--out", required=True, help="aggregate TSV (append-safe)")
    args = parser.parse_args(argv)
    args.low_score_on = args.low_score_mode == "on"

    sample = args.sample or Path(args.alignment).stem

    with tempfile.TemporaryDirectory(dir=args.tmpdir) as workdir:
        names_path = Path(workdir) / "qualifying.names"
        seen, qualifying = pass_one(args, names_path)
        if qualifying == 0:
            sys.exit(f"no read qualifies for extraction in {args.alignment}; "
                     "check --low-score-mode / --min-as against the run")
        sorted_path = pass_two(args, names_path, workdir)
        counts = classify(sorted_path, args.per_fragment)

    fragments = counts["fragments_touched"]
    orphans = counts["orphan_fragments"]
    gray = counts["lost_mate_no_AS"] + counts["lost_mate_AS_under_30"] \
        + counts["lost_mate_AS_30_49"]
    row = {
        "sample": sample,
        "low_score_mode": args.low_score_mode,
        "min_as": args.min_as,
        "primary_records_seen": seen,
        "qualifying_reads": qualifying,
        "fragments_touched": fragments,
        "both_mates_extracted": counts["both_mates_extracted"],
        "orphan_fragments": orphans,
        "orphan_fraction_of_fragments":
            f"{orphans / fragments:.6f}" if fragments else "",
        "single_end_or_mateless": counts["single_end_or_mateless"],
        "lost_mate_no_AS": counts["lost_mate_no_AS"],
        "lost_mate_AS_under_30": counts["lost_mate_AS_under_30"],
        "lost_mate_AS_30_49": counts["lost_mate_AS_30_49"],
        "lost_mate_AS_50_99": counts["lost_mate_AS_50_99"],
        "lost_mate_AS_100_plus": counts["lost_mate_AS_100_plus"],
        "lost_mate_mapq0": counts["lost_mate_mapq0"],
        "lost_mate_mapq0_fraction":
            f"{counts['lost_mate_mapq0'] / orphans:.6f}" if orphans else "",
        "lost_mate_gray_zone_fraction":
            f"{gray / orphans:.6f}" if orphans else "",
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not out_path.exists() or out_path.stat().st_size == 0
    with out_path.open("a", newline="", encoding="utf-8") as handle:
        if write_header:
            handle.write("# Tier A 4+5 mate-loss audit: fragments extraction touches,\n")
            handle.write("# how many lose a mate, and the lost mates' AS / MAPQ bins.\n")
            handle.write("# Aggregate only -- no read name, no sequence. One row per\n")
            handle.write("# sample; re-runs append.\n")
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerow(row)

    print(f"  sample {sample}: {fragments} fragments touched by extraction")
    print(f"    both mates extracted : {counts['both_mates_extracted']}")
    print(f"    orphaned (mate lost) : {orphans}"
          + (f"  ({row['orphan_fraction_of_fragments']} of touched fragments)"
             if fragments else ""))
    print(f"    single-end/mateless  : {counts['single_end_or_mateless']}")
    if orphans:
        print(f"    lost mates: no AS {counts['lost_mate_no_AS']}, "
              f"AS<30 {counts['lost_mate_AS_under_30']}, "
              f"30-49 {counts['lost_mate_AS_30_49']}, "
              f"50-99 {counts['lost_mate_AS_50_99']}, "
              f">=100 {counts['lost_mate_AS_100_plus']}; "
              f"MAPQ0 {counts['lost_mate_mapq0']}")
        print(f"    gray zone (no AS or AS<50): {row['lost_mate_gray_zone_fraction']}"
              " of lost mates")
        print()
        print("  Reading it: lost mates at AS>=100 are confident human matches that")
        print("  depletion would re-remove -- recovering them buys nothing. The case")
        print("  for building mate recovery (items 4-5) lives entirely in the gray")
        print("  zone and the MAPQ0 fraction.")
    print(f"\n  appended to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
