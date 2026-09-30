#!/usr/bin/env python3
"""WES input-alignment QC fragment: duplicate rate and capture/off-target counts.

Ludmil, 26 Sep mode memo section 3: "Report total, unmapped, host-depleted and post-fastp
reads, duplicate rate, and capture kit/off-target fraction where available. WES duplicate
behavior must be checked separately from WGS."

Modules/mode_qc.nf (wesInputQC) decodes the input alignment once, keeping primary mapped
records only (samtools view -F 2308), and tees that stream into

  * samtools flagstat -O tsv          -> --flagstat: records and 0x400-flagged records
  * samtools view -c -L <capture BED> -> --on-target: records overlapping the capture BED
                                         (only with --wes_capture_bed)

Duplicates are the submitter's own marking (flag 0x400): pksProfiler does not re-mark the
full alignment. An input nobody ran a duplicate marker on reports 0 duplicates, which is
why the fraction is read beside the header's @PG chain (provenance/input_header.tsv).

Written metrics (Sample/Metric/Value, integers or a status, as build_qc_summary.py wants):

    wes_input_primary_mapped      primary mapped records in the input alignment
    wes_input_primary_duplicates  of those, flagged 0x400
    wes_capture_on_target         of those, overlapping the capture BED (with a BED only)
    wes_capture_status            capture_bed_provided | capture_bed_not_provided
"""
import argparse
from pathlib import Path

CAPTURE_STATUSES = ("capture_bed_provided", "capture_bed_not_provided")
# Every metric this fragment can carry (tests/test_qc_denominator.py reads it).
METRICS = ("wes_input_primary_mapped", "wes_input_primary_duplicates",
           "wes_capture_on_target", "wes_capture_status")


def flagstat_counts(path):
    """{category: passed + failed} from `samtools flagstat -O tsv` output.

    Each line is <QC-passed>\t<QC-failed>\t<category>. The input is already primary
    mapped records only, so "total" is primary mapped records and "duplicates" is
    primary mapped duplicates.
    """
    counts = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        passed, failed, category = fields[0], fields[1], fields[2]
        if not (passed.isdigit() and failed.isdigit()):
            continue
        counts[category.strip()] = int(passed) + int(failed)
    return counts


def input_metrics(flagstat, on_target=None):
    """The fragment's rows, in order. `on_target` is None without a capture BED."""
    counts = flagstat_counts(flagstat)
    total = next((v for k, v in counts.items() if k.startswith("total")), None)
    duplicates = counts.get("duplicates")
    if total is None or duplicates is None:
        raise SystemExit(f"{flagstat}: no total/duplicates lines; is this samtools flagstat -O tsv output?")
    if duplicates > total:
        raise SystemExit(f"{flagstat}: {duplicates} duplicates exceed {total} records")
    rows = [("wes_input_primary_mapped", total), ("wes_input_primary_duplicates", duplicates)]
    if on_target is not None:
        if on_target > total:
            raise SystemExit(f"on-target records ({on_target}) exceed primary mapped records ({total})")
        rows.append(("wes_capture_on_target", on_target))
        rows.append(("wes_capture_status", "capture_bed_provided"))
    else:
        rows.append(("wes_capture_status", "capture_bed_not_provided"))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--flagstat", required=True, help="samtools flagstat -O tsv of the -F 2308 stream")
    ap.add_argument("--on-target", default=None,
                    help="file holding the samtools view -c -L count; omit without a capture BED")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    on_target = None
    if a.on_target:
        text = Path(a.on_target).read_text().strip()
        if not text.isdigit():
            raise SystemExit(f"{a.on_target}: expected one integer, got {text!r}")
        on_target = int(text)
    with open(a.output, "w") as out:
        out.write("Sample\tMetric\tValue\n")
        for metric, value in input_metrics(a.flagstat, on_target):
            out.write(f"{a.sample}\t{metric}\t{value}\n")


if __name__ == "__main__":
    main()
