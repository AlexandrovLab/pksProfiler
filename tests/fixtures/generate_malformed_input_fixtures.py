#!/usr/bin/env python3
"""Build a real, valid BAM and two deliberately malformed derivatives of it, for
tests/check_malformed_input_execution.sh (F21's "malformed input" scenario).

Ludmil's audit, u-3's wish list (closed piecemeal across F21): "zero, positive and
malformed inputs." Every other e2e scenario in this suite proves the pipeline computes
the right answer on good input; this one proves it fails visibly and informatively on
bad input instead of hanging, crashing without explanation, or silently producing wrong
output.

Two malformed derivatives, chosen to exercise the two distinct guard clauses
Modules/extract_reads.nf actually has, not a single generic "broken file":

  * `truncated.bam` -- valid.bam cut off partway through (a real, if boring, way for a
    file to arrive damaged: an interrupted transfer, a killed job, a full disk). Its
    leading BGZF magic bytes are untouched, so `htsfile` still reports it as BAM and
    extractReads' own format-detection branch passes it through to `samtools quickcheck
    -v` -- which is exactly the check that must catch it. Empirically confirmed
    (samtools 1.21, see this repo's own conda_envs/samtools_env.yml) to fail with
    "was missing EOF block when one should be present." and a nonzero exit, which
    `set -euo pipefail` turns into a task failure.
  * `bad_magic.bam` -- valid.bam with its first 4 bytes overwritten, so it is not a BGZF
    stream at all. `htsfile` reports "unknown data", which matches neither the `*BAM*`
    nor `*CRAM*` case in extractReads.nf's own detection branch, so the pipeline's own
    explicit message fires: "input is neither BAM nor CRAM according to htsfile: ...".

Every read is independent pseudorandom sequence (not the shipped GCF_000025745.1
reference or anything derived from real sequencing data) -- correctness of the *content*
is not what this scenario is about; only file-level integrity is.
"""
import argparse
import random
import subprocess
from pathlib import Path

READ_LEN = 100
N_READS = 200
QUAL = "I" * READ_LEN


def random_seq(rng, length):
    return "".join(rng.choices("ACGT", k=length))


def write_sam(path, seed):
    rng = random.Random(seed)
    lines = [
        "@HD\tVN:1.6\tSO:unsorted",
        "@SQ\tSN:malformed_fixture_chr1\tLN:4000",
    ]
    for i in range(N_READS):
        seq = random_seq(rng, READ_LEN)
        # flag 4 = unmapped, the exact category extractReads.nf's `samtools fastq -f 4`
        # keeps -- a valid, non-truncated version of this file is real, extractable
        # input, so what fails below is the corruption, not the fixture's design.
        lines.append(f"read{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{QUAL}")
    path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--samtools-bin", required=True)
    parser.add_argument("--seed", type=int, default=20260924)
    parser.add_argument(
        "--truncate-fraction", type=float, default=0.4,
        help="fraction of valid.bam's bytes to keep in truncated.bam")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    sam_path = args.out_dir / "valid.sam"
    valid_bam = args.out_dir / "valid.bam"
    truncated_bam = args.out_dir / "truncated.bam"
    bad_magic_bam = args.out_dir / "bad_magic.bam"

    write_sam(sam_path, args.seed)

    subprocess.run(
        [args.samtools_bin, "view", "-bS", "-o", str(valid_bam), str(sam_path)],
        check=True,
    )

    valid_bytes = valid_bam.read_bytes()
    if len(valid_bytes) < 100:
        raise SystemExit(f"valid.bam is suspiciously small ({len(valid_bytes)} bytes); "
                          "the fixture SAM may not have converted correctly")

    keep = int(len(valid_bytes) * args.truncate_fraction)
    truncated_bam.write_bytes(valid_bytes[:keep])

    corrupted = bytearray(valid_bytes)
    corrupted[0:4] = b"XXXX"
    bad_magic_bam.write_bytes(bytes(corrupted))

    print(f"valid.bam: {len(valid_bytes)} bytes, {N_READS} unmapped reads")
    print(f"truncated.bam: {keep} bytes ({args.truncate_fraction:.0%} of valid.bam)")
    print(f"bad_magic.bam: {len(corrupted)} bytes, leading magic overwritten")
    print(f"wrote fixtures to {args.out_dir}")


if __name__ == "__main__":
    main()
