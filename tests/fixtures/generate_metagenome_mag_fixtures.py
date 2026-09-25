#!/usr/bin/env python3
"""Build a synthetic metagenome fixture for tests/check_metagenome_mag_execution.sh.

This is the real-execution companion to a defect found and fixed the same session as
this file was added: build_mag_summary.py's per-sample summary table
(genomes/pks_mag_summary.tsv) requires a real CheckM2 report and a real GTDB-Tk summary
to exist for a sample before it will report ANYTHING for that sample -- including U1's
own "unbinned pool" row (Modules/pks_mag.nf, metabat2Bin's --unbinned output), which
exists specifically for samples where MetaBAT2 recovered zero real bins. Since
checkm2Predict/gtdbtkClassify only ever run on real bins (metabat2Bin.out.bins,
optional: true), a sample with zero real bins never gets either report, and
pks_mag_summary.tsv's strict join used to silently drop that sample entirely -- the
exact case U1 exists for ("common when abundance is low or binning is underpowered")
never actually surfaced anywhere in the final table. See
Modules/pks_mag.nf's stubBinlessMagQuality for the fix.

This fixture is built to reliably hit that exact scenario for real: paired reads tiled
at real, overlapping depth across the complete pks/clb island (2,193,827-2,244,594, all
19 genes) of the repo's own shipped GCF_000025745.1 (IHE3034) reference -- the same
technique tests/fixtures/generate_e2e_fixtures.py already uses, just at real assembly
depth rather than tests/check_e2e_execution.sh's tiny read-tier fixture. The assembled
contig(s) total well under MetaBAT2's 200 kb --minClsSize floor (confirmed in real
testing, 2026-09-24: an assembly this size gets --unbinned's whole-assembly pooling, not
binned into a genome), so real bins = 0 and the unbinned pool is the ONLY signal --
exactly what this fixture needs to exercise the fix.

Read names are given WITHOUT a /1 or /2 suffix, for the same reason as every other
fixture generator here: filterReads (Modules/filter_reads.nf) adds it itself.

All sequence content is public reference data already checked into this repository
(GCF_000025745.1, RefSeq, non-human). No sequencing data, patient data or any file
outside this repository is read.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from generate_e2e_fixtures import REFERENCE_FASTA, CONTIG, READ_LEN, revcomp, load_contig, write_fastq_gz  # noqa: E402

# main.nf: params.pks_start_1based / params.pks_end_1based. 50,768 bp, all 19 clb genes.
ISLAND_START_1BASED = 2193827
ISLAND_END_1BASED = 2244594

FRAGMENT_LEN = 2 * READ_LEN  # 300 bp: mate1 = first 150 bp, mate2 = revcomp(last 150 bp)
TILE_STEP = 50  # depth = FRAGMENT_LEN / TILE_STEP = 6x; << FRAGMENT_LEN so megahit's
                # k-mer graph (default k starts at 21) sees real read-to-read overlap,
                # not one disconnected 300 bp contig per fragment.


def tile_window(genome_seq, start_1based, end_1based, step):
    start0 = start_1based - 1
    end0 = end_1based
    pos = start0
    i = 0
    while pos + FRAGMENT_LEN <= end0:
        fragment = genome_seq[pos:pos + FRAGMENT_LEN]
        if "N" in fragment.upper():
            raise SystemExit(f"metagenome MAG fixture window contains an N base at {pos}; adjust coordinates")
        mate1 = fragment[:READ_LEN]
        mate2 = revcomp(fragment[READ_LEN:])
        yield (f"sim_meta_mag_{i}", mate1, mate2)
        pos += step
        i += 1


def build_metagenome_mag_positive(genome_seq, out_dir, step=TILE_STEP):
    r1, r2 = [], []
    for name, mate1, mate2 in tile_window(genome_seq, ISLAND_START_1BASED, ISLAND_END_1BASED, step):
        r1.append((name, mate1))
        r2.append((name, mate2))
    write_fastq_gz(out_dir / "metagenome_mag_R1.fastq.gz", r1)
    write_fastq_gz(out_dir / "metagenome_mag_R2.fastq.gz", r2)
    return len(r1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--step", type=int, default=TILE_STEP)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    genome_seq = load_contig(REFERENCE_FASTA, CONTIG)
    n_pairs = build_metagenome_mag_positive(genome_seq, args.out_dir, step=args.step)

    print(f"metagenome_mag pairs: {n_pairs} "
          f"({ISLAND_END_1BASED - ISLAND_START_1BASED + 1} bp window, all 19 clb genes, "
          f"well under MetaBAT2's 200 kb --minClsSize floor)")
    print(f"wrote fixtures to {args.out_dir}")


if __name__ == "__main__":
    main()
