#!/usr/bin/env python3
"""Build the synthetic inputs for tests/check_e2e_execution.sh.

Ludmil's audit, "not yet confirmed" item u-3: every existing test either pattern-
matches .nf/.py source as text, or runs `nextflow run -preview`, which wires the
DSL2 channel graph and executes zero tasks. Neither ever runs bowtie2, minimap2,
fastp or samtools on real data and checks the output -- the exact gap that let F01
(samtools fastq silently dropping a read category) through both tiers at once.

This script writes small, real, checkable fixtures for a real (non -preview) run:

  * `positive_R{1,2}.fastq.gz` -- paired reads sliced straight out of the repo's
    already-shipped GCF_000025745.1 (IHE3034) reference at four well-separated
    clbA-clbS gene bodies (see LOCI below), so a real bowtie2 alignment run
    produces real, non-zero, multi-gene clb evidence.
  * `negative_R{1,2}.fastq.gz` -- paired reads of independent pseudorandom
    sequence, unrelated to IHE3034 or to anything else in the repo, expected to
    align nowhere and produce a `negative` read-evidence tier.
  * `synthetic_hg38.fa` / `synthetic_t2t_phix.fa` -- small pseudorandom FASTA
    standing in for --hg38_db/--t2t_phix_db. These are NOT built from human data
    (there is none in this repo) and are NOT related to the E. coli reference
    above or to the negative reads, by construction (independent RNG streams) --
    the point is that neither read set should be depleted by them, so any read
    lost between filterReads and mapReads is a real regression, not an artifact
    of a host reference that accidentally matches the payload.

Every read name is deliberately given WITHOUT a /1 or /2 suffix. filterReads
(Modules/filter_reads.nf) adds the suffix itself before handing reads to fastp --
that is exactly the step u-1's mate-pair investigation is about, and generating
suffix-free input here means the real run exercises that tagging step for real,
rather than trusting pre-suffixed input to survive by accident.

All sequence content is either public reference data already checked into this
repository (GCF_000025745.1, RefSeq, non-human) or freshly generated pseudorandom
noise. No sequencing data, patient data or any file outside this repository is
read.
"""
import argparse
import gzip
import random
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_FASTA = REPO_ROOT / "indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic.fna"
CONTIG = "NC_017628.1"

# Four clbA-clbS genes, spread across the island, each with a 300 bp window that
# sits well inside the gene body (checked against ref/annotations/IHE3034.clbA-
# clbS.gff): margin from either edge so the window cannot roll into a neighbouring
# gene and split a read's featureCounts assignment across two genes.
#
# Verified by a standalone bowtie2 run against the shipped GCF_000025745.1 index
# (--very-sensitive --no-unal | samtools view -q 40): every 150 bp mate below maps
# uniquely (MAPQ >= 40) to its own window. clbK was tried first and dropped --
# its megasynthase (NRPS/PKS) domain repeats mean a 150 bp window from it maps to
# more than one place in the genome and both mates lose the MAPQ-40 filter, which
# would silently zero out that "gene" of clb evidence. clbD replaces it.
#   gene   gff span               window start  window end
#   clbS   2193827-2194339        2193850       2194149
#   clbN   2199106-2203473        2200000       2200299
#   clbD   2230049-2230918        2230100       2230399
#   clbA   2243860-2244594        2244000       2244299
LOCI = [
    ("clbS", 2193850, 300),
    ("clbN", 2200000, 300),
    ("clbD", 2230100, 300),
    ("clbA", 2244000, 300),
]

READ_LEN = 150
QUAL = "I" * READ_LEN  # uniform Phred 40; nothing here should be quality-trimmed

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def load_contig(fasta_path, contig):
    lines = []
    capture = False
    with fasta_path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                capture = line[1:].split()[0] == contig
                continue
            if capture:
                lines.append(line.strip())
    seq = "".join(lines)
    if not seq:
        raise SystemExit(f"contig {contig} not found in {fasta_path}")
    return seq


def write_fastq_gz(path, records):
    """records: iterable of (name, sequence). Quality is uniform QUAL."""
    with gzip.open(path, "wt") as handle:
        for name, seq in records:
            handle.write(f"@{name}\n{seq}\n+\n{QUAL[:len(seq)]}\n")


def random_seq(rng, length):
    return "".join(rng.choices("ACGT", k=length))


def build_positive(genome_seq, out_dir):
    r1, r2 = [], []
    for gene, start_1based, length in LOCI:
        start0 = start_1based - 1
        fragment = genome_seq[start0:start0 + length]
        if "N" in fragment.upper():
            raise SystemExit(f"fixture window for {gene} contains an N base; pick a different offset")
        if len(fragment) != length:
            raise SystemExit(f"fixture window for {gene} ran off the end of {CONTIG}")
        mate1 = fragment[:READ_LEN]
        mate2 = revcomp(fragment[READ_LEN:])
        name = f"sim_pos_{gene}"
        r1.append((name, mate1))
        r2.append((name, mate2))
    write_fastq_gz(out_dir / "positive_R1.fastq.gz", r1)
    write_fastq_gz(out_dir / "positive_R2.fastq.gz", r2)
    return len(r1)


def build_negative(out_dir, n_pairs, seed):
    rng = random.Random(seed)
    r1, r2 = [], []
    for i in range(n_pairs):
        name = f"sim_neg_{i}"
        r1.append((name, random_seq(rng, READ_LEN)))
        r2.append((name, random_seq(rng, READ_LEN)))
    write_fastq_gz(out_dir / "negative_R1.fastq.gz", r1)
    write_fastq_gz(out_dir / "negative_R2.fastq.gz", r2)
    return len(r1)


def build_synthetic_host(out_dir, name, header, length, seed):
    rng = random.Random(seed)
    seq = random_seq(rng, length)
    path = out_dir / name
    with path.open("w") as handle:
        handle.write(f">{header}\n")
        for i in range(0, len(seq), 70):
            handle.write(seq[i:i + 70] + "\n")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--negative-pairs", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260923)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    genome_seq = load_contig(REFERENCE_FASTA, CONTIG)
    n_pos = build_positive(genome_seq, args.out_dir)
    n_neg = build_negative(args.out_dir, args.negative_pairs, seed=args.seed + 1)
    build_synthetic_host(
        args.out_dir, "synthetic_hg38.fa", "synthetic_host_chr1",
        length=4000, seed=args.seed + 2,
    )
    build_synthetic_host(
        args.out_dir, "synthetic_t2t_phix.fa", "synthetic_t2t_phix_chr1",
        length=3000, seed=args.seed + 3,
    )

    print(f"positive pairs: {n_pos} (loci: {', '.join(g for g, _, _ in LOCI)})")
    print(f"negative pairs: {n_neg}")
    print(f"wrote fixtures to {args.out_dir}")


if __name__ == "__main__":
    main()
