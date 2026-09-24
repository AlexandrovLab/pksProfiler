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
  * `broad_R{1,2}.fastq.gz` / `extensive_R{1,2}.fastq.gz` -- the u-3 follow-up's
    cohort-scaling ask: two more positive fixtures, built the same way as
    `positive_*` but drawing on a wider pool of validated clb windows (see
    TIER_WINDOW_POOL below) so they clear the `broad_island` and
    `extensive_island` tier floors in scripts/classify_tumor_pks_evidence.py
    instead of just `multi_gene`.
  * `borderline_R{1,2}.fastq.gz` -- a single-locus, 2-read fixture. `reads < 5`
    fails classify_tumor_pks_evidence.py's multi_gene floor before genes or
    breadth are ever consulted, so this lands on `localized_indeterminate`
    (a real tier, not a formal positive one) regardless of which single window
    it uses -- the boundary case a 2-sample cohort has no way to exercise.
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

# u-3 follow-up (cohort aggregation beyond 2 samples): a second, independent pool of
# clb windows for the broad_island/extensive_island/localized_indeterminate fixtures
# below, so they don't have to reuse (or collide with) the exact positions LOCI
# already committed the positive/negative assertions to.
#
# Verified the same way LOCI above was: every window is a 300 bp fragment (two 150 bp
# mates, first half + revcomp(second half), no gap) at least 150 bp inside its gene's
# annotated body (ref/annotations/IHE3034.clbA-clbS.gff), with a >=100 bp gap to any
# other window in the same gene. All 12 genes here (clbA, clbC, clbD, clbF, clbG,
# clbM, clbN, clbO, clbP, clbQ, clbS) and all 38 windows across them were confirmed by
# a real standalone bowtie2 run (--very-sensitive --no-unal | samtools view -q 40,
# the exact flags Modules/pksProfiler_align.nf uses) against the shipped
# GCF_000025745.1 index: every one of the 76 mates maps uniquely, MAPQ 42, CIGAR
# 150M, at its expected coordinate. clbB, clbH, clbI, clbJ and clbK are deliberately
# absent -- the same megasynthase/NRPS domain-repeat risk documented above for clbK
# applies to all five, and none were re-checked here; clbE and clbR are absent because
# they are too short (248 bp / 212 bp) to fit a margined 300 bp window at all.
TIER_WINDOW_POOL = [
    ("clbS", 2193850), ("clbQ", 2194524),
    ("clbP", 2195239), ("clbP", 2195639), ("clbP", 2196039),
    ("clbO", 2196766), ("clbO", 2197166), ("clbO", 2197566), ("clbO", 2197966), ("clbO", 2198366),
    ("clbN", 2199256), ("clbN", 2199656), ("clbN", 2200056), ("clbN", 2200456), ("clbN", 2200856),
    ("clbN", 2201256), ("clbN", 2201656), ("clbN", 2202056), ("clbN", 2202456), ("clbN", 2202856),
    ("clbM", 2203620), ("clbM", 2204020), ("clbM", 2204420),
    ("clbL", 2205121), ("clbL", 2205521), ("clbL", 2205921),
    ("clbG", 2227522), ("clbG", 2227922),
    ("clbF", 2228787), ("clbF", 2229187),
    ("clbD", 2230199),
    ("clbC", 2231078), ("clbC", 2231478), ("clbC", 2231878), ("clbC", 2232278), ("clbC", 2232678), ("clbC", 2233078),
    ("clbA", 2244010),
]
TIER_WINDOW_LEN = 300  # READ_LEN * 2, same construction as LOCI: no gap, no overlap.

# broad_island floor (scripts/classify_tumor_pks_evidence.py): >=30 reads, >=8 genes,
# >=7.5% breadth. 17 windows (7 single-window genes + all 10 clbN windows) = 34 reads,
# 8 genes, 5100/50768 = 10.0% breadth -- comfortably inside [broad_island,
# extensive_island) on every one of the three criteria.
BROAD_WINDOWS = (
    [("clbS", 2193850), ("clbQ", 2194524), ("clbD", 2230199), ("clbA", 2244010),
     ("clbG", 2227522), ("clbF", 2228787), ("clbM", 2203620)]
    + [(gene, start) for gene, start in TIER_WINDOW_POOL if gene == "clbN"]
)

# extensive_island floor: >=100 reads, >=10 genes, >=15% breadth. All 38 pool windows
# (12 genes) give 76 reads at 11400/50768 = 22.5% breadth already clearing the breadth
# and gene floors; EXTENSIVE_DUPLICATE_WINDOWS repeats 15 of those exact windows under
# a second read name each (same genomic sequence, not a new position) to clear the
# read-count floor too, without adding any breadth a real duplicate fragment wouldn't:
# 76 + 2*15 = 106 reads.
EXTENSIVE_WINDOWS = list(TIER_WINDOW_POOL)
EXTENSIVE_DUPLICATE_WINDOWS = TIER_WINDOW_POOL[:15]

# localized_indeterminate: real tier, not a formal positive one (see
# classify_tumor_pks_evidence.py's module docstring) -- `reads < 5` fails the
# multi_gene floor before genes or breadth are even consulted, so a single window
# (2 reads, 1 gene) lands here regardless of which one it is.
BORDERLINE_WINDOWS = [("clbP", 2195239)]

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


def build_tier_sample(genome_seq, out_dir, sample_name, windows, duplicate_windows=()):
    """Like build_positive, generalized to an arbitrary window list (gene, start).

    `duplicate_windows` repeats the named (gene, start) windows under a second read
    name (suffix _dup) -- the same genomic fragment, sequenced "twice" -- to pad the
    read count for a tier floor without inventing a new, unvalidated genomic
    position. Real duplicate fragments add read depth at a position, not breadth.
    """
    r1, r2 = [], []
    seen_at = {}

    def add_pair(gene, start_1based, suffix=""):
        start0 = start_1based - 1
        fragment = genome_seq[start0:start0 + TIER_WINDOW_LEN]
        if "N" in fragment.upper():
            raise SystemExit(f"{sample_name} fixture window for {gene}@{start_1based} contains an N base")
        if len(fragment) != TIER_WINDOW_LEN:
            raise SystemExit(f"{sample_name} fixture window for {gene}@{start_1based} ran off the end of {CONTIG}")
        mate1 = fragment[:READ_LEN]
        mate2 = revcomp(fragment[READ_LEN:])
        name = f"sim_{sample_name}_{gene}_{start_1based}{suffix}"
        r1.append((name, mate1))
        r2.append((name, mate2))

    for gene, start in windows:
        if (gene, start) in seen_at:
            raise SystemExit(f"{sample_name}: duplicate window {gene}@{start} in windows list")
        seen_at[(gene, start)] = True
        add_pair(gene, start)

    for i, (gene, start) in enumerate(duplicate_windows):
        if (gene, start) not in seen_at:
            raise SystemExit(f"{sample_name}: duplicate_windows entry {gene}@{start} is not in windows")
        add_pair(gene, start, suffix=f"_dup{i}")

    write_fastq_gz(out_dir / f"{sample_name}_R1.fastq.gz", r1)
    write_fastq_gz(out_dir / f"{sample_name}_R2.fastq.gz", r2)
    genes = sorted({gene for gene, _ in windows})
    return len(r1), genes


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

    # u-3 follow-up: a slightly larger, mixed-tier cohort (broad_island,
    # extensive_island, localized_indeterminate, on top of positive's multi_gene and
    # negative's negative) so the cohort-level reducers (cohortReport,
    # masterQCSummary, masterSummary, masterTableAlign) get more than two rows to
    # aggregate, sort and de-duplicate.
    n_broad, broad_genes = build_tier_sample(genome_seq, args.out_dir, "broad", BROAD_WINDOWS)
    n_extensive, extensive_genes = build_tier_sample(
        genome_seq, args.out_dir, "extensive", EXTENSIVE_WINDOWS, EXTENSIVE_DUPLICATE_WINDOWS)
    n_borderline, borderline_genes = build_tier_sample(genome_seq, args.out_dir, "borderline", BORDERLINE_WINDOWS)

    print(f"positive pairs: {n_pos} (loci: {', '.join(g for g, _, _ in LOCI)})")
    print(f"negative pairs: {n_neg}")
    print(f"broad pairs: {n_broad} (genes: {', '.join(broad_genes)})")
    print(f"extensive pairs: {n_extensive} (genes: {', '.join(extensive_genes)})")
    print(f"borderline pairs: {n_borderline} (genes: {', '.join(borderline_genes)})")
    print(f"wrote fixtures to {args.out_dir}")


if __name__ == "__main__":
    main()
