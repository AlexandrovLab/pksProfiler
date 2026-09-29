#!/usr/bin/env python3
"""Fixtures for tests/check_mode_execution.sh: one small real input per mode.

Ludmil, 26 Sep mode memo section 6: "Add one real execution test per modality. Verify
which branches run, which outputs appear, and which tiers/assembly features are
intentionally disabled."

Everything here is public reference sequence or simulated from it:
  * generate_e2e_fixtures.py's paired pks reads (IHE3034, the repo's own reference) and its
    two synthetic host references, for the short-read DNA modes (wgs, wes, metagenome);
  * a 30 kb slice of public GRCh38 (GRCh38.d1.vd1.fa, read through its .fai) as the human
    reference for the RNA and long-read lanes, with a three-exon "transcript" placed on it
    at real GT...AG boundaries and written as a BED12 junction annotation;
  * RNA short reads: the positive pks R1 reads plus 150 bp reads spanning the exon1-exon2
    junction, which only a splice-aware human pass can remove;
  * long reads (ONT-like, 3% substitutions): 3 kb pks-island reads, 3 kb unspliced human
    reads, and 2.4 kb spliced human reads (exon1+2+3) for rna + ont.

    generate_mode_fixtures.py --out-dir DIR [--grch38 FASTA] [--region chr22:START-END]

Writes DIR/mode_fixtures.json with the expected counts check_mode_execution.sh asserts on.
"""
import argparse
import json
import random
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import generate_e2e_fixtures as e2e  # noqa: E402

GRCH38 = Path("/tscc/projects/ps-lalexandrov/shared/Reference_Genomes/GRCh38.d1.vd1/GRCh38.d1.vd1.fa")
EXON_LEN = 800
JUNCTION_READS = 20
LONG_LEN = 3000
LONG_ERROR = 0.03


def fai_slice(fasta, contig, start, end):
    """1-based inclusive slice of `contig`, read straight from the .fai offsets."""
    index = {}
    for line in Path(f"{fasta}.fai").read_text().splitlines():
        name, length, offset, line_bases, line_width = line.split("\t")[:5]
        index[name] = (int(length), int(offset), int(line_bases), int(line_width))
    length, offset, line_bases, line_width = index[contig]
    if end > length:
        raise SystemExit(f"{contig}:{start}-{end} is past the end ({length})")
    begin = offset + (start - 1) // line_bases * line_width + (start - 1) % line_bases
    rows = (end - start + 1) // line_bases + 2
    with open(fasta, "rb") as fh:
        fh.seek(begin)
        raw = fh.read(rows * line_width).decode()
    seq = raw.replace("\n", "")[: end - start + 1].upper()
    if len(seq) != end - start + 1:
        raise SystemExit("short read from the reference FASTA")
    return seq


def find(seq, frm, pattern, offset):
    """First position p >= frm where seq[p + offset: p + offset + 2] == pattern."""
    for p in range(frm, len(seq)):
        if seq[p + offset:p + offset + 2] == pattern:
            return p
    raise SystemExit(f"no {pattern} found after {frm}")


def mutate(rng, seq, rate):
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def write_fasta(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--grch38", type=Path, default=GRCH38)
    ap.add_argument("--region", default="chr22:30000001-30030000")
    ap.add_argument("--seed", type=int, default=20260929)
    a = ap.parse_args()
    out = a.out_dir
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(a.seed)

    # short-read DNA fixtures and synthetic hosts, exactly as the e2e check builds them
    subprocess.run([sys.executable, str(HERE / "generate_e2e_fixtures.py"), "--out-dir", str(out)],
                   check=True, stdout=subprocess.DEVNULL)

    contig, span = a.region.split(":")
    start, end = (int(x) for x in span.split("-"))
    human = fai_slice(a.grch38, contig, start, end)
    if human.count("N") > len(human) // 100:
        raise SystemExit(f"{a.region} is N-rich; pick another region")
    slice_name = f"{contig}_{start}_{end}"
    write_fasta(out / "human_slice.fa", slice_name, human)

    # A three-exon transcript at real splice sites: each exon ends before a GT donor and
    # the next begins after an AG acceptor.
    e1 = 2000
    d1 = find(human, e1 + EXON_LEN, "GT", 0)            # intron 1 starts at d1
    a1 = find(human, d1 + 3000, "AG", -2)               # exon 2 starts at a1
    d2 = find(human, a1 + EXON_LEN, "GT", 0)
    a2 = find(human, d2 + 3000, "AG", -2)
    exons = [(e1, d1), (a1, d2), (a2, a2 + EXON_LEN)]
    sizes = [b - s for s, b in exons]
    with open(out / "junctions.bed", "w") as fh:
        fh.write("\t".join(map(str, [slice_name, exons[0][0], exons[-1][1], "synthetic_tx", 0, "+",
                                     exons[0][0], exons[-1][1], 0, 3,
                                     ",".join(map(str, sizes)) + ",",
                                     ",".join(str(s - exons[0][0]) for s, _ in exons) + ","])) + "\n")
    exon_seq = [human[s:b] for s, b in exons]
    spliced = "".join(exon_seq)

    # RNA short: positive pks R1 reads + reads across the exon1|exon2 junction
    import gzip
    with gzip.open(out / "positive_R1.fastq.gz", "rt") as fh:
        lines = fh.read().splitlines()
    pks_short = [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 4)]
    junction = []
    for i in range(JUNCTION_READS):
        left = 50 + (i * 50) // JUNCTION_READS          # 50..99 bp on exon 1
        read = exon_seq[0][-left:] + exon_seq[1][:e2e.READ_LEN - left]
        junction.append((f"sim_human_junction_{i}", read))
    e2e.write_fastq_gz(out / "rna_short.fastq.gz", pks_short + junction)

    # Long reads
    genome = e2e.load_contig(e2e.REFERENCE_FASTA, e2e.CONTIG)
    island0 = 2193827 - 1
    pks_long = [(f"sim_pks_long_{i}", mutate(rng, genome[island0 + 2000 + 7000 * i:
                                                          island0 + 2000 + 7000 * i + LONG_LEN], LONG_ERROR))
                for i in range(6)]
    human_long = [(f"sim_human_long_{i}", mutate(rng, human[20000 + 2500 * i: 20000 + 2500 * i + LONG_LEN],
                                                 LONG_ERROR)) for i in range(3)]
    spliced_long = [(f"sim_human_spliced_long_{i}", mutate(rng, spliced, LONG_ERROR)) for i in range(2)]

    def write_long(path, records):
        with gzip.open(path, "wt") as fh:
            for name, seq in records:
                fh.write(f"@{name}\n{seq}\n+\n{'5' * len(seq)}\n")

    write_long(out / "wgs_ont.fastq.gz", pks_long + human_long)
    write_long(out / "rna_ont.fastq.gz", pks_long + spliced_long)

    expected = {
        "short_positive_reads": 2 * len(pks_short),
        "rna_short_input": len(pks_short) + len(junction),
        "rna_short_junction_reads": len(junction),
        "rna_short_pks_reads": len(pks_short),
        "long_pks_reads": len(pks_long),
        "wgs_ont_input": len(pks_long) + len(human_long),
        "wgs_ont_human_reads": len(human_long),
        "rna_ont_input": len(pks_long) + len(spliced_long),
        "rna_ont_spliced_reads": len(spliced_long),
        "exons_0based_half_open": exons,
        "slice": slice_name,
    }
    (out / "mode_fixtures.json").write_text(json.dumps(expected, indent=2) + "\n")
    print(json.dumps(expected))


if __name__ == "__main__":
    main()
