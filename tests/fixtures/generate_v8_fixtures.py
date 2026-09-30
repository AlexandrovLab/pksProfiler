#!/usr/bin/env python3
"""Tiny, real fixtures for the V8a and V8b lanes (tests/check_v8_lanes_execution.sh).

Everything here is cut from the repository's own shipped IHE3034 reference -- public
sequence, no sample data -- and mutated deterministically, so the two lanes can be run
through real minimap2 and the real decision scripts on reads whose correct answer is known
by construction:

  V8a (locus-detection specificity, pks-NEGATIVE panel)
    true reference   a 20 kb slice of the canonical clb island
    decoy panel      (a) random background sequence, and (b) a "homolog" -- a 6 kb slice
                     of the island mutated at --homolog-divergence, standing in for the
                     conserved NRPS/PKS megasynthase domains clbB/clbC/clbJ/clbK share
                     with unrelated clusters, which is the entire failure mode V8a exists
                     to catch
    reads            from the true slice (must come out confirmed_true) and from the
                     homolog (must come out confirmed_decoy). The homolog slice is taken
                     from a region of the island that is NOT in the true slice, so a
                     homolog read has one clearly better home and the expected answer is
                     not a coin flip.

  V8b (carrier attribution, pks-POSITIVE panel)
    carriers         the same island as two carriers: an E. coli copy (identical) and a
                     Klebsiella copy (mutated at --carrier-divergence)
    reads            from each copy. A read from the Klebsiella copy must be attributed to
                     the Klebsiella carrier -- the acceptance case for "which organism
                     carries it", which is a different question from V8a's.

Also writes the two per-sample inputs the V8a decision needs (read_evidence.tsv and
read_clb_gene.tsv) in the pipeline's own schemas, and the sample sheets both standalone
module entry points read.
"""
import argparse
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic.fna"
ISLAND_START_1BASED = 2193827          # params.pks_start_1based
ISLAND_END_1BASED = 2244594            # params.pks_end_1based

# The clb genes that do NOT carry the broadly conserved megasynthase domains. The V8a
# decision requires survivors to cover at least two of them, so the fixture's true reads
# are spread across three.
NON_MEGASYNTHASE = ("clbA", "clbD", "clbG")


def read_reference(path):
    sequence = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.startswith(">"):
                sequence.append(line.strip())
    return "".join(sequence)


def mutate(sequence, divergence, rng):
    """Substitute a fixed fraction of positions, deterministically."""
    bases = list(sequence)
    positions = rng.sample(range(len(bases)), int(len(bases) * divergence))
    for position in positions:
        bases[position] = rng.choice([b for b in "ACGT" if b != bases[position].upper()])
    return "".join(bases)


def write_fasta(path, records):
    with Path(path).open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i + 60] + "\n")


def write_fastq(path, reads):
    with Path(path).open("w") as handle:
        for name, sequence in reads:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def sample_reads(sequence, prefix, count, length, rng):
    reads = []
    for i in range(count):
        start = rng.randrange(0, len(sequence) - length)
        reads.append((f"{prefix}_{i}", sequence[start:start + length]))
    return reads


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--reads-per-source", type=int, default=40)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--homolog-divergence", type=float, default=0.12,
                        help="how far the NRPS/PKS look-alike is from the island: far "
                             "enough that a read from it is not the island, close enough "
                             "that it still cross-maps (default: 0.12)")
    parser.add_argument("--carrier-divergence", type=float, default=0.03,
                        help="how far the Klebsiella carrier copy is from the E. coli one "
                             "(default: 0.03)")
    parser.add_argument("--seed", type=int, default=20260930)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    genome = read_reference(REFERENCE)
    island = genome[ISLAND_START_1BASED - 1:ISLAND_END_1BASED]
    if len(island) != ISLAND_END_1BASED - ISLAND_START_1BASED + 1:
        raise SystemExit(f"island slice is {len(island)} bp; the shipped reference changed")

    true_slice = island[:20000]
    homolog_source = island[30000:36000]        # a different part of the island entirely
    homolog = mutate(homolog_source, args.homolog_divergence, rng)
    background = "".join(rng.choice("ACGT") for _ in range(30000))

    # ---- V8a panel ------------------------------------------------------------------
    write_fasta(out / "v8a_true.fna", [("true__clb_island_slice", true_slice)])
    write_fasta(out / "v8a_decoy.fna", [
        ("decoy__nrps_pks_homolog", homolog),
        ("decoy__random_background", background),
    ])
    write_fasta(out / "v8a_concatenated.fna", [
        ("true__clb_island_slice", true_slice),
        ("decoy__nrps_pks_homolog", homolog),
        ("decoy__random_background", background),
    ])
    (out / "v8a_true_contigs.txt").write_text("true__clb_island_slice\n")
    (out / "v8a_decoy_contigs.txt").write_text(
        "decoy__nrps_pks_homolog\ndecoy__random_background\n")

    # ---- V8b panel: the same island as two pks-positive carriers ---------------------
    klebsiella = mutate(true_slice, args.carrier_divergence, rng)
    write_fasta(out / "v8b_carriers.fna", [
        ("ecoli_test__clb_island", true_slice),
        ("klebsiella_test__clb_island", klebsiella),
    ])
    (out / "v8b_carrier_map.tsv").write_text(
        "contig\tcarrier_key\tcarrier_organism\tcarrier_accession\n"
        "ecoli_test__clb_island\tecoli_test\tEscherichia coli TEST\tGCF_TEST_E\n"
        "klebsiella_test__clb_island\tklebsiella_test\tKlebsiella pneumoniae TEST\tGCF_TEST_K\n")

    # ---- reads ----------------------------------------------------------------------
    true_reads = sample_reads(true_slice, "true", args.reads_per_source,
                               args.read_length, rng)
    decoy_reads = sample_reads(homolog, "decoy", args.reads_per_source,
                                args.read_length, rng)
    kleb_reads = sample_reads(klebsiella, "kleb", args.reads_per_source,
                               args.read_length, rng)

    # V8a sample: true reads plus homolog reads, the mixture the check has to separate.
    v8a_reads = true_reads + decoy_reads
    write_fastq(out / "v8a_sample.fastq", v8a_reads)
    # V8b sample: reads from the Klebsiella carrier only, so "which carrier" has one
    # correct answer and a wrong one that is 97% identical to it.
    write_fastq(out / "v8b_sample.fastq", kleb_reads)

    # ---- the two per-sample inputs the V8a decision reads ----------------------------
    # classify_tumor_pks_evidence.py's wgs schema, filled with this fixture's own numbers.
    (out / "v8a_read_evidence.tsv").write_text(
        "sample\tmode\tread_evidence\tpks_reads\tclb_genes_detected"
        "\tisland_breadth_1x\tisland_breadth_2x\tisland_breadth_3x\n"
        f"v8a_sample\twgs\tbroad_island\t{len(v8a_reads)}\t8"
        "\t0.100000\t0.080000\t0.060000\n")

    # assign_reads_to_genes.py's schema. True reads are spread across three
    # non-megasynthase clb genes so a surviving set can satisfy the decision's
    # non-megasynthase breadth requirement; homolog reads are parked on clbB and clbK,
    # the megasynthase genes they would really hit.
    lines = ["read_id\tGene\toverlap_bp"]
    for i, (name, sequence) in enumerate(true_reads):
        lines.append(f"{name}\t{NON_MEGASYNTHASE[i % len(NON_MEGASYNTHASE)]}\t{len(sequence)}")
    for i, (name, sequence) in enumerate(decoy_reads):
        lines.append(f"{name}\t{'clbB' if i % 2 else 'clbK'}\t{len(sequence)}")
    (out / "v8a_read_clb_gene.tsv").write_text("\n".join(lines) + "\n")

    # ---- sample sheets for the two standalone module entry points --------------------
    (out / "v8a_sheet.csv").write_text(
        "sample,reads_fastq,read_evidence,read_gene_assignment\n"
        f"v8a_sample,{out / 'v8a_sample.fastq'},{out / 'v8a_read_evidence.tsv'},"
        f"{out / 'v8a_read_clb_gene.tsv'}\n")
    (out / "v8b_sheet.csv").write_text(
        f"sample,reads_fastq\nv8b_sample,{out / 'v8b_sample.fastq'}\n")

    # What the assertions expect, written here so the two cannot drift apart.
    (out / "expected.tsv").write_text(
        "key\tvalue\n"
        f"true_reads\t{len(true_reads)}\n"
        f"decoy_reads\t{len(decoy_reads)}\n"
        f"carrier_reads\t{len(kleb_reads)}\n"
        "true_read_prefix\ttrue_\n"
        "decoy_read_prefix\tdecoy_\n"
        "expected_carrier\tklebsiella_test\n"
        "expected_genus\tKlebsiella\n")

    print(f"V8 fixtures in {out}")
    print(f"  true slice {len(true_slice)} bp, homolog {len(homolog)} bp "
          f"({args.homolog_divergence:.0%} divergent), background {len(background)} bp")
    print(f"  carrier copy {args.carrier_divergence:.0%} divergent")
    print(f"  reads: {len(true_reads)} true + {len(decoy_reads)} homolog (V8a), "
          f"{len(kleb_reads)} carrier (V8b)")


if __name__ == "__main__":
    main()
