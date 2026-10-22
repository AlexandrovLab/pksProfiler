#!/usr/bin/env python3
"""What is actually inside a minimap2 .mmi index: its build settings and its sequences.

rc2 review item 2.7 (Ludmil's 27 Aug ask, reopened): "Record exactly which human
assemblies/decoys/phiX are present in human-GRC-db.mmi / human-GCA-phix-db.mmi and the
build command used." The two indexes are the host-depletion step -- everything the
pipeline then counts is what they did not remove -- and until now the only thing recorded
about them anywhere was their path and their size.

A .mmi is not opaque. minimap2 writes its own build settings and its complete sequence
dictionary in the first bytes of the file, before the minimizer buckets, so both can be
read in milliseconds without loading the index, without minimap2, and without reading the
9 GB of index data behind them. From mm_idx_dump() in minimap2's mmidx.c:

    char[4]              "MMI\\x02"           the magic and format version
    uint32 x 5           w, k, b, n_seq, flag  little-endian
    n_seq times:
        uint8            length of the name (0 when the index was built with -N)
        char[len]        the name, not NUL-terminated
        uint32           the sequence length in bases

`flag` is a bitfield: 1 = MM_I_HPC (homopolymer-compressed minimizers, i.e. built with
-H), 2 = MM_I_NO_SEQ (built with no target sequence stored, which minimap2 -d does),
4 = MM_I_NO_NAME (names discarded).

What that gives us, and what it does not:

  * the contents are EXACT -- every sequence name and length, so which assembly, which
    decoys, whether phiX and EBV are in there, is not a matter of inference;
  * (w, k, flag) pins down the preset family the index was built with, because minimap2's
    presets differ in exactly those: map-ont/map-pb/default k=15 w=10, sr k=21 w=11,
    asm5/asm10/asm20 k=19 w=19/19/10, map-hifi k=19 w=19. It does NOT recover the full
    command line: the input FASTA's path, -I, -t and anything that does not change the
    index layout are not stored in the file and cannot be recovered from it.

This script prints a summary, and with --tsv the full sequence dictionary.

    python scripts/describe_minimap2_index.py <index.mmi> [--tsv out.tsv]

It is a provenance utility, not part of any pipeline lane: nothing in main.nf calls it.
"""
import argparse
import hashlib
import re
import struct
import sys
from pathlib import Path

MAGIC = b"MMI\x02"
# minimap2 mmpriv.h / minimap.h
FLAGS = {1: "MM_I_HPC (-H, homopolymer-compressed minimizers)",
         2: "MM_I_NO_SEQ (no target sequence stored; what `minimap2 -d` writes)",
         4: "MM_I_NO_NAME (sequence names discarded)"}
# (k, w) -> the presets that build an index with those settings. minimap2 2.17-2.28.
PRESETS = {
    (15, 10): "map-ont / map-pb / lr:hq / the default (no -x)",
    (21, 11): "sr (short reads)",
    (19, 19): "map-hifi / asm5 / asm10",
    (19, 10): "asm20",
    (28, 1000): "splice (ava presets differ)",
    (15, 5): "splice / splice:hq (-k15 -w5)",
}

# How a human reference's sequences are usually named. Only used to GROUP the exact
# dictionary for a reader; nothing is inferred, every name is still printed with --tsv.
CATEGORIES = (
    ("phiX", re.compile(r"^(phi ?X|NC_001422)", re.I)),
    ("EBV decoy", re.compile(r"^(chrEBV|NC_007605)", re.I)),
    ("HLA allele", re.compile(r"^HLA-", re.I)),
    ("ALT contig", re.compile(r"_alt$", re.I)),
    ("decoy (hs38d1)", re.compile(r"_decoy$|^chrUn_JTFH", re.I)),
    ("unlocalized/unplaced", re.compile(r"_random$|^chrUn|^chr\w+_[A-Z]{2}\d+", re.I)),
    ("mitochondrion", re.compile(r"^(chrM|MT)$", re.I)),
    ("primary chromosome", re.compile(r"^(chr)?([1-9]|1\d|2[0-2]|X|Y)$", re.I)),
)


def read_header(handle):
    magic = handle.read(4)
    if magic != MAGIC:
        raise SystemExit(f"not a minimap2 index (magic {magic!r}, expected {MAGIC!r})")
    w, k, b, n_seq, flag = struct.unpack("<5I", handle.read(20))
    return {"w": w, "k": k, "b": b, "n_seq": n_seq, "flag": flag}


def read_sequences(handle, n_seq):
    """[(name, length)] in index order. Reads only the dictionary, not the buckets."""
    sequences = []
    for index in range(n_seq):
        raw = handle.read(1)
        if not raw:
            raise SystemExit(f"index truncated after {index} of {n_seq} sequences")
        name_length = raw[0]
        name = handle.read(name_length).decode("utf-8", "replace") if name_length else ""
        packed = handle.read(4)
        if len(packed) != 4:
            raise SystemExit(f"index truncated in sequence {index} of {n_seq}")
        sequences.append((name, struct.unpack("<I", packed)[0]))
    return sequences


def categorize(name):
    for label, pattern in CATEGORIES:
        if pattern.search(name):
            return label
    return "other"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("index", help="the .mmi file")
    parser.add_argument("--tsv", help="write the full sequence dictionary here")
    args = parser.parse_args()

    path = Path(args.index)
    with path.open("rb") as handle:
        header = read_header(handle)
        sequences = read_sequences(handle, header["n_seq"])
        dictionary_bytes = handle.tell()
    # A real content sha256 of the part of the file that says WHAT the index is: the
    # build settings and the complete sequence dictionary. It is a few kB, so unlike a
    # checksum of the whole 9 GB index it costs nothing and can be re-verified on every
    # run if we ever want to. It changes if the index is rebuilt from a different FASTA,
    # with different settings, or with the sequences in a different order. It does NOT
    # cover the minimizer buckets, so it is a strong identity check, not an integrity
    # check: `sha256sum <index.mmi>` is the integrity check, and it reads the whole file.
    with path.open("rb") as handle:
        dictionary_sha256 = hashlib.sha256(handle.read(dictionary_bytes)).hexdigest()

    flags = [text for bit, text in sorted(FLAGS.items()) if header["flag"] & bit]
    preset = PRESETS.get((header["k"], header["w"]), "no preset known with these settings")

    groups = {}
    for name, length in sequences:
        label = categorize(name)
        count, bases = groups.get(label, (0, 0))
        groups[label] = (count + 1, bases + length)

    print(f"index              {path}")
    print(f"size_bytes         {path.stat().st_size}")
    print(f"format             MMI v2")
    print(f"minimizer k        {header['k']}")
    print(f"minimizer w        {header['w']}")
    print(f"bucket bits b      {header['b']}")
    print(f"flag               {header['flag']} "
          f"[{'; '.join(flags) if flags else 'none set'}]")
    print(f"preset family      -k {header['k']} -w {header['w']} is {preset}")
    print(f"sequences          {header['n_seq']}")
    print(f"total bases        {sum(length for _name, length in sequences)}")
    print(f"dictionary bytes   {dictionary_bytes}")
    print(f"dictionary sha256  {dictionary_sha256}")
    print()
    print(f"{'category':24s}{'sequences':>12s}{'bases':>16s}")
    for label, (count, bases) in sorted(groups.items(), key=lambda item: -item[1][1]):
        print(f"{label:24s}{count:>12d}{bases:>16d}")
    print()
    longest = sorted(sequences, key=lambda item: -item[1])[:5]
    print("five longest sequences")
    for name, length in longest:
        print(f"  {name:32s}{length:>14d}")

    if args.tsv:
        with open(args.tsv, "w") as out:
            out.write("name\tlength\tcategory\n")
            for name, length in sequences:
                out.write(f"{name}\t{length}\t{categorize(name)}\n")
        print(f"\nfull dictionary written to {args.tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
