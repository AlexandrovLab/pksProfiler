#!/usr/bin/env python3
"""E1 (Ludmil, 25 Sep): SAM fixtures that look like `bwa mem -T 0` output.

tests/check_e2e_execution.sh's seventh run. It takes the positive and negative read pairs
generate_e2e_fixtures.py already wrote (IHE3034 clb windows from the repo's own public
reference, and pseudorandom reads) and writes them as alignment records against the
synthetic host contig the way bwa mem -T 0 would report them: some reads FLAG 4, some
"mapped" to a short spurious hit with a low AS:i, some mapped with AS >= 30, one pair
mapped with no AS tag at all, plus a secondary and a supplementary record.

bwa itself cannot produce these here: -T 0 reports a read's best seed hit, and random
synthetic host sequence gives these reads no 19-mer seed to hit. So the records are
written by hand, with the fields bwa writes (FLAG, a clipped CIGAR, MAPQ 0, NM/AS/XS),
and the check runs the real samtools and pipeline on them.

    bwa_T0.sam   @PG ID:bwa PN:bwa CL:bwa mem -T 0 ...  -> auto resolves on
    star.sam     the same records under a STAR @PG      -> auto resolves off
    low_score_expected.json   per sample: extracted FASTQ records and counts

Which reads, from the four positive pairs (clbS, clbN, clbD, clbA, in LOCI order):
    clbS   both mates FLAG 4                            extracted either way
    clbN   mate 1 mapped AS 21, mate 2 FLAG 4           mate 2 either way, mate 1 only on
    clbD   both mapped, AS 25 / AS 19 (mate 2 reverse)  only on
    clbA   both mapped, AS 30 / AS 45                   never (not below 30)
and the four negative pairs: FLAG 4; mapped with no AS (counted, never added); mapped
AS 150; mapped AS 22 / 28 (on only). The secondary (clbD mate 1) and supplementary (clbN
mate 1, hard-clipped) records carry AS 5 and must never be extracted.
"""
import argparse
import gzip
import json
from pathlib import Path

HOST = "synthetic_host_chr1"   # generate_e2e_fixtures.py's synthetic_hg38.fa
HOST_LEN = 4000
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")

PG = {
    "bwa_T0": "@PG\tID:bwa\tPN:bwa\tVN:0.7.15-r1140\tCL:bwa mem -t 8 -T 0 synthetic_hg38.fa R1.fq R2.fq",
    "star": "@PG\tID:STAR\tPN:STAR\tVN:2.7.5c\tCL:STAR --runThreadN 8 --genomeDir synthetic",
}


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fastq(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    return [(lines[i][1:].split()[0], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


# Per pair: (mate1 plan, mate2 plan). A plan is None (FLAG 4) or (AS or None, reverse?).
POSITIVE_PLAN = [
    (None, None),                    # clbS
    ((21, False), None),             # clbN
    ((25, False), (19, True)),       # clbD
    ((30, False), (45, True)),       # clbA
]
NEGATIVE_PLAN = [
    (None, None),
    (("noAS", False), ("noAS", True)),
    ((150, False), (150, True)),
    ((22, False), (28, True)),
]


def records(pairs, plans, start_pos):
    """SAM lines, the reads each rule extracts, and the class counts."""
    lines, expected_on, expected_off = [], [], []
    counts = {"flag4": 0, "low": 0, "no_as": 0}
    pos = start_pos
    for (name, s1, q1), (_name2, s2, q2), plan in zip(pairs[0], pairs[1], plans):
        mates = [(s1, q1, plan[0], 0x40, "/1"), (s2, q2, plan[1], 0x80, "/2")]
        mapped = [m[2] is not None for m in mates]
        for i, (s, q, p, bit, suffix) in enumerate(mates):
            other = mates[1 - i]
            flag = 0x1 | bit
            if p is None:
                flag |= 0x4
            if not mapped[1 - i]:
                flag |= 0x8
            elif other[2][1]:
                flag |= 0x20
            if p is not None and p[1]:
                flag |= 0x10
            stored_seq, stored_qual = (revcomp(s), q[::-1]) if flag & 0x10 else (s, q)
            if p is None:
                # unmapped: placed at the mate when the mate is mapped, as bwa does
                rname, rpos = (HOST, pos) if mapped[1 - i] else ("*", 0)
                fields = [name, flag, rname, rpos, 0, "*", "=" if rname != "*" else "*",
                          rpos, 0, stored_seq, stored_qual]
                tags = ["MC:Z:21M129S"] if mapped[1 - i] else []
                counts["flag4"] += 1
                expected_on.append((name + suffix, s, q))
                expected_off.append((name + suffix, s, q))
            else:
                score = p[0]
                aligned = {150: 150, "noAS": 60}.get(score) or max(score, 19)
                cigar = f"{aligned}M{len(s) - aligned}S" if aligned < len(s) else f"{len(s)}M"
                fields = [name, flag, HOST, pos, 0, cigar, "=", pos, 0, stored_seq, stored_qual]
                tags = ["NM:i:0"]
                if score == "noAS":
                    counts["no_as"] += 1
                else:
                    tags += [f"AS:i:{score}", "XS:i:0"]
                    if score < 30:
                        counts["low"] += 1
                        expected_on.append((name + suffix, s, q))
            lines.append("\t".join(str(f) for f in fields + tags))
        pos += 40
    return lines, expected_on, expected_off, counts


def extra_records(pos_pairs):
    """A secondary (clbD mate 1) and a supplementary (clbN mate 1) record, AS 5."""
    (clbn, clbn_seq, clbn_qual) = pos_pairs[0][1]
    (clbd, clbd_seq, _q) = pos_pairs[0][2]
    return [
        "\t".join(map(str, [clbd, 0x1 | 0x40 | 0x20 | 0x100, HOST, 2000, 0, "20M130S", "=", 100, 0, "*", "*",
                            "NM:i:0", "AS:i:5"])),
        "\t".join(map(str, [clbn, 0x1 | 0x40 | 0x8 | 0x800, HOST, 3000, 0, "20M130H", "=", 100, 0,
                            clbn_seq[:20], clbn_qual[:20], "NM:i:0", "AS:i:5"])),
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--e2e-fixtures", required=True, type=Path,
                    help="generate_e2e_fixtures.py's --out-dir (positive/negative_R1/R2.fastq.gz)")
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    pos = [read_fastq(a.e2e_fixtures / f"positive_R{m}.fastq.gz") for m in (1, 2)]
    neg = [read_fastq(a.e2e_fixtures / f"negative_R{m}.fastq.gz") for m in (1, 2)]
    if len(pos[0]) != len(POSITIVE_PLAN) or len(neg[0]) != len(NEGATIVE_PLAN):
        raise SystemExit("unexpected number of positive/negative pairs; generate_e2e_fixtures.py changed")
    p_lines, p_on, p_off, p_counts = records(pos, POSITIVE_PLAN, 100)
    n_lines, n_on, n_off, n_counts = records(neg, NEGATIVE_PLAN, 1000)
    body = p_lines + n_lines + extra_records(pos)
    for name, pg in PG.items():
        head = ["@HD\tVN:1.6\tSO:unsorted", f"@SQ\tSN:{HOST}\tLN:{HOST_LEN}",
                "@RG\tID:rg1\tSM:synthetic\tPL:ILLUMINA", pg]
        (a.out_dir / f"{name}.sam").write_text("\n".join(head + body) + "\n")
    counts = {k: p_counts[k] + n_counts[k] for k in p_counts}
    pks_on, pks_off = p_on, p_off
    expected = {
        "on": {"fastq": sorted(p_on + n_on), "extracted_unmapped_reads": counts["flag4"],
               "extracted_low_score_mapped_reads": counts["low"], "no_as": counts["no_as"],
               "extracted_reads": counts["flag4"] + counts["low"], "pks_reads": len(pks_on),
               "pks_genes": sorted({r[0].split("/")[0].rsplit("_", 1)[-1] for r in pks_on})},
        "off": {"fastq": sorted(p_off + n_off), "extracted_unmapped_reads": counts["flag4"],
                "extracted_low_score_mapped_reads": 0, "no_as": None,
                "extracted_reads": counts["flag4"], "pks_reads": len(pks_off),
                "pks_genes": sorted({r[0].split("/")[0].rsplit("_", 1)[-1] for r in pks_off})},
    }
    (a.out_dir / "low_score_expected.json").write_text(json.dumps(expected, indent=1))
    print(f"wrote bwa_T0.sam / star.sam: FLAG 4 {counts['flag4']}, low-score {counts['low']}, "
          f"no AS {counts['no_as']}; pks reads extracted on {len(pks_on)}, off {len(pks_off)}")


if __name__ == "__main__":
    main()
