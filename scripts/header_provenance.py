#!/usr/bin/env python3
"""E3 (Ludmil, 25 Sep): input reference and upstream-aligner provenance from a BAM/CRAM header.

"Hash a canonicalized @SQ SN/LN/M5 list for every BAM/CRAM and capture the @PG chain
(program/version) when present. Store these in run/QC provenance." Acceptance: equivalent
BAM/CRAMs on the same reference give the same digest; different builds give different ones.

Adapted from projects/pksprofiler/validation/E3/e3_headers.py. It reads a header that
extractReads already wrote with `samtools view -H` -- header lines only, no alignment
record -- and never sees read data.

    header_provenance.py sample    --header input.header.sam --sample S --output S.header_provenance.tsv
    header_provenance.py summarize --inputs-file list.txt --table input_headers.tsv --summary input_header_summary.tsv

Per sample, one row:
  sq_fingerprint      md5 over the ordered "SN:LN" lines -- same reference, same order
  m5_set_fingerprint  md5 over the sorted, de-duplicated M5 values; NA unless every @SQ has
                      M5. Name-insensitive, so chr1 vs 1 renamings of one build still match
  n_contigs, n_contigs_with_m5, chr1_name, chr1_length, build_guess (from chr1's length)
  n_pg, pg_chain      @PG records in PP order as ID|PN|VN (CL is never recorded: command lines
                      carry paths and sample names; E1 below looks at it in memory)
  pg_programs         the same chain as PN:VN only

E1 (Ludmil, 25 Sep): the same step decides, per sample, whether extraction also takes
primary mapped records with AS:i below --extract_min_as (--extract_low_score_mapped
auto|on|off). bwa mem reports a read as unmapped only when its best local score is below
-T (default 30); GDC ran TCGA with -T 0, so reads that -T 30 would have left unmapped
(FLAG 4) are instead "mapped" with a short spurious hit and never reach `samtools fastq
-f 4`. AS < 30 among those mapped reads is exactly that set. `auto` turns it on when an
@PG record runs `bwa mem` (or `bwa-mem2 mem`), and off otherwise: the AS scale of other
aligners (STAR, minimap2, bowtie2) is not bwa's, so the threshold means nothing there.
Detecting `mem` needs the @PG CL; it is read in memory, only its program and subcommand
tokens and the -T value are looked at, and none of it is written anywhere. Four more
columns record the decision:
  aligner_bwa_mem           yes | no | unconfirmed (a bwa @PG with no CL, or not `mem`)
  bwa_min_score             the -T in that CL, "default" when not given, NA without bwa mem
  extract_low_score_mapped  on_auto_bwa_mem | off_auto_not_bwa_mem | on_requested | off_requested
  extract_min_as            the AS threshold in force (NA when off)

The cohort summary (and RUN_REPORT.txt's INPUT ALIGNMENTS section, which reads it) holds
counts only: distinct fingerprints with their build and contig count, distinct PN:VN
chains, how many inputs had low-score extraction on and off, never an ID or a path.
"""
import argparse
import collections
import csv
import hashlib
import re
import sys
from pathlib import Path

FIELDS = ["sample", "sq_fingerprint", "m5_set_fingerprint", "n_contigs", "n_contigs_with_m5",
          "chr1_name", "chr1_length", "build_guess", "n_pg", "pg_chain", "pg_programs",
          "aligner_bwa_mem", "bwa_min_score", "extract_low_score_mapped", "extract_min_as"]
CHR1_BUILD = {248956422: "GRCh38", 249250621: "GRCh37", 248387328: "T2T-CHM13"}
SUMMARY_FIELDS = ["section", "key", "detail", "n_files"]


def md5(lines):
    return hashlib.md5("\n".join(lines).encode()).hexdigest()


def tags(line):
    return dict(f.split(":", 1) for f in line.rstrip("\n").split("\t")[1:] if ":" in f)


def pg_chain(pg):
    """@PG records ordered along their PP links, roots first; unreachable ones appended."""
    by_id = {p.get("ID"): p for p in pg}
    children = collections.defaultdict(list)
    for p in pg:
        children[p.get("PP")].append(p)
    roots = [p for p in pg if p.get("PP") not in by_id] or pg[:1]
    out, seen = [], set()

    def walk(p):
        stack = [p]
        while stack:
            node = stack.pop()
            if node.get("ID") in seen:
                continue
            seen.add(node.get("ID"))
            out.append(node)
            stack.extend(reversed(children.get(node.get("ID"), [])))

    for root in roots:
        walk(root)
    out.extend(p for p in pg if p.get("ID") not in seen)
    return out


def pg_label(p):
    """PN:VN; with no PN, the ID stripped of the .N/-N suffix samtools and Picard add."""
    name = p.get("PN") or re.sub(r"([.-]\d+|-[0-9A-F]{6,})$", "", p.get("ID", ""), flags=re.I)
    return f"{name}:{p.get('VN', '')}"


# E1: the programs whose `mem` subcommand is bwa's local-alignment scoring (-T default 30).
BWA_PROGRAMS = ("bwa", "bwa-mem2")
LOW_SCORE_REQUESTS = ("auto", "on", "off")
# The extract_low_score_mode QC metric and the extract_low_score_mapped column: on/off
# and why. scripts/build_qc_summary.py accepts exactly these.
LOW_SCORE_STATES = ("on_auto_bwa_mem", "off_auto_not_bwa_mem", "on_requested", "off_requested")


def _program(token):
    return token.strip("'\"").rsplit("/", 1)[-1]


def bwa_mem_evidence(text):
    """(aligner_bwa_mem, bwa_min_score) from a header's @PG records.

    yes: some @PG's CL runs `bwa mem` / `bwa-mem2 mem` (also inside a wrapper, e.g.
    `sentieon bwa mem`); bwa_min_score is its -T, or "default". unconfirmed: a bwa @PG
    (by PN or ID) whose CL is absent or runs another subcommand (aln/sampe/samse, which
    write no AS tag anyway). no: nothing from bwa. The CL is never returned.
    """
    found, unconfirmed, scores = False, False, []
    for line in text.splitlines():
        if not line.startswith("@PG\t"):
            continue
        t = tags(line)
        tokens = t.get("CL", "").split()
        hit = next((i for i, tok in enumerate(tokens[:-1])
                    if _program(tok) in BWA_PROGRAMS and tokens[i + 1] == "mem"), None)
        if hit is not None:
            found = True
            score = "default"
            rest = tokens[hit + 2:]
            for j, tok in enumerate(rest):
                m = re.fullmatch(r"-T(-?\d+)?", tok)
                if m:
                    value = m.group(1) if m.group(1) is not None else (rest[j + 1] if j + 1 < len(rest) else "")
                    score = value if re.fullmatch(r"-?\d+", value) else "unparsed"
                    break
            scores.append(score)
            continue
        label = (t.get("PN") or re.sub(r"([.-]\d+|-[0-9A-F]{6,})$", "", t.get("ID", ""), flags=re.I)).lower()
        if label in BWA_PROGRAMS:
            unconfirmed = True
    if found:
        return "yes", ",".join(dict.fromkeys(scores))
    return ("unconfirmed" if unconfirmed else "no"), "NA"


def low_score_decision(request, bwa_mem):
    """The extract_low_score_mapped state for --extract_low_score_mapped `request`."""
    if request not in LOW_SCORE_REQUESTS:
        raise ValueError(f"--extract_low_score_mapped must be one of {', '.join(LOW_SCORE_REQUESTS)} (got {request!r})")
    if request == "on":
        return "on_requested"
    if request == "off":
        return "off_requested"
    return "on_auto_bwa_mem" if bwa_mem == "yes" else "off_auto_not_bwa_mem"


def low_score_reason(state, bwa_mem, bwa_min_score, min_as):
    """One line for the task log: what was decided and from what."""
    seen = {"yes": f"bwa mem in the @PG chain (-T {bwa_min_score})",
            "unconfirmed": "a bwa @PG without a `bwa mem` command line",
            "no": "no bwa mem in the @PG chain"}[bwa_mem]
    if state.startswith("on"):
        text = f"on ({state}): also extracting primary mapped records with AS:i < {min_as}; {seen}"
        effective = [30 if v == "default" else int(v) if re.fullmatch(r"-?\d+", v) else None
                     for v in bwa_min_score.split(",")] if bwa_mem == "yes" else [None]
        if all(t is not None and t >= min_as for t in effective):
            text += f"; bwa already left every read scoring below {min_as} unmapped, so none should be added"
        return text
    return f"off ({state}): FLAG 4 records only; {seen}"


def parse_header(text, sample, low_score_request="auto", min_as=30):
    sq, pg = [], []
    for line in text.splitlines():
        if line.startswith("@SQ\t"):
            sq.append(tags(line))
        elif line.startswith("@PG\t"):
            pg.append({k: v for k, v in tags(line).items() if k != "CL"})
    if not sq:
        raise ValueError("no @SQ lines: an unaligned or truncated header")
    m5 = [t.get("M5", "").lower() for t in sq]
    chr1 = next((t for t in sq if t.get("SN") in ("chr1", "1")), None)
    chr1_len = int(chr1["LN"]) if chr1 and chr1.get("LN", "").isdigit() else None
    chain = pg_chain(pg)
    bwa_mem, bwa_min_score = bwa_mem_evidence(text)
    state = low_score_decision(low_score_request, bwa_mem)
    return {
        "sample": sample,
        "sq_fingerprint": md5([f"{t.get('SN', '')}:{t.get('LN', '')}" for t in sq]),
        "m5_set_fingerprint": md5(sorted(set(m5))) if all(m5) else "NA",
        "n_contigs": len(sq),
        "n_contigs_with_m5": sum(1 for v in m5 if v),
        "chr1_name": chr1.get("SN") if chr1 else "NA",
        "chr1_length": chr1_len if chr1_len is not None else "NA",
        "build_guess": (CHR1_BUILD.get(chr1_len, "unknown") if chr1_len is not None else "no_chr1"),
        "n_pg": len(pg),
        "pg_chain": " > ".join(f"{p.get('ID', '')}|{p.get('PN', '')}|{p.get('VN', '')}" for p in chain) or "NA",
        "pg_programs": " > ".join(pg_label(p) for p in chain) or "NA",
        "aligner_bwa_mem": bwa_mem,
        "bwa_min_score": bwa_min_score,
        "extract_low_score_mapped": state,
        "extract_min_as": min_as if state.startswith("on") else "NA",
    }


def write_rows(path, rows, fields):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows):
    """Counts only: no sample, no @PG ID, no path."""
    out = [{"section": "files", "key": "alignments", "detail": "", "n_files": len(rows)}]
    refs = collections.Counter((r["sq_fingerprint"], r["m5_set_fingerprint"], r["build_guess"],
                                str(r["n_contigs"])) for r in rows)
    for (sq, m5, build, n), count in sorted(refs.items(), key=lambda kv: (-kv[1], kv[0])):
        out.append({"section": "reference", "key": sq,
                    "detail": f"build={build} n_contigs={n} m5_set={m5}", "n_files": count})
    chains = collections.Counter(r["pg_programs"] if r["pg_programs"] != "NA" else "(no @PG)"
                                 for r in rows)
    for chain, count in sorted(chains.items(), key=lambda kv: (-kv[1], kv[0])):
        out.append({"section": "aligner_chain", "key": chain, "detail": "", "n_files": count})
    # E1: how many inputs had low-score mapped extraction on, and why.
    states = collections.Counter(r.get("extract_low_score_mapped", "NA") for r in rows)
    for state, count in sorted(states.items()):
        out.append({"section": "extract_low_score", "key": state, "detail": "", "n_files": count})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    one = sub.add_parser("sample")
    one.add_argument("--header", required=True, type=Path)
    one.add_argument("--sample", required=True)
    one.add_argument("--output", required=True, type=Path)
    one.add_argument("--low-score-mode", default="auto", choices=LOW_SCORE_REQUESTS,
                     help="E1: --extract_low_score_mapped (default auto)")
    one.add_argument("--min-as", type=int, default=30, help="E1: --extract_min_as (default 30)")
    one.add_argument("--decision", type=Path, default=None,
                     help="E1: write 'on|off<TAB>state<TAB>reason' here for extractReads")
    many = sub.add_parser("summarize")
    many.add_argument("--inputs-file", required=True, type=Path,
                      help="one per-sample provenance TSV per line")
    many.add_argument("--table", required=True, type=Path)
    many.add_argument("--summary", required=True, type=Path)
    a = ap.parse_args()

    if a.cmd == "sample":
        if a.min_as < 1:
            sys.exit(f"ERROR: {a.sample}: --min-as must be an integer >= 1 (got {a.min_as})")
        try:
            row = parse_header(a.header.read_text(errors="replace"), a.sample,
                               a.low_score_mode, a.min_as)
        except ValueError as exc:
            sys.exit(f"ERROR: {a.sample}: {exc}")
        write_rows(a.output, [row], FIELDS)
        if a.decision:
            state = row["extract_low_score_mapped"]
            reason = low_score_reason(state, row["aligner_bwa_mem"], row["bwa_min_score"], a.min_as)
            a.decision.write_text(f"{state.split('_', 1)[0]}\t{state}\t{reason}\n")
        return

    rows = []
    for name in a.inputs_file.read_text().split():
        with open(name, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames != FIELDS:
                sys.exit(f"ERROR: {name}: unexpected header {reader.fieldnames}")
            rows.extend(reader)
    samples = [r["sample"] for r in rows]
    duplicated = sorted({s for s in samples if samples.count(s) > 1})
    if duplicated:
        sys.exit(f"ERROR: more than one provenance row for {len(duplicated)} sample(s)")
    rows.sort(key=lambda r: r["sample"])
    write_rows(a.table, rows, FIELDS)
    write_rows(a.summary, summarize(rows), SUMMARY_FIELDS)


if __name__ == "__main__":
    main()
