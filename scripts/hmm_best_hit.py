#!/usr/bin/env python3
"""Assign each read to one clb model from an nhmmscan --tblout, or to none.

F07. This replaces an awk filter with three defects:

  * `eval=$13+0` coerced a non-numeric E-value to 0, and 0 passes any stringency
    threshold, so a corrupt field became the most significant hit in the file
  * `NF < 14 { next }` dropped truncated rows without saying so
  * ties were broken by file order: when two genes matched with the same E-value and
    the same score, whichever line came first won, so reversing two rows turned
    clbA=1 into clbB=1. With --hmm_chunking the same sample can present its hits in a
    different order, which made the counts depend on how the run was executed

The rule, stated once:

  1. lowest E-value
  2. then highest bit score
  3. then longest aligned length on the read -- the same idea as featureCounts'
     --largestOverlap: the model that matched more of the read wins
  4. still identical across two different genes: the read is ambiguous and counts for
     nothing, which is what featureCounts does with a perfect overlap tie

Several hits to the same gene are not a tie; the best one stands. Ambiguity is
reported, not inferred from a gap in the counts.

nhmmscan --tblout columns, from the header of a generated file: 1 target name,
2 accession, 3 query name, 4 accession, 5 hmmfrom, 6 hmm to, 7 alifrom, 8 ali to,
9 envfrom, 10 env to, 11 modlen, 12 strand, 13 E-value, 14 score, 15 bias,
16 description.
"""
import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from clb_counts_schema import CLB_GENES, CLB_GENE_SET, validate_clb_counts

TARGET, QUERY, ALI_FROM, ALI_TO, EVALUE, SCORE = 0, 2, 6, 7, 12, 13
MINIMUM_FIELDS = 15


class TbloutError(ValueError):
    """A tblout line that cannot be trusted into a count."""


def model_name(target):
    """The clb gene a target names, with the suffix the models carry stripped."""
    return target.replace(".cds.aln", "")


def parse_row(line, line_number, source):
    fields = line.split()
    if len(fields) < MINIMUM_FIELDS:
        raise TbloutError(
            f"{source}:{line_number}: {len(fields)} fields, expected at least "
            f"{MINIMUM_FIELDS}. A truncated row is a damaged file, not a row to skip.")

    gene = model_name(fields[TARGET])
    if gene not in CLB_GENE_SET:
        raise TbloutError(
            f"{source}:{line_number}: {fields[TARGET]!r} is not a clb model. "
            "Expected clbA-clbS.")

    numbers = {}
    for label, index in (("E-value", EVALUE), ("score", SCORE),
                         ("alifrom", ALI_FROM), ("ali to", ALI_TO)):
        try:
            numbers[label] = float(fields[index])
        except ValueError:
            raise TbloutError(
                f"{source}:{line_number}: {label} is not a number: "
                f"{fields[index]!r}. Coercing it to zero would make it the most "
                "significant hit in the file.")
        if not math.isfinite(numbers[label]):
            raise TbloutError(
                f"{source}:{line_number}: {label} is not finite: {fields[index]!r}")

    # nhmmscan reports coordinates on the read, and reverse-strand hits count down.
    aligned = abs(numbers["ali to"] - numbers["alifrom"]) + 1
    return fields[QUERY], gene, numbers["E-value"], numbers["score"], aligned


def best_hits(path, max_evalue):
    """
    {read: gene} for reads with one clear winner, plus the reads that tie.

    A read is ranked on (E-value ascending, score descending, aligned length
    descending). Two different genes reaching the identical triple leave the read
    ambiguous.
    """
    best = {}          # read -> (evalue, -score, -aligned)
    winners = {}       # read -> set of genes holding that rank
    qualifying = set()

    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            read, gene, evalue, score, aligned = parse_row(line, line_number, path)
            if evalue > max_evalue:
                continue
            qualifying.add(read)

            rank = (evalue, -score, -aligned)
            current = best.get(read)
            if current is None or rank < current:
                best[read] = rank
                winners[read] = {gene}
            elif rank == current:
                winners[read].add(gene)

    assigned = {read: next(iter(genes)) for read, genes in winners.items()
                if len(genes) == 1}
    ambiguous = sorted(read for read, genes in winners.items() if len(genes) > 1)
    return assigned, ambiguous, sorted(qualifying)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tblout", required=True)
    parser.add_argument("--evalue", type=float, required=True)
    parser.add_argument("--counts-out", required=True)
    parser.add_argument("--read-ids-out", required=True)
    parser.add_argument("--ambiguous-out", default=None,
                        help="read ids that tied across genes; counted either way")
    args = parser.parse_args()

    try:
        assigned, ambiguous, qualifying = best_hits(args.tblout, args.evalue)
    except TbloutError as problem:
        raise SystemExit(f"[ERROR] {problem}")

    counts = {gene: 0 for gene in CLB_GENES}
    for gene in assigned.values():
        counts[gene] += 1

    # The same contract the cohort matrix holds every sample to (F06).
    validated = validate_clb_counts(counts.items(), args.counts_out)

    with open(args.counts_out, "w") as handle:
        handle.write("Gene\tCount\n")
        for gene in CLB_GENES:
            handle.write(f"{gene}\t{validated[gene]}\n")

    with open(args.read_ids_out, "w") as handle:
        for read in qualifying:
            handle.write(f"{read}\n")

    if args.ambiguous_out:
        with open(args.ambiguous_out, "w") as handle:
            for read in ambiguous:
                handle.write(f"{read}\n")

    print(f"{len(assigned)} reads assigned, {len(ambiguous)} ambiguous, "
          f"{len(qualifying)} with a qualifying hit")


if __name__ == "__main__":
    main()
