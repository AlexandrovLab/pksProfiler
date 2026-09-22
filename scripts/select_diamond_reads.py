#!/usr/bin/env python3
"""Select DIAMOND-supported fragments and both mates."""
import argparse, gzip
from pathlib import Path
from route_kraken_reads import iter_fastq, normalize_read_id

def load_hits(path, max_evalue, min_aa, min_query_coverage):
    selected = set()
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) != 8: raise ValueError(f"Malformed DIAMOND row {number}")
            query, _, _, length, qlen, _, evalue, _ = f
            coverage = min(1.0, int(length) * 3 / int(qlen))
            if float(evalue) <= max_evalue and int(length) >= min_aa and coverage >= min_query_coverage:
                selected.add(normalize_read_id(query))
    return selected

def main():
    p = argparse.ArgumentParser()
    for name in ("reads", "diamond", "output", "qc", "sample"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--max-evalue", type=float, default=1e-5)
    p.add_argument("--min-aa", type=int, default=25)
    p.add_argument("--min-query-coverage", type=float, default=0.5)
    a = p.parse_args(); selected = load_hits(a.diamond, a.max_evalue, a.min_aa, a.min_query_coverage)
    input_reads = rescued = 0
    with gzip.open(a.output, "wt") as out:
        for record in iter_fastq(a.reads):
            input_reads += 1
            if normalize_read_id(record[0]) in selected:
                out.writelines(record); rescued += 1
    with Path(a.qc).open("w") as out:
        out.write("Sample\tMetric\tValue\n")
        out.write(f"{a.sample}\tdiamond_input_reads\t{input_reads}\n")
        out.write(f"{a.sample}\tdiamond_rescue_fragments\t{len(selected)}\n")
        out.write(f"{a.sample}\tdiamond_rescued_reads\t{rescued}\n")
if __name__ == "__main__": main()
