#!/usr/bin/env python3
"""Route FASTQ reads independently using their exact Kraken read IDs."""
import argparse, gzip
from pathlib import Path

def normalize_read_id(value):
    """Return the exact read ID; /1 and /2 remain independent evidence units."""
    return value.split()[0].lstrip("@")

def iter_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            record = [handle.readline() for _ in range(4)]
            if not record[0]: return
            if any(x == "" for x in record): raise ValueError(f"Truncated FASTQ: {path}")
            if not record[0].startswith("@") or not record[2].startswith("+"):
                raise ValueError(f"Malformed FASTQ record: {record[0].rstrip()}")
            yield record

def read_parents(path):
    parents = {}
    with Path(path).open() as handle:
        for line in handle:
            delimiter = "|" if "|" in line else "\t"
            fields = [x.strip() for x in line.split(delimiter)]
            if len(fields) >= 2: parents[fields[0]] = fields[1]
    if not parents: raise ValueError(f"No taxonomy nodes found: {path}")
    return parents

def is_descendant(taxid, target, parents, memo):
    key = (taxid, target)
    if key in memo: return memo[key]
    seen = set()
    while taxid and taxid not in seen:
        if taxid == target:
            memo[key] = True; return True
        seen.add(taxid)
        parent = parents.get(taxid)
        if not parent or parent == taxid: break
        taxid = parent
    memo[key] = False; return False

def selected_reads(kraken_path, target, parents):
    selected, memo = set(), {}
    with Path(kraken_path).open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3: raise ValueError(f"Malformed Kraken row {number}")
            if fields[0] == "C" and is_descendant(fields[2].split()[0], target, parents, memo):
                selected.add(normalize_read_id(fields[1]))
    return selected

def route(reads, kraken, nodes, target, primary_path, non_target_path):
    selected = selected_reads(kraken, str(target), read_parents(nodes))
    counts = {"input": 0, "primary": 0, "non_target": 0}
    with gzip.open(primary_path, "wt") as primary, gzip.open(non_target_path, "wt") as other:
        for record in iter_fastq(reads):
            counts["input"] += 1
            chosen = normalize_read_id(record[0]) in selected
            (primary if chosen else other).writelines(record)
            counts["primary" if chosen else "non_target"] += 1
    return counts

def main():
    p = argparse.ArgumentParser()
    for name in ("reads", "kraken-output", "nodes", "primary", "non-target", "qc", "sample"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--target-taxid", required=True, type=int)
    a = p.parse_args()
    c = route(a.reads, a.kraken_output, a.nodes, a.target_taxid, a.primary, a.non_target)
    with Path(a.qc).open("w") as out:
        out.write("Sample\tMetric\tValue\n")
        for metric in ("input", "primary", "non_target"):
            out.write(f"{a.sample}\tprefilter_{metric}_reads\t{c[metric]}\n")
if __name__ == "__main__": main()
