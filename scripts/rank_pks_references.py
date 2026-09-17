#!/usr/bin/env python3
"""Rank PKS references from contig-to-panel minimap2 PAF alignments."""
import argparse, csv, re
from collections import defaultdict
from pathlib import Path

SUMMARY_COLUMNS = ["sample", "assembler", "best_reference_path",
 "best_reference_phylogroup", "best_reference_sequence_type",
 "alignment_coverage", "alignment_identity", "runner_up_references",
 "score_difference_from_best", "top_match_consensus_phylogroup",
 "assignment_status"]

def reference_metadata(path):
    result = {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.startswith(">"): continue
            parts = line[1:].strip().split(maxsplit=1)
            ref, header = parts[0], parts[1] if len(parts) > 1 else ""
            values = dict(re.findall(r"(?:^|\|)([^|=]+)=([^|]+)", header))
            result[ref] = {"path": header.split("|", 1)[0] or ref,
                           "phylogroup": values.get("phylogroup", "NA"),
                           "ST": values.get("ST", "NA")}
    return result

def merged_length(intervals):
    total, end = 0, -1
    for start, stop in sorted(intervals):
        if stop <= end: continue
        total += stop - max(start, end)
        end = stop
    return total

def paf_rankings(path, metadata):
    hits = defaultdict(lambda: {"length": 0, "intervals": [], "matches": 0, "block": 0})
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip(): continue
            f = line.rstrip().split("\t")
            if len(f) < 12: raise ValueError("PAF line has fewer than 12 fields")
            hit = hits[f[5]]; hit["length"] = int(f[6])
            hit["intervals"].append((int(f[7]), int(f[8])))
            hit["matches"] += int(f[9]); hit["block"] += int(f[10])
    rows = []
    for ref, hit in hits.items():
        aligned = merged_length(hit["intervals"])
        coverage = aligned / hit["length"] if hit["length"] else 0
        identity = hit["matches"] / hit["block"] if hit["block"] else 0
        meta = metadata.get(ref, {"path": ref, "phylogroup": "NA", "ST": "NA"})
        rows.append({"reference": ref, **meta, "reference_length": hit["length"],
          "aligned_reference_bases": aligned, "coverage": coverage,
          "identity": identity, "score": coverage * identity})
    return sorted(rows, key=lambda x: (-x["score"], -x["coverage"], -x["identity"], x["reference"]))

def summarize(sample, assembler, rankings, top_n=5, near_best_fraction=.98,
              min_coverage=.5, min_identity=.9, min_score_delta=.02):
    if not rankings:
        return dict(zip(SUMMARY_COLUMNS, [sample, assembler, "NA", "NA", "NA",
          "0.000000", "0.000000", "", "0.000000", "NA", "unresolved"]))
    best = rankings[0]; runner = rankings[1] if len(rankings) > 1 else None
    top = [x for x in rankings[:top_n] if x["score"] >= best["score"] * near_best_fraction]
    groups = {x["phylogroup"] for x in top if x["phylogroup"] not in ("", "NA")}
    consensus = next(iter(groups)) if len(groups) == 1 else "ambiguous" if groups else "NA"
    delta = best["score"] - (runner["score"] if runner else 0)
    other = next((x for x in rankings[1:] if x["phylogroup"] != best["phylogroup"]), None)
    cross_delta = best["score"] - (other["score"] if other else 0)
    if best["coverage"] < min_coverage or best["identity"] < min_identity: status = "unresolved"
    elif consensus == best["phylogroup"] and cross_delta >= min_score_delta: status = "confident"
    else: status = "ambiguous"
    return {"sample": sample, "assembler": assembler, "best_reference_path": best["path"],
      "best_reference_phylogroup": best["phylogroup"], "best_reference_sequence_type": best["ST"],
      "alignment_coverage": f'{best["coverage"]:.6f}', "alignment_identity": f'{best["identity"]:.6f}',
      "runner_up_references": ";".join(x["reference"] for x in rankings[1:top_n]),
      "score_difference_from_best": f"{delta:.6f}",
      "top_match_consensus_phylogroup": consensus, "assignment_status": status}

def main():
    p=argparse.ArgumentParser()
    for name in ("paf","references","summary","rankings"): p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--sample",required=True); p.add_argument("--assembler",required=True,choices=("megahit","metaspades"))
    p.add_argument("--top-n",type=int,default=5); p.add_argument("--near-best-fraction",type=float,default=.98)
    p.add_argument("--min-coverage",type=float,default=.5); p.add_argument("--min-identity",type=float,default=.9)
    p.add_argument("--min-score-delta",type=float,default=.02); a=p.parse_args()
    rows=paf_rankings(a.paf,reference_metadata(a.references))
    summary=summarize(a.sample,a.assembler,rows,a.top_n,a.near_best_fraction,a.min_coverage,a.min_identity,a.min_score_delta)
    cols=["rank","reference","path","phylogroup","ST","reference_length","aligned_reference_bases","coverage","identity","score"]
    with a.rankings.open("w",newline="") as h:
        w=csv.DictWriter(h,cols,delimiter="\t"); w.writeheader()
        for i,row in enumerate(rows,1): w.writerow({"rank":i,**row,"coverage":f'{row["coverage"]:.6f}',"identity":f'{row["identity"]:.6f}',"score":f'{row["score"]:.6f}'})
    with a.summary.open("w",newline="") as h:
        w=csv.DictWriter(h,SUMMARY_COLUMNS,delimiter="\t"); w.writeheader(); w.writerow(summary)
if __name__ == "__main__": main()
