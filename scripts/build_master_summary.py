#!/usr/bin/env python3
"""Join every per-sample result into one table.

pksProfiler writes each stage's results separately, which is right for the pipeline and
wrong for reading. This walks a finished (or partly finished) results directory and emits
one row per sample carrying whatever each stage produced, with `NA` where a stage did not
run. Stages that were never enabled simply contribute no columns.

Run it after a pipeline run:

    python3 scripts/build_master_summary.py --results results \\
        --output results/cohort/pks.master_summary.tsv

It only reads published output, so it is safe to re-run at any time and does not require
re-running the pipeline.
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

# Column groups, in the order they appear. Each group is skipped entirely when the stage
# that produces it left no files behind.
QC_COLUMNS = ["input_alignment_records", "total_primary_reads", "extracted_unmapped_reads",
              "filter_input_reads", "reads_after_fastp", "reads_after_hg38", "reads_after_t2t_phix",
              "reads_mapped_ihe3034", "num_clb_genes_align", "reads_clb_genes_align",
              "num_clb_genes_hmm", "reads_clb_genes_hmm"]
EVIDENCE_COLUMNS = ["read_evidence", "pks_reads", "clb_genes_detected", "island_breadth_1x"]
CONTIG_COLUMNS = ["final_structural_evidence", "assembler_agreement"]
MAG_COLUMNS = ["mag_bins_total", "mag_bins_pks_positive", "pks_mag_taxonomy",
               "pks_mag_completeness", "pks_mag_clb_genes"]
COMMUNITY_COLUMNS = ["community_prophages_total", "community_pks_producers",
                     "community_neighbours_assessed", "community_neighbours_with_prophage",
                     "island_nearby_integrase", "island_nearby_trna", "island_in_prophage"]
TYPING_COLUMNS = ["typing_status", "ST", "clonal_complex", "phylogroup"]

TRUE = {"true", "yes", "1", "t"}


def rows(path, delim="\t"):
    try:
        with open(path) as fh:
            yield from csv.DictReader(fh, delimiter=delim)
    except (OSError, csv.Error) as exc:
        print(f"  warning: could not read {path}: {exc}", file=sys.stderr)


def first(path):
    for row in rows(path):
        return row
    return {}


def sample_of(path, row=None):
    """The sample a published file belongs to.

    Files under by_sample/<sample>/ no longer carry a <sample>. filename prefix, because
    the directory already says which sample they are. Every published table also has a
    `sample` column; prefer it when present, since it survives a file being copied out.
    """
    if row and row.get("sample"):
        return row["sample"]
    parts = os.path.normpath(path).split(os.sep)
    if "by_sample" in parts:
        i = parts.index("by_sample")
        if i + 1 < len(parts):
            return parts[i + 1]
    return None


def read_qc(results, out):
    path = os.path.join(results, "cohort/qc/pks.qc.summary.tsv")
    if not os.path.exists(path):
        return False
    for row in rows(path):
        sample = row.get("Sample")
        if sample:
            for c in QC_COLUMNS:
                out[sample][c] = row.get(c, "NA")
    return True


def read_evidence(results, out):
    found = False
    for path in glob.glob(f"{results}/by_sample/*/read_evidence.tsv"):
        row = first(path)
        if row.get("sample"):
            found = True
            for c in EVIDENCE_COLUMNS:
                out[row["sample"]][c] = row.get(c, "NA")
    return found


def read_contigs(results, out):
    found = False
    for path in glob.glob(f"{results}/by_sample/*/contigs/final_evidence/final_pks_evidence.tsv"):
        row = first(path)
        if row.get("sample"):
            found = True
            for c in CONTIG_COLUMNS:
                out[row["sample"]][c] = row.get(c, "NA")
    return found


def read_mags(results, out):
    found = False
    for path in glob.glob(f"{results}/by_sample/*/genomes/pks_mag_summary.tsv"):
        bins = list(rows(path))
        if not bins:
            continue
        sample = sample_of(path, bins[0])
        if not sample:
            continue
        found = True
        def genes(b):
            try:
                return int(b.get("distinct_clb_genes") or b.get("clb_genes") or 0)
            except ValueError:
                return 0
        positive = [b for b in bins if genes(b) > 0]
        best = max(positive, key=genes) if positive else None
        out[sample]["mag_bins_total"] = str(len(bins))
        out[sample]["mag_bins_pks_positive"] = str(len(positive))
        out[sample]["pks_mag_taxonomy"] = (best or {}).get("classification") or (best or {}).get("host_taxonomy") or "NA"
        out[sample]["pks_mag_completeness"] = (best or {}).get("completeness", "NA")
        out[sample]["pks_mag_clb_genes"] = str(genes(best)) if best else "0"
    return found


def read_community(results, out):
    found = False
    for path in glob.glob(f"{results}/by_sample/*/community/community_prophage_inventory.tsv"):
        sample = sample_of(path)
        if not sample:
            continue
        found = True
        out[sample]["community_prophages_total"] = str(len(list(rows(path))))
    for path in glob.glob(f"{results}/by_sample/*/community/pks_community_interactions.tsv"):
        pairs = list(rows(path))
        sample = sample_of(path, pairs[0] if pairs else None)
        if not sample:
            continue
        found = True
        out[sample]["community_pks_producers"] = str(len({p.get("pks_producer_bin") for p in pairs if p.get("pks_producer_bin")}))
        out[sample]["community_neighbours_assessed"] = str(len(pairs))
        def has_prophage(p):
            try:
                return int(p.get("recipient_prophage_count") or 0) > 0
            except ValueError:
                return False
        out[sample]["community_neighbours_with_prophage"] = str(sum(1 for p in pairs if has_prophage(p)))
    for path in glob.glob(f"{results}/by_sample/*/community/pks_island_mobility.tsv"):
        mob = list(rows(path))
        sample = sample_of(path, mob[0] if mob else None)
        if not sample:
            continue
        found = True
        anyof = lambda col: "yes" if any((m.get(col) or "").strip().lower() in TRUE for m in mob) else "no"
        out[sample]["island_nearby_integrase"] = anyof("nearby_integrase")
        out[sample]["island_nearby_trna"] = anyof("nearby_trna")
    # tumour contigs report the same three flags per clb hit
    for path in glob.glob(f"{results}/by_sample/*/community/contig_pks_mobility.tsv"):
        hits = list(rows(path))
        if not hits:
            continue
        sample = sample_of(path, hits[0])
        if not sample:
            continue
        found = True
        anyof = lambda col: "yes" if any((h.get(col) or "").strip().lower() in TRUE for h in hits) else "no"
        out[sample]["island_nearby_integrase"] = anyof("has_integrase")
        out[sample]["island_nearby_trna"] = anyof("nearby_trna")
        out[sample]["island_in_prophage"] = anyof("in_prophage")
    return found


def read_typing(results, out):
    found = False
    for path in glob.glob(f"{results}/by_sample/*/strain/strain_types.tsv"):
        units = list(rows(path))
        if not units:
            continue
        sample = sample_of(path, units[0])
        if not sample:
            continue
        found = True
        typed = [u for u in units if u.get("status") == "typed"]
        pick = typed[0] if typed else units[0]
        out[sample]["typing_status"] = pick.get("status", "NA")
        for c in ("ST", "clonal_complex", "phylogroup"):
            out[sample][c] = pick.get(c, "NA")
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="a pksProfiler results directory")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    out = defaultdict(dict)
    stages = [
        ("read profiling", QC_COLUMNS, read_qc),
        ("evidence tier", EVIDENCE_COLUMNS, read_evidence),
        ("contig reassembly", CONTIG_COLUMNS, read_contigs),
        ("draft genomes", MAG_COLUMNS, read_mags),
        ("community context", COMMUNITY_COLUMNS, read_community),
        ("strain typing", TYPING_COLUMNS, read_typing),
    ]
    columns = ["sample"]
    for name, cols, fn in stages:
        present = fn(a.results, out)
        print(f"  {name:20} {'found' if present else 'not run — columns omitted'}", file=sys.stderr)
        if present:
            columns += cols

    if not out:
        sys.exit(f"No per-sample results found under {a.results}")

    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    with open(a.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", restval="NA", extrasaction="ignore")
        w.writeheader()
        for sample in sorted(out):
            w.writerow({"sample": sample, **out[sample]})
    print(f"{len(out)} samples x {len(columns)} columns -> {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
