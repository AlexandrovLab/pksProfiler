#!/usr/bin/env python3
"""Assert on the real output of tests/check_metagenome_mag_execution.sh's Nextflow run.

Checks the specific regression this session's real execution found and fixed:
Modules/pks_mag.nf's summary_input_ch used a strict join on checkm2Predict/
gtdbtkClassify output, both of which only ever run on real bins -- so a sample with
ZERO real bins (this fixture's whole point) never got a genomes/pks_mag_summary.tsv
row at all, even carrying a real, positive call from U1's own unbinned-pool detection.
stubBinlessMagQuality fixes it by giving such a sample a stand-in pair so the join
keeps it.

Exits 1 with a message naming the first failed check; exits 0 and prints one "ok" line
per check on success.
"""
import argparse
import csv
import sys
from pathlib import Path

FAILURES = []


def check(label, condition, detail=""):
    if condition:
        print(f"  ok    {label}")
    else:
        message = f"  FAIL  {label}" + (f" ({detail})" if detail else "")
        print(message)
        FAILURES.append(label)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--mag-sample", default="e2e_metagenome_mag")
    args = parser.parse_args()

    results = args.results
    sample = args.mag_sample
    by_sample = results / "by_sample" / sample

    # ---------- confirm the fixture actually landed in the scenario this test covers ----------
    mag_status_path = by_sample / "genomes" / "mag_status.tsv"
    check("mag_status.tsv exists", mag_status_path.exists())
    if mag_status_path.exists():
        row = read_tsv(mag_status_path)[0]
        check("real megahit assembly produced at least one contig",
              int(row["contig_count"]) > 0, f"got {row['contig_count']!r}")
        check("MetaBAT2 recovered ZERO real bins (fixture is sized to stay under "
              "MetaBAT2's 200 kb --minClsSize floor -- this is the scenario the fix covers, "
              "not an assumption)",
              int(row["bin_count"]) == 0, f"got {row['bin_count']!r}")

    # ---------- the regression guard itself ----------
    mag_summary_path = by_sample / "genomes" / "pks_mag_summary.tsv"
    check("genomes/pks_mag_summary.tsv exists despite zero real bins "
          "(this is exactly what stubBinlessMagQuality fixes -- before it, this file "
          "was never written at all for a zero-bin sample)",
          mag_summary_path.exists())
    if mag_summary_path.exists():
        rows = read_tsv(mag_summary_path)
        check("pks_mag_summary.tsv has at least one row", len(rows) > 0)
        unbinned_rows = [r for r in rows if r.get("unit_type") == "unbinned"]
        check("exactly one row is the unbinned-pool row (bin_id=unbinned, unit_type=unbinned)",
              len(unbinned_rows) == 1, f"got unit_types {[r.get('unit_type') for r in rows]}")
        check("no real bin rows exist (unit_type=bin) -- consistent with mag_status's bin_count=0",
              all(r.get("unit_type") != "bin" for r in rows),
              f"got {[r.get('unit_type') for r in rows]}")
        if unbinned_rows:
            row = unbinned_rows[0]
            check("unbinned row's bin_id is literally 'unbinned'",
                  row.get("bin_id") == "unbinned", f"got {row.get('bin_id')!r}")
            # locus_genes_detected, not clb_genes_detected: the unbinned pool never runs
            # through hmmsearchClb (that only runs on prokka-annotated real bins), so the
            # HMM-derived clb_genes_detected/clb_genes/best_evalue columns are correctly
            # NA for this row -- its real evidence is magBinLocusEvidence's alignment-based
            # locus_genes_detected/locus_tier/locus_breadth, which DOES run on the unbinned
            # pool by design (see locus_alignment_units_ch in Modules/pks_mag.nf).
            check("unbinned row carries real clb gene evidence (fixture floods the whole "
                  "pks island, all 19 genes) via locus_genes_detected",
                  int(row.get("locus_genes_detected") or 0) >= 10,
                  f"got {row.get('locus_genes_detected')!r}")
            check("unbinned row's locus tier clears at least multi_gene (real alignment "
                  "evidence, not just an HMM hit)",
                  row.get("locus_tier") in ("multi_gene", "broad_island", "extensive_island"),
                  f"got {row.get('locus_tier')!r}")
            # build_mag_summary.py deliberately writes "NA" (a string) for every
            # bin-specific field on the unbinned row -- not 0 or unclassified -- since
            # none of those were measured for a pool of contigs that is not a genome.
            check("unbinned row has no real CheckM2/GTDB-Tk data (stand-in, not a real bin)",
                  row.get("taxonomy") == "NA" and row.get("completeness") == "NA",
                  f"got taxonomy={row.get('taxonomy')!r} completeness={row.get('completeness')!r}")

    # ---------- no real checkm2/gtdbtk directories: neither tool ever ran ----------
    check("no real checkm2/ directory was published (no real bin ever existed for it to run on)",
          not (by_sample / "genomes" / "checkm2").exists())
    check("no real gtdbtk/ directory was published (no real bin ever existed for it to run on)",
          not (by_sample / "genomes" / "gtdbtk").exists())

    # ---------- cohort roll-up still includes this sample ----------
    master_path = results / "cohort" / "pks.master_summary.tsv"
    check("cohort/pks.master_summary.tsv exists", master_path.exists())
    if master_path.exists():
        rows = {row["sample"]: row for row in read_tsv(master_path)}
        check(f"{sample} has a master-summary row", sample in rows)
        if sample in rows:
            row = rows[sample]
            check(f"{sample}'s row carries MAG columns (mag_bins_total present, not NA)",
                  row.get("mag_bins_total", "NA") != "NA", f"got {row.get('mag_bins_total')!r}")
            check(f"{sample}'s mag_bins_total is 0 (the unbinned pool is not counted as a bin)",
                  row.get("mag_bins_total") == "0", f"got {row.get('mag_bins_total')!r}")
            # The whole point: this sample's ONLY MAG-lane signal is the unbinned pool, and
            # build_master_summary.py's own U1 columns for it (mag_unbinned_locus_tier /
            # mag_unbinned_pks_positive) can only ever be populated if read_mags finds
            # genomes/pks_mag_summary.tsv at all -- which required this session's fix.
            check(f"{sample}'s row carries the U1 unbinned-pool columns, not NA",
                  row.get("mag_unbinned_locus_tier", "NA") != "NA",
                  f"got {row.get('mag_unbinned_locus_tier')!r}")
            check(f"{sample}'s unbinned pool is reported positive",
                  row.get("mag_unbinned_pks_positive") == "yes",
                  f"got {row.get('mag_unbinned_pks_positive')!r}")

    print()
    if FAILURES:
        print(f"metagenome MAG assertions: {len(FAILURES)} FAILED: {', '.join(FAILURES)}",
              file=sys.stderr)
        return 1
    print("metagenome MAG assertions: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
