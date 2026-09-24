"""Shared parsing utilities for MAG-module scripts."""

# M1: a bin is called pks-positive by its magBinLocusEvidence tier (an alignment of
# the bin's own assembly against the canonical IHE3034 locus, scored on genes an
# alignment covers and breadth of the island it spans), not by HMM gene count alone.
# See Modules/pks_mag.nf's alignMagBinToCanonicalReference/magBinLocusEvidence and
# scripts/summarize_mag_bin_locus_evidence.py.
POSITIVE_LOCUS_TIERS = {"multi_gene", "broad_island", "extensive_island"}


def read_locus_evidence(locus_dir):
    """{bin_id: {"locus_tier": ..., "locus_genes_detected": ..., "locus_breadth": ...}}
    from a directory of `{bin_id}.locus_evidence.tsv` files (magBinLocusEvidence's
    output, one per bin). A bin absent from the directory has not been aligned.

    `bin_id` here is really "whatever unitID magBinLocusEvidence was given". Since U1,
    that directory can also carry an `unbinned.locus_evidence.tsv` -- the per-sample
    pool of contigs MetaBAT2 never placed in any bin (see Modules/pks_mag.nf's pksMAG
    workflow). This function does not distinguish the two; callers that must
    (build_mag_summary.py's unit_type column) do so explicitly.
    """
    import csv
    import glob
    import os

    result = {}
    for path in glob.glob(os.path.join(locus_dir, "*.locus_evidence.tsv")):
        bin_id = os.path.basename(path)[: -len(".locus_evidence.tsv")]
        with open(path) as fh:
            result[bin_id] = next(csv.DictReader(fh, delimiter="\t"))
    return result


def parse_hmmsearch_tblout(tblout_path, evalue_threshold=1e-5):
    """Parse an hmmsearch --tblout file, returning one dict per passing hit.

    Each dict has keys: locus_tag, clb_gene, evalue.
    Duplicate (locus_tag, clb_gene) pairs are deduplicated; the lowest
    e-value hit is kept.
    """
    best = {}  # (locus_tag, clb_gene) → lowest evalue hit dict
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 19:
                continue
            evalue = float(parts[4])
            if evalue > evalue_threshold:
                continue
            locus_tag = parts[0]
            clb_gene = parts[2]
            key = (locus_tag, clb_gene)
            if key not in best or evalue < best[key]["evalue"]:
                best[key] = {
                    "locus_tag": locus_tag,
                    "clb_gene": clb_gene,
                    "evalue": evalue,
                }
    return list(best.values())
