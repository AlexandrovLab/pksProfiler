"""Shared parsing utilities for MAG-module scripts."""

# M1: clbA, clbD, clbP and clbQ are small, single-domain tailoring genes that showed
# no cross-reactivity against the promiscuous megasynthase domains (clbB, clbC, clbH,
# clbI, clbJ, clbK, clbN, clbO) that a bare hmmsearch E-value cut lets any bacterium
# hit (v0.0.2_functional_test_20260910, ERR525841: 7 of 8 bins, including two
# Bifidobacterium bins, called pks-positive on megasynthase hits alone). A bin must
# carry at least one of these on top of the gene-count floor to be called positive.
# Kept in sync with main.nf's params.mag_specific_clb_genes default -- see
# tests/test_pks_positive_definition.py.
SPECIFIC_CLB = {"clbA", "clbD", "clbP", "clbQ"}


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
