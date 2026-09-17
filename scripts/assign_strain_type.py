#!/usr/bin/env python3
"""Join an mlst call to clonal complex and phylogroup from the reference panel.

The panel is 45,761 E. coli genomes (Lawley/Maklin) carrying ST, CC and phylogroup
labels; see ref/typing/st_phylogroup.tsv. pks-positive genomes in that panel are 92%
phylogroup B2, so a confident non-B2 pks-positive assignment is unusual and is flagged
rather than passed through silently.

No call is emitted from too little sequence: MLST needs all seven housekeeping loci, and
a partial allele profile is reported as insufficient rather than guessed at.
"""
import argparse
import csv
import re
from pathlib import Path

MLST_LOCI = 7
MISSING = {"-", "", "?", "NA"}
# Dominant phylogroup among pks+ genomes in the reference panel (4386/4760).
EXPECTED_PKS_PHYLOGROUP = "B2"


def read_lookup(path):
    with Path(path).open() as handle:
        return {row["ST"]: row for row in csv.DictReader(handle, delimiter="\t")}


ALLELE_VALUE = re.compile(r"\(([^)]*)\)")


def _allele_value(field):
    """mlst writes locus(value). Return the value, never the locus name.

    Notation: `n` exact allele, `~n` novel allele close to n, `n?` partial or
    low-identity match, `-` absent. Splitting on "(" and testing the left side
    tests the locus name, which is never a missing marker, so every locus counted
    as recovered -- including ones where nothing was found.
    """
    match = ALLELE_VALUE.search(field)
    return (match.group(1) if match else field).strip()


def parse_mlst(path):
    """mlst writes one row per input: file, scheme, ST, then allele columns."""
    text = Path(path).read_text().strip()
    if not text:
        return None, None, 0, 0, 0
    fields = text.splitlines()[0].split("\t")
    if len(fields) < 3:
        return None, None, 0, 0, 0
    scheme, st = fields[1].strip(), fields[2].strip()
    exact = inexact = missing = 0
    for field in fields[3:]:
        value = _allele_value(field)
        if not value or value in MISSING:
            missing += 1
        elif value.startswith("~") or value.endswith("?"):
            inexact += 1
        else:
            exact += 1
    return (scheme or None), (st if st not in MISSING else None), exact, inexact, missing


def classify(st, exact):
    # Only an exact allele at every locus can yield an ST. A profile of partial or
    # novel-adjacent calls is not a novel strain type, it is an untypeable assembly.
    if exact < MLST_LOCI:
        return "insufficient_loci"
    if st is None:
        return "novel_allele_combination"
    return "typed"


def main():
    p = argparse.ArgumentParser()
    for name in ("mlst", "lookup", "sample", "label", "output"):
        p.add_argument(f"--{name}", required=True)
    a = p.parse_args()

    scheme, st, exact, inexact, missing = parse_mlst(a.mlst)
    status = classify(st, exact)
    lookup = read_lookup(a.lookup)
    hit = lookup.get(st) if st else None

    record = {
        "sample": a.sample,
        "unit": a.label,
        "status": status,
        "scheme": scheme or "NA",
        "loci_called": exact,
        "loci_inexact": inexact,
        "loci_missing": missing,
        "ST": st or "NA",
        "clonal_complex": "NA",
        "phylogroup": "NA",
        "phylogroup_consistency": "NA",
        "panel_genomes": "NA",
        "panel_pks_fraction": "NA",
        "note": "",
    }

    if status != "typed":
        record["note"] = (
            f"{exact}/{MLST_LOCI} exact MLST alleles ({inexact} partial or novel, "
            f"{missing} absent) — too little sequence to type"
            if status == "insufficient_loci"
            else "all seven loci recovered but the combination is not in the scheme"
        )
    elif hit is None:
        record["note"] = f"ST{st} is not represented in the reference panel"
    else:
        record.update(
            clonal_complex=hit["clonal_complex"],
            phylogroup=hit["phylogroup"],
            phylogroup_consistency=hit["phylogroup_consistency"],
            panel_genomes=hit["n_genomes"],
            panel_pks_fraction=hit["pks_fraction"],
        )
        notes = []
        if hit["phylogroup"] not in (EXPECTED_PKS_PHYLOGROUP, "NA"):
            notes.append(
                f"phylogroup {hit['phylogroup']} is atypical for pks+ "
                f"({EXPECTED_PKS_PHYLOGROUP} in 92% of panel pks+ genomes) — verify before reporting"
            )
        if float(hit["pks_fraction"]) == 0.0:
            notes.append("no pks+ genome of this ST in the panel")
        if hit["phylogroup_consistency"] not in ("NA",) and float(hit["phylogroup_consistency"]) < 0.9:
            notes.append(f"ST{st} phylogroup is not consistent in the panel ({hit['phylogroup_consistency']})")
        record["note"] = "; ".join(notes)

    with Path(a.output).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(record), delimiter="\t")
        writer.writeheader()
        writer.writerow(record)


if __name__ == "__main__":
    main()
