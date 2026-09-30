"""Ludmil, 26 Sep mode memo section 5: separate read-level taxonomy from carrier evidence.

    Kraken/Bracken = read-level taxonomic support. Reserve carrier-genome language for
    contig/MAG evidence linking pks to a reconstructed genome.

The same goes for the DIAMOND-rescue taxonomy table, whose taxon is krakenuniq's call on
each rescued read. So wherever the docs, the cohort report or a taxonomy script talk about
those read-level outputs, "carrier" / "carries it" language must not appear unless it is
being explicitly disclaimed. Genome-bin and contig text is free to use it.

And Ludmil, 25 Sep report section 5: NM is a divergence metric relative to IHE3034; it is
never presented as specificity, uniqueness or carrier identity.
"""
import importlib.util
import re
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Read-level taxonomy outputs and the tools that make them.
READ_LEVEL = re.compile(
    r"kraken|bracken|diamond[ _-]?rescue|diamond_rescue_taxonomy|--pks_taxa|pks_community_taxa"
    r"|clb_species_support|species_support|clb_species_reported|taxonomy_status",
    re.I)
# "carrier(s)", or carry/carries/carrying the island / pks / clb / it -- not every
# "carries" (a taxDB "carries" a pseudo-taxid).
CARRIER = re.compile(r"\bcarriers?\b|\bcarr(?:y|ies|ying)\s+(?:it\b|the (?:pks )?island|"
                     r"(?:the )?pks|clb)|\bharbou?r", re.I)
# A unit may name both only to say the one is not the other.
DISCLAIMED = re.compile(
    r"not (?:a |the )?carrier|not carrier attribution|never a carrier|carrier[- ]genome language is"
    r" reserved|not the organism shown to carry|not proof that a species carries",
    re.I)
# NM / divergence text must not claim more than divergence.
DIVERGENCE = re.compile(r"\bNM\b|divergence_ihe3034|ihe3034_divergence", re.I)
OVERCLAIM = re.compile(r"\b(?:specificity|uniqueness|unique|carrier)\b", re.I)
NEGATED = re.compile(r"\bnot\b|\bnever\b|\bnothing more\b|\bno\b", re.I)


def corpus():
    """README, docs, and the workflow and script source (comments, labels, report text)."""
    files = [ROOT / "README.md", *sorted((ROOT / "docs").rglob("*.md")), ROOT / "main.nf",
             *sorted((ROOT / "Modules").glob("*.nf")), *sorted((ROOT / "scripts").glob("*.py"))]
    return [f for f in files if f.exists()]


def units(path):
    """Paragraphs, with every Markdown table row its own unit (one row, one claim)."""
    out, current = [], []
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        stripped = line.strip()
        if stripped.startswith("|"):
            if current:
                out.append(current)
                current = []
            out.append([(number, line)])
        elif not stripped or stripped in ('"""', "#", "//"):
            if current:
                out.append(current)
                current = []
        else:
            current.append((number, line))
    if current:
        out.append(current)
    return [(block[0][0], " ".join(text for _n, text in block)) for block in out]


def offenders(paths, subject, claim, allowed):
    found = []
    for path in paths:
        for line, text in units(path):
            if subject.search(text) and claim.search(text) and not allowed.search(text):
                where = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path.name
                found.append(f"{where}:{line}: {text.strip()[:160]}")
    return found


class ReadLevelTaxonomyIsNotCarrierEvidence(unittest.TestCase):
    def test_no_carrier_language_beside_read_level_taxonomy(self):
        found = offenders(corpus(), READ_LEVEL, CARRIER, DISCLAIMED)
        self.assertEqual(found, [], "\ncarrier language beside read-level taxonomy "
                                    "(say 'read-level taxonomic support'):\n" + "\n".join(found))

    def test_the_scanner_catches_what_it_is_for(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.md"
            bad.write_text("Bracken shows which species carries the island.\n\n"
                           "| `--pks_taxa` | the carrier species |\n\n"
                           "Kraken/Bracken are read-level support, not carrier attribution.\n\n"
                           "The genome bin carries clbA-clbS.\n")
            bad = bad.resolve()
            found = offenders([bad], READ_LEVEL, CARRIER, DISCLAIMED)
        self.assertEqual(len(found), 2, found)

    def test_the_docs_name_read_level_taxonomic_support(self):
        for path in ("README.md", "docs/running/taxonomy.md"):
            with self.subTest(path=path):
                self.assertIn("read-level taxonomic support", (ROOT / path).read_text())

    def test_taxonomy_output_columns_do_not_imply_a_carrier(self):
        spec = importlib.util.spec_from_file_location(
            "sdrt", ROOT / "scripts/summarize_diamond_rescue_taxonomy.py")
        sys.path.insert(0, str(ROOT / "scripts"))
        try:
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        finally:
            sys.path.pop(0)
        # `organism` read as "the organism the clb read came from"; the column is the
        # read's own krakenuniq taxon.
        self.assertIn("read_taxon", module.FIELDS)
        self.assertNotIn("organism", module.FIELDS)
        for column in module.FIELDS:
            self.assertNotRegex(column, r"carrier|host|producer|source")


    def test_the_report_still_reads_a_pre_rc1_organism_column(self):
        spec = importlib.util.spec_from_file_location("bcr_wording",
                                                      ROOT / "scripts/build_cohort_report.py")
        report = importlib.util.module_from_spec(spec)
        sys.path.insert(0, str(ROOT / "scripts"))
        try:
            spec.loader.exec_module(report)
        finally:
            sys.path.pop(0)
        with tempfile.TemporaryDirectory() as tmp:
            old = Path(tmp) / "old.tsv"
            old.write_text("sample\torganism\ttaxid\tclb_gene\tread_count\n"
                           "S1\tKlebsiella pneumoniae\t573\tclbB\t4\n")
            rows = report.diamond_rescue_reads(old)
        self.assertEqual(rows[0]["read_taxon"], "Klebsiella pneumoniae")


class DivergenceIsOnlyDivergence(unittest.TestCase):
    def test_no_specificity_or_uniqueness_claim_for_nm(self):
        paths = corpus() + [ROOT / "CHANGELOG.md"]
        found = offenders(paths, DIVERGENCE, OVERCLAIM, NEGATED)
        self.assertEqual(found, [], "\nNM/divergence text claiming more than divergence from "
                                    "IHE3034:\n" + "\n".join(found))


if __name__ == "__main__":
    unittest.main()
