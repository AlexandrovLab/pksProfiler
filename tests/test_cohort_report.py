"""One page across the cohort, derived from the pipeline's own outputs.

The pipeline writes a directory per sample; at two thousand samples that is two
thousand directories and no view of the cohort. An earlier version of this view was
built by hand for the 18 September update
(pksProfiler_analysis/mutographs_evidence_snapshots/build_mutographs_evidence.py) and
recomputed breadth with `samtools depth` per sample and re-derived the evidence class
with the v0.0.1 thresholds -- so it could disagree with the run it was describing.

This report reads and arranges; it does not derive. These tests are mostly about that
property.
"""
import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/build_cohort_report.py"
MAIN = (ROOT / "main.nf").read_text()
PLOTTING = (ROOT / "Modules/plotting.nf").read_text()
GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]

spec = importlib.util.spec_from_file_location("build_cohort_report", SCRIPT)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def build_results(directory, samples):
    """samples: (name, tier, reads, genes_detected, breadth, assembled)"""
    root = Path(directory)
    counts = {}
    for name, tier, reads, ngenes, breadth, assembled in samples:
        per = {g: 0 for g in GENES}
        for gene in GENES[:ngenes]:
            per[gene] = max(1, reads // max(ngenes, 1))
        counts[name] = per
        sample_dir = root / "by_sample" / name
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / "read_evidence.tsv").write_text(
            "sample\tread_evidence\tpks_reads\tclb_genes_detected\t"
            "island_breadth_1x\tisland_breadth_2x\tisland_breadth_3x\n"
            f"{name}\t{tier}\t{reads}\t{ngenes}\t{breadth}\t{breadth/2}\t{breadth/3}\n")
        if assembled:
            # The path pks_targeted.nf actually publishes to. The earlier fixture
            # wrote contigs/final_pks_evidence.tsv, matching the reader rather than the
            # pipeline, so the suite passed while the real report showed no contigs.
            (sample_dir / "contigs/final_evidence").mkdir(parents=True, exist_ok=True)
            (sample_dir / "contigs/final_evidence/final_pks_evidence.tsv").write_text(
                "sample\tfinal_structural_evidence\tassembler_agreement\t"
                "megahit_reference_coverage\tmetaspades_reference_coverage\n"
                f"{name}\tcomplete_island\tconcordant\t0.93\t0.91\n")
    matrix = root / "cohort/gene_counts/pks.gene.counts.align.txt"
    matrix.parent.mkdir(parents=True, exist_ok=True)
    names = list(counts)
    matrix.write_text("Gene\t" + "\t".join(names) + "\n" + "\n".join(
        gene + "\t" + "\t".join(str(counts[n][gene]) for n in names) for gene in GENES) + "\n")
    return root


COHORT = [
    ("S_ext", "extensive_island", 1400, 19, 0.61, True),
    ("S_broad", "broad_island", 180, 11, 0.12, False),
    ("S_multi", "multi_gene", 26, 5, 0.03, False),
    ("S_neg", "negative", 0, 0, 0.0, False),
]


class ItReadsRatherThanDerives(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = build_results(self.tmp.name, COHORT)

    def records(self):
        return report.collect(self.results)

    def test_the_tier_is_the_pipelines_not_recomputed(self):
        # Contradictory input on purpose: a tier that the numbers would not support.
        path = self.results / "by_sample/S_neg/read_evidence.tsv"
        path.write_text("sample\tread_evidence\tpks_reads\tclb_genes_detected\t"
                        "island_breadth_1x\tisland_breadth_2x\tisland_breadth_3x\n"
                        "S_neg\textensive_island\t0\t0\t0.0\t0.0\t0.0\n")
        tiers = {r["sample"]: r["tier"] for r in self.records()}
        self.assertEqual(tiers["S_neg"], "extensive_island",
                         "the report must report the pipeline's tier, not its own opinion")

    def test_breadth_is_taken_not_measured(self):
        # No BAM exists anywhere in this tree; the report must still have breadth.
        self.assertEqual(list(self.results.rglob("*.bam")), [])
        breadths = {r["sample"]: r["breadth_1x"] for r in self.records()}
        self.assertAlmostEqual(breadths["S_ext"], 0.61)

    def test_samples_are_ordered_by_evidence(self):
        self.assertEqual([r["sample"] for r in self.records()],
                         ["S_ext", "S_broad", "S_multi", "S_neg"])

    def test_a_sample_without_contig_analysis_is_not_invented(self):
        by_sample = {r["sample"]: r for r in self.records()}
        self.assertEqual(by_sample["S_ext"]["structural"], "complete_island")
        self.assertEqual(by_sample["S_broad"]["structural"], "")


class TheOutputs(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = build_results(self.tmp.name, COHORT)
        self.html = Path(self.tmp.name) / "report.html"
        self.tsv = Path(self.tmp.name) / "report.tsv"
        self.run = subprocess.run(
            [sys.executable, str(SCRIPT), "--results", str(self.results),
             "--output", str(self.html), "--table", str(self.tsv)],
            capture_output=True, text=True)

    def test_it_runs_and_writes_both(self):
        self.assertEqual(self.run.returncode, 0, self.run.stderr)
        self.assertTrue(self.html.is_file() and self.tsv.is_file())

    def test_the_svg_is_well_formed_and_fits_its_viewbox(self):
        text = self.html.read_text()
        svg = text[text.index("<svg"):text.index("</svg>") + 6]
        root = ET.fromstring(svg)          # raises if malformed
        width = float(root.get("width"))
        ns = "{http://www.w3.org/2000/svg}"
        rights = [float(r.get("x")) + float(r.get("width"))
                  for r in root.iter(ns + "rect") if r.get("x")]
        self.assertLessEqual(max(rights), width)

    def test_the_table_carries_every_gene_and_every_sample(self):
        with self.tsv.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), len(COHORT))
        for gene in GENES:
            self.assertIn(gene, rows[0])

    def test_it_says_when_no_sample_was_assembled(self):
        for sample, *_ in COHORT:
            evidence = (self.results / "by_sample" / sample /
                        "contigs/final_evidence/final_pks_evidence.tsv")
            if evidence.is_file():
                evidence.unlink()
        subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.results),
                        "--output", str(self.html)], capture_output=True, text=True)
        self.assertIn("did not run for any", self.html.read_text())

    def test_an_empty_results_tree_is_an_error_not_an_empty_page(self):
        empty = Path(self.tmp.name) / "empty"
        empty.mkdir()
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--results", str(empty),
             "--output", str(Path(self.tmp.name) / "x.html")],
            capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)


class Wiring(unittest.TestCase):
    def test_the_process_exists_and_is_invoked(self):
        self.assertIn("process cohortReport", PLOTTING)
        self.assertIn("cohortReport(", MAIN)

    def test_it_runs_after_the_cohort_tables(self):
        body = MAIN.split("workflow {", 1)[1]
        self.assertLess(body.index("masterQCSummary("), body.index("cohortReport("))


if __name__ == "__main__":
    unittest.main()


class HoverSurvivesASandboxedViewer(unittest.TestCase):
    """JupyterLab renders HTML files with JavaScript disabled.

    The first version put hover text only in data-tip attributes read by a script, so
    in the viewer it was opened in, nothing happened. Every hoverable thing carries a
    native <title> as well.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = build_results(self.tmp.name, COHORT)
        self.html = Path(self.tmp.name) / "report.html"
        subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.results),
                        "--output", str(self.html)], capture_output=True, text=True)
        self.text = self.html.read_text()

    def test_every_data_tip_has_a_title_beside_it(self):
        # one <title> per hoverable group or cell, so none depends on the script
        self.assertGreaterEqual(self.text.count("<title>"), self.text.count('data-tip="'))

    def test_the_row_hover_carries_the_depths_that_lost_their_bars(self):
        self.assertRegex(self.text, r"breadth [\d.]+% at ≥1×, [\d.]+% at ≥2×")

    def test_only_one_breadth_bar_is_drawn(self):
        # >=2x and >=3x are nested inside >=1x and decide nothing; three bars per row
        # is three times the ink for one number.
        self.assertEqual(self.text.count('fill="#636363"'), len(COHORT))
        self.assertNotIn('fill="#969696"', self.text)


class FiguresOpenInlineFromTheRow(unittest.TestCase):
    """A figure link sits in the row and opens the figure in place.

    Not in a new tab, and not in a list below the table: the link points at an anchor
    on this page and CSS reveals that panel. :target rather than a click handler,
    because the sandboxed viewers that strip JavaScript still follow fragments.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = build_results(self.tmp.name, COHORT)
        for sample, *_rest in COHORT[:2]:          # the two with reads
            figures = self.results / "by_sample" / sample / "figures"
            figures.mkdir(parents=True, exist_ok=True)
            (figures / "contig_validation.svg").write_text(
                '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"/>')
            (figures / "circos.pdf").write_text("pdf")
        self.html = Path(self.tmp.name) / "report.html"
        subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.results),
                        "--output", str(self.html)], capture_output=True, text=True)
        self.text = self.html.read_text()

    def test_the_link_is_in_the_row_not_a_section_below(self):
        import xml.etree.ElementTree as ET
        svg = self.text[self.text.index("<svg"):self.text.index("</svg>") + 6]
        root = ET.fromstring(svg)
        links = root.findall(".//{http://www.w3.org/2000/svg}a")
        self.assertEqual(len(links), 4)            # two samples, two figures each
        self.assertNotIn("<details>", self.text)

    def test_every_link_resolves_to_a_panel_on_this_page(self):
        import re
        links = set(re.findall(r'<a href="#(fig-[^"]+)"', self.text))
        panels = set(re.findall(r'class="figpanel" id="([^"]+)"', self.text))
        self.assertTrue(links)
        self.assertEqual(links, panels)

    def test_panels_are_hidden_until_targeted(self):
        self.assertIn(".figpanel:target", self.text)
        self.assertRegex(self.text, r"\.figpanel \{ display:none")

    def test_nothing_navigates_away(self):
        self.assertNotIn('target="_blank"', self.text)

    def test_svg_embeds_and_pdf_is_linked(self):
        self.assertIn('type="image/svg+xml"', self.text)
        self.assertIn("cannot be shown inline", self.text)

    def test_paths_are_relative_to_the_cohort_directory(self):
        self.assertIn('data="../by_sample/S_ext/figures/contig_validation.svg"', self.text)

    def test_a_sample_without_figures_says_none(self):
        self.assertIn(">none<", self.text)

    def test_opening_a_figure_requires_no_javascript(self):
        panel = self.text.split('class="figpanel"', 1)[1]
        self.assertNotIn("onclick", panel)
