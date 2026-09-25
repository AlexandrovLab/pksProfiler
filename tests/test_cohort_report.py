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


class TheReportNeverShowsAStaleSample(unittest.TestCase):
    """Ludmil, revised report finding 5.

    The report scans by_sample/ in the published tree by path, not by consuming the
    current run's channels, so a directory left over from an older run at the same
    --outdir was indistinguishable from a sample this run produced. --expected-samples
    is the filter: everything on disk still gets read, but only IDs the current run
    actually named are reported.

    His acceptance test: place a stale fake sample in outdir and delay one optional
    branch. The stale sample must never appear, and the report must not start until
    the delayed enabled branch completes. This class covers the first half; the
    Wiring class below covers the second.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = build_results(self.tmp.name, COHORT)
        # A leftover directory from an earlier run at this same --outdir, complete
        # with its own read_evidence.tsv -- exactly what a stale sample looks like.
        stale = self.results / "by_sample" / "S_stale_old_run"
        stale.mkdir(parents=True)
        (stale / "read_evidence.tsv").write_text(
            "sample\tread_evidence\tpks_reads\tclb_genes_detected\t"
            "island_breadth_1x\tisland_breadth_2x\tisland_breadth_3x\n"
            "S_stale_old_run\textensive_island\t9999\t19\t0.99\t0.9\t0.8\n")
        self.expected = Path(self.tmp.name) / "expected_samples.txt"
        self.expected.write_text("\n".join(name for name, *_ in COHORT) + "\n")

    def test_without_the_filter_the_stale_sample_leaks_through(self):
        names = {r["sample"] for r in report.collect(self.results)}
        self.assertIn("S_stale_old_run", names,
                       "fixture assumption: an unfiltered scan does pick it up")

    def test_with_the_filter_the_stale_sample_is_gone(self):
        expected = report.expected_sample_ids(self.expected)
        names = {r["sample"] for r in report.collect(self.results, expected)}
        self.assertNotIn("S_stale_old_run", names)
        self.assertEqual(names, {name for name, *_ in COHORT})

    def test_the_cli_flag_drops_it_from_the_written_table(self):
        html = Path(self.tmp.name) / "report.html"
        tsv = Path(self.tmp.name) / "report.tsv"
        run = subprocess.run(
            [sys.executable, str(SCRIPT), "--results", str(self.results),
             "--output", str(html), "--table", str(tsv),
             "--expected-samples", str(self.expected)],
            capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        with tsv.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({r["sample"] for r in rows}, {name for name, *_ in COHORT})
        self.assertNotIn("S_stale_old_run", html.read_text())

    def test_omitting_the_flag_still_works(self):
        # The flag is optional so existing manual invocations are not broken by it.
        html = Path(self.tmp.name) / "report.html"
        run = subprocess.run(
            [sys.executable, str(SCRIPT), "--results", str(self.results),
             "--output", str(html)],
            capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn("S_stale_old_run", html.read_text())


class TheCohortReportTaskWaitsForWhatItReads(unittest.TestCase):
    """The second half of finding 5's acceptance test.

    cohort_report_gate previously only tracked classifyPksReadEvidence -- the
    per-sample read-tier file -- so a run could reach cohortReport before
    pks.gene.counts.align.txt (masterTableAlign) or, worse, before targeted assembly
    published contigs/final_evidence/final_pks_evidence.tsv, since assembly is the
    slowest lane in the pipeline. Both are now mixed into the same gate.
    """

    def test_the_gene_matrix_task_feeds_the_gate(self):
        body = MAIN.split("workflow {", 1)[1]
        align_block = body.split("masterTableAlign(", 1)[1].split("\n    }", 1)[0]
        self.assertIn("cohort_report_gate = cohort_report_gate.mix(masterTableAlign.out)",
                      align_block)

    def test_targeted_assembly_feeds_the_gate_only_when_it_runs(self):
        body = MAIN.split("workflow {", 1)[1]
        assembly_block = body.split("if (tumor_targeted_assembly_b) {", 1)[1] \
                              .split("\n            }", 1)[0]
        self.assertIn("targetedPksAssembly(tumor_assembly_reads_ch, targeted_profiles_ch)",
                      assembly_block)
        self.assertIn("cohort_report_gate = cohort_report_gate.mix(", assembly_block)
        self.assertIn("targetedPksAssembly.out.evidence", assembly_block)

    def test_expected_sample_ids_is_the_same_file_masterQCSummary_already_uses(self):
        body = MAIN.split("workflow {", 1)[1]
        qc_pos = body.index("masterQCSummary(")
        self.assertIn("EXPECTED_SAMPLE_IDS", body[qc_pos:qc_pos + 300])
        report_pos = body.index("cohortReport(")
        self.assertIn("EXPECTED_SAMPLE_IDS", body[report_pos:report_pos + 150])

    def test_the_process_declares_and_passes_the_expected_samples_input(self):
        process = PLOTTING.split("process cohortReport {", 1)[1].split("\nprocess ", 1)[0]
        self.assertIn("path(expected_samples)", process)
        self.assertIn('--expected-samples "${expected_samples}"', process)


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


class RecruitmentStatsSurfaceInTheContigHoverText(unittest.TestCase):
    """recruited_fragment_ids/paired_fragments explain *why* targeted assembly
    produced nothing for a sample (no reads recruited vs. reads recruited but
    not paired), so they ride along in the same hover tooltip as the existing
    per-assembler contig counts rather than a new column of their own.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = build_results(self.tmp.name, COHORT)
        evidence = (self.results / "by_sample/S_ext/contigs/final_evidence/"
                    "final_pks_evidence.tsv")
        evidence.write_text(
            "sample\tfinal_structural_evidence\tassembler_agreement\t"
            "megahit_reference_coverage\tmetaspades_reference_coverage\t"
            "recruited_fragment_ids\tpaired_fragments\n"
            "S_ext\tcomplete_island\tconcordant\t0.93\t0.91\t9\t3\n")

    def test_collect_carries_the_recruitment_fields(self):
        record = next(r for r in report.collect(self.results) if r["sample"] == "S_ext")
        self.assertEqual(record["recruited_fragment_ids"], "9")
        self.assertEqual(record["paired_fragments"], "3")

    def test_the_hover_text_names_them(self):
        html_path = Path(self.tmp.name) / "report.html"
        subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.results),
                        "--output", str(html_path)], capture_output=True, text=True)
        self.assertIn("9 fragments recruited, 3 paired", html_path.read_text())


MAG_SUMMARY_HEADER = ("sample\tbin_id\ttaxonomy\tcompleteness\tcontamination\tgenome_size\t"
                     "contig_n50\tclb_genes_detected\tclb_genes\tbest_evalue\thas_integrase\t"
                     "has_transposase\tflanking_genes\tunexpected_taxon_flag\tlocus_tier\t"
                     "locus_genes_detected\tlocus_breadth\n")
PROPHAGE_HEADER = ("sample\tbin_id\thost_taxonomy\tprophage_id\thost_contig\tstart\tend\t"
                   "length\tvirus_score\tviral_taxonomy\thas_recA\thas_lexA\tclbS_like\t"
                   "evidence_level\n")
RESCUE_HEADER = "sample\torganism\ttaxid\tclb_gene\tread_count\n"


class GenomeBinsPanel(unittest.TestCase):
    """Item 1: magBinLocusEvidence's per-bin call, surfaced per sample.

    pks_mag_summary.tsv already carries the alignment-confirmed locus_tier -- same
    vocabulary as the read-level tier -- plus taxonomy, completeness and breadth. The
    report must show it, not recompute it.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = build_results(self.tmp.name, COHORT)
        genomes = self.results / "by_sample/S_ext/genomes"
        genomes.mkdir(parents=True, exist_ok=True)
        (genomes / "pks_mag_summary.tsv").write_text(
            MAG_SUMMARY_HEADER +
            "S_ext\tbin_001\td__Bacteria;o__Enterobacterales;s__Escherichia coli\t"
            "92.5\t1.2\t4500000\t125000\t10\tclbA,clbB\t1e-20\tTrue\tFalse\tgyrB\tFalse\t"
            "broad_island\t9\t0.09\n")
        community = self.results / "by_sample/S_ext/community"
        community.mkdir(parents=True, exist_ok=True)
        (community / "community_prophage_inventory.tsv").write_text(
            PROPHAGE_HEADER +
            "S_ext\tbin_001\tEscherichia coli\tc1|provirus_1\tc1\t100\t5000\t4900\t0.9\t"
            "Caudovirales\tTrue\tTrue\tFalse\tpredicted_provirus\n")
        self.html = Path(self.tmp.name) / "report.html"
        self.run = subprocess.run(
            [sys.executable, str(SCRIPT), "--results", str(self.results),
             "--output", str(self.html)], capture_output=True, text=True)
        self.text = self.html.read_text()

    def test_it_still_runs(self):
        self.assertEqual(self.run.returncode, 0, self.run.stderr)

    def test_collect_reads_the_bin_row_without_recomputing_it(self):
        by_sample = {r["sample"]: r for r in report.collect(self.results)}
        bins = by_sample["S_ext"]["bins"]
        self.assertEqual(len(bins), 1)
        self.assertEqual(bins[0]["locus_tier"], "broad_island")
        self.assertEqual(bins[0]["taxonomy"], "d__Bacteria;o__Enterobacterales;s__Escherichia coli")
        self.assertAlmostEqual(bins[0]["completeness"], 92.5)

    def test_prophage_regions_are_grouped_by_bin_not_reparsed_from_genomad(self):
        by_sample = {r["sample"]: r for r in report.collect(self.results)}
        self.assertEqual(by_sample["S_ext"]["bins"][0]["prophage_regions"], 1)

    def test_a_sample_with_no_mag_summary_gets_no_bins_invented(self):
        by_sample = {r["sample"]: r for r in report.collect(self.results)}
        self.assertEqual(by_sample["S_neg"]["bins"], [])

    def test_the_row_link_and_the_panel_share_an_anchor(self):
        import re
        links = set(re.findall(r'<a href="#(bins-[^"]+)"', self.text))
        panels = set(re.findall(r'class="figpanel" id="(bins-[^"]+)"', self.text))
        self.assertTrue(links)
        self.assertEqual(links, panels)

    def test_the_locus_tier_badge_reuses_the_read_tier_colour(self):
        self.assertIn(f'background:{report.TIER_FILL["broad_island"]}', self.text)
        self.assertIn(">Broad<", self.text)

    def test_prophage_context_is_a_one_line_summary_per_bin(self):
        self.assertIn("yes, 1 region", self.text)

    def test_a_bin_never_aligned_reads_not_aligned_not_zero_percent(self):
        # magBinLocusEvidence not yet having run for a bin is "NA" in the pipeline's own
        # output (mag_utils.read_locus_evidence); the report must not turn that into 0%.
        genomes = self.results / "by_sample/S_broad/genomes"
        genomes.mkdir(parents=True, exist_ok=True)
        (genomes / "pks_mag_summary.tsv").write_text(
            MAG_SUMMARY_HEADER +
            "S_broad\tbin_002\tunclassified\t50.0\t2.0\t3000000\t50000\t2\tclbA\t1e-8\t"
            "False\tFalse\t\tFalse\tNA\tNA\tNA\n")
        by_sample = {r["sample"]: r for r in report.collect(self.results)}
        self.assertEqual(by_sample["S_broad"]["bins"][0]["locus_tier"], "NA")
        html = Path(self.tmp.name) / "report2.html"
        subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.results),
                        "--output", str(html)], capture_output=True, text=True)
        self.assertIn("Not aligned", html.read_text())


class DiamondRescueTaxonomyPanel(unittest.TestCase):
    """Item 2: DIAMOND-rescued reads, joined back to their krakenPrefilter organism.

    Every row here is, by construction, a read the fast classifier's target-taxon
    routing did not keep -- see summarize_diamond_rescue_taxonomy.py's docstring.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = build_results(self.tmp.name, COHORT)
        prefilter = self.results / "by_sample/S_ext/prefilter"
        prefilter.mkdir(parents=True, exist_ok=True)
        (prefilter / "diamond_rescue_taxonomy.tsv").write_text(
            RESCUE_HEADER +
            "S_ext\tKlebsiella pneumoniae\t573\tclbB\t4\n"
            "S_ext\tKlebsiella pneumoniae\t573\tclbA\t1\n")
        self.html = Path(self.tmp.name) / "report.html"
        subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.results),
                        "--output", str(self.html)], capture_output=True, text=True)
        self.text = self.html.read_text()

    def test_collect_reads_both_rows(self):
        by_sample = {r["sample"]: r for r in report.collect(self.results)}
        rescued = by_sample["S_ext"]["diamond_rescue"]
        self.assertEqual(len(rescued), 2)
        self.assertEqual({r["clb_gene"] for r in rescued}, {"clbA", "clbB"})

    def test_a_sample_with_no_rescue_file_gets_none_invented(self):
        by_sample = {r["sample"]: r for r in report.collect(self.results)}
        self.assertEqual(by_sample["S_neg"]["diamond_rescue"], [])

    def test_the_row_link_and_the_panel_share_an_anchor(self):
        import re
        links = set(re.findall(r'<a href="#(rescue-[^"]+)"', self.text))
        panels = set(re.findall(r'class="figpanel" id="(rescue-[^"]+)"', self.text))
        self.assertTrue(links)
        self.assertEqual(links, panels)

    def test_the_panel_names_organism_gene_and_count(self):
        anchor = report.detail_anchor("rescue", "S_ext")
        panel = self.text[self.text.index(f'id="{anchor}"'):]
        self.assertIn("Klebsiella pneumoniae", panel)
        self.assertIn("clbB", panel)
        self.assertIn(">4<", panel)


class CommunitySpeciesSupportSection(unittest.TestCase):
    """Item 3: pks.clb_species_support.tsv, computed only with --pks_taxa and, until
    now, read by nothing."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = build_results(self.tmp.name, COHORT)
        self.html = Path(self.tmp.name) / "report.html"

    def _run(self):
        subprocess.run([sys.executable, str(SCRIPT), "--results", str(self.results),
                        "--output", str(self.html)], capture_output=True, text=True)
        return self.html.read_text()

    def test_absent_when_pks_taxa_never_ran(self):
        # The bullet in "Where these numbers come from" documents the path regardless;
        # the section itself -- the table and its heading -- must not appear.
        text = self._run()
        self.assertNotIn("Direct krakenuniq support for reads aligned", text)
        self.assertNotIn("Community species", text)

    def test_present_when_the_file_exists(self):
        taxonomy = self.results / "cohort/taxonomy"
        taxonomy.mkdir(parents=True, exist_ok=True)
        header = "Sample\tSpecies\tTaxID\t" + "\t".join(GENES) + "\tTotal\n"
        row = ("S_ext\tEscherichia coli\t562\t" +
               "\t".join(["3"] + ["0"] * (len(GENES) - 1)) + "\t3\n")
        (taxonomy / "pks.clb_species_support.tsv").write_text(header + row)
        text = self._run()
        self.assertIn("clb-gene support", text)
        self.assertIn("Escherichia coli", text)
        self.assertIn("<th>clbA</th>", text)

    def test_species_support_reads_the_file_verbatim(self):
        rows = report.species_support(self.results)
        self.assertEqual(rows, [])
        taxonomy = self.results / "cohort/taxonomy"
        taxonomy.mkdir(parents=True, exist_ok=True)
        header = "Sample\tSpecies\tTaxID\t" + "\t".join(GENES) + "\tTotal\n"
        row = ("S_ext\tEscherichia coli\t562\t" +
               "\t".join(["3"] + ["0"] * (len(GENES) - 1)) + "\t3\n")
        (taxonomy / "pks.clb_species_support.tsv").write_text(header + row)
        rows = report.species_support(self.results)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Species"], "Escherichia coli")
