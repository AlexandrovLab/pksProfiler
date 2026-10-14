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


# Ludmil, 2026-09-25 report P0-5b: the script reads a manifest of staged inputs, not a
# results tree. Most tests here still build a published-layout fixture -- the layout is
# the contract with the publishDir directives -- and hand the script a manifest naming
# its files, which is what main.nf's gates stage into the task. The P0-5b classes at the
# end of this file are where staged and published copies are made to disagree.
PUBLISHED_FIGURES = {name: kind for kind, name in report.FIGURE_INPUTS.items()}


def manifest_from_tree(results, manifest=None):
    """Write a manifest naming every report input present in a published-layout tree."""
    results = Path(results)
    rows = []
    for kind, rel in report.COHORT_INPUTS.items():
        if (results / rel).is_file():
            rows.append(("", kind, results / rel))
    by_sample = results / "by_sample"
    for sample_dir in sorted(by_sample.iterdir()) if by_sample.is_dir() else ():
        for kind, rel in report.SAMPLE_INPUTS.items():
            if (sample_dir / rel).is_file():
                rows.append((sample_dir.name, kind, sample_dir / rel))
        for figure in sorted((sample_dir / "figures").glob("*")):
            if figure.name in PUBLISHED_FIGURES:
                rows.append((sample_dir.name, PUBLISHED_FIGURES[figure.name], figure))
    manifest = Path(manifest) if manifest else results.parent / f"{results.name}.inputs.tsv"
    manifest.write_text("sample\tkind\tpath\n" + "".join(
        f"{sample}\t{kind}\t{path}\n" for sample, kind, path in rows))
    return manifest


def collect_tree(results, expected=None):
    return report.collect(report.read_manifest(manifest_from_tree(results)), expected)


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
        return collect_tree(self.results)

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
            [sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
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
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
                        "--output", str(self.html)], capture_output=True, text=True)
        self.assertIn("did not run for any", self.html.read_text())

    def test_an_empty_results_tree_is_an_error_not_an_empty_page(self):
        empty = Path(self.tmp.name) / "empty"
        empty.mkdir()
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(empty)),
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

    The report scanned by_sample/ in the published tree by path, not by consuming the
    current run's channels, so a directory left over from an older run at the same
    --outdir was indistinguishable from a sample this run produced. --expected-samples
    is the filter: only IDs the current run actually named are reported. Since P0-5b
    the report reads a manifest of staged inputs instead (see the P0-5b classes below);
    the filter still holds for any manifest row naming a sample outside the run.

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
        names = {r["sample"] for r in collect_tree(self.results)}
        self.assertIn("S_stale_old_run", names,
                       "fixture assumption: an unfiltered scan does pick it up")

    def test_with_the_filter_the_stale_sample_is_gone(self):
        expected = report.expected_sample_ids(self.expected)
        names = {r["sample"] for r in collect_tree(self.results, expected)}
        self.assertNotIn("S_stale_old_run", names)
        self.assertEqual(names, {name for name, *_ in COHORT})

    def test_the_cli_flag_drops_it_from_the_written_table(self):
        html = Path(self.tmp.name) / "report.html"
        tsv = Path(self.tmp.name) / "report.tsv"
        run = subprocess.run(
            [sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
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
            [sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
             "--output", str(html)],
            capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn("S_stale_old_run", html.read_text())


def mix_calls(body, gate):
    """The full text of every `<gate>.mix(...)` call, however many lines it spans."""
    calls, start = [], 0
    needle = f"{gate}.mix("
    while (start := body.find(needle, start)) >= 0:
        depth, end = 0, start + len(needle) - 1
        for end in range(start + len(needle) - 1, len(body)):
            depth += {"(": 1, ")": -1}.get(body[end], 0)
            if depth == 0:
                break
        calls.append(body[start:end + 1])
        start = end
    return calls


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
        self.assertIn("cohort_report_gate = cohort_report_gate.mix(", align_block)
        self.assertIn("masterTableAlign.out.map { f -> tuple('', 'gene_counts', f) }",
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
        call = body[report_pos:body.index("\n\n", report_pos)]
        self.assertIn("EXPECTED_SAMPLE_IDS", call)

    # Ludmil, 2026-09-25 report P0-5: every file build_cohort_report.py reads, and the
    # process whose completion publishes it. The report must wait for each of them.
    DISPLAYED_PRODUCERS = {
        "genomes/pks_mag_summary.tsv": "optional_lane_gate",
        "community/community_prophage_inventory.tsv": "optional_lane_gate",
        "prefilter/diamond_rescue_taxonomy.tsv": "diamondRescueTaxonomy",
        "cohort/taxonomy/pks.clb_species_support.tsv": "combineClbTaxonomySupport",
        "figures/circos.pdf": "plotPKS",
        "figures (island taxa)": "plotPKSTaxa",
        "figures (community taxa)": "plotCommunityTaxa",
    }

    def test_the_report_waits_on_all_three_gates(self):
        # P0-5b: the gates are now the report's staged inputs, not an ordering marker.
        body = MAIN.split("workflow {", 1)[1]
        inputs = body[body.index("COHORT_REPORT_INPUTS = sortReportInputs("):]
        inputs = inputs[:inputs.index("cohortReport(")]
        for gate in ("cohort_report_gate", "optional_lane_gate", "report_display_gate"):
            with self.subTest(gate=gate):
                self.assertIn(gate, inputs)
        call = body[body.index("cohortReport("):]
        call = call[:call.index("\n\n")]
        self.assertIn("COHORT_REPORT_INPUTS.map { rows -> rows.collect { row -> row[2] } }", call)
        self.assertIn("COHORT_REPORT_INPUTS.map { rows -> reportManifestText(rows) }", call)

    def test_every_displayed_producer_feeds_a_gate(self):
        body = MAIN.split("workflow {", 1)[1]
        mixes = mix_calls(body, "report_display_gate")
        for path, producer in self.DISPLAYED_PRODUCERS.items():
            if producer == "optional_lane_gate":
                continue
            with self.subTest(path=path):
                self.assertTrue(any(f"{producer}.out" in call for call in mixes),
                                f"{producer} publishes {path} but no gate waits for it")

    def test_the_display_gate_is_never_reassigned(self):
        body = MAIN.split("workflow {", 1)[1]
        import re
        assignments = re.findall(r"report_display_gate\s*=\s*(\S+)", body)
        self.assertTrue(assignments)
        for rhs in assignments:
            with self.subTest(rhs=rhs):
                self.assertTrue(rhs.startswith("channel.empty()") or
                                rhs.startswith("report_display_gate.mix("), rhs)

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
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
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
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
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
        record = next(r for r in collect_tree(self.results) if r["sample"] == "S_ext")
        self.assertEqual(record["recruited_fragment_ids"], "9")
        self.assertEqual(record["paired_fragments"], "3")

    def test_the_hover_text_names_them(self):
        html_path = Path(self.tmp.name) / "report.html"
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
                        "--output", str(html_path)], capture_output=True, text=True)
        self.assertIn("9 fragments recruited, 3 paired", html_path.read_text())


MAG_SUMMARY_HEADER = ("sample\tbin_id\tunit_type\ttaxonomy\tcompleteness\tcontamination\tgenome_size\t"
                     "contig_n50\tclb_genes_detected\tclb_genes\tbest_evalue\thas_integrase\t"
                     "has_transposase\tflanking_genes\tunexpected_taxon_flag\tlocus_tier\t"
                     "locus_genes_detected\tlocus_breadth\tmag_qc_status\n")
PROPHAGE_HEADER = ("sample\tbin_id\thost_taxonomy\tprophage_id\thost_contig\tstart\tend\t"
                   "length\tvirus_score\tviral_taxonomy\thas_recA\thas_lexA\tclbS_like\t"
                   "evidence_level\n")
RESCUE_HEADER = "sample\tread_taxon\ttaxid\tclb_gene\tread_count\n"


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
            "S_ext\tbin_001\tbin\td__Bacteria;o__Enterobacterales;s__Escherichia coli\t"
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
            [sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
             "--output", str(self.html)], capture_output=True, text=True)
        self.text = self.html.read_text()

    def test_it_still_runs(self):
        self.assertEqual(self.run.returncode, 0, self.run.stderr)

    def test_collect_reads_the_bin_row_without_recomputing_it(self):
        by_sample = {r["sample"]: r for r in collect_tree(self.results)}
        bins = by_sample["S_ext"]["bins"]
        self.assertEqual(len(bins), 1)
        self.assertEqual(bins[0]["locus_tier"], "broad_island")
        self.assertEqual(bins[0]["taxonomy"], "d__Bacteria;o__Enterobacterales;s__Escherichia coli")
        self.assertAlmostEqual(bins[0]["completeness"], 92.5)

    def test_prophage_regions_are_grouped_by_bin_not_reparsed_from_genomad(self):
        by_sample = {r["sample"]: r for r in collect_tree(self.results)}
        self.assertEqual(by_sample["S_ext"]["bins"][0]["prophage_regions"], 1)

    def test_a_sample_with_no_mag_summary_gets_no_bins_invented(self):
        by_sample = {r["sample"]: r for r in collect_tree(self.results)}
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
            "S_broad\tbin_002\tbin\tunclassified\t50.0\t2.0\t3000000\t50000\t2\tclbA\t1e-8\t"
            "False\tFalse\t\tFalse\tNA\tNA\tNA\n")
        by_sample = {r["sample"]: r for r in collect_tree(self.results)}
        self.assertEqual(by_sample["S_broad"]["bins"][0]["locus_tier"], "NA")
        html = Path(self.tmp.name) / "report2.html"
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
                        "--output", str(html)], capture_output=True, text=True)
        self.assertIn("Not aligned", html.read_text())


    def test_the_fixture_header_is_the_producers_header(self):
        # A hand-written schema that drifted from build_mag_summary.py is how the master
        # summary's P0-2 bug stayed green; keep this fixture tied to the producer.
        sys.path.insert(0, str(ROOT / "scripts"))
        import build_mag_summary
        self.assertEqual(MAG_SUMMARY_HEADER.rstrip("\n").split("\t"),
                         build_mag_summary.FIELDNAMES)

    def write_s_ext(self, *rows):
        (self.results / "by_sample/S_ext/genomes/pks_mag_summary.tsv").write_text(
            MAG_SUMMARY_HEADER + "".join(rows))

    UNBINNED_ROW = ("S_ext\tunbinned\tunbinned\tNA\tNA\tNA\tNA\tNA\tNA\t\tNA\tNA\tNA\t\tNA\t"
                    "broad_island\t9\t0.18\n")

    def test_the_unbinned_pool_is_not_a_genome_bin(self):
        """Ludmil, 2026-09-25 report P0-4."""
        self.write_s_ext(
            "S_ext\tbin_001\tbin\td__Bacteria;o__Enterobacterales;s__Escherichia coli\t"
            "92.5\t1.2\t4500000\t125000\t10\tclbA,clbB\t1e-20\tTrue\tFalse\tgyrB\tFalse\t"
            "negative\t0\t0.0\n", self.UNBINNED_ROW)
        record = {r["sample"]: r for r in collect_tree(self.results)}["S_ext"]
        self.assertEqual([b["bin_id"] for b in record["bins"]], ["bin_001"])
        self.assertEqual(record["unbinned"]["locus_tier"], "broad_island")
        html = Path(self.tmp.name) / "report_unbinned.html"
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
                        "--output", str(html)], capture_output=True, text=True, check=True)
        text = html.read_text()
        self.assertIn("1 bin, 0 pks+ (+ unbinned pks+)", text)
        self.assertIn("not a genome bin and not a carrier attribution", text)
        self.assertNotIn("<td>unbinned</td>", text)

    def test_an_unbinned_pool_alone_still_gets_a_panel_but_zero_bins(self):
        self.write_s_ext(self.UNBINNED_ROW)
        record = {r["sample"]: r for r in collect_tree(self.results)}["S_ext"]
        self.assertEqual(record["bins"], [])
        html = Path(self.tmp.name) / "report_unbinned_only.html"
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
                        "--output", str(html)], capture_output=True, text=True, check=True)
        self.assertIn("0 bins, 0 pks+ (+ unbinned pks+)", html.read_text())

    def test_unmeasured_completeness_reads_na_not_zero_percent(self):
        self.write_s_ext(
            "S_ext\tbin_001\tbin\tunclassified\tNA\tNA\tNA\tNA\t0\t\tNA\tFalse\tFalse\t\tFalse\t"
            "broad_island\t9\t0.09\n")
        record = {r["sample"]: r for r in collect_tree(self.results)}["S_ext"]
        self.assertIsNone(record["bins"][0]["completeness"])
        html = Path(self.tmp.name) / "report_na.html"
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
                        "--output", str(html)], capture_output=True, text=True, check=True)
        self.assertNotIn("<td>0.0%</td>", html.read_text())
        self.assertIn("<td>NA</td>", html.read_text())

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
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
                        "--output", str(self.html)], capture_output=True, text=True)
        self.text = self.html.read_text()

    def test_collect_reads_both_rows(self):
        by_sample = {r["sample"]: r for r in collect_tree(self.results)}
        rescued = by_sample["S_ext"]["diamond_rescue"]
        self.assertEqual(len(rescued), 2)
        self.assertEqual({r["clb_gene"] for r in rescued}, {"clbA", "clbB"})

    def test_a_sample_with_no_rescue_file_gets_none_invented(self):
        by_sample = {r["sample"]: r for r in collect_tree(self.results)}
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
        # 26 Sep mode memo section 5: the taxon is the read's own, not a carrier.
        self.assertIn("<th>Read-level taxon</th>", panel)
        self.assertIn("not carrier attribution", panel)
        self.assertNotIn("<th>Organism</th>", panel)
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
        subprocess.run([sys.executable, str(SCRIPT), "--inputs-manifest", str(manifest_from_tree(self.results)),
                        "--output", str(self.html)], capture_output=True, text=True)
        return self.html.read_text()

    def test_absent_when_pks_taxa_never_ran(self):
        # The bullet in "Where these numbers come from" documents the path regardless;
        # the section itself -- the table and its heading -- must not appear.
        text = self._run()
        self.assertNotIn("per-read species assignment of the reads aligned", text)
        self.assertNotIn("Read-level taxonomic support: species", text)

    def test_present_when_the_file_exists(self):
        taxonomy = self.results / "cohort/taxonomy"
        taxonomy.mkdir(parents=True, exist_ok=True)
        header = "Sample\tSpecies\tTaxID\t" + "\t".join(GENES) + "\tTotal\n"
        row = ("S_ext\tEscherichia coli\t562\t" +
               "\t".join(["3"] + ["0"] * (len(GENES) - 1)) + "\t3\n")
        (taxonomy / "pks.clb_species_support.tsv").write_text(header + row)
        text = self._run()
        self.assertIn("Read-level taxonomic support: species", text)
        self.assertIn("not carrier attribution", text)
        self.assertIn("Escherichia coli", text)
        self.assertIn("<th>clbA</th>", text)

    def test_species_support_reads_the_file_verbatim(self):
        rows = report.species_support(report.read_manifest(manifest_from_tree(self.results)))
        self.assertEqual(rows, [])
        taxonomy = self.results / "cohort/taxonomy"
        taxonomy.mkdir(parents=True, exist_ok=True)
        header = "Sample\tSpecies\tTaxID\t" + "\t".join(GENES) + "\tTotal\n"
        row = ("S_ext\tEscherichia coli\t562\t" +
               "\t".join(["3"] + ["0"] * (len(GENES) - 1)) + "\t3\n")
        (taxonomy / "pks.clb_species_support.tsv").write_text(header + row)
        rows = report.species_support(report.read_manifest(manifest_from_tree(self.results)))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Species"], "Escherichia coli")


def stage(tree, work, manifest_name="pks.cohort_report.inputs.tsv"):
    """Lay a fixture tree's report inputs out the way cohortReport's task sees them.

    One inputs/<n>/ directory per file with the file's own basename, exactly what
    `stageAs: 'inputs/*/*'` does, and a manifest whose paths are relative to it --
    reportManifestText in main.nf. Per-sample basenames therefore repeat
    (inputs/3/read_evidence.tsv, inputs/4/read_evidence.tsv, ...).
    """
    import shutil
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    scratch = manifest_from_tree(tree)
    lines = scratch.read_text().splitlines()
    out = [lines[0]]
    for index, line in enumerate(lines[1:], start=1):
        sample, kind, path = line.split("\t")
        target = work / "inputs" / str(index) / Path(path).name
        target.parent.mkdir(parents=True)
        shutil.copyfile(path, target)
        out.append(f"{sample}\t{kind}\tinputs/{index}/{Path(path).name}")
    (work / manifest_name).write_text("\n".join(out) + "\n")
    return manifest_name


class P0_5bTheReportReadsItsStagedInputsNotThePublishedTree(unittest.TestCase):
    """Ludmil, 2026-09-25 report P0-5b.

    cohortReport ran `build_cohort_report.py --results ${params.outdir}` and the script
    scanned the published tree. The task was gated on every producer, but Nextflow
    publishes asynchronously: a published copy may not exist yet when the report runs,
    and a file left at the same --outdir by an earlier run can be read in its place.

    His acceptance test: the report's values come from staged inputs, and a stale or
    missing published file cannot change them. Here the staged inputs (this run) and
    the published tree (an earlier run, laid out in the task's own working directory
    and beside it) disagree on every value the page shows.
    """

    RUN = [("S_ext", "extensive_island", 1400, 19, 0.61, True),
           ("S_neg", "negative", 0, 0, 0.0, False)]
    STALE = [("S_ext", "negative", 3, 1, 0.01, False),
             ("S_neg", "extensive_island", 9999, 19, 0.99, True),
             ("S_stale_old_run", "extensive_island", 8888, 19, 0.98, True)]

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        run_tree = build_results(root / "this_run", self.RUN)
        figures = run_tree / "by_sample/S_ext/figures"
        figures.mkdir(parents=True)
        (figures / "circos.pdf").write_text("pdf")
        self.work = root / "work"
        self.manifest = stage(run_tree, self.work)
        self.expected = self.work / "expected_samples.txt"
        self.expected.write_text("S_ext\nS_neg\n")
        # The earlier run's published tree, twice: inside the task's working directory
        # (the relative paths a tree scan would use) and as a sibling --outdir.
        self.stale_roots = [self.work, root / "results"]
        for stale_root in self.stale_roots:
            build_results(stale_root, self.STALE)
            figures = stale_root / "by_sample/S_neg/figures"
            figures.mkdir(parents=True)
            (figures / "contig_validation.svg").write_text("<svg/>")
            (stale_root / "cohort/taxonomy").mkdir(parents=True)
            (stale_root / "cohort/taxonomy/pks.clb_species_support.tsv").write_text(
                "Sample\tSpecies\tTaxID\t" + "\t".join(GENES) + "\tTotal\n"
                "S_neg\tStale organism\t1\t" + "\t".join(["1"] * len(GENES)) + "\t19\n")

    def run_report(self, tag):
        html, tsv = f"{tag}.html", f"{tag}.tsv"
        run = subprocess.run(
            [sys.executable, str(SCRIPT), "--inputs-manifest", self.manifest,
             "--output", html, "--table", tsv, "--expected-samples", self.expected.name],
            cwd=self.work, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        with (self.work / tsv).open(newline="") as handle:
            rows = {r["sample"]: r for r in csv.DictReader(handle, delimiter="\t")}
        return (self.work / html).read_text(), (self.work / tsv).read_text(), rows

    def test_the_staged_basenames_repeat_across_samples(self):
        # The collision the per-file inputs/<n>/ directories exist for.
        names = [line.split("\t")[2].rsplit("/", 1)[1] for line in
                 (self.work / self.manifest).read_text().splitlines()[1:]
                 if "\tread_evidence\t" in line]
        self.assertEqual(names, ["read_evidence.tsv", "read_evidence.tsv"])

    def test_every_value_comes_from_the_staged_inputs(self):
        html, _, rows = self.run_report("report")
        self.assertEqual(set(rows), {"S_ext", "S_neg"})
        self.assertEqual(rows["S_ext"]["read_evidence"], "extensive_island")
        self.assertEqual(rows["S_ext"]["pks_reads"], "1400")
        self.assertEqual(rows["S_ext"]["clbS"], str(1400 // 19))
        self.assertEqual(rows["S_ext"]["assembler_agreement"], "concordant")
        self.assertEqual(rows["S_neg"]["read_evidence"], "negative")
        self.assertEqual(rows["S_neg"]["pks_reads"], "0")

    def test_a_stale_published_file_the_run_did_not_stage_is_never_read(self):
        html, _, rows = self.run_report("report")
        # S_neg was not assembled this run; the stale tree says it was.
        self.assertEqual(rows["S_neg"]["final_structural_evidence"], "")
        self.assertNotIn("S_stale_old_run", html)
        self.assertNotIn("Stale organism", html)
        self.assertNotIn("Read-level taxonomic support: species", html)
        # Figures: the run made a circos plot for S_ext and nothing for S_neg.
        self.assertIn("../by_sample/S_ext/figures/circos.pdf", html)
        self.assertNotIn("../by_sample/S_neg/figures/contig_validation.svg", html)

    def test_a_missing_published_tree_changes_nothing(self):
        import shutil
        with_stale = self.run_report("with_stale")
        for stale_root in self.stale_roots:
            for top in ("by_sample", "cohort"):
                shutil.rmtree(stale_root / top)
        without = self.run_report("without")
        self.assertEqual(with_stale[0], without[0])
        self.assertEqual(with_stale[1], without[1])

    def test_an_expected_sample_with_nothing_staged_is_a_row_not_a_gap(self):
        # The run named it; it produced nothing. It must not borrow the stale tree's
        # values for itself either.
        self.expected.write_text("S_ext\nS_neg\nS_stale_old_run\n")
        _, _, rows = self.run_report("report")
        self.assertEqual(rows["S_stale_old_run"]["read_evidence"], "not_classified")
        self.assertEqual(rows["S_stale_old_run"]["pks_reads"], "0")

    def test_two_staged_files_of_one_kind_for_one_sample_is_an_error(self):
        manifest = self.work / self.manifest
        line = next(l for l in manifest.read_text().splitlines()
                    if l.startswith("S_ext\tread_evidence\t"))
        manifest.write_text(manifest.read_text() + line + "\n")
        run = subprocess.run(
            [sys.executable, str(SCRIPT), "--inputs-manifest", self.manifest,
             "--output", "x.html"], cwd=self.work, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("two staged read_evidence files for S_ext", run.stderr)

    def test_there_is_no_results_flag_left_to_pass(self):
        run = subprocess.run(
            [sys.executable, str(SCRIPT), "--results", str(self.work),
             "--output", "x.html"], cwd=self.work, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("--inputs-manifest", run.stderr)


class P0_5bWiring(unittest.TestCase):
    """The task is handed the files; the script is handed a manifest naming them."""

    def process(self, name):
        return PLOTTING.split(f"process {name} {{", 1)[1].split("\nprocess ", 1)[0]

    def test_neither_cohort_task_points_a_script_at_the_published_tree(self):
        for name in ("cohortReport", "masterSummary"):
            script = self.process(name).split('"""', 1)[1]
            with self.subTest(process=name):
                self.assertNotIn("params.outdir", script)
                self.assertNotIn("--results", script)
                self.assertIn('--inputs-manifest "${inputs_manifest}"', script)

    def test_every_input_file_gets_its_own_staging_directory(self):
        for name in ("cohortReport", "masterSummary"):
            with self.subTest(process=name):
                self.assertIn("stageAs: 'inputs/*/*'", self.process(name))
                self.assertIn("path(inputs_manifest)", self.process(name))
        # ...and the manifest names the same inputs/<n>/ directories, 1-based, in the
        # order of the one sorted list both the file list and the manifest come from.
        self.assertIn("inputs/${i + 1}/${row[2].name}", MAIN)

    def test_every_kind_main_nf_tags_is_one_a_script_reads(self):
        import re
        sys.path.insert(0, str(ROOT / "scripts"))
        import build_master_summary as summary
        tagged = set(re.findall(r"tuple\((?:sampleID|''),\s*'([a-z_]+)'", MAIN))
        # 26 Sep mode memo: the evidence kind is chosen per sample type by evidenceKind().
        helper = MAIN[MAIN.index("def evidenceKind"):]
        helper = helper[:helper.index("\n}\n")]
        tagged |= set(re.findall(r"'([a-z_]+_evidence)'", helper))
        known = (set(report.SAMPLE_INPUTS) | set(report.COHORT_INPUTS)
                 | set(report.FIGURE_INPUTS) | set(summary.PUBLISHED)
                 | set(summary.COHORT_PUBLISHED))
        self.assertTrue(tagged)
        self.assertEqual(tagged - known, set(), "main.nf stages a kind nothing reads")
        report_kinds = (set(report.SAMPLE_INPUTS) | set(report.COHORT_INPUTS)
                        | set(report.FIGURE_INPUTS))
        self.assertEqual(report_kinds - tagged, set(),
                         "the report reads a kind main.nf never stages")
        # contig_mobility has no producer today (see build_master_summary.PUBLISHED).
        self.assertEqual(set(summary.PUBLISHED) - tagged, {"contig_mobility"})
        self.assertIn("qc_summary", tagged)

