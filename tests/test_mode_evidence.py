"""Mode-specific evidence schemas, and the mode field everywhere they are joined.

Ludmil, 26 Sep mode memo:
  section 3 -- WES is off-target DNA evidence only: no WGS tier, its own summary.
  section 4 -- RNA is expressed clb sequence: `clb_rna_reads`, `rna_clb_genes_detected`,
               never genomic positivity fields; breadth is expression breadth.
  section 6 -- every per-sample evidence file, the master summary and the cohort report
               carry `mode`; a field name must not imply genomic presence for WES/RNA.

The measurements are the same three in every mode (clb-assigned reads, clb genes with a
read, island breadth); only their names and what they may be read as differ. These tests
hold the producer (classify_tumor_pks_evidence.py), the module that runs it, main.nf's
report-manifest kind, and both consumers to one vocabulary.
"""
import csv
import importlib.util
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
MAIN = (ROOT / "main.nf").read_text()
TARGETED = (ROOT / "Modules/pks_targeted.nf").read_text()
PLOTTING = (ROOT / "Modules/plotting.nf").read_text()
GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(SCRIPTS))
    spec.loader.exec_module(module)
    return module


classify = load("classify_tumor_pks_evidence", SCRIPTS / "classify_tumor_pks_evidence.py")
report = load("build_cohort_report", SCRIPTS / "build_cohort_report.py")
summary = load("build_master_summary", SCRIPTS / "build_master_summary.py")

FROZEN = {"multi_gene": (5, 3, .01), "broad_island": (30, 8, .075),
          "extensive_island": (100, 10, .15)}


def write_inputs(directory, reads=500, genes=12, depth=5, length=100):
    d = Path(directory)
    (d / "qc.tsv").write_text(f"Sample\tMetric\tValue\nS\treads_clb_genes_align\t{reads}\n")
    (d / "counts.tsv").write_text(
        "# featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tbam\n"
        + "".join(f"{g}\tc\t1\t2\t+\t2\t50\n" for g in GENES[:genes]))
    (d / "depth.tsv").write_text("".join(f"c\t{i}\t{depth}\n" for i in range(1, length + 1)))
    return d


def run_classify(directory, *extra):
    d = Path(directory)
    out = d / "out.tsv"
    subprocess.run([sys.executable, str(SCRIPTS / "classify_tumor_pks_evidence.py"),
                    "--sample", "S", "--qc", str(d / "qc.tsv"), "--counts", str(d / "counts.tsv"),
                    "--depth", str(d / "depth.tsv"), "--island-length", "100",
                    "--output", str(out), *extra], check=True)
    with out.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames, rows[0]


class Schemas(unittest.TestCase):
    def test_rna_schema_is_the_memo_s(self):
        self.assertEqual(classify.SCHEMAS["rna"],
                         ["sample", "mode", "clb_rna_reads", "rna_clb_genes_detected",
                          "rna_expression_breadth_1x", "rna_expression_breadth_2x",
                          "rna_expression_breadth_3x", "rna_library_type",
                          # memo section 4 RNA QC (docs/running/wes_rna.md)
                          "clb_transcripts_detected", "rna_residual_human_fraction",
                          "rna_residual_checked_reads"])

    def test_wes_schema(self):
        self.assertEqual(classify.SCHEMAS["wes"],
                         ["sample", "mode", "clb_wes_reads", "wes_clb_genes_detected",
                          "wes_island_breadth_1x", "wes_island_breadth_2x",
                          "wes_island_breadth_3x"])

    def test_every_schema_carries_mode_second(self):
        for mode, fields in classify.SCHEMAS.items():
            with self.subTest(mode=mode):
                self.assertEqual(fields[:2], ["sample", "mode"])

    def test_wes_and_rna_share_no_field_with_the_genomic_schema(self):
        genomic = set(classify.SCHEMAS["wgs"]) - {"sample", "mode"}
        for mode in ("wes", "rna"):
            with self.subTest(mode=mode):
                self.assertEqual(set(classify.SCHEMAS[mode]) & genomic, set())
                self.assertNotIn("read_evidence", classify.SCHEMAS[mode])

    def test_modes_are_the_four_sample_types(self):
        self.assertEqual(set(classify.MODES), {"wgs", "metagenome", "wes", "rna"})


class ProducerOutput(unittest.TestCase):
    def test_wes_writes_evidence_without_a_tier(self):
        with tempfile.TemporaryDirectory() as d:
            fields, row = run_classify(write_inputs(d), "--mode", "wes")
        self.assertEqual(fields, classify.SCHEMAS["wes"])
        self.assertEqual(row["mode"], "wes")
        self.assertEqual(row["clb_wes_reads"], "500")
        self.assertEqual(row["wes_clb_genes_detected"], "12")
        self.assertEqual(row["wes_island_breadth_1x"], "1.000000")

    def test_rna_writes_expression_evidence_and_the_library_type(self):
        with tempfile.TemporaryDirectory() as d:
            fields, row = run_classify(write_inputs(d, depth=2), "--mode", "rna",
                                       "--rna-library-type", "polyA")
        self.assertEqual(fields, classify.SCHEMAS["rna"])
        self.assertEqual(row["mode"], "rna")
        self.assertEqual(row["clb_rna_reads"], "500")
        self.assertEqual(row["rna_clb_genes_detected"], "12")
        self.assertEqual(row["rna_expression_breadth_2x"], "1.000000")
        self.assertEqual(row["rna_expression_breadth_3x"], "0.000000")
        self.assertEqual(row["rna_library_type"], "polyA")

    def test_rna_library_type_defaults_to_unknown(self):
        with tempfile.TemporaryDirectory() as d:
            _, row = run_classify(write_inputs(d), "--mode", "rna")
        self.assertEqual(row["rna_library_type"], "unknown")

    def test_an_unknown_library_type_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            write_inputs(d)
            with self.assertRaises(subprocess.CalledProcessError):
                run_classify(d, "--mode", "rna", "--rna-library-type", "ribo")

    def test_a_wes_or_rna_sample_with_tier_level_evidence_still_gets_no_tier(self):
        # Evidence that would be extensive_island under WGS: no tier in any other mode.
        for mode in ("wes", "rna", "metagenome"):
            row = classify.evidence_row(mode, "S", 500, 19, (0.9, 0.8, 0.7), FROZEN)
            with self.subTest(mode=mode):
                self.assertNotIn("extensive_island", row.values())
        self.assertEqual(classify.evidence_row("wgs", "S", 500, 19, (0.9, 0.8, 0.7),
                                               FROZEN)["read_evidence"], "extensive_island")

    def test_wgs_is_the_default_mode(self):
        with tempfile.TemporaryDirectory() as d:
            _, row = run_classify(write_inputs(d))
        self.assertEqual((row["mode"], row["read_evidence"]), ("wgs", "extensive_island"))


class ModuleAndManifestWiring(unittest.TestCase):
    def block(self, name):
        start = TARGETED.index(f"process {name}")
        return TARGETED[start:TARGETED.index("\n}\n", start)]

    def test_the_module_runs_the_sample_type_as_the_mode(self):
        block = self.block("classifyPksReadEvidence")
        self.assertIn('"--mode ${params.sample_type}"', block)
        self.assertIn('--rna-library-type ${params.rna_library_type}', block)

    def test_each_mode_publishes_its_own_file_name(self):
        block = self.block("classifyPksReadEvidence")
        expr = "${params.sample_type.toString() in ['wes', 'rna'] ? params.sample_type : 'read'}_evidence.tsv"
        self.assertEqual(block.count(expr), 2, "the output declaration and the script must agree")
        for mode, kind in classify.EVIDENCE_KIND.items():
            with self.subTest(mode=mode):
                prefix = mode if mode in ("wes", "rna") else "read"
                self.assertEqual(f"{prefix}_evidence", kind)

    def test_main_nf_tags_the_same_kind_the_producer_names(self):
        helper = MAIN[MAIN.index("def evidenceKind"):]
        helper = helper[:helper.index("\n}\n")]
        mapping = dict(re.findall(r"(\w+): '(\w+)'", helper))
        self.assertEqual(mapping, {m: k for m, k in classify.EVIDENCE_KIND.items()
                                   if k != "read_evidence"})
        self.assertIn("'read_evidence')", helper)
        self.assertIn("tuple(sampleID, evidence_kind, evidenceFile)", MAIN)

    def test_both_consumers_read_every_evidence_kind(self):
        for kind in set(classify.EVIDENCE_KIND.values()):
            with self.subTest(kind=kind):
                self.assertEqual(report.SAMPLE_INPUTS[kind], f"{kind}.tsv")
                self.assertEqual(summary.PUBLISHED[kind], f"{kind}.tsv")

    def test_both_consumers_are_given_the_sample_type(self):
        for process in ("cohortReport", "masterSummary"):
            start = PLOTTING.index(f"process {process}")
            block = PLOTTING[start:PLOTTING.index("\n}\n", start)]
            with self.subTest(process=process):
                self.assertIn('--sample-type "${params.sample_type}"', block)

    def test_the_two_consumers_label_modes_identically(self):
        for mode, tier in (("wgs", "extensive_island"), ("wgs", "untiered"),
                           ("metagenome", "untiered"), ("wes", None), ("rna", None),
                           ("NA", None)):
            with self.subTest(mode=mode, tier=tier):
                self.assertEqual(report.evidence_label(mode, tier),
                                 summary.evidence_label(mode, tier))

    def test_labels_say_not_tiered_for_wes_and_rna(self):
        self.assertEqual(summary.evidence_label("rna"), "RNA expression evidence (not tiered)")
        self.assertEqual(summary.evidence_label("wes"), "WES off-target DNA evidence (not tiered)")


def stage(directory, rows):
    """Write each (sample, kind, header, values) as its own staged input + manifest."""
    d = Path(directory)
    lines = ["sample\tkind\tpath"]
    for n, (sample, kind, header, values) in enumerate(rows, start=1):
        target = d / "inputs" / str(n) / f"{kind}.tsv"
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            for v in values:
                w.writerow(v)
        lines.append(f"{sample}\t{kind}\tinputs/{n}/{kind}.tsv")
    manifest = d / "manifest.tsv"
    manifest.write_text("\n".join(lines) + "\n")
    return manifest


def mixed_cohort(directory):
    """One sample per mode, each with evidence that would be a high WGS tier."""
    matrix_header = ["Gene", "W", "M", "X", "R"]
    matrix = [[g, "10", "10", "10", "10"] for g in GENES]
    qc_header = ["Sample", "reads_clb_genes_align"]
    return stage(directory, [
        ("W", "read_evidence", classify.SCHEMAS["wgs"],
         [["W", "wgs", "extensive_island", "190", "19", "0.9", "0.8", "0.7"]]),
        ("M", "read_evidence", classify.SCHEMAS["metagenome"],
         [["M", "metagenome", "untiered", "190", "19", "0.9", "0.8", "0.7"]]),
        ("X", "wes_evidence", classify.SCHEMAS["wes"],
         [["X", "wes", "190", "19", "0.9", "0.8", "0.7"]]),
        ("R", "rna_evidence", classify.SCHEMAS["rna"],
         [["R", "rna", "190", "19", "0.9", "0.8", "0.7", "total"]]),
        ("", "gene_counts", matrix_header, matrix),
        ("", "qc_summary", qc_header, [[s, "190"] for s in "WMXR"]),
    ])


class CohortReportModes(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.manifest = mixed_cohort(self.dir)

    def build(self, *extra):
        html, table = self.dir / "r.html", self.dir / "r.tsv"
        run = subprocess.run([sys.executable, str(SCRIPTS / "build_cohort_report.py"),
                              "--inputs-manifest", str(self.manifest), "--output", str(html),
                              "--table", str(table), *extra],
                             capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        with table.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = {r["sample"]: r for r in reader}
        return html.read_text(), reader.fieldnames, rows

    def test_the_table_has_mode_and_label_right_after_sample(self):
        _, fields, _ = self.build()
        self.assertEqual(fields[:3], ["sample", "mode", "evidence_label"])

    def test_every_row_carries_its_mode(self):
        _, _, rows = self.build()
        self.assertEqual({s: r["mode"] for s, r in rows.items()},
                         {"W": "wgs", "M": "metagenome", "X": "wes", "R": "rna"})

    def test_rna_and_wes_are_labelled_as_their_own_evidence_never_a_tier(self):
        page, _, rows = self.build()
        self.assertEqual(rows["R"]["evidence_label"], "RNA expression evidence (not tiered)")
        self.assertEqual(rows["X"]["evidence_label"], "WES off-target DNA evidence (not tiered)")
        for sample in ("R", "X", "M"):
            with self.subTest(sample=sample):
                self.assertNotIn(rows[sample]["read_evidence"], report.TIERS)
        self.assertEqual(rows["W"]["read_evidence"], "extensive_island")
        self.assertIn("RNA expression evidence (not tiered)", page)
        self.assertIn("WES off-target DNA evidence (not tiered)", page)

    def test_the_measurements_come_from_the_mode_s_own_columns(self):
        _, _, rows = self.build()
        for sample in "WMXR":
            with self.subTest(sample=sample):
                self.assertEqual(rows[sample]["pks_reads"], "190")
                self.assertEqual(rows[sample]["genes_detected"], "19")
                self.assertEqual(rows[sample]["island_breadth_1x"], "0.900000")

    def test_only_the_wgs_row_is_counted_in_a_tier_tile(self):
        page, _, _ = self.build()
        tiles = re.findall(r'<div class="n">(\d+)</div><div class="l">([^<]+)</div>', page)
        counted = {label: int(n) for n, label in tiles}
        self.assertEqual(counted["Extensive"], 1)
        self.assertEqual(counted["RNA expression evidence (not tiered)"], 1)
        self.assertEqual(counted["WES off-target DNA evidence (not tiered)"], 1)
        self.assertEqual(counted["Metagenome read evidence (not tiered)"], 1)

    def test_a_sample_with_no_evidence_takes_the_run_s_sample_type(self):
        (self.dir / "expected.txt").write_text("W\nM\nX\nR\nLOST\n")
        _, _, rows = self.build("--expected-samples", str(self.dir / "expected.txt"),
                                "--sample-type", "rna")
        self.assertEqual(rows["LOST"]["mode"], "rna")
        self.assertEqual(rows["LOST"]["read_evidence"], "not_classified")
        self.assertEqual(rows["X"]["mode"], "wes")

    def test_two_modes_for_one_sample_is_an_error(self):
        manifest = stage(self.dir / "dup", [
            ("S", "read_evidence", classify.SCHEMAS["wgs"],
             [["S", "wgs", "negative", "0", "0", "0", "0", "0"]]),
            ("S", "rna_evidence", classify.SCHEMAS["rna"],
             [["S", "rna", "1", "1", "0.1", "0", "0", "unknown"]])])
        run = subprocess.run([sys.executable, str(SCRIPTS / "build_cohort_report.py"),
                              "--inputs-manifest", str(manifest),
                              "--output", str(self.dir / "dup.html")],
                             capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("more than one mode", run.stderr)


class MasterSummaryModes(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.manifest = mixed_cohort(self.dir)

    def build(self, *extra):
        out = self.dir / "summary.tsv"
        run = subprocess.run([sys.executable, str(SCRIPTS / "build_master_summary.py"),
                              "--inputs-manifest", str(self.manifest), "--output", str(out),
                              *extra], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        with out.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            return reader.fieldnames, {r["sample"]: r for r in reader}

    def test_mode_and_label_follow_sample(self):
        fields, _ = self.build()
        self.assertEqual(fields[:3], ["sample", "mode", "evidence_label"])

    def test_wes_and_rna_have_their_own_columns_and_no_tier(self):
        fields, rows = self.build()
        for column in summary.WES_EVIDENCE_COLUMNS + summary.RNA_EVIDENCE_COLUMNS:
            self.assertIn(column, fields)
        self.assertEqual(rows["R"]["clb_rna_reads"], "190")
        self.assertEqual(rows["R"]["rna_library_type"], "total")
        self.assertEqual(rows["R"]["read_evidence"], "NA")
        self.assertEqual(rows["R"]["evidence_label"], "RNA expression evidence (not tiered)")
        self.assertEqual(rows["X"]["clb_wes_reads"], "190")
        self.assertEqual(rows["X"]["read_evidence"], "NA")
        self.assertEqual(rows["X"]["evidence_label"], "WES off-target DNA evidence (not tiered)")
        self.assertEqual(rows["W"]["clb_rna_reads"], "NA")
        self.assertEqual(rows["W"]["read_evidence"], "extensive_island")
        self.assertEqual(rows["M"]["read_evidence"], "untiered")

    def test_mode_columns_are_omitted_only_when_nothing_could_say(self):
        manifest = stage(self.dir / "qc_only", [("", "qc_summary", ["Sample", "reads_clb_genes_align"],
                                                 [["S", "3"]])])
        out = self.dir / "qc_only.tsv"
        for extra, expected in (([], False), (["--sample-type", "wes"], True)):
            subprocess.run([sys.executable, str(SCRIPTS / "build_master_summary.py"),
                            "--inputs-manifest", str(manifest), "--output", str(out), *extra],
                           check=True, capture_output=True)
            with out.open() as fh:
                row = next(csv.DictReader(fh, delimiter="\t"))
            with self.subTest(extra=extra):
                self.assertEqual("mode" in row, expected)
                if expected:
                    self.assertEqual(row["mode"], "wes")
                    self.assertEqual(row["evidence_label"],
                                     "WES off-target DNA evidence (not tiered)")


if __name__ == "__main__":
    unittest.main()
