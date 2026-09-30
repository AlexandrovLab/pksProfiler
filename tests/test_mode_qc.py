"""WES and RNA QC (Ludmil, 26 Sep mode memo sections 3 and 4).

Section 3: "Report total, unmapped, host-depleted and post-fastp reads, duplicate rate,
and capture kit/off-target fraction where available. WES duplicate behavior must be
checked separately from WGS."
Section 4: "total, unmapped, host-depleted, and residual-human read counts" and
transcription-specific outputs such as `clb_transcripts_detected`.

Definitions are in docs/running/wes_rna.md. These tests hold the three helpers
(wes_input_qc.py, pks_duplicate_rate.py, rna_residual_human.py), the evidence classifier's
RNA columns, the mode-aware QC table, the master summary / cohort report, and main.nf's
wiring to them. Synthetic inputs only.
"""
import csv
import gzip
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
MODE_QC = (ROOT / "Modules/mode_qc.nf").read_text()
TARGETED = (ROOT / "Modules/pks_targeted.nf").read_text()
GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]


def load(name):
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


wes = load("wes_input_qc")
dup = load("pks_duplicate_rate")
residual = load("rna_residual_human")
qc = load("build_qc_summary")
classify = load("classify_tumor_pks_evidence")
report = load("build_cohort_report")

# `samtools flagstat -O tsv` of 20 primary mapped records, 5 flagged 0x400 (samtools 1.21).
FLAGSTAT = """20\t0\ttotal (QC-passed reads + QC-failed reads)
20\t0\tprimary
0\t0\tsecondary
0\t0\tsupplementary
4\t1\tduplicates
4\t1\tprimary duplicates
20\t0\tmapped
100.00%\tN/A\tmapped %
20\t0\tprimary mapped
100.00%\tN/A\tprimary mapped %
0\t0\tpaired in sequencing
"""


def sam(name, flag, pos, cigar, contig="NC_017628.1"):
    return f"{name}\t{flag}\t{contig}\t{pos}\t42\t{cigar}\t*\t0\t0\t*\t*\n"


def process_block(text, name):
    start = text.index(f"process {name} {{")
    nxt = re.search(r"\nprocess \w+ \{", text[start + 1:])
    return text[start:] if nxt is None else text[start:start + 1 + nxt.start()]


class WesInputQC(unittest.TestCase):
    def write(self, d, text=FLAGSTAT):
        path = Path(d) / "flagstat.tsv"
        path.write_text(text)
        return path

    def test_counts_qc_failed_records_too(self):
        with tempfile.TemporaryDirectory() as d:
            rows = dict(wes.input_metrics(self.write(d)))
        self.assertEqual(rows["wes_input_primary_mapped"], 20)
        self.assertEqual(rows["wes_input_primary_duplicates"], 5)
        self.assertEqual(rows["wes_capture_status"], "capture_bed_not_provided")
        self.assertNotIn("wes_capture_on_target", rows)

    def test_with_a_capture_bed(self):
        with tempfile.TemporaryDirectory() as d:
            rows = dict(wes.input_metrics(self.write(d), on_target=12))
        self.assertEqual(rows["wes_capture_on_target"], 12)
        self.assertEqual(rows["wes_capture_status"], "capture_bed_provided")

    def test_impossible_counts_are_refused(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(SystemExit):
                wes.input_metrics(self.write(d), on_target=21)
            with self.assertRaises(SystemExit):
                wes.input_metrics(self.write(d, FLAGSTAT.replace("4\t1\tduplicates", "30\t0\tduplicates")))
            with self.assertRaises(SystemExit):
                wes.input_metrics(self.write(d, "not flagstat\n"))

    def test_cli_writes_a_qc_fragment(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "on").write_text("7\n")
            out = Path(d) / "out.tsv"
            subprocess.run([sys.executable, str(SCRIPTS / "wes_input_qc.py"), "--sample", "S",
                            "--flagstat", str(self.write(d)), "--on-target", str(Path(d) / "on"),
                            "--output", str(out)], check=True)
            lines = out.read_text().splitlines()
        self.assertEqual(lines[0], "Sample\tMetric\tValue")
        self.assertIn("S\twes_capture_on_target\t7", lines)
        self.assertEqual({line.split("\t")[1] for line in lines[1:]}, set(wes.METRICS))


class PksDuplicateRate(unittest.TestCase):
    def test_unclipped_five_prime(self):
        self.assertEqual(dup.five_prime(0, 100, "100M"), 100)
        self.assertEqual(dup.five_prime(0, 100, "5S95M"), 95)
        self.assertEqual(dup.five_prime(16, 100, "100M"), 199)
        self.assertEqual(dup.five_prime(16, 100, "90M2D8M10S"), 209)

    def test_position_collapse(self):
        lines = ["@HD\tVN:1.6\n",
                 sam("a", 0, 100, "150M"), sam("a_dup", 0, 100, "150M"),
                 # soft-clipped copy of the same fragment start: still a duplicate
                 sam("a_clip", 0, 103, "3S147M"),
                 # same position, other strand: not a duplicate
                 sam("b", 16, 100, "150M"),
                 sam("c", 16, 200, "150M"), sam("c_dup", 16, 210, "10D140M"),
                 # secondary / supplementary / unmapped are never counted
                 sam("a", 256, 100, "150M"), sam("a", 2048, 100, "150M"), sam("u", 4, 0, "*")]
        self.assertEqual(dup.duplicate_counts(lines), (6, 2))

    def test_no_reads(self):
        self.assertEqual(dup.duplicate_counts([]), (0, 0))

    def test_cli(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out.tsv"
            subprocess.run([sys.executable, str(SCRIPTS / "pks_duplicate_rate.py"), "--sample", "S",
                            "--output", str(out)], input=sam("a", 0, 1, "10M") * 3, text=True, check=True)
            self.assertEqual(out.read_text().splitlines()[1:],
                             ["S\twes_pks_island_reads\t3", "S\twes_pks_island_duplicates\t2"])


class RnaResidualHuman(unittest.TestCase):
    def fastq(self, d, n):
        path = Path(d) / "reads.fastq.gz"
        with gzip.open(path, "wt") as fh:
            for i in range(n):
                fh.write(f"@r{i}\nACGT\n+\nIIII\n")
        return path

    def names(self, path):
        return [line[1:].strip() for i, line in enumerate(Path(path).read_text().splitlines()) if i % 4 == 0]

    def test_step(self):
        self.assertEqual(residual.step_for(0, 200000), 1)
        self.assertEqual(residual.step_for(200000, 200000), 1)
        self.assertEqual(residual.step_for(200001, 200000), 2)
        self.assertEqual(residual.step_for(1000000, 200000), 5)

    def test_systematic_sample_respects_the_cap(self):
        with tempfile.TemporaryDirectory() as d:
            total, sampled = residual.subsample(self.fastq(d, 10), 3, Path(d) / "s.fq")
            names = self.names(Path(d) / "s.fq")
        self.assertEqual((total, sampled), (10, 3))
        self.assertEqual(names, ["r0", "r4", "r8"])   # every 4th read, spread across the file

    def test_below_the_cap_everything_is_checked(self):
        with tempfile.TemporaryDirectory() as d:
            total, sampled = residual.subsample(self.fastq(d, 4), 200000, Path(d) / "s.fq")
        self.assertEqual((total, sampled), (4, 4))

    def test_empty_input(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(residual.subsample(self.fastq(d, 0), 5, Path(d) / "s.fq"), (0, 0))

    def test_mapped_reads_counts_primary_mapped_names(self):
        lines = ["@SQ\tSN:c\tLN:10\n", "a\t0\tc\t1\t0\t4M\n", "a\t2048\tc\t1\t0\t4M\n",
                 "b\t4\t*\t0\t0\t*\n", "d\t256\tc\t1\t0\t4M\n", "e\t16\tc\t1\t0\t4M\n"]
        self.assertEqual(residual.mapped_reads(lines), 2)

    def test_cli_round_trip(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            subprocess.run([sys.executable, str(SCRIPTS / "rna_residual_human.py"), "subsample",
                            "--input", str(self.fastq(d, 6)), "--cap", "4", "--output", str(d / "s.fq"),
                            "--stats", str(d / "stats")], check=True)
            subprocess.run([sys.executable, str(SCRIPTS / "rna_residual_human.py"), "count",
                            "--sample", "S", "--stats", str(d / "stats"), "--output", str(d / "qc.tsv")],
                           input="r0\t0\tc\t1\t0\t4M\n", text=True, check=True)
            lines = (d / "qc.tsv").read_text().splitlines()
        self.assertEqual(lines[1:], ["S\trna_residual_input_reads\t6",
                                     "S\trna_residual_checked_reads\t3",
                                     "S\trna_residual_human_reads\t1"])

    def test_more_mapped_than_checked_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "stats").write_text("1\t1\n")
            result = subprocess.run([sys.executable, str(SCRIPTS / "rna_residual_human.py"), "count",
                                     "--sample", "S", "--stats", str(d / "stats"),
                                     "--output", str(d / "qc.tsv")],
                                    input="a\t0\tc\t1\t0\t4M\nb\t0\tc\t1\t0\t4M\n", text=True,
                                    capture_output=True)
        self.assertNotEqual(result.returncode, 0)


class RnaEvidenceColumns(unittest.TestCase):
    def run_classify(self, d, counts, *extra):
        d = Path(d)
        (d / "qc.tsv").write_text(f"Sample\tMetric\tValue\nS\treads_clb_genes_align\t{sum(counts)}\n")
        (d / "counts.tsv").write_text(
            "# featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tbam\n"
            + "".join(f"{g}\tc\t1\t2\t+\t2\t{n}\n" for g, n in zip(GENES, counts)))
        (d / "depth.tsv").write_text("".join(f"c\t{i}\t1\n" for i in range(1, 11)))
        subprocess.run([sys.executable, str(SCRIPTS / "classify_tumor_pks_evidence.py"),
                        "--sample", "S", "--qc", str(d / "qc.tsv"), "--counts", str(d / "counts.tsv"),
                        "--depth", str(d / "depth.tsv"), "--island-length", "10", "--mode", "rna",
                        "--output", str(d / "out.tsv"), *extra], check=True)
        with (d / "out.tsv").open() as fh:
            return next(csv.DictReader(fh, delimiter="\t"))

    def test_clb_transcripts_detected_is_genes_at_the_read_floor(self):
        counts = [5, 3, 2, 1, 0, 9]
        with tempfile.TemporaryDirectory() as d:
            row = self.run_classify(d, counts)
            self.assertEqual(row["rna_clb_genes_detected"], "5")       # >= 1 read
            self.assertEqual(row["clb_transcripts_detected"], "3")     # >= 3 reads (default)
            row = self.run_classify(d, counts, "--rna-transcript-min-reads", "6")
            self.assertEqual(row["clb_transcripts_detected"], "1")

    def test_residual_human_fraction_from_the_fragment(self):
        with tempfile.TemporaryDirectory() as d:
            frag = Path(d) / "residual.tsv"
            frag.write_text("Sample\tMetric\tValue\nS\trna_residual_input_reads\t900\n"
                            "S\trna_residual_checked_reads\t400\nS\trna_residual_human_reads\t10\n")
            row = self.run_classify(d, [1], "--residual-human-qc", str(frag))
        self.assertEqual(row["rna_residual_human_fraction"], "0.025000")
        self.assertEqual(row["rna_residual_checked_reads"], "400")

    def test_na_without_the_check_or_with_nothing_left(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(self.run_classify(d, [1])["rna_residual_human_fraction"], "NA")
            frag = Path(d) / "residual.tsv"
            frag.write_text("Sample\tMetric\tValue\nS\trna_residual_input_reads\t0\n"
                            "S\trna_residual_checked_reads\t0\nS\trna_residual_human_reads\t0\n")
            row = self.run_classify(d, [1], "--residual-human-qc", str(frag))
        self.assertEqual(row["rna_residual_human_fraction"], "NA")
        self.assertEqual(row["rna_residual_checked_reads"], "0")

    def test_no_other_mode_gains_the_rna_columns(self):
        for mode in ("wgs", "metagenome", "wes"):
            with self.subTest(mode=mode):
                self.assertNotIn("clb_transcripts_detected", classify.SCHEMAS[mode])
                self.assertNotIn("rna_residual_human_fraction", classify.SCHEMAS[mode])


def base_metrics(**extra):
    values = dict(filter_input_reads=100, reads_after_fastp=100, reads_after_hg38=90,
                  reads_after_t2t_phix=80)
    values.update(extra)
    return values


class ModeAwareQcTable(unittest.TestCase):
    def write(self, metrics, stages):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "qc.tsv"
            qc.write_summary(metrics, out, stages=stages)
            with out.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                return reader.fieldnames, {r["Sample"]: r for r in reader}

    def test_wgs_table_is_unchanged(self):
        fields, _ = self.write({"S": base_metrics()}, ["read_qc", "host_depletion"])
        self.assertEqual(fields, qc.OUTPUT_COLUMNS)

    def test_wes_columns_only_when_wes_stages_are_expected(self):
        metrics = {"S": base_metrics(wes_input_primary_mapped=200, wes_input_primary_duplicates=50,
                                     wes_capture_on_target=150, wes_capture_status="capture_bed_provided",
                                     wes_pks_island_reads=10, wes_pks_island_duplicates=2)}
        fields, rows = self.write(metrics, ["wes_input_qc", "read_qc", "host_depletion", "wes_pks_duplicates"])
        self.assertEqual(fields[fields.index("status") - len(qc.WES_QC_COLUMNS):fields.index("status")],
                         qc.WES_QC_COLUMNS)
        self.assertNotIn("rna_residual_human_fraction", fields)
        row = rows["S"]
        self.assertEqual(row["wes_input_duplicate_fraction"], "0.250000")
        self.assertEqual(row["wes_capture_on_target_fraction"], "0.750000")
        self.assertEqual(row["wes_capture_off_target_fraction"], "0.250000")
        self.assertEqual(row["wes_pks_island_duplicate_fraction"], "0.200000")
        self.assertEqual(row["terminal_state"], "completed_untiered")

    def test_capture_bed_not_provided(self):
        metrics = {"S": base_metrics(wes_input_primary_mapped=200, wes_input_primary_duplicates=0,
                                     wes_capture_status="capture_bed_not_provided",
                                     wes_pks_island_reads=0, wes_pks_island_duplicates=0)}
        _, rows = self.write(metrics, ["wes_input_qc", "read_qc", "host_depletion", "wes_pks_duplicates"])
        self.assertEqual(rows["S"]["wes_capture_off_target_fraction"], "capture_bed_not_provided")
        self.assertEqual(rows["S"]["wes_input_duplicate_fraction"], "0.000000")
        self.assertEqual(rows["S"]["wes_pks_island_duplicate_fraction"], "NA")   # no island reads

    def test_fastq_wes_has_no_input_qc(self):
        metrics = {"S": base_metrics(wes_pks_island_reads=4, wes_pks_island_duplicates=0)}
        fields, rows = self.write(metrics, ["read_qc", "host_depletion", "wes_pks_duplicates"])
        self.assertIn("wes_input_duplicate_fraction", fields)
        self.assertEqual(rows["S"]["wes_input_duplicate_fraction"], "NA")
        self.assertEqual(rows["S"]["wes_capture_off_target_fraction"], "NA")

    def test_missing_wes_input_qc_is_a_failed_stage(self):
        _, rows = self.write({"S": base_metrics()}, ["wes_input_qc", "read_qc", "host_depletion"])
        self.assertEqual(rows["S"]["terminal_state"], "failed_at_wes_input_qc")

    def test_impossible_wes_counts(self):
        for bad in (dict(wes_input_primary_mapped=10, wes_input_primary_duplicates=11),
                    dict(wes_input_primary_mapped=10, wes_capture_on_target=11),
                    dict(wes_pks_island_reads=1, wes_pks_island_duplicates=2)):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                qc.output_row("S", base_metrics(**bad))

    def test_rna_columns(self):
        metrics = {"S": base_metrics(rna_residual_input_reads=80, rna_residual_checked_reads=80,
                                     rna_residual_human_reads=4)}
        fields, rows = self.write(metrics, ["read_qc", "host_depletion", "rna_residual_human"])
        self.assertIn("rna_residual_human_fraction", fields)
        self.assertNotIn("wes_input_duplicate_fraction", fields)
        self.assertEqual(rows["S"]["rna_residual_human_fraction"], "0.050000")

    def test_the_residual_check_must_sample_the_final_depleted_reads(self):
        with self.assertRaises(ValueError):
            qc.output_row("S", base_metrics(rna_residual_input_reads=81, rna_residual_checked_reads=10,
                                            rna_residual_human_reads=0))
        with self.assertRaises(ValueError):
            qc.output_row("S", base_metrics(rna_residual_input_reads=80, rna_residual_checked_reads=10,
                                            rna_residual_human_reads=11))
        qc.output_row("S", base_metrics(reads_after_pangenome=70, rna_residual_input_reads=70,
                                        rna_residual_checked_reads=70, rna_residual_human_reads=0))

    def test_an_unknown_capture_status_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            frag = Path(d) / "f.tsv"
            frag.write_text("Sample\tMetric\tValue\nS\twes_capture_status\tmaybe\n")
            with self.assertRaises(ValueError):
                qc.load_fragments([frag])

    def test_inferred_stages_add_the_mode_columns(self):
        metrics = {"S": base_metrics(rna_residual_input_reads=80, rna_residual_checked_reads=8,
                                     rna_residual_human_reads=0)}
        self.assertIn("rna_residual_human", qc.infer_stages(metrics))


class Consumers(unittest.TestCase):
    def test_master_summary_carries_the_mode_qc(self):
        summary = load("build_master_summary")
        self.assertEqual(summary.WES_QC_COLUMNS, qc.WES_QC_COLUMNS)
        self.assertIn("clb_transcripts_detected", summary.RNA_EVIDENCE_COLUMNS)
        self.assertIn("rna_residual_human_fraction", summary.RNA_EVIDENCE_COLUMNS)
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "qc.tsv").write_text("Sample\t" + "\t".join(qc.WES_QC_COLUMNS) + "\n"
                                      "S\t200\t0.250000\tcapture_bed_not_provided\tcapture_bed_not_provided\t10\t0.200000\n")
            (d / "wes.tsv").write_text("sample\tmode\tclb_wes_reads\twes_clb_genes_detected\twes_island_breadth_1x\n"
                                       "S\twes\t10\t3\t0.010000\n")
            (d / "manifest.tsv").write_text("sample\tkind\tpath\n\tqc_summary\tqc.tsv\nS\twes_evidence\twes.tsv\n")
            subprocess.run([sys.executable, str(SCRIPTS / "build_master_summary.py"), "--inputs-manifest",
                            str(d / "manifest.tsv"), "--output", str(d / "out.tsv")], check=True,
                           capture_output=True)
            with (d / "out.tsv").open() as fh:
                row = next(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(row["wes_input_duplicate_fraction"], "0.250000")
        self.assertEqual(row["wes_capture_off_target_fraction"], "capture_bed_not_provided")

    def test_master_summary_without_wes_qc_has_no_wes_columns(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "qc.tsv").write_text("Sample\treads_after_hg38\nS\t10\n")
            (d / "manifest.tsv").write_text("sample\tkind\tpath\n\tqc_summary\tqc.tsv\n")
            subprocess.run([sys.executable, str(SCRIPTS / "build_master_summary.py"), "--inputs-manifest",
                            str(d / "manifest.tsv"), "--output", str(d / "out.tsv")], check=True,
                           capture_output=True)
            header = (d / "out.tsv").read_text().splitlines()[0].split("\t")
        self.assertNotIn("wes_input_duplicate_fraction", header)

    def test_cohort_report_rna_rows(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            fields = classify.SCHEMAS["rna"]
            row = classify.evidence_row("rna", "R", 12, 4, (0.1, 0.05, 0.0), {}, transcripts=2,
                                        residual=("0.010000", 100))
            (d / "rna.tsv").write_text("\t".join(fields) + "\n" + "\t".join(str(row[f]) for f in fields) + "\n")
            (d / "manifest.tsv").write_text("sample\tkind\tpath\nR\trna_evidence\trna.tsv\n")
            subprocess.run([sys.executable, str(SCRIPTS / "build_cohort_report.py"), "--inputs-manifest",
                            str(d / "manifest.tsv"), "--output", str(d / "r.html"),
                            "--table", str(d / "t.tsv")], check=True, capture_output=True)
            with (d / "t.tsv").open() as fh:
                table = next(csv.DictReader(fh, delimiter="\t"))
            page = (d / "r.html").read_text()
        self.assertEqual(table["clb_transcripts_detected"], "2")
        self.assertEqual(table["rna_residual_human_fraction"], "0.010000")
        self.assertIn("2 clb transcripts detected", page)


class Wiring(unittest.TestCase):
    def test_params(self):
        self.assertIn("params.wes_capture_bed   = null", MAIN)
        self.assertIn("params.rna_transcript_min_reads        = 3", MAIN)
        self.assertIn("params.rna_residual_human_sample_reads = 200000", MAIN)
        self.assertIn('params.rna_residual_human_args         = "-n 2 -m 20 -s 30"', MAIN)
        self.assertRegex(MAIN, r"wes_capture_bed\s*:\s*Provenance\.code\(params\.wes_capture_bed\)")

    def test_wes_qc_runs_for_wes_only(self):
        self.assertRegex(MAIN, r'if \(params\.sample_type == "wes"\) \{\s*wesInputQC\(sample_sheet\)')
        self.assertRegex(MAIN, r'if \(params\.sample_type == "wes"\) \{\s*wesPksDuplicates\(ALIGN_OUT\.profile')
        self.assertIn("wesInputQC.out.qc.map", MAIN)
        self.assertIn("wesPksDuplicates.out.qc.map", MAIN)
        # wesInputQC sits inside the alignment-input branch: a FASTQ sheet has no alignment.
        alignment_branch = MAIN[MAIN.index("    if (alignment_input) {"):MAIN.index('} else if (input_data_type == "fastq")')]
        self.assertIn("wesInputQC(sample_sheet)", alignment_branch)

    def test_rna_residual_runs_for_rna_and_reaches_the_evidence(self):
        self.assertRegex(MAIN, r'if \(params\.sample_type == "rna"\) \{\s*rnaResidualHuman\(MAPPED_READS\)')
        self.assertIn('params.sample_type == "rna" ? rna_residual_qc', MAIN)
        self.assertIn(".join(evidence_mode_qc, by: 0)", MAIN)
        block = process_block(TARGETED, "classifyPksReadEvidence")
        self.assertIn("path(mode_qc, stageAs: 'mode_qc/*')", block)
        self.assertIn("--residual-human-qc ${mode_qc}", block)
        self.assertIn("--rna-transcript-min-reads ${params.rna_transcript_min_reads}", block)

    def test_residual_check_uses_the_human_pass_plus_relaxed_floors(self):
        block = process_block(MODE_QC, "rnaResidualHuman")
        self.assertIn("HostDepletion.humanPass('rna'", block)
        self.assertIn("params.rna_residual_human_args", block)
        self.assertIn("HostDepletion.isLong(params.read_type) ? ''", block)
        self.assertIn("--cap ${params.rna_residual_human_sample_reads}", block)

    def test_wes_input_qc_counts_primary_mapped_records_once(self):
        block = process_block(MODE_QC, "wesInputQC")
        self.assertEqual(block.count("samtools view -@"), 2)   # BED / no-BED branches, one decode each
        self.assertIn("-u -F 2308", block)
        self.assertIn('samtools view -c -L "\\$CAPTURE_BED" on_target.fifo', block)
        self.assertIn("samtools flagstat -O tsv -", block)

    def test_pks_duplicates_read_the_island_only(self):
        block = process_block(MODE_QC, "wesPksDuplicates")
        self.assertIn('"${params.pks_contig}:${params.pks_start_1based}-${params.pks_end_1based}"', block)
        self.assertIn("samtools view -F 2308", block)

    def test_run_report_and_preflight_know_the_new_params(self):
        provenance = (SCRIPTS / "finalize_provenance.py").read_text()
        for name in ("wes_capture_bed", "rna_transcript_min_reads",
                     "rna_residual_human_sample_reads", "rna_residual_human_args"):
            with self.subTest(name=name):
                self.assertIn(f'"{name}"', provenance)
        self.assertIn('"--wes_capture_bed"', (SCRIPTS / "preflight.py").read_text())


if __name__ == "__main__":
    unittest.main()
