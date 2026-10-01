"""L4 (Ludmil, 1 Oct): QC/provenance of E1's AS<30 rescue at every stage.

Extraction already recorded extract_low_score_mode, extracted_low_score_mapped_reads and
extract_low_score_no_as_records, and input_header.tsv the aligner and -T. What was missing:
how many of the rescued reads survive each later stage. The reads are tagged by name, the
least invasive way that leaves every output unchanged when the rule is off:

  extractReads           writes <sample>.low_score_rescued.names (the FASTQ names, /1 /2
                         from the mate bits) -- an optional output, absent when nothing
                         was rescued
  mapReads               counts those names on the streams it already reads: its input
                         (after fastp), the hg38 pass's output, the T2T+phiX pass's output
  lowScoreRescueProfile  records in the profiling BAM, and featureCounts-assigned to clb
                         with the pipeline's own annotation filter and arguments

QC summary columns low_score_rescued_* (only when some sample has them), RUN_REPORT's
LOW-SCORE RESCUE (L4) section. The live check is the bwa -T 0 run of
tests/check_e2e_execution.sh (tests/fixtures/assert_low_score_outputs.py).
"""
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAP = (ROOT / "Modules/map_reads.nf").read_text()
EXTRACT = (ROOT / "Modules/extract_reads.nf").read_text()
ALIGN = (ROOT / "Modules/pksProfiler_align.nf").read_text()
MAIN = (ROOT / "main.nf").read_text()
sys.path.insert(0, str(ROOT / "scripts"))
import build_qc_summary as bqs  # noqa: E402
import finalize_provenance as fp  # noqa: E402

PROFILE = EXTRACT.split("process lowScoreRescueProfile", 1)[1]


def counter_awk(stage):
    """mapReads' counter for one pass, as bash receives it."""
    start = MAP.index(f"    awk -v names=\"${{rescue}}\" -v hits={stage}.rescued")
    end = MAP.index(f"< {stage}.count.fifo", start)
    text = MAP[start:end].replace("\\$", "$").replace("\\\n", "\n")
    return text.split("'", 1)[1].rsplit("'", 1)[0]


def fastq(names):
    return "".join(f"@{n} extra\nACGT\n+\nIIII\n" for n in names)


def count(stage, stream_names, rescued):
    with tempfile.TemporaryDirectory() as tmp:
        names = ""
        if rescued is not None:
            (Path(tmp) / "names").write_text("".join(f"{n}\n" for n in rescued))
            names = "names"
        result = subprocess.run(["awk", "-v", f"names={names}", "-v", f"hits={stage}.rescued",
                                 counter_awk(stage)], input=fastq(stream_names), cwd=tmp,
                                capture_output=True, text=True, check=True)
        hits = Path(tmp) / f"{stage}.rescued"
        return result.stdout, hits.read_text() if hits.exists() else None


class MapReadsCountsByName(unittest.TestCase):
    def test_hits_are_exact_names_mate_suffix_included(self):
        for stage in ("hg38", "t2t"):
            with self.subTest(stage=stage):
                total, hits = count(stage, ["a/1", "a/2", "b/1", "c", "d/2"], ["a/1", "c", "zz/1"])
                self.assertEqual(total, "5\n")
                self.assertEqual(hits, "2\n")    # a/1 and c; a/2 is the other mate

    def test_without_names_the_counter_prints_exactly_what_it_did(self):
        total, hits = count("hg38", ["a/1", "a/2", "b"], None)
        self.assertEqual(total, "3\n")
        self.assertIsNone(hits)

    def test_the_rescued_metrics_are_written_only_with_names(self):
        block = MAP.split('if [[ -n "${rescue}" ]]; then\n        printf', 1)[1]
        for metric in ("low_score_rescued_after_fastp", "low_score_rescued_after_hg38",
                       "low_score_rescued_after_t2t_phix"):
            self.assertIn(metric, block.split("\n    fi", 1)[0])

    def test_the_reads_themselves_are_untouched(self):
        # The minimap2 | samtools fastq | tee | bgzip pipes carry no new stage.
        for stage_output in ('> "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz"',
                             '> "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"'):
            pipe = MAP.split(stage_output, 1)[0].rsplit("minimap2 -2", 1)[1]
            self.assertNotIn("rescue", pipe)


class ExtractionListsTheNames(unittest.TestCase):
    def test_optional_output_from_the_on_branch_only(self):
        self.assertIn('path("${sampleID}.low_score_rescued.names"), optional: true, emit: rescue_names',
                      EXTRACT)
        on = EXTRACT.split('if [[ "\\$LOW_SCORE" == "on" ]]; then', 1)[1].split("\n    else\n", 1)[0]
        self.assertIn('-v names="${sampleID}.low_score_rescued.names"', on)
        self.assertIn('> names', on)
        off = EXTRACT.split("\n    else\n", 1)[1].split("\n    fi\n", 1)[0]
        self.assertNotIn("low_score_rescued", off)


class ProfiledAndAssignedUseThePipelinesOwnRule(unittest.TestCase):
    def test_same_annotation_filter_and_featurecounts_arguments(self):
        def gff_filter(text):
            return re.search(r"awk -F '\\t' '\n\s*BEGIN \{ OFS=\"\\t\" \}\n(.*?)\n", text).group(1).strip()

        def fc_args(text):
            block = text.split("featureCounts \\", 1)[1].split("\n\n", 1)[0]
            return [line.strip() for line in block.splitlines()
                    if line.strip().startswith("-") and not line.strip().startswith("-o")]
        self.assertEqual(gff_filter(PROFILE), gff_filter(ALIGN))
        self.assertEqual(fc_args(PROFILE), fc_args(ALIGN))
        self.assertIn("featurecounts_long = HostDepletion.isLong(params.read_type) ? \"-L\" : \"\"", PROFILE)

    def test_only_the_rescued_reads_are_counted(self):
        self.assertIn('samtools view -b -N "${rescue_names}" -o rescued.bam "${bam}"', PROFILE)


class Wiring(unittest.TestCase):
    def test_main_passes_the_names_and_profiles_them(self):
        body = MAIN.split("\nworkflow {", 1)[1]
        self.assertIn("RESCUE_NAMES = EXTRACT_OUT.rescue_names", body)
        self.assertIn("mapReads(withRescueNames(FILTERED_UNMAPPED_READS, RESCUE_NAMES))", body)
        self.assertIn("lowScoreRescueProfile(RESCUE_NAMES.join(ALIGN_OUT.profile", body)
        self.assertIn("lowScoreRescueProfile.out.qc", body)
        self.assertIn("def withRescueNames(reads, names)", MAIN)
        self.assertIn("path(rescue_names, stageAs: 'rescue/*')", MAP)


class QcSummaryAndRunReport(unittest.TestCase):
    BASE = {"extracted_unmapped_reads": 10, "extracted_low_score_mapped_reads": 5,
            "extracted_reads": 15, "filter_input_reads": 15, "reads_after_fastp": 14,
            "reads_after_hg38": 12, "reads_after_t2t_phix": 12}
    RESCUE = {"low_score_rescued_after_fastp": 5, "low_score_rescued_after_hg38": 4,
              "low_score_rescued_after_t2t_phix": 4, "low_score_rescued_profiled": 3,
              "low_score_rescued_clb_assigned": 3}

    def test_columns_only_when_some_sample_rescued(self):
        self.assertFalse(set(bqs.LOW_SCORE_RESCUE_COLUMNS) & set(bqs.table_columns([], {"S": self.BASE})))
        columns = bqs.table_columns([], {"S": dict(self.BASE, **self.RESCUE), "T": self.BASE})
        at = columns.index("extract_low_score_mode") + 1
        self.assertEqual(columns[at:at + 5], bqs.LOW_SCORE_RESCUE_COLUMNS)
        self.assertEqual(bqs.output_row("T", self.BASE)["low_score_rescued_after_hg38"], "NA")

    def test_a_rescued_read_can_only_be_lost(self):
        for bad in (dict(self.RESCUE, low_score_rescued_after_fastp=6),        # > extracted
                    dict(self.RESCUE, low_score_rescued_after_t2t_phix=5),     # > after hg38
                    dict(self.RESCUE, low_score_rescued_clb_assigned=4)):      # > profiled
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                bqs.output_row("S", dict(self.BASE, **bad))

    def test_run_report_sums_the_stages(self):
        with tempfile.TemporaryDirectory() as tmp:
            qc = Path(tmp) / "qc"
            qc.mkdir()
            columns = ["Sample", "extracted_low_score_mapped_reads"] + bqs.LOW_SCORE_RESCUE_COLUMNS
            rows = ["\t".join(columns), "\t".join(["A", "5", "5", "4", "4", "3", "3"]),
                    "\t".join(["B", "2", "2", "2", "1", "1", "0"]),
                    "\t".join(["C", "0", "NA", "NA", "NA", "NA", "NA"])]
            (qc / "pks.qc.summary.tsv").write_text("\n".join(rows) + "\n")
            rollup = fp.qc_rollup(tmp)
        self.assertEqual(rollup["low_score_rescue"]["samples"], 2)
        self.assertEqual(rollup["low_score_rescue"]["totals"]["low_score_rescued_after_hg38"], 6)
        text = "\n".join(fp.render_qc_rollup(rollup))
        self.assertIn("LOW-SCORE RESCUE (L4)", text)
        self.assertRegex(text, r"assigned to a clb gene\s+3")
        self.assertRegex(text, r"extracted \(AS below --extract_min_as\)\s+7")

    def test_documented(self):
        doc = (ROOT / "docs/output.md").read_text()
        for column in bqs.LOW_SCORE_RESCUE_COLUMNS:
            self.assertIn(column.replace("low_score_rescued_after_hg38", "_after_hg38")
                          .replace("low_score_rescued_after_t2t_phix", "_after_t2t_phix"), doc)


if __name__ == "__main__":
    unittest.main()
