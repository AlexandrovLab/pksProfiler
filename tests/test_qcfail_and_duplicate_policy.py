"""v3 item 10: what extraction does with QCFAIL (0x200) and DUPLICATE (0x400).

The two flags need opposite treatment and both are easy to get wrong later:

  0x400 DUPLICATE is retained ON PURPOSE. The submitter's markdup assigned it in human
  coordinate space, where these reads -- unmapped to human, or weakly mapped -- had no
  usable coordinates, so it says nothing about duplication among them. Duplicates are
  re-marked after IHE3034 alignment instead, by compute_fragment_support_variants.py.
  Adding 0x400 to the extraction mask would apply a human-coordinate verdict to a
  microbial question and would change every published result. The review asks for that
  to be documented so nobody "fixes" it; these tests make the documentation enforceable.

  0x200 QCFAIL is retained too, but only for want of a reason: it has no standard
  meaning across pipelines. It is now counted (extracted_qcfail_reads) so a policy can
  be chosen from what cohorts actually show.

The counting arithmetic is checked by running the module's own awk over a synthetic SAM,
rather than by re-implementing it here -- a second implementation could agree with
itself and disagree with the pipeline.
"""
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXTRACT_MODULE = (ROOT / "Modules" / "extract_reads.nf").read_text()
# Only the extractReads process. The file also holds lowScoreRescueProfile, which counts
# MAPPED records in an already-rescued BAM with `-F 4` -- a legitimate use of -F that has
# nothing to do with what extraction keeps from the submitter's alignment.
EXTRACT = EXTRACT_MODULE[EXTRACT_MODULE.index("process extractReads {"):
                         EXTRACT_MODULE.index("process inputHeaderSummary {")]
QC_SUMMARY = (ROOT / "scripts" / "build_qc_summary.py").read_text()
MASTER = (ROOT / "scripts" / "build_master_summary.py").read_text()

HAVE_AWK = shutil.which("awk") is not None

# FLAG bits, so the tests read as the policy rather than as magic numbers.
UNMAPPED = 4
SECONDARY = 0x100
QCFAIL = 0x200
DUPLICATE = 0x400
SUPPLEMENTARY = 0x800


class ExtractionMask(unittest.TestCase):
    """The mask must stay secondary+supplementary only."""

    def test_the_mask_is_2304_and_nothing_more(self):
        masks = set(re.findall(r"-F (\d+)", EXTRACT))
        self.assertEqual(masks, {"2304"},
                         "extraction filters with -F 2304 (0x100 secondary + 0x800 "
                         f"supplementary) and nothing else; found {sorted(masks)}")

    def test_duplicates_are_not_excluded(self):
        # 2304 + 0x400 = 3328; 2304 + 0x200 = 2816; both = 3840. Each would silently
        # change what reaches profiling for every cohort.
        for forbidden, why in ((3328, "duplicates"), (2816, "QCFAIL"),
                               (3840, "duplicates and QCFAIL")):
            self.assertNotIn(f"-F {forbidden}", EXTRACT,
                             f"-F {forbidden} would exclude {why} at extraction. That is "
                             "not a cleanup: see the policy comment in extract_reads.nf "
                             "and v3 review item 10 before changing it.")

    def test_the_duplicate_decision_is_written_down_where_the_mask_is(self):
        # The review's actual ask: document the retention so it is not "fixed" later.
        # It has to live at the mask, not in a distant file nobody reads first.
        head = EXTRACT[:EXTRACT.index("mkfifo extract.count.fifo")]
        self.assertIn("0x400", head)
        self.assertIn("compute_fragment_support_variants.py", head,
                      "the comment should say where duplicates ARE re-marked, or the "
                      "retention reads like an oversight")
        self.assertIn("0x200", head)


@unittest.skipUnless(HAVE_AWK, "needs awk")
class QcfailCounting(unittest.TestCase):
    """Run the module's own awk over a synthetic SAM and check its tallies."""

    def extract_awk(self):
        """The classifier awk program as the module really writes it.

        The module body is a Groovy string, so a shell backslash is doubled there;
        undouble it to get what the task script actually runs.
        """
        start = EXTRACT.index("/^@/ { print; next }")
        end = EXTRACT.index("> counts }'", start) + len("> counts }")
        return EXTRACT[start:end].replace("\\\\", "\\").replace("\\$", "$")

    def run_awk(self, records, min_as=30):
        program = self.extract_awk()
        with tempfile.TemporaryDirectory() as tmp:
            sam = Path(tmp) / "in.sam"
            sam.write_text("".join(records))
            counts = Path(tmp) / "classes"
            names = Path(tmp) / "names"
            subprocess.run(
                ["awk", "-F", "\t", "-v", f"min_as={min_as}",
                 "-v", f"counts={counts}", "-v", f"names={names}", program, str(sam)],
                check=True, capture_output=True, text=True)
            return [int(v) for v in counts.read_text().split()]

    def record(self, name, flag, score=None):
        fields = [name, str(flag), "chr1", "100", "60", "10M",
                  "*", "0", "0", "ACGTACGTAC", "IIIIIIIIII"]
        if score is not None:
            fields.append(f"AS:i:{score}")
        return "\t".join(fields) + "\n"

    def test_qcfail_is_counted_among_extracted_reads(self):
        rows = [
            self.record("a", UNMAPPED),                        # extracted, clean
            self.record("b", UNMAPPED | QCFAIL),               # extracted, flagged
            self.record("c", 0, score=5),                      # extracted (low AS), clean
            self.record("d", QCFAIL, score=5),                 # extracted (low AS), flagged
        ]
        unmapped, low, no_as, other, qcf = self.run_awk(rows)
        self.assertEqual((unmapped, low), (2, 2))
        self.assertEqual(qcf, 2, "both flagged reads reached the FASTQ and must be counted")

    def test_a_duplicate_is_extracted_and_is_not_a_qcfail(self):
        rows = [self.record("a", UNMAPPED | DUPLICATE)]
        unmapped, low, no_as, other, qcf = self.run_awk(rows)
        self.assertEqual(unmapped, 1, "a duplicate-flagged read is still extracted")
        self.assertEqual(qcf, 0, "0x400 is not 0x200")

    def test_records_that_are_not_extracted_do_not_raise_the_count(self):
        # A high-AS mapped read is dropped by the awk (`other`), and a mapped read with
        # no AS tag is counted but never written. Neither reaches profiling, so neither
        # belongs in a count of extracted QCFAIL reads.
        rows = [
            self.record("keep", UNMAPPED | QCFAIL),
            self.record("highscore", QCFAIL, score=500),
            self.record("noas", QCFAIL),
        ]
        unmapped, low, no_as, other, qcf = self.run_awk(rows)
        self.assertEqual((unmapped, low, no_as, other), (1, 0, 1, 1))
        self.assertEqual(qcf, 1,
                         "only the extracted read counts; the dropped ones would "
                         "overstate how much flagged data reached profiling")

    def test_counting_does_not_change_what_is_extracted(self):
        # The whole point: visibility, not filtering. Every record that was extracted
        # before is extracted now, flagged or not.
        rows = [self.record(f"r{i}", UNMAPPED | (QCFAIL if i % 2 else 0)) for i in range(6)]
        unmapped, low, no_as, other, qcf = self.run_awk(rows)
        self.assertEqual(unmapped, 6, "every unmapped record is extracted regardless of 0x200")
        self.assertEqual(qcf, 3)


class MetricReaches(unittest.TestCase):

    def test_the_metric_is_known_to_the_qc_summary(self):
        self.assertIn('"extracted_qcfail_reads"', QC_SUMMARY)

    def test_the_metric_reaches_the_master_summary(self):
        self.assertIn('"extracted_qcfail_reads"', MASTER)

    def test_extraction_writes_it(self):
        self.assertIn("extracted_qcfail_reads", EXTRACT)

    def test_the_rescue_off_branch_says_why_it_cannot_measure(self):
        # It is blank there, which must read as a documented limit rather than a bug:
        # that branch pipes into samtools fastq, which emits no FLAGs to count.
        self.assertIn("QCFAIL_EXTRACTED=\"\"", EXTRACT)
        self.assertIn("samtools view -c -f 0x204 -F 2304", EXTRACT,
                      "the off branch should name the command that measures it on demand")


if __name__ == "__main__":
    unittest.main()
