"""What a consumer reads must be where a process publishes it.

build_cohort_report.py read by_sample/<sample>/contigs/final_pks_evidence.tsv while
pks_targeted.nf publishes it to contigs/final_evidence/. Nothing errored: the reader
returned no rows, `assembled` became False, and the report showed no contig columns for
samples whose assembly had in fact succeeded. The only symptom was an absence, noticed by
eye on 2026-09-22.

The suite did not catch it because the fixture in test_cohort_report.py wrote to the
reader's path rather than the pipeline's -- the same mistake as M4, where a CheckM2
fixture was written to match the reader and hid a real column-name error for as long as
it existed.

This is the second path-contract defect in this work; cohortReport wrote its report one
directory below its own output declaration. A publishDir and the code that reads from it
are a contract with nothing enforcing it, so it is asserted here.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = (ROOT / "scripts/build_cohort_report.py").read_text()
TARGETED = (ROOT / "Modules/pks_targeted.nf").read_text()
FIXTURE = (ROOT / "tests/test_cohort_report.py").read_text()


class ReadersAgreeWithPublishDir(unittest.TestCase):
    def test_final_evidence_is_read_where_it_is_published(self):
        published = re.findall(r'publishDir\s*\{\s*"\$\{params\.sample_dir\}/\$\{sampleID\}/'
                               r'(contigs/[a-z_]+)"', TARGETED)
        self.assertIn("contigs/final_evidence", published,
                      "pks_targeted.nf no longer publishes to contigs/final_evidence; "
                      "update the reader and this test together")
        self.assertIn("contigs/final_evidence/final_pks_evidence.tsv", REPORT,
                      "build_cohort_report.py must read the published path")
        self.assertNotIn('"contigs/final_pks_evidence.tsv"', REPORT)

    def test_the_fixture_uses_the_published_layout(self):
        # A fixture written to match the reader cannot catch a reader that is wrong.
        self.assertIn("contigs/final_evidence/final_pks_evidence.tsv", FIXTURE)

    def test_every_sample_dir_path_the_report_reads_is_published_somewhere(self):
        modules = "\n".join(p.read_text() for p in (ROOT / "Modules").glob("*.nf"))
        reads = set(re.findall(r'by_sample / sample / "([^"]+)"', REPORT))
        missing = []
        for rel in sorted(reads):
            leaf = rel.rsplit("/", 1)[-1]
            directory = rel.rsplit("/", 1)[0] if "/" in rel else ""
            if directory and directory not in modules:
                missing.append(f"{rel} (no publishDir mentions {directory})")
            elif leaf.split(".", 1)[-1] not in modules and leaf not in modules:
                missing.append(f"{rel} (no process emits {leaf})")
        self.assertEqual(missing, [], "\n" + "\n".join(missing))


if __name__ == "__main__":
    unittest.main()
