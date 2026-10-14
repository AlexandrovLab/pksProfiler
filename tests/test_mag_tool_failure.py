"""rc2: a CheckM2/GTDB-Tk task that failed and was ignored must not crash the run or
drop the sample.

The v1.0.0-rc1 Mayo metagenome run (--enable_mags, ~37 h) stopped with
`Invalid method invocation 'call' with arguments: [<sample>, null,
[bin.18.locus_evidence.tsv, ...]]`: several gtdbtkClassify tasks exited 1 under
errorStrategy 'ignore', never emitted a summary, and the remainder:true join that
attaches per-bin locus evidence then produced a left-orphaned 3-element row the next
6-parameter closure could not destructure.

Execution half of the guard: tests/check_mag_channel_regression.sh (real Nextflow,
no conda). This file covers build_mag_summary.py's handling of the stand-in tables
stubFailedMagQuality writes, and the static wiring in Modules/pks_mag.nf.
"""
import csv
import importlib.util
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKS_MAG = (ROOT / "Modules/pks_mag.nf").read_text()
sys.path.insert(0, str(ROOT / "scripts"))


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CHECKM2_OK = ("Name\tCompleteness\tContamination\tGenome_Size\tContig_N50\n"
              "bin.18\t90.0\t1.0\t4600000\t50000\n"
              "bin.19\t80.0\t2.0\t3000000\t20000\n")
# stubFailedMagQuality's CheckM2 stand-in: header plus Status, no rows (so rc1's parser,
# and any reader that ignores Status, sees every bin as unmeasured).
CHECKM2_FAILED = "Name\tCompleteness\tContamination\tGenome_Size\tContig_N50\tStatus\n"
GTDBTK_OK = ("user_genome\tclassification\n"
             "bin.18\td__Bacteria;o__Enterobacterales;g__Escherichia\n"
             "bin.19\td__Bacteria;o__Bacteroidales\n")
GTDBTK_FAILED = ("user_genome\tclassification\n"
                 "bin.18\tgtdbtk_failed\n"
                 "bin.19\tgtdbtk_failed\n")
LOCUS_HEADER = "sample\tbin_id\tlocus_tier\tlocus_genes_detected\tlocus_genes\tlocus_breadth\n"


class BuildMagSummaryWithFailedTools(unittest.TestCase):
    def run_summary(self, checkm2_text, gtdbtk_text):
        mod = load("build_mag_summary")
        with tempfile.TemporaryDirectory() as tmp:
            checkm2 = os.path.join(tmp, "checkm2.tsv")
            gtdbtk = os.path.join(tmp, "gtdbtk.tsv")
            Path(checkm2).write_text(checkm2_text)
            Path(gtdbtk).write_text(gtdbtk_text)
            for d in ("tblout", "context", "locus"):
                os.makedirs(os.path.join(tmp, d))
            Path(tmp, "locus", "bin.18.locus_evidence.tsv").write_text(
                LOCUS_HEADER + "S1\tbin.18\tbroad_island\t12\tclbA\t0.6\n")
            Path(tmp, "locus", "bin.19.locus_evidence.tsv").write_text(
                LOCUS_HEADER + "S1\tbin.19\tnone\t0\t\t0.0\n")
            out = os.path.join(tmp, "summary.tsv")
            argv = sys.argv
            sys.argv = ["build_mag_summary.py", "--checkm2", checkm2, "--gtdbtk", gtdbtk,
                        "--tblout_dir", os.path.join(tmp, "tblout"),
                        "--context_dir", os.path.join(tmp, "context"),
                        "--locus_dir", os.path.join(tmp, "locus"),
                        "--sample", "S1", "--out", out]
            try:
                mod.main()
            finally:
                sys.argv = argv
            with open(out) as fh:
                return {r["bin_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    def test_both_tools_ok(self):
        rows = self.run_summary(CHECKM2_OK, GTDBTK_OK)
        self.assertEqual(rows["bin.18"]["mag_qc_status"], "ok")
        self.assertEqual(rows["bin.18"]["completeness"], "90.0")

    def test_gtdbtk_failed_is_named_per_bin_and_positivity_unchanged(self):
        rows = self.run_summary(CHECKM2_OK, GTDBTK_FAILED)
        self.assertEqual(set(rows), {"bin.18", "bin.19"})
        for row in rows.values():
            self.assertEqual(row["mag_qc_status"], "gtdbtk_failed")
            self.assertEqual(row["taxonomy"], "gtdbtk_failed")
            # unknown taxonomy is not an anomalous one
            self.assertEqual(row["unexpected_taxon_flag"], "False")
        # positivity is the locus tier, untouched by the taxonomy failure
        self.assertEqual(rows["bin.18"]["locus_tier"], "broad_island")
        self.assertEqual(rows["bin.19"]["locus_tier"], "none")

    def test_checkm2_failed_reads_as_na_not_zero(self):
        rows = self.run_summary(CHECKM2_FAILED, GTDBTK_OK)
        self.assertEqual(rows["bin.18"]["completeness"], "NA")
        self.assertEqual(rows["bin.18"]["genome_size"], "NA")
        self.assertEqual(rows["bin.18"]["mag_qc_status"], "checkm2_failed")

    def test_a_header_only_real_report_is_not_a_failure(self):
        # CheckM2 header with no rows and no Status column: unmeasured, not failed.
        rows = self.run_summary("Name\tCompleteness\tContamination\tGenome_Size\tContig_N50\n", GTDBTK_OK)
        self.assertEqual(rows["bin.18"]["completeness"], "NA")
        self.assertEqual(rows["bin.18"]["mag_qc_status"], "ok")

    def test_both_failed(self):
        rows = self.run_summary(CHECKM2_FAILED, GTDBTK_FAILED)
        self.assertEqual(rows["bin.19"]["mag_qc_status"], "checkm2_failed;gtdbtk_failed")


class PksMagWiring(unittest.TestCase):
    def workflow(self, name):
        start = PKS_MAG.index(f"workflow {name} {{")
        nxt = re.search(r"\nworkflow \w+ \{", PKS_MAG[start + 1:])
        return PKS_MAG[start:] if nxt is None else PKS_MAG[start:start + 1 + nxt.start()]

    def test_stub_process_exists_unpublished_and_without_conda(self):
        start = PKS_MAG.index("process stubFailedMagQuality {")
        block = PKS_MAG[start:PKS_MAG.index("\n}\n", start)]
        self.assertNotIn("publishDir", block)
        self.assertNotIn("conda", block)
        self.assertIn("Contig_N50\\\\tStatus", block)
        self.assertIn("gtdbtk_failed", block)

    def test_pksmag_routes_both_join_chains_through_the_guarded_subworkflows(self):
        wf = self.workflow("pksMAG")
        self.assertIn("magQualityFallback(metabat2Bin.out.bins, checkm2Predict.out.report, "
                      "gtdbtkClassify.out.summary)", wf)
        self.assertIn(".mix(magQualityFallback.out.checkm2_stub)", wf)
        self.assertIn(".mix(magQualityFallback.out.gtdbtk_stub)", wf)
        self.assertIn("summary_input_ch = magSummaryInput.out.rows", wf)
        self.assertIn("community_input_ch = magCommunityInput.out.rows", wf)
        # the rc1 crash site: a multi-parameter closure right after a remainder join
        self.assertNotIn("locus_evidence_per_sample_ch, by: 0, remainder: true)\n        .map { sampleID", wf)

    def test_closures_after_remainder_joins_take_one_parameter(self):
        for name in ("magSummaryInput", "magCommunityInput", "magQualityFallback"):
            wf = self.workflow(name)
            tail = wf[wf.index("remainder: true"):]
            for closure in re.findall(r"\.(?:map|filter) \{([^}]*?)->", tail):
                with self.subTest(workflow=name, closure=closure):
                    self.assertEqual(closure.strip().count(","), 0, closure)


if __name__ == "__main__":
    unittest.main()
