"""Cohort aggregation passes a list, not one command-line argument per sample.

F02, from Ludmil's 2026-09-19 report. Every cohort-wide step built its argument list
in Groovy -- one filename per sample interpolated into .command.sh -- which overflows
the OS argument limit at cohort scale and fails as "Argument list too long", after the
per-sample work is already done.

Each of our own aggregators now takes one filename naming the inputs. combine_mpa.py
is KrakenTools' and has no such option, so that step builds its list with find and
checks the size against ARG_MAX before calling, which turns the overflow into a
message that says what happened.
"""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLOTTING = (ROOT / "Modules/plotting.nf").read_text()
TAXA = (ROOT / "Modules/pks_taxa.nf").read_text()
MAIN = (ROOT / "main.nf").read_text()
GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]

# process name -> (module text, what it must use instead of an argument list)
COHORT_PROCESSES = {
    "masterTableAlign": (PLOTTING, "--manifest"),
    "masterTableHMM": (PLOTTING, "--manifest"),
    "masterQCSummary": (PLOTTING, "--inputs-file"),
    "combineClbTaxonomySupport": (TAXA, "--species-files-from"),
    "process_bracken": (TAXA, "find . -maxdepth 1"),
}


def block_for(process, text):
    self_block = text.split(f"process {process} {{", 1)
    assert len(self_block) == 2, f"process {process} not found"
    return self_block[1].split("\n}", 1)[0]


def run(args, cwd):
    return subprocess.run([sys.executable, *[str(a) for a in args]],
                          cwd=cwd, capture_output=True, text=True)


class NoPerSampleArguments(unittest.TestCase):
    def test_no_cohort_process_joins_its_inputs_into_the_command(self):
        for process, (text, _) in COHORT_PROCESSES.items():
            with self.subTest(process=process):
                self.assertNotIn("join(' ')", block_for(process, text))

    def test_each_cohort_process_reads_a_list_instead(self):
        for process, (text, expected) in COHORT_PROCESSES.items():
            with self.subTest(process=process):
                self.assertIn(expected, block_for(process, text))

    def test_the_lists_are_built_next_to_the_files_in_main(self):
        for channel in ("ALIGN_MERGE", "HMM_MERGE", "QC_MERGE", "CLB_SUPPORT_MERGE"):
            with self.subTest(channel=channel):
                self.assertIn(f"{channel}.rows.collectFile", MAIN)

    def test_the_krakentools_step_checks_the_limit_before_calling(self):
        block = block_for("process_bracken", TAXA)
        self.assertIn("getconf ARG_MAX", block)
        self.assertIn("combine_mpa.py", block)


class ListFileInputs(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.d = Path(self.dir.name)
        self.addCleanup(self.dir.cleanup)

    def test_qc_summary_reads_a_list_file(self):
        for sample, value in (("S1", 10), ("S2", 20)):
            (self.d / f"{sample}.qc.tsv").write_text(
                f"Sample\tMetric\tValue\n{sample}\textracted_unmapped_reads\t{value}\n")
        (self.d / "qc.list").write_text("S1.qc.tsv\nS2.qc.tsv\n")
        (self.d / "expected.txt").write_text("S1\nS2\n")
        out = self.d / "qc.summary.tsv"
        result = run([ROOT / "scripts/build_qc_summary.py", "--inputs-file", "qc.list",
                      "--expected-samples", "expected.txt", "--output", out], self.d)
        self.assertEqual(result.returncode, 0, result.stderr)
        samples = [line.split("\t")[0] for line in out.read_text().splitlines()[1:]]
        self.assertEqual(samples, ["S1", "S2"])

    def test_qc_summary_refuses_a_list_naming_an_absent_file(self):
        (self.d / "qc.list").write_text("not_staged.qc.tsv\n")
        out = self.d / "qc.summary.tsv"
        result = run([ROOT / "scripts/build_qc_summary.py", "--inputs-file", "qc.list",
                      "--output", out], self.d)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(out.exists())

    def test_qc_summary_refuses_both_input_forms_at_once(self):
        (self.d / "qc.list").write_text("")
        result = run([ROOT / "scripts/build_qc_summary.py", "--inputs", "a.tsv",
                      "--inputs-file", "qc.list", "--output", "out.tsv"], self.d)
        self.assertNotEqual(result.returncode, 0)

    def test_species_support_reads_a_list_file(self):
        columns = ["Species", "TaxID", *GENES, "Total"]
        for sample, value in (("S1", 1), ("S2", 2)):
            row = ["Escherichia coli", "562", *[str(value)] * len(GENES), str(value * len(GENES))]
            (self.d / f"{sample}.clb_species_support.tsv").write_text(
                "\t".join(columns) + "\n" + "\t".join(row) + "\n")
        (self.d / "clb.list").write_text(
            "S1.clb_species_support.tsv\nS2.clb_species_support.tsv\n")
        out = self.d / "combined.tsv"
        result = run([ROOT / "scripts/combine_clb_species_support.py",
                      "--species-files-from", "clb.list", "--species-output", out], self.d)
        self.assertEqual(result.returncode, 0, result.stderr)
        samples = [line.split("\t")[0] for line in out.read_text().splitlines()[1:]]
        self.assertEqual(samples, ["S1", "S2"])

    def test_species_support_requires_exactly_one_input_form(self):
        result = run([ROOT / "scripts/combine_clb_species_support.py",
                      "--species-output", "out.tsv"], self.d)
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
