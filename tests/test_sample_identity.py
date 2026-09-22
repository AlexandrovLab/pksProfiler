"""Sample identifiers come from the sample sheet and survive to the master table.

F03, from Ludmil's 2026-09-19 report. Identifiers were rebuilt from filenames with
chained unanchored str.replace(), so `case.counts.txt` and `case.txt.counts.txt` --
produced by the two legal, distinct sheet identifiers `case` and `case.txt` -- both
reduced to `case`, pd.concat wrote two columns of that name, and the merge exited 0.

The pipeline now carries the identifier alongside the file in a manifest. These tests
cover both mergers on both routes, plus the static wiring that keeps the identifier
from being dropped in main.nf again.
"""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MERGE = ROOT / "scripts/mergeGeneCounts.py"
HMM = ROOT / "scripts/build_hmm_matrix.py"
GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]


def write_counts(path, value):
    rows = ["Geneid\tChr\tStart\tEnd\tStrand\tLength\tX"]
    rows += [f"{gene}\tctg\t1\t10\t+\t10\t{value}" for gene in GENES]
    Path(path).write_text("\n".join(rows) + "\n")


def write_hmm_counts(path, value):
    rows = ["Gene\tCount"] + [f"{gene}\t{value}" for gene in GENES]
    Path(path).write_text("\n".join(rows) + "\n")


def run(args, cwd):
    return subprocess.run([sys.executable, *[str(a) for a in args]],
                          cwd=cwd, capture_output=True, text=True)


class MergeGeneCountsIdentity(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.d = Path(self.dir.name)
        self.addCleanup(self.dir.cleanup)
        # the two files the colliding pair of sheet identifiers produces
        write_counts(self.d / "case.counts.txt", 1)
        write_counts(self.d / "case.txt.counts.txt", 2)

    def manifest(self, rows, name="m.tsv"):
        path = self.d / name
        path.write_text("sample_id\tfile\n" + "".join(f"{s}\t{f}\n" for s, f in rows))
        return path

    def test_case_and_case_txt_stay_two_columns(self):
        m = self.manifest([("case", "case.counts.txt"), ("case.txt", "case.txt.counts.txt")])
        out = self.d / "merged.tsv"
        result = run([MERGE, "--manifest", m, "--output", out], self.d)
        self.assertEqual(result.returncode, 0, result.stderr)
        header = out.read_text().splitlines()[0].split("\t")
        self.assertEqual(header, ["Gene", "case", "case.txt"])

    def test_a_genuine_duplicate_is_refused_and_nothing_is_written(self):
        m = self.manifest([("case", "case.counts.txt"), ("case", "case.txt.counts.txt")])
        out = self.d / "merged.tsv"
        result = run([MERGE, "--manifest", m, "--output", out], self.d)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Duplicate sample identifier", result.stderr)
        self.assertFalse(out.exists())

    def test_column_order_does_not_depend_on_row_order(self):
        forward = self.manifest([("case", "case.counts.txt"), ("case.txt", "case.txt.counts.txt")], "f.tsv")
        reverse = self.manifest([("case.txt", "case.txt.counts.txt"), ("case", "case.counts.txt")], "r.tsv")
        outputs = []
        for manifest, name in ((forward, "a.tsv"), (reverse, "b.tsv")):
            out = self.d / name
            result = run([MERGE, "--manifest", manifest, "--output", out], self.d)
            self.assertEqual(result.returncode, 0, result.stderr)
            outputs.append(out.read_text())
        self.assertEqual(outputs[0], outputs[1])

    def test_a_manifest_naming_an_absent_file_is_refused(self):
        m = self.manifest([("case", "not_here.counts.txt")])
        out = self.d / "merged.tsv"
        result = run([MERGE, "--manifest", m, "--output", out], self.d)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(out.exists())

    def test_legacy_filenames_strip_one_anchored_suffix(self):
        # Without a manifest the name still comes off the filename, but once and from
        # the end -- so the two files remain two samples.
        out = self.d / "merged.tsv"
        result = run([MERGE, "case.counts.txt", "case.txt.counts.txt", out], self.d)
        self.assertEqual(result.returncode, 0, result.stderr)
        header = out.read_text().splitlines()[0].split("\t")
        self.assertEqual(header, ["Gene", "case", "case.txt"])


class HmmMatrixIdentity(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.d = Path(self.dir.name)
        self.addCleanup(self.dir.cleanup)
        write_hmm_counts(self.d / "case.hmm_counts.tsv", 1)
        write_hmm_counts(self.d / "case.txt.hmm_counts.tsv", 2)

    def test_manifest_identifiers_are_used_verbatim(self):
        m = self.d / "m.tsv"
        m.write_text("sample_id\tfile\ncase\tcase.hmm_counts.tsv\ncase.txt\tcase.txt.hmm_counts.tsv\n")
        out = self.d / "matrix.tsv"
        result = run([HMM, "--manifest", m, "--out", out], self.d)
        self.assertEqual(result.returncode, 0, result.stderr)
        header = out.read_text().splitlines()[0].split("\t")
        self.assertEqual(header[1:], ["case", "case.txt"])

    def test_duplicate_identifiers_are_refused(self):
        m = self.d / "m.tsv"
        m.write_text("sample_id\tfile\ncase\tcase.hmm_counts.tsv\ncase\tcase.txt.hmm_counts.tsv\n")
        out = self.d / "matrix.tsv"
        result = run([HMM, "--manifest", m, "--out", out], self.d)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Duplicate sample identifier", result.stderr)
        self.assertFalse(out.exists())

    def test_manifest_and_inputs_are_mutually_exclusive(self):
        out = self.d / "matrix.tsv"
        result = run([HMM, "--out", out], self.d)
        self.assertNotEqual(result.returncode, 0)


class WiringTests(unittest.TestCase):
    """Static: the identifier must not be dropped on the way to the merge again."""

    MAIN = (ROOT / "main.nf").read_text()
    PLOTTING = (ROOT / "Modules/plotting.nf").read_text()

    def test_sample_ids_are_trimmed_where_the_tuples_are_built(self):
        self.assertNotIn("tuple(row.patient,", self.MAIN)
        self.assertGreaterEqual(self.MAIN.count("row.patient.toString().trim()"), 2)

    def test_the_merge_gets_a_manifest_not_a_bare_file_list(self):
        for channel in ("ALIGN_MERGE", "HMM_MERGE"):
            self.assertIn(f"{channel}.rows.collectFile", self.MAIN)
        self.assertNotIn("--strip-suffix", self.PLOTTING)
        self.assertEqual(self.PLOTTING.count("--manifest"), 2)

    def test_neither_master_table_interpolates_one_argument_per_sample(self):
        # F02: a filename per sample on the command line overflows the argument limit.
        for process in ("masterTableAlign", "masterTableHMM"):
            block = self.PLOTTING.split(f"process {process} {{", 1)[1].split("\n}", 1)[0]
            self.assertNotIn("${inputs}", block)
            self.assertIn("path(manifest)", block)


if __name__ == "__main__":
    unittest.main()
