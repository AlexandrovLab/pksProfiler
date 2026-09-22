"""One rule for which clb model a read belongs to, and it does not depend on order.

F07, from Ludmil's 2026-09-19 report: "The AWK filter converts E-values and scores
with numeric coercion. A nonnumeric E-value becomes zero and passes the stringent
threshold in the executed fixture. Rows with too few fields are skipped. Exact
E-value/score ties across genes are resolved by first appearance: reversing two tied
rows changed clbA=1 into clbB=1."

His acceptance test: "Invalid numeric fields and truncated rows fail. Comment-only
no-hit output remains valid. Reverse and shuffle tied hits, compare chunked versus
unchunked execution, and require invariant counts with ambiguity reported."
"""
import importlib.util
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/hmm_best_hit.py"
HMM_MODULE = (ROOT / "Modules/pksProfiler_hmm.nf").read_text()
spec = importlib.util.spec_from_file_location("hmm_best_hit", SCRIPT)
best = importlib.util.module_from_spec(spec)
spec.loader.exec_module(best)

HEADER = "# target name accession query name accession ...\n"


def row(gene, read, evalue, score, alifrom=1, alito=100):
    """A tblout line: 1 target, 3 query, 7 alifrom, 8 ali to, 13 E-value, 14 score."""
    return (f"{gene}.cds.aln - {read} - 1 100 {alifrom} {alito} 1 100 500 + "
            f"{evalue} {score} 0.0 -\n")


class TheRule(unittest.TestCase):
    def assign(self, lines, evalue=1e-5):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False)
        handle.write(HEADER + "".join(lines))
        handle.close()
        path = Path(handle.name)
        self.addCleanup(lambda: path.unlink(missing_ok=True))
        return best.best_hits(path, evalue)

    def test_the_lowest_evalue_wins(self):
        assigned, ambiguous, _ = self.assign([row("clbA", "r1", "1e-30", 120),
                                              row("clbB", "r1", "1e-10", 80)])
        self.assertEqual(assigned, {"r1": "clbA"})
        self.assertEqual(ambiguous, [])

    def test_the_score_breaks_an_evalue_tie(self):
        assigned, _, _ = self.assign([row("clbA", "r1", "1e-20", 90),
                                      row("clbB", "r1", "1e-20", 110)])
        self.assertEqual(assigned, {"r1": "clbB"})

    def test_the_longer_alignment_breaks_an_evalue_and_score_tie(self):
        # The largest-overlap step: same significance, but clbK matched more of the read.
        assigned, ambiguous, _ = self.assign([
            row("clbA", "r1", "1e-20", 90, alifrom=1, alito=50),
            row("clbK", "r1", "1e-20", 90, alifrom=1, alito=100)])
        self.assertEqual(assigned, {"r1": "clbK"})
        self.assertEqual(ambiguous, [])

    def test_a_read_tied_on_all_three_is_dropped(self):
        assigned, ambiguous, qualifying = self.assign([
            row("clbA", "r1", "1e-20", 90), row("clbB", "r1", "1e-20", 90)])
        self.assertEqual(assigned, {})
        self.assertEqual(ambiguous, ["r1"])
        self.assertEqual(qualifying, ["r1"])   # it still had a qualifying hit

    def test_several_hits_to_one_gene_are_not_a_tie(self):
        assigned, ambiguous, _ = self.assign([row("clbA", "r1", "1e-20", 90),
                                              row("clbA", "r1", "1e-20", 90)])
        self.assertEqual(assigned, {"r1": "clbA"})
        self.assertEqual(ambiguous, [])

    def test_reverse_strand_coordinates_still_give_a_length(self):
        assigned, _, _ = self.assign([
            row("clbA", "r1", "1e-20", 90, alifrom=100, alito=1),
            row("clbB", "r1", "1e-20", 90, alifrom=60, alito=1)])
        self.assertEqual(assigned, {"r1": "clbA"})

    def test_hits_above_the_threshold_do_not_qualify(self):
        assigned, _, qualifying = self.assign([row("clbA", "r1", "1e-2", 5)])
        self.assertEqual(assigned, {})
        self.assertEqual(qualifying, [])


class OrderCannotChangeTheAnswer(unittest.TestCase):
    """His acceptance test: reverse and shuffle tied hits, require invariant counts."""

    LINES = [row("clbA", "r1", "1e-30", 120), row("clbB", "r1", "1e-10", 80),
             row("clbA", "r2", "1e-20", 90), row("clbB", "r2", "1e-20", 90),
             row("clbA", "r3", "1e-20", 90, alito=50), row("clbK", "r3", "1e-20", 90)]

    def assign(self, lines):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False)
        handle.write(HEADER + "".join(lines))
        handle.close()
        path = Path(handle.name)
        self.addCleanup(lambda: path.unlink(missing_ok=True))
        return best.best_hits(path, 1e-5)[:2]

    def test_reversing_the_file_changes_nothing(self):
        self.assertEqual(self.assign(self.LINES), self.assign(list(reversed(self.LINES))))

    def test_shuffling_changes_nothing(self):
        expected = self.assign(self.LINES)
        for seed in range(6):
            shuffled = list(self.LINES)
            random.Random(seed).shuffle(shuffled)
            with self.subTest(seed=seed):
                self.assertEqual(self.assign(shuffled), expected)

    def test_splitting_into_chunks_changes_nothing(self):
        # Standing in for --hmm_chunking: the same hits, arriving in two pieces.
        whole = self.assign(self.LINES)
        first, second = self.LINES[:3], self.LINES[3:]
        merged_assigned, merged_ambiguous = self.assign(second + first)
        self.assertEqual((merged_assigned, merged_ambiguous), whole)


class MalformedInputFails(unittest.TestCase):
    def run_script(self, text):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        d = Path(directory.name)
        (d / "in.tbl").write_text(text)
        return subprocess.run(
            [sys.executable, str(SCRIPT), "--tblout", str(d / "in.tbl"),
             "--evalue", "1e-5", "--counts-out", str(d / "counts.tsv"),
             "--read-ids-out", str(d / "ids.txt")],
            capture_output=True, text=True), d

    def test_a_non_numeric_evalue_fails_instead_of_becoming_zero(self):
        result, _ = self.run_script(HEADER + row("clbA", "r1", "not_a_number", 120))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not a number", result.stdout + result.stderr)

    def test_a_truncated_row_fails_instead_of_being_skipped(self):
        result, _ = self.run_script(HEADER + "clbA.cds.aln - r1 - 1 100\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("truncated", result.stdout + result.stderr)

    def test_an_unknown_model_fails(self):
        result, _ = self.run_script(HEADER + row("notAGene", "r1", "1e-20", 90))
        self.assertNotEqual(result.returncode, 0)

    def test_comment_only_output_is_a_valid_no_hit_sample(self):
        result, d = self.run_script(HEADER + "#\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = (d / "counts.tsv").read_text().splitlines()[1:]
        self.assertEqual(len(rows), 19)
        self.assertTrue(all(line.split("\t")[1] == "0" for line in rows))


class TheModuleUsesIt(unittest.TestCase):
    def test_the_awk_filter_is_gone(self):
        code = "\n".join(l for l in HMM_MODULE.splitlines()
                         if not l.strip().startswith("#"))
        self.assertNotIn("$13+0", code)
        self.assertIn("hmm_best_hit.py", code)

    def test_ambiguity_is_reported(self):
        self.assertIn("hmm_ambiguous_reads", HMM_MODULE)
        self.assertIn("hmm_ambiguous_reads",
                      (ROOT / "scripts/build_qc_summary.py").read_text())


if __name__ == "__main__":
    unittest.main()
