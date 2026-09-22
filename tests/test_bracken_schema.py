"""A Bracken report the plotter does not understand is refused, not guessed at.

F11, from Ludmil's 2026-09-19 report, verbatim: "The plotting parser falls back to
another numeric column when the expected abundance column is absent and catches parse
exceptions by returning an empty frame. A file containing only name and taxonomy_id
was accepted by the plotting parser as 562 reads for Escherichia coli because 562 was
the taxonomy identifier."

That is the case these tests are built around: not a crash, a plausible-looking
abundance chart made of taxonomy ids.
"""
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "plot_bracken_taxa", ROOT / "scripts/plot_bracken_taxa.py")
plot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plot)


class TheSchemaIsRequired(unittest.TestCase):
    def report(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        handle.write(text)
        handle.close()
        path = Path(handle.name)
        self.addCleanup(lambda: path.unlink(missing_ok=True))
        return path

    def test_the_taxonomy_id_is_no_longer_plotted_as_an_abundance(self):
        # Ludmil's file, exactly: E. coli would have been drawn at 562 reads.
        path = self.report("name\ttaxonomy_id\nEscherichia coli\t562\n")
        with self.assertRaises(ValueError) as caught:
            plot.read_bracken_report(path)
        message = str(caught.exception)
        self.assertIn("missing new_est_reads", message)
        self.assertIn("Refusing to guess", message)

    def test_the_message_says_what_the_file_did_contain(self):
        path = self.report("name\ttaxonomy_id\nEscherichia coli\t562\n")
        with self.assertRaises(ValueError) as caught:
            plot.read_bracken_report(path)
        self.assertIn("taxonomy_id", str(caught.exception))

    def test_a_real_report_still_reads(self):
        path = self.report("name\ttaxonomy_id\tnew_est_reads\n"
                           "Escherichia coli\t562\t1840\n"
                           "Bacteroides fragilis\t817\t95\n")
        frame = plot.read_bracken_report(path)
        self.assertEqual(list(frame["reads"]), [1840, 95])
        self.assertEqual(list(frame["taxon"]), ["Escherichia coli", "Bacteroides fragilis"])

    def test_column_case_does_not_matter(self):
        path = self.report("Name\tNew_Est_Reads\nEscherichia coli\t7\n")
        self.assertEqual(list(plot.read_bracken_report(path)["reads"]), [7])

    def test_a_non_numeric_abundance_column_is_refused(self):
        path = self.report("name\tnew_est_reads\nEscherichia coli\tlots\n")
        with self.assertRaises(ValueError) as caught:
            plot.read_bracken_report(path)
        self.assertIn("holds no numbers", str(caught.exception))

    def test_missing_name_is_refused_rather_than_taking_the_first_column(self):
        path = self.report("taxonomy_id\tnew_est_reads\n562\t1840\n")
        with self.assertRaises(ValueError):
            plot.read_bracken_report(path)


class GenuinelyEmptyIsStillEmpty(unittest.TestCase):
    """The workflow writes zero-byte reports for samples with no hits, so those must
    stay a quiet 'no taxa' rather than becoming an error."""

    def report(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        handle.write(text)
        handle.close()
        path = Path(handle.name)
        self.addCleanup(lambda: path.unlink(missing_ok=True))
        return path

    def test_a_zero_byte_report_reads_as_no_taxa(self):
        self.assertTrue(plot.read_bracken_report(self.report("")).empty)

    def test_a_header_only_report_reads_as_no_taxa(self):
        path = self.report("name\ttaxonomy_id\tnew_est_reads\n")
        self.assertTrue(plot.read_bracken_report(path).empty)

    def test_a_missing_file_reads_as_no_taxa(self):
        self.assertTrue(plot.read_bracken_report(Path("/no/such/report.tsv")).empty)

    def test_zero_read_rows_are_dropped(self):
        path = self.report("name\tnew_est_reads\nEscherichia coli\t0\nBacteroides\t5\n")
        self.assertEqual(list(plot.read_bracken_report(path)["taxon"]), ["Bacteroides"])


if __name__ == "__main__":
    unittest.main()
