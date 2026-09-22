"""The taxDB column order is checked, not assumed.

F04 reported that this parser swaps rank and scientific name. It does not, for the
database in use: across all 2,518,002 rows of krakenUniq_8_8_2023's taxDB, field 4
holds a recognised rank in 97.9% and field 3 in none. Applying the fix as written
would have created the bug it describes -- every organism named "species" or
"no rank", and any rank filter matching nothing.

The legitimate half of the finding is that nothing checked. The parser assumed the
order, so a build that differed would have been mis-parsed in silence. It now asserts
the order at load and says which way round the file looks.
"""
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "build_clb_species_matrix", ROOT / "scripts/build_clb_species_matrix.py")
matrix = importlib.util.module_from_spec(spec)
spec.loader.exec_module(matrix)


class ColumnOrderIsAsserted(unittest.TestCase):
    CORRECT = [["562", "561", "Escherichia coli", "species"],
               ["2", "131567", "Bacteria", "superkingdom"],
               ["561", "543", "Escherichia", "genus"]]

    def test_the_real_order_is_accepted(self):
        matrix.assert_taxdb_column_order(self.CORRECT, "taxDB")

    def test_a_file_with_rank_in_field_three_is_refused(self):
        swapped = [[row[0], row[1], row[3], row[2]] for row in self.CORRECT]
        with self.assertRaises(ValueError) as caught:
            matrix.assert_taxdb_column_order(swapped, "taxDB")
        self.assertIn("field 3, not field 4", str(caught.exception))

    def test_a_file_that_is_neither_is_refused(self):
        junk = [["1", "1", "alpha", "beta"], ["2", "1", "gamma", "delta"]]
        with self.assertRaises(ValueError) as caught:
            matrix.assert_taxdb_column_order(junk, "taxDB")
        self.assertIn("does not look like a KrakenUniq taxDB", str(caught.exception))

    def test_an_unfamiliar_rank_does_not_fail_the_file(self):
        # The vocabulary is for telling ranks from names, not for validating ranks.
        rows = self.CORRECT + [["99", "1", "Some organism", "pathovar"]]
        matrix.assert_taxdb_column_order(rows, "taxDB")

    def test_short_rows_are_ignored_rather_than_crashing(self):
        matrix.assert_taxdb_column_order(self.CORRECT + [["1"], []], "taxDB")

    def test_the_parser_still_reads_name_from_three_and_rank_from_four(self):
        # The point of the assertion is to keep this true, not to change it.
        source = (ROOT / "scripts/build_clb_species_matrix.py").read_text()
        self.assertIn("scientific_name = fields[2].strip()", source)
        self.assertIn("rank = fields[3].strip()", source)


if __name__ == "__main__":
    unittest.main()


class TheDatabasesOwnTaxonomyWins(unittest.TestCase):
    """F05: taxDB is the tree that explains this database's labels.

    Ludmil's reasoning, verbatim: "When taxonomy/nodes.dmp and names.dmp exist,
    load_taxonomy uses them without incorporating additional taxDB nodes. KrakenUniq
    can add assembly/sequence pseudo-taxids specific to a database build. A fixture
    containing a valid pseudo-node in taxDB and its species in NCBI files yielded no
    species ancestor because the pseudo-node was never loaded."

    So the fixtures below differ on exactly that: taxDB carries a pseudo-taxid that
    the NCBI dump does not. Which tree was loaded is observable from whether that node
    is present, rather than from which code path happened to run.
    """

    # 9000001 is the kind of id KrakenUniq invents for an assembly; NCBI has no such node
    TAXDB = ("562\t561\tEscherichia coli\tspecies\n"
             "561\t543\tEscherichia\tgenus\n"
             "9000001\t562\tEscherichia coli assembly GCF_000005845\tsequence\n")
    NODES = "562\t|\t561\t|\tspecies\t|\n561\t|\t543\t|\tgenus\t|\n"
    NAMES = ("562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n"
             "561\t|\tEscherichia\t|\t\t|\tscientific name\t|\n")

    def build(self, taxdb=True, ncbi=True):
        import tempfile
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        db = Path(directory.name)
        if taxdb:
            (db / "taxDB").write_text(self.TAXDB)
        if ncbi:
            (db / "taxonomy").mkdir()
            (db / "taxonomy/nodes.dmp").write_text(self.NODES)
            (db / "taxonomy/names.dmp").write_text(self.NAMES)
        return db

    def test_with_both_present_the_pseudo_taxid_survives(self):
        parents, _ranks, _names = matrix.load_taxonomy(self.build())
        self.assertIn("9000001", parents,
                      "the NCBI dump won, so the database's own taxids are missing")

    def test_the_ignored_dump_is_announced(self):
        import contextlib, io
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            matrix.load_taxonomy(self.build())
        self.assertIn("is present and ignored", stderr.getvalue())

    def test_falling_back_to_ncbi_warns_what_it_costs(self):
        import contextlib, io
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            parents, _ranks, _names = matrix.load_taxonomy(self.build(taxdb=False))
        self.assertIn("will not resolve to a species", stderr.getvalue())
        self.assertNotIn("9000001", parents)      # the point of the warning

    def test_taxdb_alone_is_silent(self):
        import contextlib, io
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            matrix.load_taxonomy(self.build(ncbi=False))
        self.assertEqual(stderr.getvalue(), "")

    def test_neither_is_still_an_error(self):
        with self.assertRaises(FileNotFoundError):
            matrix.load_taxonomy(self.build(taxdb=False, ncbi=False))
