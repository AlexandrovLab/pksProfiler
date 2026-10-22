"""The host-depletion references are described, and the describer is correct.

rc2 review item 2.7 (Ludmil's 27 Aug ask, reopened): "Document host-depletion reference
contents. Record exactly which human assemblies/decoys/phiX are present in
human-GRC-db.mmi / human-GCA-phix-db.mmi and the build command used."

Everything the pipeline reports is what those two indexes did not remove, and until now
the only recorded fact about either was a path and a size.

Two halves here:

  * scripts/describe_minimap2_index.py is what the documentation was read from, so it is
    tested against synthetic .mmi files this file builds byte by byte -- a correct one,
    the flag combinations, and truncated ones. If the parser is wrong, every number in
    docs/host_depletion_references.md is wrong, so this is the part that matters.
  * the documentation and conf/databases.tsv must actually carry the findings, including
    the uncomfortable one (the index named ...phix... contains no phiX) and the honest
    one (the original build command could not be recovered). A doc that quietly dropped
    either would close the item without answering it.

No real index is opened here: the tests run anywhere, with no databases present.
"""
import importlib.util
import io
import struct
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/describe_minimap2_index.py"
DOC = ROOT / "docs/host_depletion_references.md"
MANIFEST = ROOT / "conf/databases.tsv"

spec = importlib.util.spec_from_file_location("describe_minimap2_index", SCRIPT)
describe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(describe)


def build_index(sequences, k=21, w=11, b=14, flag=0, trailing=b"\x00" * 64):
    """A synthetic .mmi: real header and dictionary, junk where the buckets would be."""
    out = bytearray(b"MMI\x02")
    out += struct.pack("<5I", w, k, b, len(sequences), flag)
    for name, length in sequences:
        encoded = name.encode()
        out += bytes([len(encoded)]) + encoded + struct.pack("<I", length)
    return bytes(out) + trailing


class TheParserReadsWhatMinimap2Wrote(unittest.TestCase):
    """mm_idx_dump(): magic, then w/k/b/n_seq/flag, then (len, name, length) per sequence."""

    def parse(self, blob):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.mmi"
            path.write_bytes(blob)
            with path.open("rb") as handle:
                header = describe.read_header(handle)
                sequences = describe.read_sequences(handle, header["n_seq"])
            return header, sequences

    def test_header_fields_are_read_in_minimap2_s_order(self):
        # w comes before k in the file; swapping them would silently mislabel the preset.
        header, _ = self.parse(build_index([("chr1", 1000)], k=19, w=19, b=12, flag=3))
        self.assertEqual((header["k"], header["w"], header["b"], header["flag"]),
                         (19, 19, 12, 3))
        self.assertEqual(header["n_seq"], 1)

    def test_names_and_lengths_round_trip(self):
        pairs = [("NC_000001.11", 248956422), ("chrM", 16569), ("phiX174", 5386)]
        _, sequences = self.parse(build_index(pairs))
        self.assertEqual(sequences, pairs)

    def test_a_long_dictionary_round_trips(self):
        pairs = [(f"NW_{i:09d}.1", 1000 + i) for i in range(500)]
        header, sequences = self.parse(build_index(pairs))
        self.assertEqual(header["n_seq"], 500)
        self.assertEqual(sequences, pairs)

    def test_a_nameless_index_is_read_as_empty_names_not_as_a_crash(self):
        _, sequences = self.parse(build_index([("", 1000), ("", 2000)], flag=4))
        self.assertEqual(sequences, [("", 1000), ("", 2000)])

    def test_a_file_that_is_not_an_index_is_refused(self):
        with self.assertRaises(SystemExit):
            self.parse(b"BAM\x01" + b"\x00" * 100)

    def test_a_truncated_dictionary_is_refused_rather_than_half_read(self):
        blob = build_index([("chr1", 1000), ("chr2", 2000)], trailing=b"")
        with self.assertRaises(SystemExit):
            self.parse(blob[:-3])

    def test_a_dictionary_claiming_more_sequences_than_it_holds_is_refused(self):
        body = build_index([("chr1", 1000)], trailing=b"")
        lied = body[:16] + struct.pack("<I", 9) + body[20:]
        with self.assertRaises(SystemExit):
            self.parse(lied)


class TheReportedFactsAreTheOnesTheFileCarries(unittest.TestCase):
    def report(self, blob):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.mmi"
            path.write_bytes(blob)
            buffer = io.StringIO()
            argv = sys.argv
            sys.argv = ["describe_minimap2_index.py", str(path)]
            try:
                with redirect_stdout(buffer):
                    describe.main()
            finally:
                sys.argv = argv
            return buffer.getvalue()

    def test_sr_settings_are_named_as_sr(self):
        text = self.report(build_index([("chr1", 10)], k=21, w=11))
        self.assertIn("sr (short reads)", text)

    def test_unknown_settings_are_not_guessed_at(self):
        text = self.report(build_index([("chr1", 10)], k=7, w=3))
        self.assertIn("no preset known", text)

    def test_the_hpc_flag_is_named(self):
        text = self.report(build_index([("chr1", 10)], flag=1))
        self.assertIn("MM_I_HPC", text)

    def test_counts_and_totals_are_the_dictionary_s(self):
        text = self.report(build_index([("a", 100), ("b", 250)]))
        self.assertIn("sequences          2", text)
        self.assertIn("total bases        350", text)

    def test_the_dictionary_digest_covers_the_dictionary_and_not_the_buckets(self):
        """Two indexes with the same dictionary and different buckets digest alike."""
        pairs = [("chr1", 100)]
        one = self.report(build_index(pairs, trailing=b"\x01" * 64))
        two = self.report(build_index(pairs, trailing=b"\x02" * 999))
        digest = lambda text: [l for l in text.splitlines() if "dictionary sha256" in l][0]
        self.assertEqual(digest(one), digest(two))
        # ... and a different dictionary does not.
        other = self.report(build_index([("chr1", 101)], trailing=b"\x01" * 64))
        self.assertNotEqual(digest(one), digest(other))

    def test_it_runs_as_a_program_and_writes_the_full_dictionary(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.mmi"
            path.write_bytes(build_index([("NC_000001.11", 100), ("phiX174", 5386)]))
            out = Path(tmp) / "d.tsv"
            result = subprocess.run([sys.executable, "-I", str(SCRIPT), str(path),
                                     "--tsv", str(out)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = out.read_text().splitlines()
        self.assertEqual(rows[0], "name\tlength\tcategory")
        self.assertEqual(rows[1].split("\t")[:2], ["NC_000001.11", "100"])
        self.assertEqual(rows[2].split("\t")[:2], ["phiX174", "5386"])


class TheDocumentationAnswersTheItem(unittest.TestCase):
    TEXT = DOC.read_text() if DOC.exists() else ""

    def test_the_document_exists_and_is_linked_from_the_readme(self):
        self.assertTrue(DOC.exists(), "docs/host_depletion_references.md is missing")
        self.assertIn("docs/host_depletion_references.md", (ROOT / "README.md").read_text())

    def test_it_names_both_indexes_and_both_parameters(self):
        for token in ("human-GRC-db.mmi", "human-GCA-phix-db.mmi",
                      "--hg38_db", "--t2t_phix_db"):
            with self.subTest(token=token):
                self.assertIn(token, self.TEXT)

    def test_it_records_the_assemblies_by_accession(self):
        # The point of the item: which assembly, not "a human reference".
        self.assertIn("GCF_000001405.40", self.TEXT)   # GRCh38.p14, RefSeq
        self.assertIn("GCA_009914755.4", self.TEXT)    # T2T-CHM13v2.0
        self.assertIn("3298430636", self.TEXT.replace(",", ""))
        self.assertIn("3117292070", self.TEXT.replace(",", ""))

    def test_it_answers_the_decoy_and_phix_question_explicitly(self):
        for token in ("decoy", "EBV", "HLA", "phiX"):
            with self.subTest(token=token):
                self.assertIn(token, self.TEXT)

    def test_it_states_that_the_phix_index_contains_no_phix(self):
        # Found while answering 2.7, and the reason the item could not just be closed
        # with "T2T + phiX": the file name claims something the index does not carry.
        self.assertIn("no phiX", self.TEXT)
        self.assertIn("NC_001422", self.TEXT)

    def test_it_says_plainly_that_the_build_command_was_not_recovered(self):
        self.assertIn("could not be recovered", self.TEXT)

    def test_it_records_what_is_verifiable_instead(self):
        for token in ("minimap2 -x sr -d", "k=21", "w=11", "scripts/describe_minimap2_index.py"):
            with self.subTest(token=token):
                self.assertIn(token, self.TEXT)


class TheManifestCarriesTheSameFacts(unittest.TestCase):
    ROWS = {row.split("\t")[0]: row.split("\t")
            for row in MANIFEST.read_text().splitlines() if not row.startswith("#")}

    def test_both_host_depletion_indexes_are_rows(self):
        self.assertIn("hg38_db", self.ROWS)
        self.assertIn("t2t_phix_db", self.ROWS)

    def test_the_release_column_names_the_assembly_not_just_a_build_date(self):
        self.assertIn("GCF_000001405.40", self.ROWS["hg38_db"][2])
        self.assertIn("GCA_009914755.4", self.ROWS["t2t_phix_db"][2])

    def test_the_description_points_at_the_document(self):
        for name in ("hg38_db", "t2t_phix_db"):
            with self.subTest(database=name):
                self.assertIn("docs/host_depletion_references.md", self.ROWS[name][1])

    def test_the_phix_row_carries_the_warning_where_a_reader_will_meet_it(self):
        self.assertIn("NO phiX", self.ROWS["t2t_phix_db"][1])

    def test_a_key_file_digest_replaced_the_placeholder(self):
        # rc3 froze both rows with key_file "-": there was no small release-defining file
        # to checksum. The index's own header+dictionary is exactly that file.
        for name in ("hg38_db", "t2t_phix_db"):
            row = self.ROWS[name]
            with self.subTest(database=name):
                self.assertNotEqual(row[8], "-")
                self.assertRegex(row[9], r"^[0-9a-f]{64}$")

    def test_the_source_column_says_the_command_is_a_reconstruction(self):
        for name in ("hg38_db", "t2t_phix_db"):
            with self.subTest(database=name):
                self.assertIn("command not recorded", self.ROWS[name][11])
                self.assertIn("reconstruction", self.ROWS[name][11])


if __name__ == "__main__":
    unittest.main()
