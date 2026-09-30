"""scripts/build_specificity_panels.py: the island screen that decides panel membership.

The V8a and V8b panels have opposite membership tests, and the reviewer's split only means
something if they are enforced from sequence rather than taken from the accession list:

  V8a decoy   must NOT carry the clb island -- a pks-carrying "decoy" absorbs genuine pks
              reads and turns locus specificity into a coin flip
  V8b carrier MUST carry it -- a pks-negative "carrier" can never win a read

Nothing here downloads anything. The screen is driven through real minimap2 PAF text and
real FASTA/GenBank fixtures written in the test, so the arithmetic (merged interval
coverage, not summed), the accept/reject direction and the MIBiG colibactin exclusion are
all exercised without network access.
"""
import importlib.util
import subprocess
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_specificity_panels",
                                               ROOT / "scripts/build_specificity_panels.py")
bp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bp)


def paf(intervals, query_len=1000):
    """minimap2 -c PAF: query name/len/start/end are columns 1-4."""
    return "\n".join(
        f"island\t{query_len}\t{start}\t{end}\t+\ttarget\t5000000\t0\t{end - start}"
        f"\t{end - start}\t{end - start}\t60"
        for start, end in intervals) + "\n"


class IslandCoverageIsAMergedUnionNotASum(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.island = Path(self.tmp.name) / "island.fna"
        self.island.write_text(">island\n" + "A" * 1000 + "\n")
        self.genome = Path(self.tmp.name) / "genome.fna"
        self.genome.write_text(">target\n" + "C" * 5000 + "\n")

    def tearDown(self):
        self.tmp.cleanup()

    def coverage(self, intervals):
        with mock.patch.object(bp.subprocess, "run",
                                return_value=subprocess.CompletedProcess(
                                    [], 0, stdout=paf(intervals), stderr="")):
            return bp.island_coverage("minimap2", self.island, self.genome, "asm20", 1)

    def test_one_alignment_covering_half_the_island(self):
        fraction, covered, total = self.coverage([(0, 500)])
        self.assertAlmostEqual(fraction, 0.5)
        self.assertEqual((covered, total), (500, 1000))

    def test_overlapping_alignments_are_merged_not_double_counted(self):
        # Summed these are 800 bp of a 1000 bp island; merged they are 600.
        fraction, covered, _total = self.coverage([(0, 400), (200, 600)])
        self.assertEqual(covered, 600)
        self.assertAlmostEqual(fraction, 0.6)

    def test_a_contained_alignment_adds_nothing(self):
        _fraction, covered, _total = self.coverage([(0, 800), (100, 200)])
        self.assertEqual(covered, 800)

    def test_no_alignment_is_zero_coverage_not_an_error(self):
        fraction, covered, _total = self.coverage([])
        self.assertEqual((fraction, covered), (0.0, 0))


class PanelMembershipIsScreenedNotAssumed(unittest.TestCase):
    """The accept/reject direction, read straight out of the source: the two thresholds
    exist, they point opposite ways, and a failure is recorded rather than dropped."""

    SOURCE = (ROOT / "scripts/build_specificity_panels.py").read_text()

    def test_a_decoy_is_accepted_only_below_its_ceiling(self):
        self.assertIn("ok = coverage < args.max_decoy_island_coverage", self.SOURCE)

    def test_a_carrier_is_accepted_only_above_its_floor(self):
        self.assertIn("ok = coverage >= args.min_carrier_island_coverage", self.SOURCE)

    def test_a_rejected_candidate_is_written_to_the_screen_table_before_the_skip(self):
        build = self.SOURCE.split("def build(", 1)[1]
        append_pos = build.index('screen_rows.append(')
        skip_pos = build.index("if not ok:\n            continue")
        self.assertLess(append_pos, skip_pos,
                         "a rejected candidate must still appear in the screen table")

    def test_organism_identity_is_confirmed_against_the_assemblys_own_report(self):
        self.assertIn("identity_ok = bool(genus_expected) and organism.startswith(genus_expected)",
                       self.SOURCE)
        self.assertIn("rejected_organism_mismatch", self.SOURCE)

    def test_an_empty_panel_is_an_error_not_an_empty_index(self):
        self.assertIn("every candidate failed the island screen; no panel was built",
                       self.SOURCE)


class ContigsAreTaggedWithTheirMember(unittest.TestCase):
    def test_every_contig_name_carries_its_member_key(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "in.fna"
            src.write_text(">NC_000913.3 E. coli K-12\nACGT\n>plasmid_1 pX\nTTTT\n")
            dest = Path(tmp) / "out.fna"
            n_contigs, n_bases, contigs = bp.write_tagged(src, dest, "ecoli_k12_mg1655")
            self.assertEqual(n_contigs, 2)
            self.assertEqual(n_bases, 8)
            self.assertEqual(contigs, ["ecoli_k12_mg1655__NC_000913.3",
                                        "ecoli_k12_mg1655__plasmid_1"])
            # The original header is kept as the description, so provenance survives.
            self.assertIn(">ecoli_k12_mg1655__NC_000913.3 NC_000913.3 E. coli K-12",
                           dest.read_text())


class MibigClustersExcludeColibactin(unittest.TestCase):
    """The one cluster that can never be a decoy for colibactin reads is colibactin's."""

    def gbk(self, name, text, sequence):
        lines = [f"LOCUS       {name}    {len(sequence)} bp", text, "ORIGIN"]
        for i in range(0, len(sequence), 60):
            lines.append(f"     {i + 1} " + sequence[i:i + 60].lower())
        lines.append("//")
        return "\n".join(lines) + "\n"

    def test_a_colibactin_cluster_is_excluded_and_an_nrps_cluster_is_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            members = {
                "BGC0000984.gbk": self.gbk("BGC0000984",
                                            "  /product=\"colibactin\"\n  /biosyn=\"NRPS\"",
                                            "ACGT" * 2000),
                "BGC0000001.gbk": self.gbk("BGC0000001",
                                            "  /product=\"bacillaene\"\n  /biosyn=\"transAT PKS\"",
                                            "GGCC" * 2000),
                "BGC0000002.gbk": self.gbk("BGC0000002",
                                            "  /product=\"a ribosomal peptide\"",
                                            "TTAA" * 2000),
                # Under the 5 kb floor: a fragment is not a cluster.
                "BGC0000003.gbk": self.gbk("BGC0000003", "  /biosyn=\"NRPS\"", "ACGT" * 100),
            }
            for name, text in members.items():
                (tmp / name).write_text(text)
            tarball = tmp / "mibig.tar.gz"
            with tarfile.open(tarball, "w:gz") as archive:
                for name in members:
                    archive.add(tmp / name, arcname=f"mibig_gbk_4.0/{name}")

            kept, seen, skipped = bp.mibig_clusters(tarball, tmp / "out", 100)
            names = [bgc for bgc, _path, _length in kept]
            self.assertEqual(seen, 4)
            self.assertEqual(skipped, 1)
            self.assertIn("BGC0000001", names)
            self.assertNotIn("BGC0000984", names)      # colibactin
            self.assertNotIn("BGC0000002", names)      # not NRPS/PKS
            self.assertNotIn("BGC0000003", names)      # too short

    def test_the_whole_archive_is_scanned_even_once_the_cap_is_reached(self):
        # Otherwise "0 colibactin clusters excluded" would mean "the scan stopped early",
        # and a colibactin cluster sitting late in the archive would go unreported.
        source = (ROOT / "scripts/build_specificity_panels.py").read_text()
        body = source.split("def mibig_clusters(", 1)[1].split("\ndef ", 1)[0]
        self.assertNotIn("break", body.split("for entry in archive:", 1)[1])


class PanelListsAreComplete(unittest.TestCase):
    """The accession lists are inputs to this script, and each row must say why it is
    there: a panel nobody can justify is a panel nobody can review."""

    def rows(self, name):
        text = (ROOT / "ref/specificity_panels" / name).read_text()
        lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
        header = lines[0].split("\t")
        return header, [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]

    def test_the_v8a_panel_rows_are_well_formed_and_justified(self):
        header, rows = self.rows("v8a_locus_decoy_panel.tsv")
        self.assertEqual(header, ["key", "accession", "category", "expected_organism",
                                   "rationale"])
        self.assertGreaterEqual(len(rows), 20)
        for row in rows:
            self.assertRegex(row["accession"], r"^GCF_\d+\.\d+$")
            self.assertTrue(row["rationale"].strip(), f"{row['key']} has no rationale")

    def test_the_v8b_panel_rows_are_well_formed_and_justified(self):
        header, rows = self.rows("v8b_carrier_panel.tsv")
        self.assertEqual(header, ["key", "accession", "genus_group", "expected_organism",
                                   "rationale"])
        for row in rows:
            self.assertRegex(row["accession"], r"^GCF_\d+\.\d+$")
            self.assertTrue(row["rationale"].strip(), f"{row['key']} has no rationale")

    def test_the_v8b_panel_carries_more_than_one_strain_per_named_genus(self):
        # "Klebsiella carries it" must not rest on one assembly.
        _header, rows = self.rows("v8b_carrier_panel.tsv")
        counts = {}
        for row in rows:
            counts[row["genus_group"]] = counts.get(row["genus_group"], 0) + 1
        for genus in ("escherichia", "klebsiella", "citrobacter"):
            self.assertGreaterEqual(counts.get(genus, 0), 2,
                                     f"{genus} has fewer than two candidate strains")

    def test_no_key_is_reused_across_the_two_panels_in_the_same_role(self):
        _h, decoys = self.rows("v8a_locus_decoy_panel.tsv")
        _h, carriers = self.rows("v8b_carrier_panel.tsv")
        shared = {r["key"] for r in decoys} & {r["key"] for r in carriers}
        # K. aerogenes KCTC 2190 is deliberately in both lists as a CANDIDATE: the screen
        # decides which panel it ends up in. Any other overlap would be an accident.
        self.assertLessEqual(shared, {"klebsiella_aerogenes_kctc2190"})


if __name__ == "__main__":
    unittest.main()
