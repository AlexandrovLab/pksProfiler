"""V8b carrier attribution: which pks-positive organism carries the island?

Ludmil, 25 Sep report, section 4: "V8b. Carrier attribution. Separately compare genuine
pks reads against alternative pks-positive carriers/strains (pks+ E. coli, Klebsiella,
Citrobacter, Enterobacter, etc.). This asks which organism is the likely carrier and
should not be conflated with locus-detection specificity."

These cases are built from hand-written SAM records with explicit AS tags, so they
exercise the attribution LOGIC, not minimap2. The panel is pks-positive throughout, so
"this read fits Klebsiella better than IHE3034" is a statement about the organism; there
is no true/decoy axis here at all, and tests/test_v8_lane_separation.py holds the two
lanes' vocabularies apart.
"""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("carrier_attribution",
                                               ROOT / "scripts/carrier_attribution.py")
ca = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ca)

CARRIERS = {
    "ecoli_ihe3034__NC_017628.1": {"carrier": "ecoli_ihe3034",
                                    "organism": "Escherichia coli IHE3034",
                                    "accession": "GCF_000025745.1", "genus": "Escherichia"},
    "ecoli_nissle_1917__NZ_CP007799.1": {"carrier": "ecoli_nissle_1917",
                                          "organism": "Escherichia coli Nissle 1917",
                                          "accession": "GCF_000714595.1", "genus": "Escherichia"},
    "klebsiella_pneumoniae_sgh10__NZ_CP025080.1": {"carrier": "klebsiella_pneumoniae_sgh10",
                                                    "organism": "Klebsiella pneumoniae SGH10",
                                                    "accession": "GCF_002813595.1", "genus": "Klebsiella"},
    "citrobacter_koseri_atccbaa895__NC_009792.1": {"carrier": "citrobacter_koseri_atccbaa895",
                                                    "organism": "Citrobacter koseri ATCC BAA-895",
                                                    "accession": "GCF_000018045.1", "genus": "Citrobacter"},
}


class Args:
    min_delta_as = 1.0
    tie_margin = 0.0
    min_attributed_reads = 5
    min_top_carrier_fraction = 0.8


def sam_record(read_id, contig, as_score, flag=0):
    return (read_id, flag, contig, as_score)


def records_for(reads):
    """reads: {read_id: {contig: AS}} -> the flat record stream the script parses."""
    out = []
    for read_id, hits in reads.items():
        if not hits:
            out.append((read_id, ca.FLAG_UNMAPPED, "*", None))
            continue
        for i, (contig, score) in enumerate(hits.items()):
            out.append(sam_record(read_id, contig, score, flag=0 if i == 0 else 0x100))
    return out


def attribute(reads, args=None):
    args = args or Args()
    return ca.attribute_reads(records_for(reads), CARRIERS, args.min_delta_as, args.tie_margin)


class PerReadAttribution(unittest.TestCase):
    def test_a_clear_winner_is_attributed_to_that_strain(self):
        rows = attribute({"r1": {"ecoli_ihe3034__NC_017628.1": 100,
                                  "klebsiella_pneumoniae_sgh10__NZ_CP025080.1": 80}})
        self.assertEqual(rows[0]["attribution"], ca.ATTRIBUTED)
        self.assertEqual(rows[0]["best_carrier"], "ecoli_ihe3034")
        self.assertEqual(rows[0]["delta_as"], 20)
        self.assertEqual(rows[0]["best_genus"], "Escherichia")

    def test_a_klebsiella_read_is_attributed_to_klebsiella_not_to_the_pipeline_reference(self):
        # The acceptance case: the panel's job is to let a non-IHE3034 carrier win.
        rows = attribute({"r1": {"ecoli_ihe3034__NC_017628.1": 70,
                                  "klebsiella_pneumoniae_sgh10__NZ_CP025080.1": 96}})
        self.assertEqual(rows[0]["attribution"], ca.ATTRIBUTED)
        self.assertEqual(rows[0]["best_carrier"], "klebsiella_pneumoniae_sgh10")
        self.assertEqual(rows[0]["best_genus"], "Klebsiella")
        self.assertEqual(rows[0]["runner_up_carrier"], "ecoli_ihe3034")

    def test_an_exact_tie_is_ambiguous_never_force_assigned(self):
        rows = attribute({"r1": {"ecoli_ihe3034__NC_017628.1": 90,
                                  "klebsiella_pneumoniae_sgh10__NZ_CP025080.1": 90}})
        self.assertEqual(rows[0]["attribution"], ca.AMBIGUOUS)
        self.assertEqual(rows[0]["carriers_within_margin"], 2)

    def test_a_margin_below_min_delta_as_is_ambiguous(self):
        # One AS point apart with --min-delta-as 2: too close to name a strain.
        class Strict(Args):
            min_delta_as = 2.0
        rows = attribute({"r1": {"ecoli_ihe3034__NC_017628.1": 91,
                                  "ecoli_nissle_1917__NZ_CP007799.1": 90}}, Strict())
        self.assertEqual(rows[0]["attribution"], ca.AMBIGUOUS)

    def test_a_read_matching_no_panel_contig_is_unmapped_not_attributed(self):
        rows = attribute({"r1": {}})
        self.assertEqual(rows[0]["attribution"], ca.UNMAPPED)
        self.assertEqual(rows[0]["best_carrier"], "")

    def test_alignments_to_contigs_outside_the_map_are_ignored_not_guessed(self):
        records = [("r1", 0, "some_other_contig", 200),
                   ("r1", 0x100, "ecoli_ihe3034__NC_017628.1", 90)]
        rows = ca.attribute_reads(records, CARRIERS, 1.0, 0.0)
        self.assertEqual(rows[0]["best_carrier"], "ecoli_ihe3034")


class SampleLevelCall(unittest.TestCase):
    def summarize(self, reads, args=None):
        args = args or Args()
        rows = attribute(reads, args)
        return ca.summarize("S1", rows, CARRIERS, args)

    def test_a_dominant_carrier_is_attributed(self):
        reads = {f"r{i}": {"klebsiella_pneumoniae_sgh10__NZ_CP025080.1": 100,
                            "ecoli_ihe3034__NC_017628.1": 60} for i in range(8)}
        summary = self.summarize(reads)
        self.assertEqual(summary["carrier_attribution"], ca.CARRIER_ATTRIBUTED)
        self.assertEqual(summary["top_carrier"], "klebsiella_pneumoniae_sgh10")
        self.assertEqual(summary["top_genus"], "Klebsiella")
        self.assertTrue(summary["strain_level_resolved"])

    def test_two_strains_of_one_genus_resolve_the_genus_but_not_the_strain(self):
        # Every read ties between two E. coli strains: the genus is supported, the strain
        # is not, and the call says exactly that rather than picking one.
        reads = {f"r{i}": {"ecoli_ihe3034__NC_017628.1": 100,
                            "ecoli_nissle_1917__NZ_CP007799.1": 100} for i in range(8)}
        summary = self.summarize(reads)
        self.assertEqual(summary["carrier_attribution"], ca.CARRIER_GENUS_ONLY)
        self.assertFalse(summary["strain_level_resolved"])
        self.assertTrue(summary["genus_level_resolved"])
        self.assertEqual(summary["top_genus"], "Escherichia")

    def test_reads_split_across_genera_are_ambiguous(self):
        reads = {}
        for i in range(4):
            reads[f"e{i}"] = {"ecoli_ihe3034__NC_017628.1": 100,
                              "klebsiella_pneumoniae_sgh10__NZ_CP025080.1": 60}
        for i in range(4):
            reads[f"k{i}"] = {"klebsiella_pneumoniae_sgh10__NZ_CP025080.1": 100,
                              "ecoli_ihe3034__NC_017628.1": 60}
        summary = self.summarize(reads)
        self.assertEqual(summary["carrier_attribution"], ca.CARRIER_AMBIGUOUS)
        self.assertLess(float(summary["top_carrier_fraction"]), 0.8)

    def test_too_few_reads_is_undetermined_whatever_the_point_estimate(self):
        reads = {"r1": {"klebsiella_pneumoniae_sgh10__NZ_CP025080.1": 100,
                         "ecoli_ihe3034__NC_017628.1": 50}}
        summary = self.summarize(reads)
        self.assertEqual(summary["carrier_attribution"], ca.CARRIER_UNDETERMINED)

    def test_the_panel_size_is_reported_so_a_call_names_what_it_competed_against(self):
        reads = {f"r{i}": {"ecoli_ihe3034__NC_017628.1": 100} for i in range(8)}
        self.assertEqual(self.summarize(reads)["panel_carriers"], 4)

    def test_per_carrier_counts_are_kept_not_collapsed_to_the_winner(self):
        reads = {f"r{i}": {"ecoli_ihe3034__NC_017628.1": 100,
                            "klebsiella_pneumoniae_sgh10__NZ_CP025080.1": 60} for i in range(6)}
        reads["x1"] = {"citrobacter_koseri_atccbaa895__NC_009792.1": 100,
                       "ecoli_ihe3034__NC_017628.1": 60}
        summary = self.summarize(reads)
        self.assertIn("ecoli_ihe3034:6", summary["carrier_read_counts"])
        self.assertIn("citrobacter_koseri_atccbaa895:1", summary["carrier_read_counts"])


class VocabularyIsV8bOnly(unittest.TestCase):
    """Ludmil: carrier attribution "should not be conflated with locus-detection
    specificity". The script must not emit any of V8a's words."""

    SOURCE = (ROOT / "scripts/carrier_attribution.py").read_text()

    def test_no_v8a_verdict_name_appears_in_the_output_vocabulary(self):
        for banned in ("confirmed_true", "confirmed_decoy", "SURVIVES",
                        "NOT_DISTINGUISHABLE", "UNDERPOWERED"):
            self.assertNotIn(f'"{banned}"', self.SOURCE)

    def test_the_summary_fields_are_all_carrier_named(self):
        for field in ("carrier_attribution", "top_carrier", "top_genus"):
            self.assertIn(field, ca.SUMMARY_FIELDS)
        self.assertNotIn("verdict", ca.SUMMARY_FIELDS)
        self.assertNotIn("retention", ca.SUMMARY_FIELDS)


class CommandLine(unittest.TestCase):
    def test_it_writes_both_tables_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "map.tsv").write_text(
                "contig\tcarrier_key\tcarrier_organism\tcarrier_accession\n"
                + "".join(f"{c}\t{e['carrier']}\t{e['organism']}\t{e['accession']}\n"
                          for c, e in CARRIERS.items()))
            sam = tmp / "s.sam"
            lines = ["@SQ\tSN:ecoli_ihe3034__NC_017628.1\tLN:5108383"]
            for i in range(8):
                lines.append(f"r{i}\t0\tklebsiella_pneumoniae_sgh10__NZ_CP025080.1\t1\t60\t100M"
                             f"\t*\t0\t0\t*\t*\tAS:i:100\tNM:i:0")
                lines.append(f"r{i}\t256\tecoli_ihe3034__NC_017628.1\t1\t0\t100M"
                             f"\t*\t0\t0\t*\t*\tAS:i:60\tNM:i:8")
            sam.write_text("\n".join(lines) + "\n")
            ca.main(["--sample", "S1", "--sam", str(sam),
                     "--carrier-map", str(tmp / "map.tsv"),
                     "--min-attributed-reads", "5",
                     "--output-summary", str(tmp / "sum.tsv"),
                     "--output-per-read", str(tmp / "reads.tsv")])
            with (tmp / "sum.tsv").open() as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["carrier_attribution"], ca.CARRIER_ATTRIBUTED)
            self.assertEqual(summary["top_carrier"], "klebsiella_pneumoniae_sgh10")
            with (tmp / "reads.tsv").open() as handle:
                per_read = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(per_read), 8)


if __name__ == "__main__":
    unittest.main()
