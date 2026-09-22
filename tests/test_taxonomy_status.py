"""A sample with no species says so, in the table that lists every sample.

F10, from Ludmil's 2026-09-19 report: "The workflow writes zero-byte MPA files for
no-hit/low-support cases... A separate guard sums rank-level reads and runs Bracken
when their total is at least two, although `-t 2` is applied per taxon. Two species
with one read each pass the guard but leave no eligible species. Header-only
species-support inputs produce no rows in the long table, which is acceptable only
with an explicit companion sample-status table."

The companion table already exists -- pks.qc.summary.tsv has a row per sample and a
status column -- so the taxonomy outcome goes there rather than into a new file, and
the species table keeps only samples with species, which is what it is for.
"""
import csv
import importlib.util
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAXA = (ROOT / "Modules/pks_taxa.nf").read_text()
MAIN = (ROOT / "main.nf").read_text()
spec = importlib.util.spec_from_file_location(
    "build_qc_summary", ROOT / "scripts/build_qc_summary.py")
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)

# the awk the module ships, with Nextflow's escaping resolved
GUARD_AWK = (r'$8 == "species" && $2 ~ /^[0-9]+$/ && $2+0 > best { best = $2+0 } '
             r'END { print best+0 }')


def kraken_report(rows):
    """(reads, rank) per line, in the column positions the guard reads."""
    return "".join(f"pct\t{reads}\tx\tx\tx\tx\tx\t{rank}\n" for reads, rank in rows)


class TheGuardAsksBrackensQuestion(unittest.TestCase):
    def best(self, text):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as handle:
            handle.write(text)
            path = handle.name
        try:
            out = subprocess.run(["awk", "-F", "\t", GUARD_AWK, path],
                                 capture_output=True, text=True, check=True)
            return int(out.stdout.strip())
        finally:
            Path(path).unlink()

    def test_two_species_with_one_read_each_no_longer_qualify(self):
        # The case from the report: summing gave 2 and ran Bracken, which then
        # reported neither species because -t 2 is per taxon.
        self.assertEqual(self.best(kraken_report([(1, "species"), (1, "species")])), 1)

    def test_one_species_with_two_reads_qualifies(self):
        self.assertEqual(self.best(kraken_report([(2, "species")])), 2)

    def test_the_largest_taxon_decides_not_the_total(self):
        self.assertEqual(self.best(kraken_report([(1, "species"), (7, "species"),
                                                  (1, "species")])), 7)

    def test_other_ranks_are_ignored(self):
        self.assertEqual(self.best(kraken_report([(9, "genus"), (1, "species")])), 1)

    def test_no_species_rows_is_zero(self):
        self.assertEqual(self.best(kraken_report([(9, "genus")])), 0)

    def test_the_module_no_longer_sums_the_rank(self):
        guard = TAXA.split("GENUS_READS=", 1)[1][:400]
        self.assertIn("best", guard)
        self.assertNotIn("sum +=", guard)


class TheOutcomeReachesTheQcTable(unittest.TestCase):
    def summarise(self, metrics, sample="S1"):
        with tempfile.TemporaryDirectory() as directory:
            fragment = Path(directory) / f"{sample}.qc.tsv"
            lines = ["Sample\tMetric\tValue"]
            lines += [f"{sample}\t{metric}\t{value}" for metric, value in metrics]
            fragment.write_text("\n".join(lines) + "\n")
            output = Path(directory) / "summary.tsv"
            qc.write_summary(qc.load_fragments([fragment]), output, [sample])
            with output.open(newline="") as handle:
                return list(csv.DictReader(handle, delimiter="\t"))[0]

    BASE = [("filter_input_reads", 400), ("reads_after_fastp", 360),
            ("reads_after_hg38", 300), ("reads_after_t2t_phix", 250)]

    def test_a_sample_with_no_island_reads_says_so(self):
        row = self.summarise(self.BASE + [("taxonomy_status", "no_pks_reads"),
                                          ("clb_species_reported", 0)])
        self.assertEqual(row["taxonomy_status"], "no_pks_reads")
        self.assertEqual(row["clb_species_reported"], "0")

    def test_a_sample_below_the_threshold_is_distinguishable(self):
        row = self.summarise(self.BASE + [("taxonomy_status", "below_rank_threshold"),
                                          ("clb_species_reported", 0)])
        self.assertEqual(row["taxonomy_status"], "below_rank_threshold")

    def test_a_sample_with_species_reports_how_many(self):
        row = self.summarise(self.BASE + [("taxonomy_status", "species_identified"),
                                          ("clb_species_reported", 3)])
        self.assertEqual(row["clb_species_reported"], "3")

    def test_taxonomy_not_run_is_na_rather_than_a_verdict(self):
        row = self.summarise(self.BASE)
        self.assertEqual(row["taxonomy_status"], "NA")
        self.assertEqual(row["clb_species_reported"], "NA")


class Wiring(unittest.TestCase):
    def test_the_taxonomy_lane_emits_a_qc_fragment(self):
        self.assertIn("taxonomy.qc.tsv", TAXA)
        self.assertIn("emit: qc", TAXA)
        for metric in ("taxonomy_status", "clb_species_reported"):
            self.assertIn(metric, TAXA)

    def test_the_fragment_reaches_the_summary(self):
        self.assertIn("Bracken.out.qc", MAIN)

    def test_every_status_the_module_writes_is_one_of_the_documented_four(self):
        written = set(re.findall(r'STATUS="(\w+)"', TAXA))
        written |= set(re.findall(r'taxonomy_status\\\\t(\w+)', TAXA))
        self.assertTrue(written)
        self.assertTrue(written <= {"species_identified", "below_rank_threshold",
                                    "no_species_identified", "no_pks_reads"}, written)


if __name__ == "__main__":
    unittest.main()


class TextMetricsAreValidatedNotWavedThrough(unittest.TestCase):
    """The QC loader required every value to be a number. A status is not a number,
    but it should not become a free-text field either."""

    def load(self, metric, value):
        with tempfile.TemporaryDirectory() as directory:
            fragment = Path(directory) / "S1.qc.tsv"
            fragment.write_text(f"Sample\tMetric\tValue\nS1\t{metric}\t{value}\n")
            return qc.load_fragments([fragment])

    def test_a_known_status_is_accepted(self):
        loaded = self.load("taxonomy_status", "no_pks_reads")
        self.assertEqual(loaded["S1"]["taxonomy_status"], "no_pks_reads")

    def test_an_unknown_status_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            self.load("taxonomy_status", "probably_fine")
        self.assertIn("Unknown taxonomy_status", str(caught.exception))

    def test_a_count_metric_must_still_be_a_number(self):
        with self.assertRaises(ValueError):
            self.load("reads_after_fastp", "lots")

    def test_the_status_stays_text_rather_than_becoming_an_int(self):
        loaded = self.load("taxonomy_status", "species_identified")
        self.assertIsInstance(loaded["S1"]["taxonomy_status"], str)
