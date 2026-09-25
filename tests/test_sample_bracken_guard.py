"""sampleBracken must not crash on a sparse community-taxonomy sample.

Ludmil, revised report finding 4: sampleBracken called `bracken -t 2`
unconditionally. Bracken itself fails when no taxon at the requested rank
reaches its read threshold, so a thin sample turned into a task failure
instead of an explicit below-threshold/no-call state. pks_taxa.nf's own
Bracken process already guards this exact way -- the same pre-check, on the
same read-count definition Bracken itself uses (the largest single-taxon
count at that rank, not a sum across taxa, per F10) -- so sampleBracken now
mirrors it rather than inventing a second rule.

His acceptance test: 0, 1 and 2 reads assigned to a single genus/species --
the first two must complete as valid below-threshold outputs, the third must
invoke Bracken successfully.
"""
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREFILTER = (ROOT / "Modules/pks_prefilter.nf").read_text()
TAXA = (ROOT / "Modules/pks_taxa.nf").read_text()

# The awk exactly as both modules ship it (identical in each), Nextflow's
# escaping resolved. Column 8 is the rank name, column 2 the read count --
# the krakenuniq report format both processes consume.
GENUS_AWK = r"""$8 == "genus" && $2 ~ /^[0-9]+$/ && $2+0 > best {
        best = $2+0
      }
      END {
        print best+0
      }"""
SPECIES_AWK = r"""$8 == "species" && $2 ~ /^[0-9]+$/ && $2+0 > best {
        best = $2+0
      }
      END {
        print best+0
      }"""


def kraken_report(rows):
    """rows: (rank, reads) tuples. Real krakenuniq reports carry more columns;
    only 2 (reads) and 8 (rank) matter to this awk."""
    lines = []
    for rank, reads in rows:
        fields = ["0"] * 9
        fields[1] = str(reads)
        fields[7] = rank
        lines.append("\t".join(fields))
    return "\n".join(lines) + "\n"


def sample_bracken_script():
    """The sampleBracken process block. Splits on a top-level `process ` keyword
    (always at column 0 after a newline) rather than the bare word, which also
    appears in this file's own prose comments."""
    after = PREFILTER.split("\nprocess sampleBracken")[1]
    return after.split("\nprocess ")[0]


def run_awk(awk, text):
    with tempfile.NamedTemporaryFile("w", suffix=".report", delete=False) as handle:
        handle.write(text)
        path = handle.name
    try:
        result = subprocess.run(["awk", "-F", "\t", awk, path],
                                capture_output=True, text=True, check=True)
        return int(result.stdout.strip())
    finally:
        Path(path).unlink()


class BestSingleTaxonReadCount(unittest.TestCase):
    """F10's rule: the largest count on one taxon, not a sum across taxa."""

    def test_no_genus_rows_is_zero(self):
        self.assertEqual(run_awk(GENUS_AWK, kraken_report([("species", 5)])), 0)

    def test_one_read_does_not_reach_two(self):
        self.assertEqual(run_awk(GENUS_AWK, kraken_report([("genus", 1)])), 1)

    def test_exactly_two_reaches_the_threshold(self):
        self.assertEqual(run_awk(GENUS_AWK, kraken_report([("genus", 2)])), 2)

    def test_two_singletons_do_not_sum_to_the_threshold(self):
        # The F10 bug this awk exists to avoid: summing 1+1 across two genera
        # would wrongly clear -t 2, when neither genus alone does.
        text = kraken_report([("genus", 1), ("genus", 1)])
        self.assertEqual(run_awk(GENUS_AWK, text), 1)

    def test_species_rank_is_read_independently_of_genus(self):
        text = kraken_report([("genus", 10), ("species", 1)])
        self.assertEqual(run_awk(GENUS_AWK, text), 10)
        self.assertEqual(run_awk(SPECIES_AWK, text), 1)


class SampleBrackenMirrorsTheEstablishedGuard(unittest.TestCase):
    def test_the_unconditional_call_is_gone(self):
        # No bracken invocation outside the per-level if/continue guard.
        script = sample_bracken_script()
        self.assertIn('if [[ "\\$LVL_READS" -lt 2 ]]', script)
        self.assertIn("continue", script)
        self.assertIn("bracken -d", script)

    def test_below_threshold_still_writes_valid_empty_outputs(self):
        script = sample_bracken_script()
        self.assertIn(': > "\\$bracken_output"', script)
        self.assertIn(': > "\\$bracken_kraken_report"', script)

    def test_uses_the_same_awk_as_the_taxonomy_lane_bracken(self):
        # One rule for "does this rank clear Bracken's threshold," not two.
        prefilter_script = sample_bracken_script()
        self.assertIn('$8 == "genus"', prefilter_script)
        self.assertIn('$8 == "species"', prefilter_script)
        self.assertIn('$2+0 > best', prefilter_script)
        # And the taxonomy lane's own Bracken process still has its version.
        self.assertIn('$8 == "genus"', TAXA)
        self.assertIn('$8 == "species"', TAXA)

    def test_threshold_matches_bracken_s_own_dash_t(self):
        script = sample_bracken_script()
        self.assertIn("-t 2", script)
        self.assertIn("-lt 2", script)


if __name__ == "__main__":
    unittest.main()
