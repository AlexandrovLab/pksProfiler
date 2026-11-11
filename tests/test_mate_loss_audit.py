"""Tier A 4+5: the mate-loss audit.

The audit replays the extraction predicate from Modules/extract_reads.nf over an
alignment and counts fragments that lose a mate. These tests pin the predicate to the
module's semantics (unmapped always; low-AS only when rescue is on; mapped-no-AS never)
and run the whole two-pass audit over a synthetic SAM covering every fragment fate, so
the numbers it reports on real cohorts rest on verified accounting rather than on the
script agreeing with itself.

The synthetic SAM is invented here; no sample data is involved.
"""
import csv
import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "a45_mate_loss_audit.py"

_spec = importlib.util.spec_from_file_location("a45_mate_loss_audit", SCRIPT)
audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(audit)

HAVE_SAMTOOLS = shutil.which("samtools") is not None

HEADER = "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:10000\n"


def rec(name, flag, mapq=0, cigar="*", rname="*", pos=0, score=None):
    fields = [name, str(flag), rname, str(pos), str(mapq), cigar,
              "*", "0", "0", "ACGTACGTAC", "IIIIIIIIII"]
    if score is not None:
        fields.append(f"AS:i:{score}")
    return "\t".join(fields) + "\n"


def mapped(name, flag, score, mapq=60):
    return rec(name, flag, mapq=mapq, cigar="10M", rname="chr1", pos=100, score=score)


# Every fragment fate the audit must classify:
#   F1  both mates unmapped                         -> both_mates_extracted
#   F2  R1 unmapped, R2 mapped AS 140               -> orphan, lost mate AS_100_plus
#   F3  R1 mapped AS 25, R2 mapped AS 35            -> orphan only when rescue is ON,
#                                                      lost mate in the 30-49 gray zone
#   F4  R1 unmapped, R2 mapped, NO AS tag, MAPQ 0   -> orphan, lost_mate_no_AS + mapq0
#   F5  unpaired unmapped read                      -> single_end_or_mateless
#   F6  both mates confidently mapped               -> untouched by extraction
SAM = HEADER + "".join([
    rec("F1", 77), rec("F1", 141),
    rec("F2", 69), mapped("F2", 137, 140),
    mapped("F3", 65, 25), mapped("F3", 129, 35),
    rec("F4", 69), mapped("F4", 137, None, mapq=0),
    rec("F5", 4),
    mapped("F6", 65, 200), mapped("F6", 129, 150),
])


class ExtractionPredicate(unittest.TestCase):
    """qualifies() must say exactly what Modules/extract_reads.nf does."""

    def test_unmapped_always_qualifies(self):
        for mode in (True, False):
            self.assertTrue(audit.qualifies(4, None, mode, 30))
            self.assertTrue(audit.qualifies(77, 500, mode, 30))

    def test_low_as_qualifies_only_when_rescue_is_on(self):
        self.assertTrue(audit.qualifies(0, 25, True, 30))
        self.assertFalse(audit.qualifies(0, 25, False, 30))

    def test_the_threshold_is_strict(self):
        # The module's awk uses score < min_as, so AS exactly at the threshold stays out.
        self.assertFalse(audit.qualifies(0, 30, True, 30))

    def test_mapped_without_as_never_qualifies(self):
        # The rescue awk counts no-AS mapped records and drops them (no_as++; next).
        self.assertFalse(audit.qualifies(0, None, True, 30))
        self.assertFalse(audit.qualifies(0, None, False, 30))

    def test_segment_bits(self):
        self.assertEqual(audit.segment_of(69), 1)
        self.assertEqual(audit.segment_of(137), 2)
        self.assertEqual(audit.segment_of(4), 0)      # unpaired: neither bit
        self.assertEqual(audit.segment_of(192), 0)    # both bits: not a clean mate


@unittest.skipUnless(HAVE_SAMTOOLS, "needs samtools to decode the synthetic SAM")
class EndToEnd(unittest.TestCase):

    def run_audit(self, mode):
        with tempfile.TemporaryDirectory() as workdir:
            sam = Path(workdir) / "synthetic.sam"
            sam.write_text(SAM)
            out = Path(workdir) / "audit.tsv"
            audit.main(["--alignment", str(sam), "--sample", "SYNTH",
                        "--low-score-mode", mode, "--min-as", "30",
                        "--tmpdir", workdir, "--out", str(out)])
            with out.open() as handle:
                rows = list(csv.DictReader(
                    (line for line in handle if not line.startswith("#")),
                    delimiter="\t"))
            return rows[0]

    def test_rescue_on_counts_every_fate(self):
        row = self.run_audit("on")
        self.assertEqual(row["fragments_touched"], "5")       # F1 F2 F3 F4 F5
        self.assertEqual(row["both_mates_extracted"], "1")    # F1
        self.assertEqual(row["orphan_fragments"], "3")        # F2 F3 F4
        self.assertEqual(row["single_end_or_mateless"], "1")  # F5
        self.assertEqual(row["lost_mate_AS_100_plus"], "1")   # F2's mate, AS 140
        self.assertEqual(row["lost_mate_AS_30_49"], "1")      # F3's mate, AS 35
        self.assertEqual(row["lost_mate_no_AS"], "1")         # F4's mate
        self.assertEqual(row["lost_mate_mapq0"], "1")         # F4's mate
        # Gray zone = no AS or AS<50: F3 and F4 of the three orphans.
        self.assertEqual(row["lost_mate_gray_zone_fraction"], f"{2 / 3:.6f}")

    def test_rescue_off_never_touches_the_low_as_fragment(self):
        row = self.run_audit("off")
        self.assertEqual(row["fragments_touched"], "4",
                         "F3 has no unmapped read, so with rescue off extraction "
                         "never sees it and it must not appear in the census")
        self.assertEqual(row["orphan_fragments"], "2")        # F2 F4 only
        self.assertEqual(row["lost_mate_AS_30_49"], "0")

    def test_untouched_fragments_stay_out_of_the_denominator(self):
        # F6 qualifies nowhere in either mode; orphan_fraction is over touched
        # fragments, not over the whole library.
        row = self.run_audit("on")
        self.assertEqual(row["orphan_fraction_of_fragments"], f"{3 / 5:.6f}")


if __name__ == "__main__":
    unittest.main()
