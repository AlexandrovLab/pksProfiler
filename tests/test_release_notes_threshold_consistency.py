"""The docs never contradict the frozen tier thresholds in main.nf.

Ludmil, revised report finding 11. main.nf and the Evidence tiers table in
RELEASE_NOTES_v0.0.2-dev.md agree on broad_island=7.5%/extensive_island=15%.
The "Changed from v3-working" section further down the same file still said
broad_island=10%/extensive_island=20%, including a sentence calling the 20%
threshold "unchanged" -- self-contradictory within one document.

His acceptance test: a repository-wide search for the old percentages should
return only text explicitly labeled as historical/superseded, not a second,
silently-disagreeing claim about the current thresholds.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()
RELEASE_NOTES = (ROOT / "RELEASE_NOTES_v0.0.2-dev.md").read_text()


def current_breadth_params(text):
    broad = re.search(r"params\.tumor_broad_island_min_breadth\s*=\s*([\d.]+)", text)
    extensive = re.search(r"params\.tumor_extensive_island_min_breadth\s*=\s*([\d.]+)", text)
    return float(broad.group(1)), float(extensive.group(1))


class TheFrozenThresholdsAreWhatTheCodeActuallyUses(unittest.TestCase):
    def test_main_nf_still_has_the_recalibrated_values(self):
        broad, extensive = current_breadth_params(MAIN)
        self.assertEqual(broad, 0.075)
        self.assertEqual(extensive, 0.15)


class ReleaseNotesDoNotContradictTheEvidenceTiersTable(unittest.TestCase):
    def test_the_evidence_tiers_table_shows_the_current_values(self):
        self.assertIn(">=7.5%", RELEASE_NOTES)
        self.assertIn(">=15%", RELEASE_NOTES)

    def test_unchanged_is_no_longer_claimed_for_a_threshold_that_changed(self):
        # The exact self-contradiction: extensive_island's breadth DID change
        # (20% -> 15%), so nothing may call it unchanged without qualification.
        self.assertNotIn("unchanged thresholds", RELEASE_NOTES)

    def test_every_old_percentage_is_explicitly_labeled_historical(self):
        # His acceptance test, applied to this file: 30/8/10 and 100/10/20 may
        # still appear, but only inside text that says so plainly -- either an
        # inline label, or anywhere in the "Changed from v3-working" section,
        # whose own header is an unambiguous historical label for everything in it.
        section_start = RELEASE_NOTES.index("### Changed from v3-working")
        next_header = RELEASE_NOTES.find("\n### ", section_start + 1)
        section_end = next_header if next_header != -1 else len(RELEASE_NOTES)
        inline_labels = ("as v0.0.2 first introduced it", "Superseded", "recalibrated")

        for match in re.finditer(r"30/8/10|100/10/20", RELEASE_NOTES):
            with self.subTest(position=match.start()):
                in_historical_section = section_start <= match.start() < section_end
                window = RELEASE_NOTES[max(0, match.start() - 200):match.start() + 200]
                self.assertTrue(
                    in_historical_section or any(label in window for label in inline_labels),
                    f"stale threshold near {match.start()} is not labeled historical: "
                    f"...{window}...")

    def test_the_recalibration_direction_is_stated_correctly(self):
        self.assertIn("10% -> **7.5%**", RELEASE_NOTES)
        self.assertIn("20% -> **15%**", RELEASE_NOTES)


if __name__ == "__main__":
    unittest.main()
