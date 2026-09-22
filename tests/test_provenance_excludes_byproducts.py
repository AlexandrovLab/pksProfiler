"""The code digest must not depend on build byproducts.

Provenance.code() hashed every file under scripts/ recursively, including __pycache__ --
which is gitignored, so the task hash depended on files that are not source. Fifteen
processes carry that digest, so running any helper script rewrote a .pyc and invalidated
all of them: a -resume after executing one script re-ran the whole cohort.

Observed 2026-09-22. One regenerated .pyc (build_qc_summary.cpython-310.pyc, rewritten by
running the script) turned a cohort-only resume into a full re-run, including a 6-hour
nhmmscan.

It also made the digest non-reproducible between machines -- different interpreter
versions and optimisation levels emit different bytecode -- which undercuts F16, whose
whole purpose is an equivalence record two people can compare.

Verified end to end against the real code path: rewriting a .pyc leaves the recorded
digest identical, adding a .py changes it.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = (ROOT / "lib/Provenance.groovy").read_text()


class CodeDigestSkipsByproducts(unittest.TestCase):
    def test_the_recursive_walk_filters(self):
        walk = re.search(r"eachFileRecurse\s*\{[^}]*\}", PROVENANCE)
        self.assertIsNotNone(walk, "the recursive walk is gone; update this test")
        self.assertIn("generated(f)", walk.group(0),
                      "code() must skip build byproducts when walking a directory")

    def test_the_exclusions_cover_what_bit_us(self):
        for fragment in ("__pycache__", ".pyc", ".ipynb_checkpoints"):
            self.assertIn(fragment, PROVENANCE,
                          f"{fragment} must be excluded from the code digest")

    def test_data_digest_is_untouched(self):
        # data() fingerprints large references by name/size/mtime and has no such walk;
        # this pins that the fix did not leak into it.
        self.assertIn("static String data(", PROVENANCE)

    def test_the_reason_is_recorded_in_the_source(self):
        # A future reader must not "tidy up" the filter without knowing what it cost.
        self.assertIn("gitignored", PROVENANCE)
