"""`-resume` across the rc1 boundary is documented, and the mechanism it describes is real.

rc2 review item 2.8: "Document rc1-to-rc2 resume behaviour. Because the scripts digest
changed, -resume will not transparently reuse most rc1 tasks. State this explicitly so
production reruns do not silently span mixed code states."

The risk is silent: a results tree half-written by rc1 and half by v1.0.0, with one
RUN_REPORT.txt describing only the later session, is a cohort nobody can attribute to a
version string. So the document has to say the uncomfortable part -- plan for a full
re-run -- and has to tell a reader how to check a tree they already have.

Two things are tested. First, that docs/running/resume.md exists, is reachable, and makes
the claims the item asks for. Second, that the two mechanisms it blames are still in the
code: a document explaining an invalidation path that no longer exists would be worse
than none. The arithmetic in the document (53 of 60 processes) is a measurement against
the v1.0.0-rc1 tag and is deliberately not recomputed here -- CI checkouts are shallow
and a test that silently passes when the tag is absent would prove nothing.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs/running/resume.md"
TEXT = DOC.read_text() if DOC.exists() else ""
MAIN = (ROOT / "main.nf").read_text()
MODULES = {p.name: p.read_text() for p in (ROOT / "Modules").glob("*.nf")}


class TheDocumentExistsAndIsReachable(unittest.TestCase):
    def test_it_exists(self):
        self.assertTrue(DOC.exists(), "docs/running/resume.md is missing")

    def test_the_readme_links_it(self):
        self.assertIn("docs/running/resume.md", (ROOT / "README.md").read_text())

    def test_the_changelog_says_it_too(self):
        # Someone upgrading reads the changelog, not the docs tree.
        changelog = (ROOT / "CHANGELOG.md").read_text()
        self.assertIn("resume", changelog.lower())
        self.assertIn("docs/running/resume.md", changelog)


class ItMakesTheClaimsTheItemAsksFor(unittest.TestCase):
    def test_it_names_both_invalidation_causes(self):
        self.assertIn("dep_digest", TEXT)
        self.assertIn("conda", TEXT.lower())

    def test_it_states_the_conclusion_rather_than_implying_it(self):
        self.assertIn("full re-run", TEXT)

    def test_it_says_which_processes_can_still_be_reused(self):
        # The honest answer is "one, and only on FASTQ input".
        self.assertIn("filterReads", TEXT)
        self.assertIn("FASTQ", TEXT)

    def test_it_tells_a_reader_how_to_check_a_tree_they_already_have(self):
        for token in ("runs/provenance", "code.commit", "trace.txt", "CACHED"):
            with self.subTest(token=token):
                self.assertIn(token, TEXT)

    def test_it_says_what_does_not_invalidate_a_task(self):
        # Otherwise every resource tweak gets treated as a code change.
        self.assertIn("errorStrategy", TEXT)
        self.assertIn("publishDir", TEXT)

    def test_it_names_the_execution_check_that_proves_the_mechanism(self):
        self.assertIn("check_resume_reference_change.sh", TEXT)
        self.assertTrue((ROOT / "tests/check_resume_reference_change.sh").exists())


class TheMechanismItDescribesStillExists(unittest.TestCase):
    def test_the_scripts_digest_is_still_computed_and_still_interpolated(self):
        self.assertIn("scripts             : Provenance.code(params.scripts)", MAIN)
        carriers = [name for name, text in MODULES.items()
                    if "dep_digest?.scripts" in text]
        self.assertGreaterEqual(len(carriers), 8,
                                "almost no module carries the scripts digest any more: "
                                "docs/running/resume.md describes a mechanism that is gone")

    def test_most_processes_still_declare_a_conda_environment(self):
        declared = sum(len(re.findall(r"\n\s*conda\s", text)) for text in MODULES.values())
        self.assertGreaterEqual(declared, 40)

    def test_every_declared_environment_file_exists(self):
        names = set()
        for text in MODULES.values():
            names |= set(re.findall(r"conda_envs/(\w+\.yml)", text))
        for name in sorted(re.findall(r"conda_envs/(\w+\.yml)", MAIN)):
            names.add(name)
        self.assertTrue(names)
        for name in sorted(names):
            with self.subTest(env=name):
                self.assertTrue((ROOT / "conda_envs" / name).exists())


if __name__ == "__main__":
    unittest.main()
