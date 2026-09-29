"""One version, everywhere it is written (GAMEPLAN v1.0.0-rc1, "Release mechanics").

"nextflow.config, CITATION.cff, CHANGELOG.md, git SHA ... all agree ... Add a test that
fails if any of these disagree." The manifest in nextflow.config is the source of truth;
everything that prints a version must either read it or state the same string:

  * CITATION.cff `version:`
  * the top release heading of CHANGELOG.md (`## v<version> — ...`)
  * the README version badge
  * the workflow figure's legend (workflow_logo/make_figure.py reads the manifest, and the
    committed v2.svg must carry what it would draw)
  * main.nf: `--help` and the launch record take `workflow.manifest.version`, and no
    literal "pksProfiler <version>" string is hard-coded anywhere in main.nf or lib/
  * RUN_REPORT.txt (scripts/finalize_provenance.py) prints the recorded manifest version
    and the git commit, and hard-codes no version of its own

Bumping the version is then one edit to nextflow.config plus the three files this test
names, and a release cannot go out with two different numbers on it.
"""
import importlib.util
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = (ROOT / "nextflow.config").read_text()
MAIN = (ROOT / "main.nf").read_text()

MANIFEST = re.search(r"manifest\s*\{(.*?)\n\}", CONFIG, re.S).group(1)
VERSION = re.search(r"^\s*version\s*=\s*'([^']+)'", MANIFEST, re.M).group(1)

# A version-looking literal: 0.0.2, v0.0.2dev, 1.0.0-rc1 ...
VERSION_LITERAL = re.compile(r"\bv?\d+\.\d+\.\d+(?:[-.]?(?:dev|rc\d*|alpha\d*|beta\d*))?\b")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class VersionConsistency(unittest.TestCase):
    def test_the_manifest_declares_a_version(self):
        self.assertRegex(VERSION, r"^\d+\.\d+\.\d+")

    def test_citation_cff(self):
        cff = (ROOT / "CITATION.cff").read_text()
        match = re.search(r"^version:\s*['\"]?([^'\"\s]+)", cff, re.M)
        self.assertIsNotNone(match, "CITATION.cff has no version field")
        self.assertEqual(match.group(1), VERSION)

    def test_top_changelog_release_heading(self):
        changelog = (ROOT / "CHANGELOG.md").read_text()
        headings = re.findall(r"^## v(\S+)", changelog, re.M)
        self.assertTrue(headings, "CHANGELOG.md has no '## v<version>' release heading")
        self.assertEqual(headings[0], VERSION,
                         "the newest CHANGELOG release heading must be the manifest version")

    def test_readme_badge(self):
        readme = (ROOT / "README.md").read_text()
        badges = re.findall(r"badge/version-([^-)]+(?:--[^-)]+)*)-", readme)
        self.assertTrue(badges, "README has no version badge")
        for badge in badges:
            # shields.io escapes a literal '-' as '--'
            self.assertEqual(badge.replace("--", "-"), VERSION)

    def test_workflow_figure_reads_the_manifest_and_the_svg_agrees(self):
        figure = load("make_figure", ROOT / "workflow_logo/make_figure.py")
        self.assertEqual(figure.VERSION, VERSION)
        source = (ROOT / "workflow_logo/make_figure.py").read_text()
        self.assertNotRegex(source, r'"pksProfiler v\d')
        svg = (ROOT / "workflow_logo/v2.svg").read_text()
        labels = re.findall(r"pksProfiler v([^<\s]+)", svg)
        self.assertEqual(labels, [VERSION],
                         "workflow_logo/v2.svg is stale: rerun python3 workflow_logo/make_figure.py")

    def test_help_and_launch_record_take_the_manifest_version(self):
        self.assertIn("pksProfiler ${workflow.manifest.version} -- ${workflow.manifest.description}", MAIN)
        self.assertIn("manifest_version: workflow.manifest.version", MAIN)

    def test_no_version_literal_is_hard_coded_in_the_workflow(self):
        sources = {"main.nf": MAIN}
        sources.update({p.name: p.read_text() for p in (ROOT / "lib").glob("*.groovy")})
        for name, text in sources.items():
            code = "\n".join(line for line in text.splitlines()
                             if not line.strip().startswith(("//", "*", "/*")))
            with self.subTest(file=name):
                self.assertNotRegex(code, r"pksProfiler\s+v?\d+\.\d+",
                                    f"{name} hard-codes a pksProfiler version")

    def test_run_report_prints_the_recorded_version_and_commit(self):
        source = (ROOT / "scripts/finalize_provenance.py").read_text()
        self.assertIn("f\"pksProfiler {code.get('manifest_version', '?')}\"", source)
        self.assertIn("{code.get('commit')} ({code.get('commit_source')})", source)
        strings = re.findall(r"(?:f?\"[^\"\n]*\"|f?'[^'\n]*')", source)
        for literal in strings:
            with self.subTest(literal=literal):
                self.assertNotRegex(literal, VERSION_LITERAL,
                                    "finalize_provenance.py hard-codes a version string")

    def test_run_report_header_is_the_manifest_version(self):
        provenance = load("finalize_provenance", ROOT / "scripts/finalize_provenance.py")
        text = provenance.render({
            "status": "succeeded",
            "code": {"manifest_version": VERSION, "commit": "0123abcd" * 5,
                     "commit_source": "git", "dirty": False},
            "invocation": {"params": {}},
            "run": {"session_id": "S"},
        })
        self.assertIn(f"pksProfiler {VERSION}", text)
        self.assertIn("0123abcd" * 5, text)
        other = {v.lstrip("v") for v in VERSION_LITERAL.findall(text)} - {VERSION}
        self.assertFalse(other, f"RUN_REPORT prints another version: {other}")


if __name__ == "__main__":
    unittest.main()
