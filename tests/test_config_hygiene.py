"""nextflow.config must declare a real manifest and no profile it cannot honour.

Two failure modes this guards. A missing manifest means the pipeline cannot state its own
version, cannot enforce a minimum Nextflow release, and cannot be run reproducibly with
`-r <tag>`. A container profile with no `container` directive anywhere means
`-profile singularity` silently runs on the host instead of failing.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = (ROOT / "nextflow.config").read_text()
MODULES = list((ROOT / "Modules").glob("*.nf")) + [ROOT / "main.nf"]

CONTAINER_ENGINES = ("docker", "singularity", "podman", "shifter", "charliecloud")


def profiles():
    block = CONFIG[CONFIG.index("profiles {"):]
    return set(re.findall(r"^\s{4}(\w+)\s*[{]", block, re.M))


class ManifestTests(unittest.TestCase):
    def test_manifest_exists_with_the_required_keys(self):
        self.assertIn("manifest {", CONFIG)
        for key in ("name", "version", "nextflowVersion", "homePage", "mainScript"):
            self.assertRegex(CONFIG, rf"{key}\s*=", f"manifest is missing {key}")

    def test_nextflow_version_is_a_floor_not_a_comment(self):
        m = re.search(r"nextflowVersion\s*=\s*'([^']+)'", CONFIG)
        self.assertIsNotNone(m)
        self.assertTrue(m.group(1).startswith(">="), "nextflowVersion must be enforced, e.g. '>=24.10.0'")

    def test_version_agrees_with_the_changelog_and_citation_file(self):
        version = re.search(r"version\s*=\s*'([^']+)'", CONFIG).group(1)
        self.assertIn(version, (ROOT / "CHANGELOG.md").read_text())
        self.assertIn(version, (ROOT / "CITATION.cff").read_text())


class ProfileTests(unittest.TestCase):
    def test_no_container_profile_without_container_directives(self):
        declares_container = any(
            re.search(r"^\s*container\s+", p.read_text(), re.M) for p in MODULES
        )
        offered = profiles() & set(CONTAINER_ENGINES)
        if not declares_container:
            self.assertEqual(offered, set(),
                             f"profiles {sorted(offered)} offered but no process declares a container")

    def test_every_scheduler_profile_has_a_config_file(self):
        for name in profiles():
            m = re.search(rf"{name}\s*{{\s*includeConfig '([^']+)'", CONFIG)
            if m:
                self.assertTrue((ROOT / m.group(1)).exists(), f"{name} -> missing {m.group(1)}")


class TemplateCruftTests(unittest.TestCase):
    def test_no_unfilled_template_placeholders(self):
        for marker in ("{{", "TODO nf-core"):
            self.assertNotIn(marker, CONFIG, f"nextflow.config still contains {marker!r}")

    def test_no_commented_out_profiles(self):
        # A commented-out profile reads as a feature that exists.
        self.assertNotRegex(CONFIG, r"^\s*//\s*\w+\s*{\s*includeConfig")


if __name__ == "__main__":
    unittest.main()


class CondaEnvironmentWinsOnPath(unittest.TestCase):
    """29 Sep: activation swaps the conda base bin for the env's in place, so a tool
    earlier on the launching shell's PATH (a user minimap2 2.24) shadowed the pinned env.
    beforeScript puts the base bin first, so the env bin is first after activation."""

    def test_before_script_puts_the_conda_base_first(self):
        config = (Path(__file__).resolve().parents[1] / "nextflow.config").read_text()
        self.assertIn('if command -v conda >/dev/null 2>&1; then '
                      'export PATH="$(conda info --base)/bin:$PATH"; fi\'', config)

    def test_before_script_disables_the_user_site_directory(self):
        """
        T6 (14 Oct 2026): the same trap one level down, for Python packages.

        Python puts ~/.local/lib/pythonX.Y/site-packages AHEAD of the environment's own
        site-packages, so one `pip install --user` silently beats the pinned conda build
        in every environment and every task. On TSCC ~/.local/lib/python3.12 held
        pandas 1.5.3 and ~/.local/lib/python3.9 held matplotlib 3.8.2, while the
        environments pinned pandas 2.3.0 and matplotlib 3.9.4 -- and the tasks imported
        the ~/.local ones. conda-meta (what RUN_REPORT's SOFTWARE section reads) said the
        pinned versions were installed, truthfully, and they still were not used.

        Without this export, every version pin in conda_envs/ is advisory.
        """
        config = (Path(__file__).resolve().parents[1] / "nextflow.config").read_text()
        self.assertRegex(config, r"process\.beforeScript\s*=\s*'export PYTHONNOUSERSITE=1;")
