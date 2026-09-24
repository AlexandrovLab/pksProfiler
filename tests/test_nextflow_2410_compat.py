"""main.nf actually compiles on the Nextflow version the manifest declares as its floor.

Ludmil's email, tracker item m-1: "Github shows your tests failing so you should
probably fix that." Investigation found this had nothing to do with the
pandas/numpy environment issue seen in local test runs all session (CI's unit-tests
job never installs pandas, so that failure mode cannot occur there) -- CI's actual,
separate "workflow-validation" job has been failing since before this session
started, because `nextflow.config` declares `nextflowVersion = '>=24.10.0'` and CI
pins exactly that version, but `nextflow run main.nf --help` never actually compiled
on 24.10.0. Reproduced directly (NXF_VER=24.10.0) against the pre-fix baseline
commit (870e4fa) and fixed in three places, all variants of the same underlying
Nextflow 24.10.0 quirk: a `def <newVar> = <ImplicitBinding>...` declaration, where
the initializer's first token is a reserved/implicit binding name (`log`, `Channel`),
fails to compile with "Variable `X` already defined in the process scope". Nextflow
26.x (used throughout the rest of this session) has no such restriction, which is
exactly why this went unnoticed until checked against the literal declared floor.

    1. nextflow.config referenced params.conda_cache_dir at config-evaluation time,
       before main.nf's own `params.conda_cache_dir = null` declaration had run --
       "Unknown config attribute". Fixed by declaring it in nextflow.config's own
       params {} block too, the same pattern already used for `outdir`.
    2. `Preflight.run(..., log)` passed the bare implicit `log` binding as a
       positional argument. Fixed by wrapping it in a closure whose body only ever
       uses the syntactically-safe `log.info(...)`/`log.warn(...)`/`log.error(...)`
       method-call form.
    3. `def strain_typing_summary_ch = Channel.empty()` (capital C) in two places in
       pks_mag.nf. Fixed by using the already-proven-safe lowercase `channel.empty()`
       alias -- Nextflow documents both spellings as equivalent.

These tests guard against reintroducing any of the three, plus one direct
compile-and-run check for anyone with Nextflow available.
"""
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()
NEXTFLOW_CONFIG = (ROOT / "nextflow.config").read_text()
PREFLIGHT_GROOVY = (ROOT / "lib/Preflight.groovy").read_text()
MODULE_FILES = {p.name: p.read_text() for p in (ROOT / "Modules").glob("*.nf")}

import re


def code_only(text):
    return "\n".join(line for line in text.splitlines()
                      if not line.strip().startswith(("#", "//")))


MAIN_CODE = code_only(MAIN)
MODULE_CODE = {name: code_only(text) for name, text in MODULE_FILES.items()}


class NoBareImplicitBindingIsAssignedToANewLocal(unittest.TestCase):
    """The exact shape that fails to compile on Nextflow 24.10.0."""

    def test_no_def_assigns_a_new_local_from_capitalized_channel(self):
        pattern = re.compile(r"def\s+\w+\s*=\s*Channel\.")
        for name, text in {"main.nf": MAIN_CODE, **MODULE_CODE}.items():
            with self.subTest(file=name):
                self.assertNotRegex(text, pattern,
                                    f"{name} declares a new local from Channel.something() "
                                    "(capitalized) -- use lowercase channel. instead")

    def test_no_def_assigns_a_new_local_from_bare_log(self):
        pattern = re.compile(r"def\s+\w+\s*=\s*log\b")
        for name, text in {"main.nf": MAIN_CODE, **MODULE_CODE}.items():
            with self.subTest(file=name):
                self.assertNotRegex(text, pattern,
                                    f"{name} assigns the bare `log` binding to a new local")


class PreflightRunNeverReceivesTheBareLogBinding(unittest.TestCase):
    def test_the_call_site_passes_a_closure_not_bare_log(self):
        call_index = MAIN_CODE.index('Preflight.run("${params.scripts}/preflight.py"')
        call_site = MAIN_CODE[call_index:call_index + 200]
        self.assertNotRegex(call_site, r",\s*log\)",
                            "Preflight.run's last argument must not be the bare `log` binding")
        self.assertIn("preflight_log", call_site)

    def test_the_closure_only_uses_the_safe_method_call_form(self):
        closure = MAIN_CODE[MAIN_CODE.index("def preflight_log ="):
                            MAIN_CODE.index("def preflight_status")]
        for safe_call in ("log.error(", "log.warn(", "log.info("):
            with self.subTest(call=safe_call):
                self.assertIn(safe_call, closure)

    def test_preflight_groovy_accepts_a_closure_not_an_object_named_log(self):
        signature = PREFLIGHT_GROOVY[PREFLIGHT_GROOVY.index("static int run("):
                                     PREFLIGHT_GROOVY.index(")", PREFLIGHT_GROOVY.index("static int run(")) + 1]
        self.assertIn("Closure logSink", signature)
        self.assertNotIn("Object log", signature)


class NextflowConfigDeclaresEveryParamItReferencesAtConfigTime(unittest.TestCase):
    def test_conda_cache_dir_is_declared_in_nextflow_configs_own_params_block(self):
        params_block = NEXTFLOW_CONFIG[NEXTFLOW_CONFIG.index("params {"):
                                       NEXTFLOW_CONFIG.index("\n}", NEXTFLOW_CONFIG.index("params {"))]
        self.assertIn("conda_cache_dir", params_block)

    def test_it_is_declared_before_conda_cachedir_reads_it(self):
        params_decl = NEXTFLOW_CONFIG.index("conda_cache_dir = null")
        conda_cachedir_use = NEXTFLOW_CONFIG.index("conda.cacheDir = params.conda_cache_dir")
        self.assertLess(params_decl, conda_cachedir_use)


@unittest.skipUnless(shutil.which("nextflow"), "nextflow not on PATH")
class TheWorkflowActuallyCompilesOnTheDeclaredFloor(unittest.TestCase):
    """Direct proof, not just a structural proxy -- only runs if nextflow is available.

    Does not pin NXF_VER itself (that would force a slow first-time download in
    environments without 24.10.0 already cached); running under whatever
    `nextflow` resolves to on PATH still catches a regression on newer versions,
    and the three checks above are what specifically pins the 24.10.0 contract.
    """

    def test_help_compiles_and_runs(self):
        result = subprocess.run(["nextflow", "run", str(ROOT / "main.nf"), "--help"],
                                capture_output=True, text=True, cwd=ROOT)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("pksProfiler", result.stdout)


if __name__ == "__main__":
    unittest.main()
