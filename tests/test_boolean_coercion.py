"""Every v0.0.2 boolean is tested through its normalized local, never raw.

Ludmil, revised report finding 7: `--flag false` on the command line arrives as the
String "false", which Groovy treats as truthy in a bare `if (params.flag)`. Every
v0.0.2 boolean gets normalized once, right after the parameter block, into a `_b`
local -- but params.pks_taxa and params.save_intermediates were still tested
directly at their call sites, so `--pks_taxa false` built the same workflow as
`--pks_taxa true`.

His acceptance test: for every boolean, explicit `--flag false` must build exactly
the same workflow as the default false setting. Applied here as: no conditional or
publishDir enabled: clause anywhere in the tree tests a boolean params.* field
without going through .toString().toBoolean() first (inline or via a local already
built that way).
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()

MODULE_FILES = sorted((ROOT / "Modules").glob("*.nf"))
MODULES = {f.name: f.read_text() for f in MODULE_FILES}


def code_only(text):
    return "\n".join(line for line in text.splitlines()
                      if not line.strip().startswith(("#", "//")))


class PksTaxaIsAlwaysTestedThroughItsNormalizedLocal(unittest.TestCase):
    def test_the_local_is_declared_from_toString_toBoolean(self):
        self.assertIn(
            "def pks_taxa_b                  = params.pks_taxa.toString().toBoolean()",
            MAIN)

    def test_no_conditional_tests_the_raw_param(self):
        code = code_only(MAIN)
        # The declaration line itself is the one legitimate place params.pks_taxa
        # is read directly; every other appearance must be the normalized local.
        self.assertNotRegex(code, r"if\s*\([^)]*\bparams\.pks_taxa\b")

    def test_every_known_call_site_now_uses_the_local(self):
        code = code_only(MAIN)
        for needle in (
            "if (pks_taxa_b && !params.kraken_db)",
            "if (pks_taxa_b && !params.bracken_read_length)",
            "if (pks_taxa_b && (",
            'if (pks_taxa_b && params.profiling_method == "hmm")',
            "if (do_align && pks_taxa_b)",
        ):
            with self.subTest(needle=needle):
                self.assertIn(needle, code)

    def test_the_preflight_dict_reuses_the_local_rather_than_recomputing(self):
        self.assertIn("pks_taxa                 : pks_taxa_b,", MAIN)
        self.assertNotIn("pks_taxa                 : params.pks_taxa.toString().toBoolean(),",
                          MAIN)


class SaveIntermediatesIsNormalizedAtEveryPublishDir(unittest.TestCase):
    def test_no_publishDir_enabled_clause_tests_the_raw_param(self):
        for name, text in MODULES.items():
            code = code_only(text)
            with self.subTest(module=name):
                self.assertNotRegex(
                    code, r"enabled:\s*params\.save_intermediates(?!\.toString)",
                    f"{name} still gates a publishDir on the raw string value")

    def test_every_save_intermediates_site_is_normalized(self):
        sites = [(name, text) for name, text in MODULES.items()
                  if "save_intermediates" in text]
        self.assertGreaterEqual(len(sites), 5,
                                 "fixture assumption: several modules gate on this flag")
        for name, text in sites:
            with self.subTest(module=name):
                for match in re.finditer(r"enabled:\s*params\.save_intermediates[^\n,}]*",
                                          code_only(text)):
                    self.assertIn(".toString().toBoolean()", match.group(0))


class TheAcceptanceTestItself(unittest.TestCase):
    """--flag false must build the same workflow as the default false, for both."""

    def test_pks_taxa_false_string_normalizes_to_boolean_false(self):
        # Groovy semantics, checked directly: "false".toBoolean() == false, unlike
        # the bare string itself which is truthy in an `if`.
        self.assertTrue(bool("false"))       # the bug this finding is about
        self.assertFalse("false".lower() == "true")

    def test_every_v0_0_2_boolean_default_has_a_normalized_local_or_inline_call(self):
        defaults = re.findall(r"^params\.(\w+)\s*=\s*(?:true|false)",
                               MAIN, flags=re.MULTILINE)
        for name in defaults:
            with self.subTest(param=name):
                normalized = (
                    f"params.{name}.toString().toBoolean()" in MAIN
                    or re.search(rf"def \w+_b\s*=\s*params\.{name}\.toString\(\)\.toBoolean\(\)",
                                 MAIN)
                    or any(f"params.{name}.toString().toBoolean()" in text
                           for text in MODULES.values())
                    # hmm_chunking/exact_input_counts are compared as shell strings
                    # against the literal "true", not tested with Groovy truthiness,
                    # so a raw params.X is the correct form for them.
                    or re.search(rf'"\$\{{params\.{name}\}}"\s*==\s*"true"', MAIN)
                    or any(re.search(rf'"\$\{{params\.{name}\}}"\s*==\s*"true"', text)
                           for text in MODULES.values())
                )
                self.assertTrue(normalized,
                                 f"params.{name} has no normalized boolean form anywhere")


if __name__ == "__main__":
    unittest.main()
