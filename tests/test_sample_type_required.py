"""--sample_type has no default and "auto" is no longer an accepted value.

Ludmil, revised report finding 9: "auto" was accepted by both the default and
the --sample_type validation, but nothing ever resolved it to a biological lane
-- every routing decision in main.nf tests for a concrete value ("metagenome",
"tumor_wgs"), so a run with --sample_type auto (or none at all, since it was the
default) matched none of them. Read-tier classification, contig analysis, and
MAG/taxonomy all silently never ran; alignment and the gene-count table still
did, so the run looked complete rather than failed.

His acceptance test: every auto mode must either resolve deterministically before
task submission and log the resolution, or fail loudly before any task runs.
Applied here as the second branch, his own first preference: no default,
required, fails at the same "Missing required parameter" check as --sample,
--hg38_db, and --t2t_phix_db, before Preflight.run and before any task is
submitted.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()


def code_only(text):
    return "\n".join(line for line in text.splitlines()
                      if not line.strip().startswith(("#", "//")))


class SampleTypeHasNoDefault(unittest.TestCase):
    def test_the_default_is_null_not_auto(self):
        self.assertIn("params.sample_type = null", MAIN)
        self.assertNotRegex(MAIN, r'params\.sample_type\s*=\s*"auto"')

    def test_auto_is_no_longer_an_accepted_value(self):
        code = code_only(MAIN)
        m = re.search(r"def valid_sample_types = (\[[^\]]*\])", code)
        self.assertIsNotNone(m, "valid_sample_types list not found")
        values = eval(m.group(1))
        self.assertNotIn("auto", values)
        self.assertEqual(set(values), {"metagenome", "tumor_wgs", "tumor_wes", "tumor_rna"})


class MissingSampleTypeFailsBeforeAnyTaskRuns(unittest.TestCase):
    def test_it_is_checked_alongside_the_other_required_parameters(self):
        block = MAIN[MAIN.index("Required input validation"):
                      MAIN.index("Required input validation") + 1200]
        self.assertIn("if (!params.sample_type)", block)
        self.assertIn("Missing required parameter: --sample_type", block)

    def test_the_check_happens_before_the_first_process_is_invoked(self):
        # Preflight.run itself is a plain script call, not a submitted Nextflow
        # task; the first real task is extractReads(...), further down.
        body = MAIN.split("workflow {", 1)[1]
        sample_type_check = body.index("if (!params.sample_type)")
        first_task = body.index("extractReads(")
        self.assertLess(sample_type_check, first_task,
                         "sample_type must be validated before any task is submitted")

    def test_the_error_message_names_every_supported_value(self):
        self.assertIn(
            "Missing required parameter: --sample_type. Supported: "
            "metagenome, tumor_wgs, tumor_wes, tumor_rna",
            MAIN)


class DownstreamRoutingIsUnaffectedByTheRemovalOfAuto(unittest.TestCase):
    """auto never matched any of these; removing it changes nothing here."""

    def test_prefilter_mode_auto_resolution_still_only_checks_metagenome(self):
        code = code_only(MAIN)
        self.assertIn('params.sample_type == "metagenome" && params.kraken_db', code)

    def test_mag_gating_is_unchanged(self):
        code = code_only(MAIN)
        self.assertIn('enable_mags_b && params.sample_type != "metagenome"', code)


if __name__ == "__main__":
    unittest.main()
