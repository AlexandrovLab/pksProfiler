"""V8a and V8b must stay two checks, not one.

Ludmil, 25 Sep report, section 4:

  V8a. Locus-detection specificity ... "This asks whether the sequence really belongs to
       the pks/clb locus."
  V8b. Carrier attribution ... "This asks which organism is the likely carrier and should
       not be conflated with locus-detection specificity."

The two lanes share an input (the retained pks-island reads) and an aligner, which is
exactly why they are easy to collapse back into one "specificity" number. This file is the
guard: separate panels, separate scripts, separate output directories, separate columns,
separate verdict vocabularies, and neither lane's panel able to serve the other's purpose.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
V8A_MODULE = (ROOT / "Modules/pks_specificity.nf").read_text()
V8B_MODULE = (ROOT / "Modules/pks_carrier.nf").read_text()
V8A_SCRIPT = (ROOT / "scripts/specificity_decision.py").read_text()
V8B_SCRIPT = (ROOT / "scripts/carrier_attribution.py").read_text()
MAIN = (ROOT / "main.nf").read_text()
BUILDER = (ROOT / "scripts/build_specificity_panels.py").read_text()
V8A_PANEL = (ROOT / "ref/specificity_panels/v8a_locus_decoy_panel.tsv").read_text()
V8B_PANEL = (ROOT / "ref/specificity_panels/v8b_carrier_panel.tsv").read_text()

V8A_VERDICTS = ("SURVIVES", "NOT_DISTINGUISHABLE_FROM_CROSS_MAPPING", "UNDERPOWERED")
V8B_VERDICTS = ("CARRIER_ATTRIBUTED", "CARRIER_GENUS_ONLY", "CARRIER_AMBIGUOUS",
                "CARRIER_UNDETERMINED")


class TwoPanelsWithOppositeMembershipTests(unittest.TestCase):
    def test_the_two_panel_lists_exist_and_are_different_files(self):
        self.assertTrue((ROOT / "ref/specificity_panels/v8a_locus_decoy_panel.tsv").is_file())
        self.assertTrue((ROOT / "ref/specificity_panels/v8b_carrier_panel.tsv").is_file())

    def test_the_decoy_panel_covers_every_category_the_reviewer_named(self):
        for category in ("enterobacterales_pks_negative", "nrps_pks_homolog_genome",
                          "gut_background", "oral_background"):
            self.assertIn(category, V8A_PANEL, f"V8a panel has no {category} member")
        # GRCh38 and UniVec are not accessions, so they are builder flags, not TSV rows.
        self.assertIn("--grch38", BUILDER)
        self.assertIn("--univec", BUILDER)
        self.assertIn("UniVec", BUILDER)

    def test_the_carrier_panel_covers_every_genus_the_reviewer_named(self):
        for genus in ("escherichia", "klebsiella", "citrobacter", "enterobacter"):
            self.assertIn(genus, V8B_PANEL, f"V8b panel has no {genus} candidate")

    def test_the_two_panels_are_screened_in_opposite_directions(self):
        self.assertIn("--max-decoy-island-coverage", BUILDER)
        self.assertIn("--min-carrier-island-coverage", BUILDER)
        self.assertIn("rejected_pks_positive", BUILDER)
        self.assertIn("rejected_not_pks_positive", BUILDER)

    def test_the_builder_refuses_v8a_only_sources_on_a_v8b_build(self):
        # UniVec, GRCh38 and MIBiG clusters are pks-negative by construction, so putting
        # them in a carrier panel would be a category error, not a configuration choice.
        self.assertIn("belong to the V8a locus panel", BUILDER)


def nf_code_only(text):
    """A .nf file's executable lines: '//' comments and in-script '#' comments dropped."""
    return "\n".join(line for line in text.splitlines()
                      if not line.strip().startswith(("//", "#")))


def python_code_only(text):
    """Executable lines only: the module docstring and '#' comments are stripped.

    Each script points at the other one by name in its docstring -- that cross-reference
    is the whole point of keeping them apart, and it must not read as a leak. What must
    not cross is a verdict the script can actually EMIT.
    """
    body = text.split('"""', 2)[2] if text.count('"""') >= 2 else text
    lines = []
    in_doc = False
    for line in body.splitlines():
        stripped = line.strip()
        if stripped.count('"""') == 1:
            in_doc = not in_doc
            continue
        if in_doc or stripped.startswith("#"):
            continue
        lines.append(line)
    return "\n".join(lines)


class VocabulariesDoNotLeak(unittest.TestCase):
    V8A_CODE = python_code_only(V8A_SCRIPT)
    V8B_CODE = python_code_only(V8B_SCRIPT)

    def test_no_v8b_verdict_can_be_emitted_by_the_v8a_decision(self):
        for verdict in V8B_VERDICTS:
            self.assertNotIn(verdict, self.V8A_CODE)

    def test_no_v8a_verdict_can_be_emitted_by_the_v8b_attribution(self):
        for verdict in V8A_VERDICTS:
            self.assertNotIn(verdict, self.V8B_CODE)

    def test_each_script_still_points_at_the_other_one(self):
        # The separation is documented, not silent: a reader of either script is told
        # which question it does NOT answer and where that one lives.
        self.assertIn("carrier_attribution.py", V8A_SCRIPT)
        self.assertIn("specificity_decision.py", V8B_SCRIPT)

    def test_the_v8a_lane_never_names_a_carrier_in_its_executable_part(self):
        # The header comment names V8b on purpose; nothing the lane RUNS may.
        code = nf_code_only(V8A_MODULE)
        self.assertNotIn("carrier", code.lower())

    def test_the_v8b_lane_never_claims_specificity_in_its_executable_part(self):
        # build_specificity_panels.py builds BOTH panels, so its filename is not a claim
        # about this lane; every other mention of the word would be.
        code = nf_code_only(V8B_MODULE).replace("build_specificity_panels.py", "")
        code = code.replace("ref/specificity_panels/", "")
        self.assertNotIn("specificity", code.lower())

    def test_the_v8b_script_says_a_lost_read_is_not_evidence_against_the_locus(self):
        self.assertIn("is not evidence", V8B_SCRIPT)
        self.assertIn("evidence about the organism", V8B_SCRIPT)


class SeparateOutputsAndColumns(unittest.TestCase):
    REPORT = (ROOT / "scripts/build_cohort_report.py").read_text()
    SUMMARY = (ROOT / "scripts/build_master_summary.py").read_text()

    def test_each_lane_publishes_under_its_own_directory(self):
        self.assertIn("/locus_specificity\"", V8A_MODULE)
        self.assertIn("/carrier_attribution\"", V8B_MODULE)

    def test_the_two_manifest_kinds_are_distinct(self):
        self.assertIn("tuple(sampleID, 'locus_specificity', summary)", MAIN)
        self.assertIn("tuple(sampleID, 'carrier_attribution', summary)", MAIN)

    def test_the_master_summary_keeps_two_column_groups(self):
        self.assertIn("V8A_COLUMNS", self.SUMMARY)
        self.assertIn("V8B_COLUMNS", self.SUMMARY)
        v8a = re.search(r"V8A_COLUMNS = \[(.*?)\]", self.SUMMARY, re.S).group(1)
        v8b = re.search(r"V8B_COLUMNS = \[(.*?)\]", self.SUMMARY, re.S).group(1)
        self.assertNotIn("carrier", v8a)
        self.assertNotIn("specificity", v8b)

    def test_the_cohort_report_renders_two_sections(self):
        self.assertIn("def locus_specificity_section(", self.REPORT)
        self.assertIn("def carrier_attribution_section(", self.REPORT)
        self.assertIn("locus_specificity_section(records)", self.REPORT)
        self.assertIn("carrier_attribution_section(records)", self.REPORT)


class BothLanesAreOptionalAndOffByDefault(unittest.TestCase):
    def test_both_flags_default_to_false(self):
        self.assertRegex(MAIN, r"params\.pks_locus_specificity\s*=\s*false")
        self.assertRegex(MAIN, r"params\.pks_carrier_attribution\s*=\s*false")

    def test_both_flags_go_through_a_normalized_boolean_local(self):
        for local in ("pks_locus_specificity_b", "pks_carrier_attribution_b"):
            self.assertRegex(MAIN, rf"def {local}\s*=\s*params\.\w+\.toString\(\)\.toBoolean\(\)")
        code = "\n".join(line for line in MAIN.splitlines()
                          if not line.strip().startswith("//"))
        self.assertNotRegex(code, r"if\s*\([^)]*\bparams\.pks_locus_specificity\b")
        self.assertNotRegex(code, r"if\s*\([^)]*\bparams\.pks_carrier_attribution\b")

    def test_neither_lane_can_run_without_alignment_profiling(self):
        self.assertIn("--pks_locus_specificity / --pks_carrier_attribution need alignment profiling",
                       MAIN)

    def test_the_island_reads_are_extracted_once_for_all_three_consumers(self):
        # extractPksIslandReads can be invoked only once; taxonomy, V8a and V8b share it.
        self.assertIn("if (do_align && (pks_taxa_b || pks_locus_specificity_b || pks_carrier_attribution_b))",
                       MAIN)
        self.assertEqual(MAIN.count("extractPksIslandReads("), 1)

    def test_both_panels_are_tracked_dependencies(self):
        self.assertIn("pks_specificity_panel : Provenance.data(", MAIN)
        self.assertIn("pks_carrier_panel     : Provenance.data(", MAIN)
        self.assertIn("${params.dep_digest?.pks_carrier_panel}", V8B_MODULE)

    def test_the_run_report_records_which_check_ran_against_which_panel(self):
        for key in ("pks_locus_specificity    : pks_locus_specificity_b",
                     "pks_carrier_attribution  : pks_carrier_attribution_b",
                     "pks_carrier_panel_index            : params.pks_carrier_panel_index"):
            self.assertIn(key, MAIN)


if __name__ == "__main__":
    unittest.main()
