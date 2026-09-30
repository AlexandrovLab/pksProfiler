"""Matched normals are processed identically to tumours -- by construction.

Ludmil, 26 Sep mode memo: "Process matched normals identically. Matched normal/blood
provides a built-in specificity control and should use the same extraction, depletion,
profiling, and QC."

How a normal is represented: it is not. A sample sheet row is `patient` plus its input
(`alignment`/`bam`, or `fastq1`[,`fastq2`]); a tumour and its matched normal are two rows
with two identifiers, in the same run and so the same --sample_type lane. main.nf turns
each row into `tuple(sampleID, files)` before any process sees it, and the identifier is
used for naming only. So the proof that a normal is treated like its tumour is that
nothing CAN tell them apart:

  1. the sheet is read through those five columns and no other -- a `sample_role`,
     `tumor_normal` or similar column is ignored, not acted on;
  2. what reaches every process is (identifier, files), with no role;
  3. no workflow, process or helper-script code inspects an identifier's value, other
     than the validations enumerated below;
  4. the `tumor_*` parameters are run-wide (they apply to every wgs sample, tumour or
     normal), and no parameter selects samples by role;
  5. every place where two samples of one run can take different paths is enumerated
     below, and each is decided by the sample's own data under the same rule.

A new divergence -- a role column, a role parameter, a branch on a sample's name, a new
per-sample filter -- fails here until it is added to the enumerations with its reason.
The execution half runs one synthetic BAM as a tumour and as its normal in one run and
compares task scripts, tool records and outputs: tests/check_mode_execution.sh, mode
normal_parity (tests/fixtures/assert_normal_parity.py).
"""
import re
import stat
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()
PREFLIGHT = (ROOT / "scripts/preflight.py").read_text()
MODULES = {p.name: p.read_text() for p in sorted((ROOT / "Modules").glob("*.nf"))}
LIBS = {p.name: p.read_text() for p in sorted((ROOT / "lib").glob("*.groovy"))}
SCRIPTS = {p.name: p.read_text() for p in sorted((ROOT / "scripts").glob("*.py"))}
MODE_SCRIPT = (ROOT / "tests/check_mode_execution.sh").read_text()
PARITY = ROOT / "tests/fixtures/assert_normal_parity.py"

SHEET_COLUMNS = {"patient", "alignment", "bam", "fastq1", "fastq2"}

# Names a tumour/normal role would plausibly be given, as a column or a parameter.
ROLE_WORDS = r"(?:normal|tumou?r|matched|role|germline|blood|control|tissue)"

# 4. Parameters whose NAME says tumour. All are run-wide: read once, applied to every
# sample of a --sample_type wgs run whatever it is. Named for the lane they were
# calibrated on (Ludmil's tier thresholds from TCGA tumour breadth), not a sample role.
RUN_WIDE_TUMOR_PARAMS = {
    "tumor_targeted_assembly",         # the wgs lane's targeted assembly on/off switch
    "tumor_contig_tiers",              # which evidence tiers are eligible for it
    "tumor_contig_min_aligned_bp",     # contig-validation threshold, same for every sample
    *(f"tumor_{tier}_min_{field}"      # the read-evidence tier thresholds
      for tier in ("multi_gene", "broad_island", "extensive_island")
      for field in ("pks_reads", "clb_genes", "breadth")),
}

# 3. The only code that looks at an identifier's value: validating it.
IDENTIFIER_INSPECTIONS = {
    ("main.nf", r"sample_id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/"),     # allowed characters
    ("pks_mag.nf", r"sampleID != NO_REAL_BINS_SENTINEL"),              # drop the join sentinel
}

# 5. Every per-sample `.filter` in the graph: where two samples of one run can take
# different paths. Each decides from run-wide flags or from the sample's own data, by
# the same rule for a tumour and its normal.
PER_SAMPLE_FILTERS = {
    # main.nf: MAG assembly reads -- run-wide flags only.
    ("main.nf", 'params.sample_type == "metagenome" && enable_mags_b'),
    # main.nf: targeted assembly on samples whose OWN read-evidence tier is in
    # --tumor_contig_tiers. A normal that reaches the tier is assembled like a tumour;
    # one that does not is not, exactly as for a tumour below the tier.
    ("main.nf", "tier in contig_tiers"),
    # pks_targeted.nf: splits the aligned contigs by assembler, not by sample.
    ("pks_targeted.nf", "assembler == 'megahit'"),
    ("pks_targeted.nf", "assembler == 'metaspades'"),
    # pks_mag.nf (metagenome only): contig counts, bins, the join sentinel.
    ("pks_mag.nf", "contig_count == 0"),
    ("pks_mag.nf", "contig_count > 0"),
    ("pks_mag.nf", 'unitID != "unbinned"'),
    ("pks_mag.nf", 'unitID == "unbinned"'),
    ("pks_mag.nf", "gene_count"),
    ("pks_mag.nf", "has_bins == null"),
    ("pks_mag.nf", "sampleID != NO_REAL_BINS_SENTINEL"),
}

# 5, continued: in-process decisions that differ between two samples only through the
# sample's own input. E1 (Ludmil, 25 Sep): --extract_low_score_mapped auto extracts
# low-AS mapped reads when the input's own @PG chain runs bwa mem; a tumour and normal
# aligned by different pipelines get different extraction, by the same rule, and each
# sample's decision is logged and in by_sample/<s>/provenance/input_header.tsv.
DATA_DRIVEN_EXTRACTION = "--low-score-mode \"${params.extract_low_score_mapped}\""


def code_only(text):
    return "\n".join(line for line in text.splitlines()
                     if not line.strip().startswith(("//", "#", "*", "/*")))


def closure_body(text, start):
    """The text of a `{ ... }` closure whose opening brace ends just before `start`."""
    depth = 1
    for index in range(start, len(text)):
        depth += {"{": 1, "}": -1}.get(text[index], 0)
        if depth == 0:
            return " ".join(text[start:index].split())
    raise AssertionError("unbalanced closure")


def workflow_block(text):
    return text.split("\nworkflow {", 1)[1]


class TheSheetCarriesNoRole(unittest.TestCase):
    def test_main_reads_only_the_input_columns(self):
        named = set(re.findall(r"\brow\.(\w+)", code_only(MAIN)))
        self.assertEqual(named, SHEET_COLUMNS)
        # row[column] only ever indexes the required columns (patient, fastq1).
        self.assertIn('def required_sample_columns = alignment_input ? ["patient"] : ["patient", "fastq1"]',
                      MAIN)

    def test_preflight_reads_only_the_input_columns(self):
        named = set(re.findall(r'row\.get\("(\w+)"', PREFLIGHT))
        looped = set(re.findall(r'for column in \(([^)]*)\)', PREFLIGHT)[0].replace('"', "").replace(" ", "")
                     .split(","))
        self.assertLessEqual(named | looped, SHEET_COLUMNS)

    def test_no_role_column_is_read_anywhere(self):
        for name, text in {"main.nf": MAIN, "preflight.py": PREFLIGHT, **MODULES, **LIBS}.items():
            with self.subTest(file=name):
                self.assertNotRegex(code_only(text), rf"\brow(?:\.|\.get\(\"|\[\"){ROLE_WORDS}\w*")

    def test_each_row_becomes_identifier_plus_files_and_nothing_else(self):
        self.assertIn("tuple(row.patient.toString().trim(), file(alignment, checkIfExists: true))", MAIN)
        self.assertRegex(MAIN, r"tuple\(\s*row\.patient\.toString\(\)\.trim\(\),\s*fastq_files\s*\)")


class TheIdentifierIsOnlyAName(unittest.TestCase):
    def test_groovy_never_branches_on_an_identifiers_value(self):
        pattern = re.compile(r"(?:sampleID|sample_id|sampleId)\s*(?:==~|=~|==|!=|\bin\b)\s*[^\n]*"
                             r"|(?:sampleID|sample_id)\.(?:contains|startsWith|endsWith|matches|find|split)\(")
        found = set()
        for name, text in {"main.nf": MAIN, **MODULES, **LIBS}.items():
            for m in pattern.finditer(code_only(text)):
                found.add((name, m.group(0).strip()))
        unexplained = {f for f in found
                       if not any(f[0] == n and f[1].startswith(p) for n, p in IDENTIFIER_INSPECTIONS)}
        self.assertEqual(unexplained, set())

    def test_task_scripts_never_branch_on_the_sample_name(self):
        # Bash tests or case statements with the identifier as an operand, not a filename.
        operand = re.compile(r'(?:\[\[|\btest\b|\bcase\b)[^\n]*(?<![\w.])"?\$\{sampleID\}"?(?![\w.])\s*'
                             r'(?:==|!=|=~|\bin\b)')
        for name, text in MODULES.items():
            with self.subTest(module=name):
                self.assertIsNone(operand.search(text))

    def test_helper_scripts_never_branch_on_the_sample_name(self):
        pattern = re.compile(r"args\.sample\s*(?:==|!=|\bin\b)|args\.sample\.(?:startswith|endswith|split|"
                             r"lower|upper|find)|re\.\w+\([^)]*args\.sample")
        for name, text in SCRIPTS.items():
            with self.subTest(script=name):
                self.assertIsNone(pattern.search(text))


class NoParameterSelectsByRole(unittest.TestCase):
    def test_role_named_parameters_are_the_enumerated_run_wide_ones(self):
        declared = set(re.findall(rf"^params\.(\w*{ROLE_WORDS}\w*)\s*=", MAIN, re.M | re.I))
        self.assertEqual(declared, RUN_WIDE_TUMOR_PARAMS)

    def test_the_tumor_parameters_are_never_combined_with_a_sample_identifier(self):
        for name, text in {"main.nf": MAIN, **MODULES}.items():
            for line in code_only(text).splitlines():
                if re.search(r"params\.tumor_", line):
                    with self.subTest(file=name, line=line.strip()[:80]):
                        self.assertNotRegex(line, r"(?:sampleID|sample_id)\s*(?:==|!=|==~|=~)")

    def test_sample_type_is_one_value_for_the_whole_run(self):
        self.assertIn("params.sample_type = null", MAIN)
        self.assertNotRegex(code_only(MAIN), r"params\.sample_type\s*\[")


class EveryPerSampleDivergenceIsEnumerated(unittest.TestCase):
    def test_every_filter_is_data_driven_and_listed(self):
        found = set()
        for name, text in {"main.nf": workflow_block(MAIN), **MODULES}.items():
            for m in re.finditer(r"\.filter\s*\{", text):
                found.add((name, closure_body(text, m.end())))
        unexplained = {f for f in found if not any(f[0] == n and needle in f[1]
                                                   for n, needle in PER_SAMPLE_FILTERS)}
        self.assertEqual(unexplained, set(), "a new per-sample filter: can a normal now diverge?")
        # ... and a listed one that no longer exists is removed from the list.
        stale = {(n, needle) for n, needle in PER_SAMPLE_FILTERS
                 if not any(f[0] == n and needle in f[1] for f in found)}
        self.assertEqual(stale, set())

    def test_the_multiline_tier_filter_is_the_evidence_gate(self):
        gate = workflow_block(MAIN).split("tumor_eligible_profiles_ch = classifyPksReadEvidence.out.evidence", 1)[1]
        gate = gate.split("tumor_eligible_ids_ch", 1)[0]
        self.assertIn("tier in contig_tiers", gate)
        self.assertNotRegex(gate, r"sampleID\s*(?:==|!=|=~)")

    def test_no_process_is_skipped_by_a_when_block(self):
        for name, text in MODULES.items():
            with self.subTest(module=name):
                self.assertNotRegex(code_only(text), r"^\s*when:", )

    def test_the_data_driven_extraction_rule_is_the_same_flag_for_every_sample(self):
        extract = MODULES["extract_reads.nf"]
        self.assertIn(DATA_DRIVEN_EXTRACTION, extract)


class TheExecutionComparisonExists(unittest.TestCase):
    def test_check_mode_execution_runs_a_tumour_and_its_normal_together(self):
        code = code_only(MODE_SCRIPT)
        self.assertIn("normal_parity", re.search(r'MODES=\$\{MODES:-"([^"]+)"\}', MODE_SCRIPT).group(1))
        self.assertIn("parity_tumour \"$fixtures/wes_input.bam\" parity_normal \"$fixtures/wes_input.bam\"",
                      code)
        self.assertIn("run_mode normal_parity --sample_type wgs --profiling_method bowtie2", code)
        self.assertIn("--save_intermediates true", code.split("run_mode normal_parity", 1)[1])
        self.assertIn("assert_normal_parity.py", code)

    def test_the_comparison_covers_the_four_stages_scripts_and_outputs(self):
        self.assertTrue(PARITY.stat().st_mode & stat.S_IXUSR)
        text = PARITY.read_text()
        for needle in ("extractReads", "filterReads", "mapReads", "pksProfilerAlign",
                       "classifyPksReadEvidence", ".command.sh", ".tool_versions.tsv",
                       "depletion.qc.tsv", "counts.txt", "read_evidence.tsv",
                       "host_depleted.fastq.gz", "pks.qc.summary.tsv", "pks.master_summary.tsv",
                       "the normal ran exactly the tumour's processes"):
            with self.subTest(needle=needle):
                self.assertIn(needle, text)


if __name__ == "__main__":
    unittest.main()
