"""The run records the prefilter mode it RAN under, not the one it was asked for.

Ludmil, rc2 review section 2 (reproducibility blockers), first item:

    Record the resolved prefilter mode, not only the requested value. RUN_REPORT
    currently records the requested mode (for example, auto) rather than the
    resolved off/balanced state. Persist the resolved value in run provenance and
    sample-level evidence/QC.

`--prefilter_mode auto` is the default, and `auto` is not a mode anything runs in:
it becomes `balanced` for a metagenome with a --kraken_db and `off` otherwise. Up
to rc3 the only recorded value was the request, so RUN_REPORT.txt, the provenance
JSON and every per-sample QC row said `auto` -- and whether the pks counts beside
them came from the complete host-depleted stream or from a krakenuniq+DIAMOND
narrowed subset of it was not recoverable from any output. Two runs with
byte-identical recorded parameters could have profiled different read sets.

The fix has one moving part: params.prefilter_mode_resolved, resolved once at
script level in main.nf. From there it reaches

  * the workflow body, which builds the channel graph from it (so the recorded
    value cannot drift from the graph that ran);
  * runs/provenance/<session>.json, through the record's params dump, and
    RUN_REPORT.txt, through finalize_provenance.PARAMETER_GROUPS;
  * every per-sample QC fragment, through the two profiling lanes -- which are
    exactly what the prefilter narrows, and one of which always runs -- and from
    there the prefilter_mode column of cohort/cohort_qc_summary.tsv.

It cannot be resolved inside the workflow: Nextflow 26 ignores assignments to
params made there ("`params.x` is defined multiple times -- Assignments following
the first are ignored"), so processes would see null. That is asserted here as a
structural rule, because it is silent when it is got wrong.
"""
import csv
import importlib.util
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()
WORKFLOW = MAIN.split("workflow {", 1)[1]
SCRIPT_LEVEL = MAIN.split("workflow {", 1)[0]
MODULES = {p.name: p.read_text() for p in (ROOT / "Modules").glob("*.nf")}

spec = importlib.util.spec_from_file_location(
    "build_qc_summary", ROOT / "scripts/build_qc_summary.py")
build_qc_summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build_qc_summary)

sys.path.insert(0, str(ROOT / "scripts"))
import finalize_provenance  # noqa: E402

RESOLVED = "params.prefilter_mode_resolved"
# The profiling lanes: what the prefilter narrows, and the only stages whose counts
# depend on it. One of the two always runs (--profiling_method bowtie2 | hmm | both).
PROFILING_LANES = ("pksProfiler_align.nf", "pksProfiler_hmm.nf")


def resolution_block():
    """The script-level definition of params.prefilter_mode_resolved."""
    start = SCRIPT_LEVEL.index(f"{RESOLVED} =")
    return SCRIPT_LEVEL[start:].split("\n\n", 1)[0]


def resolve(prefilter_mode, sample_type, kraken_db):
    """Groovy's rule, re-stated here so the cases below read as cases."""
    if prefilter_mode != "auto":
        return prefilter_mode
    return "balanced" if (sample_type == "metagenome" and kraken_db) else "off"


class ResolvedOnceAtScriptLevel(unittest.TestCase):
    def test_it_is_defined_outside_the_workflow(self):
        # Inside the workflow Nextflow 26 drops the assignment and every task sees null.
        self.assertIn(f"{RESOLVED} =", SCRIPT_LEVEL)
        self.assertNotIn(f"{RESOLVED} =", WORKFLOW)

    def test_it_is_defined_after_the_parameters_it_reads(self):
        for name in ("params.prefilter_mode =", "params.sample_type =", "params.kraken_db"):
            with self.subTest(parameter=name):
                self.assertLess(SCRIPT_LEVEL.index(name), SCRIPT_LEVEL.index(f"{RESOLVED} ="),
                                f"{name} is declared after the resolution reads it")

    def test_the_rule_is_the_documented_one(self):
        block = resolution_block()
        self.assertIn('params.prefilter_mode.toString() != "auto"', block)
        self.assertIn('params.sample_type == "metagenome" && params.kraken_db', block)
        self.assertIn('"balanced"', block)
        self.assertIn('"off"', block)

    def test_an_explicit_mode_resolves_to_itself_so_validation_still_sees_a_typo(self):
        # The resolution must not rewrite an unknown mode to "off": validatePrefilterParams
        # is what reports `Unknown --prefilter_mode`, and it needs the original.
        self.assertIn('exit 1, "Unknown --prefilter_mode: ${params.prefilter_mode}', MAIN)
        self.assertEqual(resolve("bogus", "wgs", None), "bogus")

    def test_the_four_cases(self):
        self.assertEqual(resolve("auto", "metagenome", "/db/kraken"), "balanced")
        self.assertEqual(resolve("auto", "metagenome", None), "off")
        self.assertEqual(resolve("auto", "wgs", "/db/kraken"), "off")
        self.assertEqual(resolve("auto", "rna", None), "off")
        self.assertEqual(resolve("off", "metagenome", "/db/kraken"), "off")
        self.assertEqual(resolve("balanced", "metagenome", "/db/kraken"), "balanced")

    def test_the_workflow_builds_its_graph_from_the_resolved_value(self):
        self.assertIn(f"def prefilter_mode = {RESOLVED}.toString()", WORKFLOW)
        # ... and does not resolve a second time, which could drift from the record.
        self.assertNotIn('prefilter_mode = "balanced"', WORKFLOW)
        self.assertNotIn('prefilter_mode = "off"', WORKFLOW)
        self.assertIn('if (prefilter_mode == "balanced") {', WORKFLOW)

    def test_the_auto_resolution_is_still_logged_at_launch(self):
        self.assertIn('log.info "prefilter_mode auto -> ${prefilter_mode}', WORKFLOW)


class ItReachesTheRunProvenance(unittest.TestCase):
    def test_the_launch_record_dumps_every_param_including_this_one(self):
        # RunRecord.start's invocation.params is params.findAll{...}, so a script-level
        # param is in the JSON by construction. Assert the dump is still unfiltered
        # apart from dep_digest, which is the only deliberate exclusion.
        self.assertIn("params.findAll { key, _value -> key != 'dep_digest' }", MAIN)

    def test_run_report_prints_the_requested_and_the_resolved_mode(self):
        names = [n for title, members in finalize_provenance.PARAMETER_GROUPS
                 if title == "PROFILING" for n in members]
        self.assertIn("prefilter_mode", names)
        self.assertIn("prefilter_mode_resolved", names)
        self.assertEqual(names.index("prefilter_mode") + 1,
                         names.index("prefilter_mode_resolved"),
                         "the resolved value must be printed next to the request")

    def test_the_reported_name_is_a_real_pipeline_parameter(self):
        # The same guard tests/test_provenance_record.py applies to the whole group:
        # a name that is not a params.* is skipped exactly like an unset one, so a typo
        # would silently drop the resolved mode back out of the report (P0-1).
        declared = set(re.findall(r"params\.(\w+)\s*=", MAIN))
        self.assertIn("prefilter_mode_resolved", declared)

    def test_a_report_built_from_a_record_shows_both(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            runs, cohort = root / "runs", root / "cohort"
            (runs / "provenance").mkdir(parents=True)
            (cohort / "qc").mkdir(parents=True)
            (cohort / "qc" / "pks.qc.summary.tsv").write_text("Sample\tstatus\nS1\tcomplete\n")
            record = {
                "record_version": 1, "status": "started",
                "code": {"commit": "abc", "dirty": False, "modified_files": []},
                "invocation": {"command_line": "nextflow run main.nf",
                               "params": {"sample_type": "metagenome",
                                          "profiling_method": "bowtie2",
                                          "prefilter_mode": "auto",
                                          "prefilter_mode_resolved": "balanced"}},
                "dependencies": {}, "environment": {},
                "run": {"session_id": "S", "started": "2026-10-22T00:00:00+00:00"},
                "outputs": {}, "notes": [],
            }
            (runs / "provenance" / "S.json").write_text(json.dumps(record) + "\n")
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/finalize_provenance.py"),
                 "--outdir", str(root), "--runs-dir", str(runs), "--cohort-dir", str(cohort),
                 "--session", "S", "--success", "true", "--exit-status", "0",
                 "--duration", "1m", "--error-message", ""],
                capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = (root / "RUN_REPORT.txt").read_text()
        self.assertRegex(report, r"prefilter_mode\s+auto")
        self.assertRegex(report, r"prefilter_mode_resolved\s+balanced")


class ItReachesEverySampleQcRow(unittest.TestCase):
    def test_both_profiling_lanes_emit_the_metric_from_the_resolved_param(self):
        for name in PROFILING_LANES:
            with self.subTest(module=name):
                text = MODULES[name]
                self.assertIn("prefilter_mode", text)
                self.assertIn(f'"${{{RESOLVED}}}"', text,
                              f"{name} must write the resolved param, not params.prefilter_mode")

    def test_every_qc_fragment_the_lane_writes_carries_it(self):
        """The HMM lane has a zero-read early-exit branch with its own fragment."""
        for name in PROFILING_LANES:
            text = MODULES[name]
            headers = len(re.findall(r'printf "Sample\\*tMetric\\*tValue', text))
            emits = len(re.findall(r'printf "%s\\*tprefilter_mode\\*t%s', text))
            with self.subTest(module=name):
                self.assertGreater(headers, 0)
                self.assertEqual(emits, headers,
                                 f"{name} writes {headers} QC fragment(s) but emits the "
                                 f"resolved prefilter mode in {emits} of them")

    def test_no_lane_emits_the_requested_value_by_mistake(self):
        for name in PROFILING_LANES:
            with self.subTest(module=name):
                self.assertNotIn("${params.prefilter_mode}", MODULES[name])

    def test_it_is_a_column_of_the_cohort_qc_table(self):
        self.assertIn("prefilter_mode", build_qc_summary.OUTPUT_COLUMNS)

    def test_only_resolved_values_are_accepted(self):
        self.assertEqual(build_qc_summary.TEXT_METRICS["prefilter_mode"], {"off", "balanced"})

    def summarise(self, metrics, sample="S1"):
        with tempfile.TemporaryDirectory() as directory:
            fragment = Path(directory) / f"{sample}.qc.tsv"
            lines = ["Sample\tMetric\tValue"]
            lines += [f"{sample}\t{metric}\t{value}" for metric, value in metrics]
            fragment.write_text("\n".join(lines) + "\n")
            output = Path(directory) / "summary.tsv"
            build_qc_summary.write_summary(
                build_qc_summary.load_fragments([fragment]), output, [sample])
            with output.open(newline="") as handle:
                return list(csv.DictReader(handle, delimiter="\t"))[0]

    BASE = [("filter_input_reads", 400), ("reads_after_fastp", 360),
            ("reads_after_hg38", 300), ("reads_after_t2t_phix", 250)]

    def test_the_resolved_mode_reaches_the_row(self):
        for mode in ("off", "balanced"):
            with self.subTest(mode=mode):
                row = self.summarise(self.BASE + [("prefilter_mode", mode)])
                self.assertEqual(row["prefilter_mode"], mode)

    def test_a_sample_whose_profiling_never_ran_says_na_rather_than_guessing(self):
        row = self.summarise(self.BASE)
        self.assertEqual(row["prefilter_mode"], "NA")

    def test_the_unresolved_request_is_refused(self):
        # If params.prefilter_mode_resolved ever stopped resolving, the fragment would
        # carry `auto` -- a value no run can be reproduced from. Fail, do not record it.
        with self.assertRaises(ValueError):
            self.summarise(self.BASE + [("prefilter_mode", "auto")])

    def test_the_two_lanes_may_both_write_it_as_long_as_they_agree(self):
        """--profiling_method both: two fragments, one sample, same resolved mode."""
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for lane, mode in (("align", "balanced"), ("hmm", "balanced")):
                path = Path(directory) / f"S1.{lane}.qc.tsv"
                path.write_text("Sample\tMetric\tValue\nS1\tprefilter_mode\t%s\n" % mode)
                paths.append(path)
            metrics = build_qc_summary.load_fragments(paths)
            self.assertEqual(metrics["S1"]["prefilter_mode"], "balanced")

    def test_two_lanes_disagreeing_is_an_error_not_a_coin_toss(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for lane, mode in (("align", "balanced"), ("hmm", "off")):
                path = Path(directory) / f"S1.{lane}.qc.tsv"
                path.write_text("Sample\tMetric\tValue\nS1\tprefilter_mode\t%s\n" % mode)
                paths.append(path)
            with self.assertRaises(ValueError):
                build_qc_summary.load_fragments(paths)


class ItIsDocumented(unittest.TestCase):
    def test_the_prefilter_doc_says_where_the_resolved_value_is_recorded(self):
        doc = (ROOT / "docs/prefilter.md").read_text()
        self.assertIn("prefilter_mode_resolved", doc)
        self.assertIn("cohort_qc_summary.tsv", doc)


if __name__ == "__main__":
    unittest.main()


class CohortExtractionModeUniformity(unittest.TestCase):
    """Ludmil's v3 review, Tier A item 2.

    extract_low_score_mode is written per sample and nothing compared it across the
    cohort, so a run could hold rescued and FLAG4-only samples side by side in silence.
    Those samples saw different read pools; a prevalence difference between them can be
    extraction rather than biology. C1 found this live in TCGA WES: 1,287 rescued and 17
    not, because the 17 are bwa sampe.
    """

    SCRIPT = ROOT / "scripts" / "build_qc_summary.py"

    def _qc(self, tmp, modes):
        path = Path(tmp) / "pks.qc.summary.tsv"
        with path.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "extract_low_score_mode", "terminal_state"])
            for i, m in enumerate(modes):
                w.writerow([f"S{i:03d}", m, "completed"])
        return path

    def _check(self, modes, required):
        import importlib.util
        spec = importlib.util.spec_from_file_location("bqs", self.SCRIPT)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with tempfile.TemporaryDirectory() as tmp:
            path = self._qc(tmp, modes)
            return mod.check_uniform_low_score_mode(path, required)

    def test_a_uniform_cohort_passes_quietly(self):
        self.assertIsNone(self._check(["on_auto_bwa_mem"] * 40, required=True))

    def test_auto_and_requested_on_are_not_a_mix(self):
        # Both rescued, so equally sensitive; which switch said so does not matter.
        self.assertIsNone(self._check(["on_auto_bwa_mem"] * 30 + ["on_requested"] * 10,
                                      required=True))

    def test_a_mixed_cohort_is_refused_when_required(self):
        with self.assertRaises(SystemExit) as caught:
            self._check(["on_auto_bwa_mem"] * 1287 + ["off_auto_no_as_tag"] * 17,
                        required=True)
        message = str(caught.exception)
        self.assertIn("1287", message)
        self.assertIn("17", message)
        self.assertIn("--extract_low_score_mapped", message,
                      "the error must say how to resolve it")

    def test_a_mixed_cohort_only_warns_by_default(self):
        # A mix is legitimate when the inputs genuinely differ, so the default carries on.
        self.assertIsNone(self._check(["on_auto_bwa_mem"] * 1287 + ["off_auto_no_as_tag"] * 17,
                                      required=False))

    def test_a_cohort_with_no_extraction_column_is_not_an_error(self):
        # FASTQ-only cohorts never extract, so there is nothing to compare.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pks.qc.summary.tsv"
            with path.open("w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["sample", "terminal_state"])
                w.writerow(["S000", "completed"])
            import importlib.util
            spec = importlib.util.spec_from_file_location("bqs", self.SCRIPT)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            self.assertIsNone(mod.check_uniform_low_score_mode(path, True))

    def test_the_workflow_declares_the_parameter_and_passes_it_through(self):
        main_nf = (ROOT / "main.nf").read_text()
        self.assertIn("params.require_low_score_mode", main_nf)
        plotting = (ROOT / "Modules" / "plotting.nf").read_text()
        self.assertIn("--require-uniform-low-score-mode", plotting)

    def test_the_run_report_warns_on_a_mixed_cohort(self):
        prov = (ROOT / "scripts" / "finalize_provenance.py").read_text()
        self.assertIn("MIXED EXTRACTION MODE IN ONE COHORT", prov)
        self.assertIn("extraction_modes", prov)
