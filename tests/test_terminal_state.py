"""Every expected sample ends in exactly one terminal state (GAMEPLAN v1.0.0-rc1).

"Extend QC so every sample ends as exactly one of: completed positive / completed
negative / failed at a named stage / excluded for a named reason." The cohort QC summary
(build_qc_summary.py) carries it as `terminal_state`, and the master summary copies it:

    completed_positive   tiered short-read wgs, a formal tier (multi_gene and up)
    completed_negative   tiered short-read wgs, 0 reads or localized_indeterminate
    completed_untiered   finished, but the evidence carries no call (wes, rna,
                         metagenome, long-read wgs, or no evidence step at all)
    failed_at_<stage>    the first expected stage, in pipeline order, with no QC
    excluded_<reason>    named in --excluded-samples
"""
import csv
import importlib.util
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/build_qc_summary.py"
MAIN = (ROOT / "main.nf").read_text()
STATE = re.compile(r"^(completed_(positive|negative|untiered)|failed_at_[a-z_]+|excluded_[a-z0-9_]+)$")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


qc = load("build_qc_summary", SCRIPT)
classify = load("classify_tumor_pks_evidence", ROOT / "scripts/classify_tumor_pks_evidence.py")

ALIGNMENT_WGS = ["extraction", "read_qc", "host_depletion", "alignment_profiling", "evidence"]
COMPLETE = {"extracted_unmapped_reads": 10, "filter_input_reads": 10, "reads_after_fastp": 9,
            "reads_after_hg38": 8, "reads_after_t2t_phix": 7,
            "reads_clb_genes_align": 3, "num_clb_genes_align": 2}


class TerminalStateRule(unittest.TestCase):
    def state(self, values, stages=ALIGNMENT_WGS, excluded=None):
        return qc.terminal_state(values, stages, excluded)

    def test_completed_states_come_from_the_evidence_outcome(self):
        for outcome, expected in (("positive", "completed_positive"),
                                  ("negative", "completed_negative"),
                                  ("untiered", "completed_untiered")):
            with self.subTest(outcome=outcome):
                self.assertEqual(self.state({**COMPLETE, "evidence_outcome": outcome}), expected)

    def test_failed_stage_is_the_first_missing_in_pipeline_order(self):
        values = {**COMPLETE, "evidence_outcome": "positive"}
        for stage, metrics in (("extraction", ["extracted_unmapped_reads"]),
                               ("read_qc", ["reads_after_fastp"]),
                               ("host_depletion", ["reads_after_t2t_phix"]),
                               ("alignment_profiling", ["reads_clb_genes_align"]),
                               ("evidence", ["evidence_outcome"])):
            with self.subTest(stage=stage):
                broken = {k: v for k, v in values.items() if k not in metrics}
                self.assertEqual(self.state(broken), f"failed_at_{stage}")
        # a failed depletion leaves everything after it missing too: named once, upstream
        upstream = {k: v for k, v in COMPLETE.items() if k in ("extracted_unmapped_reads",
                                                                "filter_input_reads",
                                                                "reads_after_fastp")}
        self.assertEqual(self.state(upstream), "failed_at_host_depletion")

    def test_a_sample_with_no_qc_failed_at_the_first_stage(self):
        self.assertEqual(self.state({}), "failed_at_extraction")
        self.assertEqual(self.state({}, stages=["read_qc", "host_depletion"]), "failed_at_read_qc")

    def test_exclusion_wins(self):
        self.assertEqual(self.state({**COMPLETE, "evidence_outcome": "positive"},
                                    excluded="consent_withdrawn"), "excluded_consent_withdrawn")

    def test_no_evidence_step_is_completed_untiered(self):
        hmm_only = ["read_qc", "host_depletion", "hmm_profiling"]
        values = {"filter_input_reads": 5, "reads_after_fastp": 5, "reads_after_hg38": 5,
                  "reads_after_t2t_phix": 5, "reads_clb_genes_hmm": 1, "num_clb_genes_hmm": 1}
        self.assertEqual(self.state(values, hmm_only), "completed_untiered")

    def test_stage_names_are_validated_and_ordered(self):
        self.assertEqual(qc.parse_stages("evidence,read_qc"), ["read_qc", "evidence"])
        with self.assertRaises(ValueError):
            qc.parse_stages("read_qc,bogus")


class EvidenceOutcome(unittest.TestCase):
    """classify_tumor_pks_evidence.py decides positive/negative/untiered, once."""

    def outcome(self, mode, reads, genes, breadth, no_tiers=False):
        thresholds = {"multi_gene": (5, 3, .01), "broad_island": (30, 8, .075),
                      "extensive_island": (100, 10, .15)}
        row = classify.evidence_row(mode, "S", reads, genes, (breadth, 0, 0), thresholds,
                                    no_tiers=no_tiers)
        return classify.evidence_outcome(row)

    def test_formal_tiers_are_positive_and_everything_below_negative(self):
        self.assertEqual(self.outcome("wgs", 10, 4, .02), "positive")
        self.assertEqual(self.outcome("wgs", 200, 12, .5), "positive")
        self.assertEqual(self.outcome("wgs", 2, 1, .001), "negative")   # localized_indeterminate
        self.assertEqual(self.outcome("wgs", 0, 0, 0), "negative")

    def test_untiered_modes(self):
        for mode in ("metagenome", "wes", "rna"):
            with self.subTest(mode=mode):
                self.assertEqual(self.outcome(mode, 500, 19, .9), "untiered")
        self.assertEqual(self.outcome("wgs", 500, 19, .9, no_tiers=True), "untiered")


class EndToEnd(unittest.TestCase):
    def test_every_expected_sample_gets_exactly_one_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            fragments = []

            def frag(name, rows):
                path = d / name
                path.write_text("Sample\tMetric\tValue\n" + "".join(f"{s}\t{m}\t{v}\n" for s, m, v in rows))
                fragments.append(path)

            full = ["POS", "NEG", "NOEVID"]
            frag("extract.tsv", [(s, "extracted_unmapped_reads", 10) for s in full + ["DEPFAIL"]])
            frag("filter.tsv", [(s, m, v) for s in full + ["DEPFAIL"]
                                for m, v in (("filter_input_reads", 10), ("reads_after_fastp", 9))])
            frag("deplete.tsv", [(s, m, v) for s in full
                                 for m, v in (("reads_after_hg38", 8), ("reads_after_t2t_phix", 7))])
            frag("align.tsv", [(s, m, v) for s in full
                               for m, v in (("reads_clb_genes_align", 3), ("num_clb_genes_align", 2))])
            frag("evidence.tsv", [("POS", "evidence_outcome", "positive"),
                                  ("NEG", "evidence_outcome", "negative")])
            (d / "list.txt").write_text("".join(f"{p}\n" for p in fragments))
            (d / "expected.txt").write_text("POS\nNEG\nNOEVID\nDEPFAIL\nNEVER\nGONE\n")
            (d / "excluded.tsv").write_text("GONE\tfailed_identity_check\n")
            out = d / "qc.tsv"
            result = subprocess.run([sys.executable, str(SCRIPT), "--inputs-file", str(d / "list.txt"),
                                     "--expected-samples", str(d / "expected.txt"),
                                     "--expected-stages", ",".join(ALIGNMENT_WGS),
                                     "--excluded-samples", str(d / "excluded.tsv"),
                                     "--output", str(out)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with out.open() as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            states = {r["Sample"]: r["terminal_state"] for r in rows}
            self.assertEqual(len(rows), len(states), "one row per sample")
            self.assertEqual(states, {
                "POS": "completed_positive", "NEG": "completed_negative",
                "NOEVID": "failed_at_evidence", "DEPFAIL": "failed_at_host_depletion",
                "NEVER": "failed_at_extraction", "GONE": "excluded_failed_identity_check"})
            for state in states.values():
                self.assertRegex(state, STATE)
            self.assertEqual(rows[0].keys().__iter__().__next__(), "Sample")
            self.assertEqual(list(rows[0])[-1], "terminal_state")

    def test_a_malformed_exclusion_list_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.tsv"
            path.write_text("S1\tNot A Reason\n")
            with self.assertRaises(ValueError):
                qc.read_excluded(path)


class Wiring(unittest.TestCase):
    def test_classify_emits_the_outcome_fragment_and_it_reaches_the_qc_summary(self):
        targeted = (ROOT / "Modules/pks_targeted.nf").read_text()
        self.assertIn('--qc-fragment "${sampleID}.evidence.qc.tsv"', targeted)
        self.assertIn('path("${sampleID}.evidence.qc.tsv"), emit: qc', targeted)
        self.assertIn('pattern: "*_evidence.tsv"', targeted)
        self.assertIn("classifyPksReadEvidence.out.qc.map { _sampleID, qc_file -> qc_file }", MAIN)

    def test_main_passes_the_expected_stages_in_order(self):
        body = MAIN[MAIN.index("def qc_stages"):MAIN.index("masterQCSummary(")]
        order = [s for s in qc.STAGE_ORDER if f"'{s}'" in body]
        self.assertEqual(order, list(qc.STAGE_ORDER), "main.nf names every stage")
        positions = [body.index(f"'{s}'") for s in qc.STAGE_ORDER]
        self.assertEqual(positions, sorted(positions), "main.nf lists them in pipeline order")
        call = MAIN[MAIN.index("masterQCSummary("):]
        self.assertIn("qc_stages.join(',')", call[:call.index("\n    )")])
        plotting = (ROOT / "Modules/plotting.nf").read_text()
        self.assertIn('--expected-stages "${expected_stages}"', plotting)

    def test_master_summary_carries_it_after_the_mode(self):
        source = (ROOT / "scripts/build_master_summary.py").read_text()
        self.assertIn('TERMINAL_COLUMNS = ["terminal_state"]', source)
        main_block = source[source.index("columns += MODE_COLUMNS"):]
        self.assertLess(main_block.index("TERMINAL_COLUMNS"), main_block.index("for name, cols, fn in stages"))

    def test_documented(self):
        docs = (ROOT / "docs/output.md").read_text()
        for word in ("terminal_state", "completed_positive", "completed_negative",
                     "completed_untiered", "failed_at_", "excluded_"):
            self.assertIn(word, docs)


if __name__ == "__main__":
    unittest.main()
