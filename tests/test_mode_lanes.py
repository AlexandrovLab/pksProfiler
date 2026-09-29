"""Each sample type runs a real, deliberate lane -- and only its own features.

Ludmil, 26 Sep mode memo section 6: "Each of wgs, wes, rna and metagenome must map to a
real execution branch and dedicated contract test." Section 7: "no accidental WGS-tier
reuse". These are the static half of that contract (the channel graph for each mode is
built for real by tests/check_workflow_construction.sh): which branches a mode enters, and
which WGS-only features -- tiers, targeted assembly, MAGs -- it is kept out of.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()
BODY = MAIN.split("workflow {", 1)[1]


def code_only(text):
    return "\n".join(line for line in text.splitlines()
                      if not line.strip().startswith(("#", "//")))


CODE = code_only(BODY)


def accepted():
    m = re.search(r"def valid_sample_types = (\[[^\]]*\])", CODE)
    return eval(m.group(1))


def block_after(marker):
    """The brace-balanced block that opens on the line containing `marker`."""
    start = CODE.index(marker)
    open_at = CODE.index("{", start)
    depth = 0
    for i in range(open_at, len(CODE)):
        depth += {"{": 1, "}": -1}.get(CODE[i], 0)
        if depth == 0:
            return CODE[open_at:i + 1]
    raise AssertionError(f"unbalanced block after {marker}")


EVIDENCE_GATE = re.search(r'if \(params\.sample_type in (\[[^\]]*\])\) \{\s*read_evidence_input_ch',
                          CODE)


class EveryModeIsRouted(unittest.TestCase):
    def test_the_evidence_branch_takes_every_accepted_mode(self):
        self.assertIsNotNone(EVIDENCE_GATE, "evidence branch not found")
        self.assertEqual(set(eval(EVIDENCE_GATE.group(1))), set(accepted()))

    def test_the_evidence_branch_runs_the_classifier(self):
        start = EVIDENCE_GATE.start()
        self.assertIn("classifyPksReadEvidence(read_evidence_input_ch)",
                      CODE[start:start + 800])


class OnlyWgsIsTieredAndAssembled(unittest.TestCase):
    def test_targeted_assembly_is_inside_the_wgs_branch_only(self):
        wgs = block_after('if (params.sample_type == "wgs") {')
        self.assertIn("targetedPksAssembly(", wgs)
        self.assertEqual(CODE.count("targetedPksAssembly("), 1)
        self.assertIn("tier in contig_tiers", wgs)
        self.assertEqual(CODE.count("tier in contig_tiers"), 1)

    def test_mags_are_metagenome_only(self):
        self.assertIn('params.sample_type == "metagenome" && enable_mags_b', CODE)
        self.assertIn('enable_mags_b && params.sample_type != "metagenome"', CODE)


class WesLane(unittest.TestCase):
    """Section 3: standard depletion, profiling, wes_evidence.tsv; no tier, no assembly."""

    def test_wes_is_accepted(self):
        self.assertIn("wes", accepted())

    def test_wes_is_refused_hmm_only_profiling(self):
        self.assertIn('params.sample_type == "wes" && params.profiling_method == "hmm"', CODE)

    def test_wes_writes_its_own_evidence_kind(self):
        helper = MAIN[MAIN.index("def evidenceKind"):]
        self.assertIn("wes: 'wes_evidence'", helper[:helper.index("\n}\n")])

    def test_the_default_on_assembly_flag_is_inert_not_refused_for_wes(self):
        self.assertIn("--tumor_targeted_assembly does not apply to --sample_type", MAIN)

    def test_wes_qc_todo_is_documented(self):
        doc = (ROOT / "docs/running/wes_rna.md").read_text()
        self.assertIn("TODO", doc)
        self.assertIn("duplicate", doc)


if __name__ == "__main__":
    unittest.main()
