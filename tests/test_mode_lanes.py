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
        self.assertIn('params.sample_type in ["wes", "rna"] && params.profiling_method == "hmm"', CODE)

    def test_wes_writes_its_own_evidence_kind(self):
        helper = MAIN[MAIN.index("def evidenceKind"):]
        self.assertIn("wes: 'wes_evidence'", helper[:helper.index("\n}\n")])

    def test_the_default_on_assembly_flag_is_inert_not_refused_for_wes(self):
        self.assertIn("--tumor_targeted_assembly does not apply to --sample_type", MAIN)

    def test_wes_qc_todo_is_documented(self):
        doc = (ROOT / "docs/running/wes_rna.md").read_text()
        self.assertIn("TODO", doc)
        self.assertIn("duplicate", doc)


MAP = (ROOT / "Modules/map_reads.nf").read_text()
PLAN = (ROOT / "lib/HostDepletion.groovy").read_text()


def groovy_map(name):
    m = re.search(rf"{name}\s*=\s*\[(.*?)\]", PLAN, re.S)
    return dict(re.findall(r"(\w+)\s*:\s*'([^']*)'", m.group(1)))


def load_script(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RnaLane(unittest.TestCase):
    """Section 4: splice-aware human depletion, profiling as WGS, rna_evidence.tsv."""

    def test_rna_is_accepted(self):
        self.assertIn("rna", accepted())

    def test_only_rna_s_human_pass_is_spliced(self):
        presets = groovy_map("HUMAN_PRESET")
        self.assertEqual(presets, {"wgs": "sr", "wes": "sr", "metagenome": "sr",
                                   "rna": "splice:sr"})
        self.assertEqual(set(presets), set(accepted()))

    def test_t2t_and_pangenome_passes_stay_sr(self):
        self.assertIn("static final String RESIDUAL_PRESET = 'sr'", PLAN)
        for plan in ("t2tPass", "pangenomePass"):
            body = PLAN[PLAN.index(f"static Map {plan}("):]
            self.assertIn("preset: RESIDUAL_PRESET", body[:body.index("\n    }\n")])

    def test_the_human_pass_uses_the_splice_index_and_junctions_for_rna_only(self):
        body = PLAN[PLAN.index("static Map humanPass("):]
        body = body[:body.index("\n    }\n")]
        self.assertIn("if (rna && rnaJunctionsBed) args << \"-j", body)
        self.assertIn("index : (rna ? rnaSpliceDb : hg38Db)", body)

    def test_map_reads_builds_every_call_from_the_plan(self):
        self.assertIn("HostDepletion.humanPass(params.sample_type, params.minimap2_frag_mode,\n"
                      "                                        params.hg38_db, params.rna_splice_db,\n"
                      "                                        params.rna_junctions_bed)", MAP)
        self.assertNotRegex(MAP, r"minimap2 -2 -ax sr\b", "a hard-coded -ax sr call bypasses the plan")
        self.assertIn('"${human.index}"', MAP)

    def test_qc_metric_names_are_unchanged(self):
        for metric in ("reads_after_hg38", "reads_after_t2t_phix", "reads_after_pangenome"):
            self.assertIn(f"\\t{metric}\\t", MAP)

    def test_hg38_db_is_not_required_for_rna_and_the_splice_db_is(self):
        self.assertIn("static String humanIndexParam(Object sampleType)", PLAN)
        self.assertIn("sampleType?.toString() == 'rna' ? 'rna_splice_db' : 'hg38_db'", PLAN)
        self.assertIn("def human_index_param = HostDepletion.humanIndexParam(params.sample_type)", CODE)
        self.assertIn("if (!params[human_index_param])", CODE)
        self.assertNotIn("if (!params.hg38_db)", CODE)
        pre = load_script("preflight")
        self.assertEqual(pre.human_index_param({"sample_type": "rna"}), "rna_splice_db")
        for mode in ("wgs", "wes", "metagenome"):
            self.assertEqual(pre.human_index_param({"sample_type": mode}), "hg38_db")

    def test_main_nf_hands_preflight_the_rna_references(self):
        block = MAIN[MAIN.index("def preflight_config = ["):]
        block = block[:block.index("]\n")]
        for key in ("rna_splice_db", "rna_junctions_bed", "sample_type"):
            self.assertIn(f"{key}", block)

    def test_library_type_is_validated_and_defaults_to_unknown(self):
        self.assertIn('params.rna_library_type  = "unknown"', MAIN)
        self.assertIn('params.rna_library_type.toString() in ["polyA", "total", "unknown"]', CODE)

    def test_prefilter_auto_resolves_to_off_for_rna(self):
        # auto -> balanced only for metagenomes (with a kraken db).
        resolve = CODE[CODE.index('if (prefilter_mode == "auto") {'):]
        resolve = resolve[:resolve.index('log.info "prefilter_mode auto')]
        self.assertIn('params.sample_type == "metagenome" && params.kraken_db', resolve)
        self.assertIn('prefilter_mode = "off"', resolve)

    def test_rna_writes_its_own_evidence_kind_with_the_library_type(self):
        helper = MAIN[MAIN.index("def evidenceKind"):]
        self.assertIn("rna: 'rna_evidence'", helper[:helper.index("\n}\n")])
        targeted = (ROOT / "Modules/pks_targeted.nf").read_text()
        self.assertIn("--rna-library-type ${params.rna_library_type}", targeted)

    def test_run_report_records_the_rna_parameters(self):
        prov = load_script("finalize_provenance")
        names = {n for _group, members in prov.PARAMETER_GROUPS for n in members}
        for name in ("sample_type", "rna_splice_db", "rna_junctions_bed", "rna_library_type"):
            self.assertIn(name, names)

    def test_rna_references_invalidate_the_depletion_task(self):
        self.assertRegex(MAIN, r"rna_splice_db\s*:\s*Provenance\.data\(params\.rna_splice_db\)")
        self.assertRegex(MAIN, r"rna_junctions_bed\s*:\s*Provenance\.code\(params\.rna_junctions_bed\)")


if __name__ == "__main__":
    unittest.main()
