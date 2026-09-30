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

    def test_wes_and_rna_qc_are_documented(self):
        # The memo section 3/4 QC TODOs are implemented (Modules/mode_qc.nf); the doc now
        # defines each column instead of listing them as open.
        doc = (ROOT / "docs/running/wes_rna.md").read_text()
        self.assertNotIn("TODO", doc)
        for column in ("wes_input_duplicate_fraction", "wes_capture_off_target_fraction",
                       "capture_bed_not_provided", "wes_pks_island_duplicate_fraction",
                       "clb_transcripts_detected", "rna_residual_human_fraction"):
            self.assertIn(column, doc)


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
        self.assertEqual(groovy_map("DNA_PRESET")["short"], "sr")
        self.assertEqual(groovy_map("RNA_PRESET")["short"], "splice:sr")
        body = PLAN[PLAN.index("static Map humanPass("):]
        body = body[:body.index("\n    }\n")]
        self.assertIn("preset: (rna ? RNA_PRESET : DNA_PRESET)[rt]", body)
        self.assertIn("boolean rna = sampleType?.toString() == 'rna'", body)

    def test_t2t_and_pangenome_passes_are_never_spliced(self):
        for plan in ("t2tPass", "pangenomePass"):
            body = PLAN[PLAN.index(f"static Map {plan}("):]
            body = body[:body.index("\n    }\n")]
            with self.subTest(plan=plan):
                self.assertIn("preset: DNA_PRESET[rt]", body)
                self.assertNotIn("RNA_PRESET", body)

    def test_the_human_pass_uses_the_splice_index_and_junctions_for_rna_only(self):
        body = PLAN[PLAN.index("static Map humanPass("):]
        body = body[:body.index("\n    }\n")]
        self.assertIn("if (rna && rnaJunctionsBed) args << \"-j", body)
        self.assertIn("rna ? (lr ? rnaSpliceDbLong : rnaSpliceDb) : (lr ? hg38DbLong : hg38Db)", body)

    def test_map_reads_builds_every_call_from_the_plan(self):
        self.assertIn("HostDepletion.humanPass(params.sample_type, params.read_type, params.minimap2_frag_mode,\n"
                      "                                        params.hg38_db, params.hg38_db_long,\n"
                      "                                        params.rna_splice_db, params.rna_splice_db_long,\n"
                      "                                        params.rna_junctions_bed)", MAP)
        self.assertNotRegex(MAP, r"minimap2 -2 -ax sr\b", "a hard-coded -ax sr call bypasses the plan")
        self.assertIn('"${human.index}"', MAP)

    def test_qc_metric_names_are_unchanged(self):
        for metric in ("reads_after_hg38", "reads_after_t2t_phix", "reads_after_pangenome"):
            self.assertIn(f"\\t{metric}\\t", MAP)

    def test_hg38_db_is_not_required_for_rna_and_the_splice_db_is(self):
        self.assertIn("static String humanIndexParam(Object sampleType, Object readType)", PLAN)
        self.assertIn("sampleType?.toString() == 'rna' ? 'rna_splice_db' : 'hg38_db'", PLAN)
        self.assertIn("def human_index_param = HostDepletion.humanIndexParam(params.sample_type, params.read_type)", CODE)
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


class LongReadLane(unittest.TestCase):
    """Owner decision, 29 Sep: --read_type short|ont|hifi; long reads for wgs and rna only.

    The preset table itself is executed by tests/check_depletion_plan.sh; this pins the
    source and the wiring around it.
    """

    FILTER = (ROOT / "Modules/filter_reads.nf").read_text()
    ALIGN = (ROOT / "Modules/pksProfiler_align.nf").read_text()
    TARGETED = (ROOT / "Modules/pks_targeted.nf").read_text()

    def test_read_types_and_the_modes_that_take_long_reads(self):
        self.assertIn("READ_TYPES = ['short', 'ont', 'hifi']", PLAN)
        self.assertIn("LONG_READ_SAMPLE_TYPES = ['wgs', 'rna']", PLAN)
        self.assertIn('params.read_type             = "short"', MAIN)

    def test_preset_per_read_type(self):
        self.assertEqual(groovy_map("DNA_PRESET"), {"short": "sr", "ont": "map-ont", "hifi": "map-hifi"})
        self.assertEqual(groovy_map("RNA_PRESET"),
                         {"short": "splice:sr", "ont": "splice", "hifi": "splice:hq"})

    def test_read_type_is_validated_and_refused_outside_wgs_and_rna(self):
        self.assertIn("if (!(read_type in HostDepletion.READ_TYPES))", CODE)
        self.assertIn("long-read input is supported for wgs and rna only", MAIN)
        body = MAIN.split("workflow {", 1)[1]
        self.assertLess(body.index("long-read input is supported"), body.index("extractReads("))

    # rc1 pass, item 1: switched off (the flag defaults on), refused only when the
    # command line asked for it; resolved before preflight and provenance see it.
    def test_long_reads_switch_targeted_assembly_off_unless_explicitly_requested(self):
        body = MAIN.split("workflow {", 1)[1]
        block = body[body.index("def long_read_input = HostDepletion.isLong(params.read_type)"):]
        block = block[:block.index("tumor_targeted_assembly_b = false") + 40]
        self.assertIn("Preflight.explicitlyTrue(workflow.commandLine, 'tumor_targeted_assembly')", block)
        self.assertIn("if (long_read_input && tumor_targeted_assembly_b)", block)
        self.assertIn("exit 1", block)
        self.assertIn("log.info", block)
        self.assertIn("tumor_targeted_assembly_b = false", block)
        self.assertLess(body.index("tumor_targeted_assembly_b = false"),
                        body.index("tumor_targeted_assembly  : tumor_targeted_assembly_b"))
        self.assertNotIn("Pass --tumor_targeted_assembly false", MAIN)
        helper = (ROOT / "lib/Preflight.groovy").read_text()
        self.assertIn("static boolean explicitlyTrue(Object commandLine, String name)", helper)

    def test_long_reads_refuse_bracken_taxonomy(self):
        self.assertIn("if (long_reads && (pks_taxa_b || pks_community_taxa_b))", CODE)
        self.assertIn("Bracken's read-length model (--bracken_read_length) is a short-read concept", MAIN)

    def test_long_read_wgs_is_untiered(self):
        self.assertIn("params.sample_type.toString() == 'wgs' && HostDepletion.isLong(params.read_type) ? ' --no-tiers' : ''",
                      self.TARGETED)

    def test_long_reads_need_the_long_indexes(self):
        self.assertIn('return isLong(readType) ? "${base}_long".toString() : base', PLAN)
        self.assertIn("return isLong(readType) ? 't2t_phix_db_long' : 't2t_phix_db'", PLAN)
        self.assertIn("def t2t_index_param = HostDepletion.t2tIndexParam(params.read_type)", CODE)
        self.assertNotIn("if (!params.t2t_phix_db)", CODE)
        pre = load_script("preflight")
        cases = {("wgs", "short"): ("hg38_db", "t2t_phix_db"),
                 ("wgs", "ont"): ("hg38_db_long", "t2t_phix_db_long"),
                 ("rna", "short"): ("rna_splice_db", "t2t_phix_db"),
                 ("rna", "hifi"): ("rna_splice_db_long", "t2t_phix_db_long")}
        for (mode, rt), expected in cases.items():
            config = {"sample_type": mode, "read_type": rt}
            with self.subTest(mode=mode, read_type=rt):
                self.assertEqual((pre.human_index_param(config), pre.t2t_index_param(config)), expected)

    def test_frag_is_never_passed_to_long_reads(self):
        body = PLAN[PLAN.index("static String fragArg("):]
        self.assertIn("!isLong(readType)", body[:body.index("\n    }\n")])

    def test_long_reads_are_filtered_with_fastplong_keeping_the_metric_names(self):
        self.assertIn("HostDepletion.isLong(params.read_type) ? params.fastplong_env : params.fastp_env", self.FILTER)
        self.assertIn("-l ${params.long_read_min_length}", self.FILTER)
        for metric in ("filter_input_reads", "reads_after_fastp"):
            self.assertIn(f"\\t{metric}\\t", self.FILTER)
        env = (ROOT / "conda_envs/fastplong_env.yml").read_text()
        self.assertIn("fastplong=", env)
        self.assertIn("params.long_read_min_length  = 200", MAIN)

    def test_long_reads_are_profiled_with_minimap2_not_bowtie2(self):
        self.assertIn("minimap2 -ax ${long_preset} --secondary=no", self.ALIGN)
        self.assertIn('"${params.pks_reference_fasta}" "${reads}"', self.ALIGN)
        self.assertIn("-q ${params.pks_align_min_mapq}", self.ALIGN.split("minimap2 -ax ${long_preset}", 1)[1][:400])
        self.assertIn('def featurecounts_long = long_read ? "-L" : ""', self.ALIGN)
        self.assertIn("--largestOverlap ${featurecounts_long}", self.ALIGN)
        self.assertIn("minimap2=2.31", (ROOT / "conda_envs/pks_align_env.yml").read_text())
        self.assertIn("static String profilingPreset(Object readType)", PLAN)

    def test_long_reads_refuse_a_mate_file(self):
        self.assertIn("if (long_reads && row.fastq2?.toString()?.trim())", CODE)

    def test_run_report_records_read_type_and_the_long_indexes(self):
        prov = load_script("finalize_provenance")
        names = {n for _group, members in prov.PARAMETER_GROUPS for n in members}
        for name in ("read_type", "long_read_min_length", "hg38_db_long",
                     "t2t_phix_db_long", "rna_splice_db_long"):
            self.assertIn(name, names)
        self.assertIn("minimap2 -ax map-ont", prov.profiling_aligner({"read_type": "ont", "pks_align_min_mapq": "40"}))
        self.assertIn("minimap2 -ax map-hifi", prov.profiling_aligner({"read_type": "hifi"}))
        self.assertTrue(prov.profiling_aligner({"read_type": "short"}).startswith("bowtie2"))
        self.assertTrue(prov.profiling_aligner({}).startswith("bowtie2"))

    def test_long_indexes_carry_digests(self):
        for key in ("hg38_db_long", "t2t_phix_db_long", "rna_splice_db_long"):
            self.assertRegex(MAIN, rf"{key}\s*:\s*Provenance\.data\(params\.{key}\)")


if __name__ == "__main__":
    unittest.main()
