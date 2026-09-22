"""Profiling may narrow the read stream; assembly never may.

The prefilter (krakenuniq + DIAMOND rescue) produces PROFILING_READS, which feeds alignment
and HMM profiling only. Both assembly paths -- metagenome MAG reconstruction and tumour
targeted assembly -- draw on MAPPED_READS, the complete host-depleted stream, in every lane
and regardless of --prefilter_mode. Binning needs whole-community coverage, and targeted
assembly needs full depth over the island; handing either a prefiltered subset would be a
silent scientific regression.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()
BODY = MAIN.split("workflow {", 1)[1]


def assignments_of(channel):
    """Lines where `channel` is the source being consumed or piped."""
    return [line.strip() for line in BODY.splitlines() if channel in line and not line.strip().startswith("//")]


class ReadStreamInvariantTests(unittest.TestCase):
    def test_both_assembly_paths_read_the_complete_stream(self):
        self.assertRegex(BODY, r"MAPPED_READS\s*\n\s*\.filter\s*\{[^}]*enable_mags_b")
        self.assertIn("tumor_assembly_reads_ch = MAPPED_READS", BODY)

    def test_no_assembly_path_reads_the_narrowed_stream(self):
        for forbidden in ("pksMAG(PROFILING_READS)",
                          "tumor_assembly_reads_ch = PROFILING_READS",
                          "tumorWGS(PROFILING_READS)",
                          "MAG_ASSEMBLY_READS = PROFILING_READS"):
            self.assertNotIn(forbidden, BODY)

    def test_profiling_reads_only_feed_profiling(self):
        consumers = re.findall(r"(\w+)\(PROFILING_READS\)", BODY)
        self.assertEqual(sorted(consumers), ["pksProfilerAlign", "pksProfilerHMM"])

    def test_mag_branch_precedes_the_prefilter(self):
        # If the branch moved below the prefilter block it would silently inherit the
        # narrowed stream even though the code still says MAPPED_READS.
        self.assertLess(BODY.index("MAG_ASSEMBLY_READS"), BODY.index("krakenPrefilter(MAPPED_READS)"))

    def test_prefilter_off_is_a_passthrough(self):
        self.assertIn("MAPPED_READS.set { PROFILING_READS }", BODY)

    def test_auto_resolves_to_balanced_only_for_metagenomes_with_a_kraken_db(self):
        self.assertIn('if (prefilter_mode == "auto")', BODY)
        block = BODY[BODY.index('if (prefilter_mode == "auto")'):]
        block = block[:block.index("log.info")]
        self.assertIn('params.sample_type == "metagenome" && params.kraken_db', block)
        self.assertIn('prefilter_mode = "balanced"', block)
        self.assertIn('prefilter_mode = "off"', block)

    def test_prefilter_is_not_gated_on_sample_type(self):
        # The prefilter applies to tumour and metagenome alike; only assembly is exempt.
        # Take the branch that actually runs it, not the earlier validation block.
        call = BODY.index("krakenPrefilter(MAPPED_READS)")
        branch = BODY.rindex('if (prefilter_mode == "balanced") {', 0, call)
        self.assertNotIn("sample_type", BODY[branch:call])


class ReadEvidenceScopeTests(unittest.TestCase):
    def test_evidence_is_computed_for_tumour_and_metagenome(self):
        self.assertIn('if (params.sample_type in ["tumor_wgs", "metagenome"])', BODY)
        self.assertIn("classifyPksReadEvidence(read_evidence_input_ch)", BODY)

    def test_only_the_tumour_lane_gates_on_the_tier(self):
        # For metagenomes the tier is an annotation; the contig-count gate is the real gate.
        gate = "tier in contig_tiers"
        self.assertIn(gate, BODY)
        tumour_block = BODY[BODY.index('if (params.sample_type == "tumor_wgs") {'):]
        self.assertIn(gate, tumour_block)
        # pksMAG on the metagenome branch is not tier-conditioned.
        mag_call = BODY[BODY.index("if (enable_mags_b) {"):]
        self.assertNotIn(gate, mag_call[:mag_call.index("}")])


if __name__ == "__main__":
    unittest.main()
