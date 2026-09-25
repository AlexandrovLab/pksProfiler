#!/usr/bin/env python3
"""Static integration contracts for the optional MAG branch."""

import re
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MAIN = (REPO / "main.nf").read_text()
MAG = (REPO / "Modules" / "pks_mag.nf").read_text()
BASE_CONFIG = (REPO / "conf" / "base.config").read_text()


class MagRoutingTests(unittest.TestCase):
    def test_mag_branch_consumes_preserved_reads_not_profile_candidates(self):
        self.assertIn(".set { MAG_ASSEMBLY_READS }", MAIN)
        self.assertIn("pksMAG(MAG_ASSEMBLY_READS)", MAIN)
        self.assertNotIn("pksMAG(PROFILING_READS)", MAIN)

    def test_mag_routing_requires_explicit_metagenome_opt_in(self):
        self.assertRegex(
            MAIN,
            re.compile(
                r'MAPPED_READS\s*\.filter\s*\{.*?'
                r'params\.sample_type\s*==\s*"metagenome"\s*&&\s*'
                r'enable_mags_b.*?\.set\s*\{\s*MAG_ASSEMBLY_READS\s*\}',
                re.DOTALL,
            ),
        )

    def test_mag_branch_does_not_consume_prefilter_outputs(self):
        mag_call = MAIN.index("pksMAG(MAG_ASSEMBLY_READS)")
        routing = MAIN.rfind("MAPPED_READS", 0, mag_call)
        self.assertGreaterEqual(routing, 0)
        self.assertNotIn("KRAKEN_PREFILTER_OUT", MAIN[routing:mag_call])
        self.assertNotIn("DIAMOND_RESCUE_OUT", MAIN[routing:mag_call])


    def test_tumor_wgs_targeted_validation_is_evidence_gated(self):
        for parameter in ('tumor_multi_gene_min_clb_genes', 'tumor_multi_gene_min_pks_reads', 'tumor_multi_gene_min_breadth', 'tumor_broad_island_min_clb_genes', 'tumor_broad_island_min_pks_reads', 'tumor_broad_island_min_breadth', 'tumor_extensive_island_min_clb_genes', 'tumor_extensive_island_min_pks_reads', 'tumor_extensive_island_min_breadth'):
            self.assertIn(parameter, MAIN)
        self.assertIn('tier in contig_tiers', MAIN)
        self.assertIn("tumor_assembly_reads_ch = MAPPED_READS", MAIN)
        self.assertIn("targetedPksAssembly(tumor_assembly_reads_ch, targeted_profiles_ch)", MAIN)


class MagModuleContractTests(unittest.TestCase):
    def test_single_and_paired_megahit_modes_are_present(self):
        self.assertIn('"-1 ${reads_list[0]} -2 ${reads_list[1]}"', MAG)
        self.assertIn('"-r ${reads_list[0]}"', MAG)

    def test_invalid_read_cardinality_is_rejected(self):
        self.assertEqual(MAG.count("reads_list.size() in [1, 2]"), 2)

    def test_empty_bin_output_is_explicitly_optional(self):
        self.assertRegex(
            MAG,
            r'path\("bins/bin\.\[0-9\]\*\.fa"\), emit: bins, optional: true',
        )
        self.assertIn('emit: status', MAG)
        self.assertIn('contig_count', MAG)
        self.assertIn('bin_count', MAG)

    def test_unbinned_contigs_are_pooled_and_explicitly_optional(self):
        # U1: real bins (numeric suffix) and the unbinned pool are two distinct,
        # non-overlapping glob patterns -- the unbinned file must never be counted
        # as a bin, which is exactly the false positivity this pipeline argues
        # against elsewhere (the --tumor_enable_mags removal, same commit series).
        self.assertRegex(
            MAG,
            r'path\("bins/bin\.unbinned\.fa"\), emit: unbinned, optional: true',
        )
        self.assertIn("--unbinned", MAG)

    def test_sample_status_distinguishes_successful_mag_outcomes(self):
        for status in ("no_contigs", "contigs_no_bins", "bins_no_pks", "pks_positive_bins"):
            self.assertIn(status, MAG)
        self.assertIn("magSampleStatus(mag_status_input_ch)", MAG)
        self.assertIn("no_contig_counts_ch", MAG)
        self.assertIn("reads.join(nonempty_contigs_ch", MAG)

    def test_no_pks_context_is_retained_for_summary(self):
        self.assertIn("remainder: true", MAG)
        self.assertIn("contexts ?: []", MAG)

    def test_prophages_are_called_for_all_bins_not_only_pks_positive_bins(self):
        self.assertIn("genomadProphages(bins_flat_ch)", MAG)
        self.assertNotIn("genomadProphages(pks_pos", MAG)
        self.assertIn("communityProphageSummary(community_input_ch)", MAG)

    def test_community_outputs_are_explicitly_hypothesis_generating(self):
        script = (REPO / "scripts" / "build_community_prophage.py").read_text()
        self.assertIn("co-occurrence_only_not_evidence_of_induction", script)
        self.assertIn("recipient_clbS_like", script)
        self.assertIn("recipient_has_recA", script)
        self.assertIn("recipient_has_lexA", script)

    def test_shell_processes_enable_strict_failure_handling(self):
        process_count = len(re.findall(r"^process\s+\w+\s*\{", MAG, re.MULTILINE))
        strict_count = MAG.count("set -euo pipefail")
        self.assertEqual(strict_count, process_count)

    def test_all_mag_resource_labels_are_configured(self):
        for label in ("mag_assembly", "mag_binning", "mag_gtdbtk", "mag_prophage", "mag_hmm"):
            self.assertIn(f"withLabel:{label}", BASE_CONFIG)


if __name__ == "__main__":
    unittest.main()
