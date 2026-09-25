import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

class TargetedAssemblyContractTests(unittest.TestCase):
    def test_targeted_branch_is_default_for_positive_tiers(self):
        main = (ROOT / "main.nf").read_text()
        self.assertRegex(main, r"params\.tumor_targeted_assembly\s*=\s*true")
        self.assertRegex(main, r'params\.tumor_contig_tiers\s*=\s*"broad_island,extensive_island"')
        self.assertIn('tier in contig_tiers', main)
        self.assertIn("targetedPksAssembly(tumor_assembly_reads_ch, targeted_profiles_ch)", main)

    def test_both_assemblers_preserve_pairs(self):
        module = (ROOT / "Modules/pks_targeted.nf").read_text()
        self.assertIn('args+=(-1 "${r1}" -2 "${r2}")', module)
        self.assertIn("targetedPksMegahit(targetedPksRecruit.out.reads)", module)
        self.assertIn("targetedPksSpades(targetedPksRecruit.out.reads)", module)

    def test_recruitment_uses_identity_not_mapq(self):
        script_path = ROOT / "scripts/recover_recruited_mates.py"
        spec = importlib.util.spec_from_file_location("recover", script_path)
        recover = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(recover)
        with tempfile.TemporaryDirectory() as directory:
            sam = Path(directory) / "reads.sam"
            sam.write_text("r/1\t0\tref\t1\t0\t90M\t*\t0\t0\t" + "A"*90 + "\t" + "I"*90 + "\tNM:i:2\n")
            ids, alignments = recover.recruited_ids(sam, 60, .90)
            self.assertEqual(ids, {"r"})
            self.assertEqual(alignments, 1)

    def test_both_assemblers_are_aligned_to_canonical_reference(self):
        module = (ROOT / "Modules/pks_targeted.nf").read_text()
        self.assertIn("alignTargetedContigsToCanonicalReference", module)
        self.assertNotIn("alignTargetedContigsToPksPanel", module)
        self.assertIn("tuple(sampleID, 'megahit', contigs)", module)
        self.assertIn("tuple(sampleID, 'metaspades', contigs)", module)
        self.assertIn("minimap2 -x asm10 -c", module)

    def test_structural_calls_use_best_assembler_coverage(self):
        script_path = ROOT / "scripts/summarize_tumor_pks_contigs.py"
        spec = importlib.util.spec_from_file_location("summarize", script_path)
        summarize = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(summarize)
        self.assertEqual(summarize.structural_call(.95, .92), ("complete_island", "concordant"))
        self.assertEqual(summarize.structural_call(.95, .10), ("complete_island", "discordant"))
        self.assertEqual(summarize.structural_call(.80, .10), ("near_complete_island", "discordant"))
        self.assertEqual(summarize.structural_call(0, 0), ("read_evidence_only", "no_contig_support"))

    def test_fragmented_alignments_use_reference_interval_union(self):
        script_path = ROOT / "scripts/summarize_tumor_pks_contigs.py"
        spec = importlib.util.spec_from_file_location("summarize_union", script_path)
        summarize = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(summarize)
        hits = [{"start": 0, "end": 60}, {"start": 40, "end": 100}]
        self.assertEqual(summarize.union_length(hits, 0, 100), 100)


class RecruitmentStatsReachFinalEvidence(unittest.TestCase):
    """targetedPksRecruit.out.stats (pks_recruitment.tsv) used to have zero
    downstream consumers, so a sample where targeted assembly produced nothing
    gave no way to tell "no reads recruited" apart from "reads recruited but
    couldn't be paired". Both now reach final_pks_evidence.tsv.
    """

    def test_recruitment_stats_are_joined_into_the_summarize_call(self):
        module = (ROOT / "Modules/pks_targeted.nf").read_text()
        self.assertIn(
            "combined_ch = megahit_aligned_ch.join(metaspades_aligned_ch, by: 0)"
            ".join(profiles, by: 0).join(targetedPksRecruit.out.stats, by: 0)",
            module)
        process = module.split("process summarizeTargetedPksEvidence {", 1)[1]
        self.assertIn("path(recruitment_stats)", process.split("\n    script:", 1)[0])
        self.assertIn('--recruitment-stats "${recruitment_stats}"', process)

    def test_recruitment_helper_parses_the_metric_value_tsv(self):
        script_path = ROOT / "scripts/summarize_tumor_pks_contigs.py"
        spec = importlib.util.spec_from_file_location("summarize_recruit", script_path)
        summarize = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(summarize)
        with tempfile.TemporaryDirectory() as directory:
            stats = Path(directory) / "pks_recruitment.tsv"
            stats.write_text(
                "metric\tvalue\ncredible_alignments\t42\nrecruited_fragment_ids\t7\n"
                "fragment_ids_found\t7\npaired_fragments\t3\nsingleton_reads\t1\n"
                "min_aligned_bases\t60\nmin_identity\t0.9\n")
            parsed = summarize.recruitment(stats)
        self.assertEqual(parsed["recruited_fragment_ids"], "7")
        self.assertEqual(parsed["paired_fragments"], "3")

    def test_end_to_end_merges_recruitment_fields_into_the_output_row(self):
        script_path = ROOT / "scripts/summarize_tumor_pks_contigs.py"
        spec = importlib.util.spec_from_file_location("summarize_e2e", script_path)
        summarize = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(summarize)
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            megahit_contigs = d / "megahit.fa"; megahit_contigs.write_text("")
            metaspades_contigs = d / "metaspades.fa"; metaspades_contigs.write_text("")
            megahit_paf = d / "megahit.paf"; megahit_paf.write_text("")
            metaspades_paf = d / "metaspades.paf"; metaspades_paf.write_text("")
            raw_depth = d / "raw_depth.tsv"; raw_depth.write_text("")
            gff = d / "ref.gff"; gff.write_text("##gff-version 3\n")
            read_evidence = d / "read_evidence.tsv"
            read_evidence.write_text(
                "sample\tread_evidence\tpks_reads\tclb_genes_detected\t"
                "island_breadth_1x\tisland_breadth_2x\tisland_breadth_3x\n"
                "S1\tbroad_island\t120\t9\t0.20\t0.10\t0.05\n")
            recruitment_stats = d / "pks_recruitment.tsv"
            recruitment_stats.write_text(
                "metric\tvalue\ncredible_alignments\t15\nrecruited_fragment_ids\t9\n"
                "fragment_ids_found\t4\npaired_fragments\t0\nsingleton_reads\t4\n"
                "min_aligned_bases\t60\nmin_identity\t0.9\n")
            output_tsv = d / "final_pks_evidence.tsv"
            output_svg = d / "contig_validation.svg"
            old_argv = sys.argv
            try:
                sys.argv = [
                    "summarize_tumor_pks_contigs.py",
                    "--sample", "S1",
                    "--megahit-contigs", str(megahit_contigs),
                    "--metaspades-contigs", str(metaspades_contigs),
                    "--megahit-paf", str(megahit_paf),
                    "--metaspades-paf", str(metaspades_paf),
                    "--raw-depth", str(raw_depth),
                    "--read-evidence", str(read_evidence),
                    "--recruitment-stats", str(recruitment_stats),
                    "--gff", str(gff),
                    "--contig", "c1",
                    "--region-start", "0", "--region-end", "1000",
                    "--island-start", "0", "--island-end", "1000",
                    "--output-tsv", str(output_tsv),
                    "--output-svg", str(output_svg),
                ]
                summarize.main()
            finally:
                sys.argv = old_argv
            with output_tsv.open() as handle:
                row = list(csv.DictReader(handle, delimiter="\t"))[0]
        # Reads were recruited (9 fragment ids) but none of them paired up --
        # exactly the "recruited but couldn't pair" case the audit flagged.
        self.assertEqual(row["recruited_fragment_ids"], "9")
        self.assertEqual(row["paired_fragments"], "0")
        # Still carries the read-tier fields, merged the same way as before.
        self.assertEqual(row["read_evidence"], "broad_island")
        self.assertEqual(row["final_structural_evidence"], "read_evidence_only")


if __name__ == "__main__":
    unittest.main()
