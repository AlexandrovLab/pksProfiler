import importlib.util
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

if __name__ == "__main__":
    unittest.main()
