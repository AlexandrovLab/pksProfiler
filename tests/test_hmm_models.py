"""The shipped clb profile HMMs are the exact_v1 population build.

exact_v1 deduplicates at 100% identity; the superseded build clustered at 99%, which on
an island this conserved collapsed most genes to one or two sequences. See
ref/hmm/PROVENANCE.md.
"""
import hashlib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HMM_DIR = ROOT / "ref/hmm"
MAIN = (ROOT / "main.nf").read_text()

PROTEIN = HMM_DIR / "clb_population_protein_exact_v1.hmm"
DNA = HMM_DIR / "clb_population_dna_exact_v1.hmm"
GENES = [f"clb{c}" for c in "ABCDEFGHIJKLMNOPQRS"]

# From pksProfiler_reference_build/models/exact_v1/SHA256SUMS.
SHA256 = {
    PROTEIN.name: "6d68ce28736c40408c3d368fd64567f04262012bca0cf00b27a4bcfca1b0272c",
    DNA.name: "95242c9f4ef86030157c8c0c89817915a2fc5df25420d1f8bb6c1e574132627c",
}


def headers(path, field):
    return [line.split(None, 1)[1].strip() for line in path.read_text().splitlines()
            if line.startswith(field)]


class ShippedModelTests(unittest.TestCase):
    def test_defaults_point_at_exact_v1(self):
        self.assertIn('params.hmm_model        = "${projectDir}/ref/hmm/clb_population_dna_exact_v1.hmm"', MAIN)
        self.assertIn('params.clb_protein_hmm         = "${projectDir}/ref/hmm/clb_population_protein_exact_v1.hmm"', MAIN)

    def test_both_models_are_present(self):
        self.assertTrue(PROTEIN.exists(), PROTEIN)
        self.assertTrue(DNA.exists(), DNA)

    def test_checksums_match_the_reference_build(self):
        for path in (PROTEIN, DNA):
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(digest, SHA256[path.name], f"{path.name} does not match the reference build")

    def test_all_nineteen_genes_in_both_alphabets(self):
        self.assertEqual(headers(PROTEIN, "NAME"), GENES)
        self.assertEqual(headers(DNA, "NAME"), GENES)
        self.assertEqual(set(headers(PROTEIN, "ALPH")), {"amino"})
        self.assertEqual(set(headers(DNA, "ALPH")), {"DNA"})

    def test_no_gene_is_built_from_a_single_sequence(self):
        # The failure mode of the superseded 99%-clustered build: five protein genes and
        # eight DNA genes had NSEQ 1, so the profile could model no variation at all.
        for path in (PROTEIN, DNA):
            counts = {name: int(n) for name, n in zip(headers(path, "NAME"), headers(path, "NSEQ"))}
            singletons = sorted(g for g, n in counts.items() if n < 2)
            self.assertEqual(singletons, [], f"{path.name} has single-sequence models: {singletons}")

    def test_nucleotide_model_is_pressed_but_protein_need_not_be(self):
        # nhmmscan reads the pressed index; hmmsearch reads the plain .hmm.
        for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
            self.assertTrue(Path(str(DNA) + suffix).exists(), f"missing {DNA.name}{suffix}")

    def test_v001_benchmark_model_is_retained(self):
        # Kept immutable for reproducing v0.0.1; selectable via --hmm_model.
        self.assertTrue((HMM_DIR / "clb_all_dna.hmm").exists())
        self.assertTrue((HMM_DIR / "clb_all_dna.hmm.h3i").exists())

    def test_unknown_provenance_protein_model_is_gone(self):
        # The model carried over from v3-working matched no known build; it was replaced.
        self.assertFalse((HMM_DIR / "clb_all_protein.hmm").exists())

    def test_pressed_index_guard_exists(self):
        self.assertIn("is not pressed", MAIN)


if __name__ == "__main__":
    unittest.main()
