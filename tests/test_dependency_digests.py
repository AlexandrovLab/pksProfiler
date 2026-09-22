"""Every task that uses a reference, model or helper carries its digest.

F15, from Ludmil's 2026-09-19 report. Nextflow hashes a task from its declared inputs
and its rendered script. A reference or helper named only as a path string inside the
script is neither, so rebuilding an index in place or editing a helper leaves the hash
unchanged and `-resume` reuses tasks that no longer match the current code.

main.nf computes a digest per dependency (lib/Provenance.groovy) and each task script
carries the digests of the dependencies it actually uses -- only those, so a database
change does not invalidate extraction. This test is the part that has to hold as
processes are added: use a tracked dependency, carry its digest.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()
PROVENANCE = ROOT / "lib/Provenance.groovy"

# params name as written in a script -> key in params.dep_digest
TRACKED = {
    "scripts": "scripts", "hmm_model": "hmm_model", "clb_protein_hmm": "clb_protein_hmm",
    "pks_genome_annotation": "pks_annotation", "pks_cytoband": "pks_cytoband",
    "adapters": "adapters", "hg38_db": "hg38_db", "t2t_phix_db": "t2t_phix_db",
    "pangenome_db": "pangenome_db", "pks_genome": "pks_genome",
    "pks_recruit_index": "pks_recruit_index", "pks_reference_fasta": "pks_reference_fasta",
    "cram_reference": "cram_reference", "kraken_db": "kraken_db",
    "genomad_db": "genomad_db", "gtdbtk_db": "gtdbtk_db", "checkm2_db": "checkm2_db",
}


def process_blocks():
    """(module name, process name, block text) for every process in Modules/."""
    for module in sorted((ROOT / "Modules").glob("*.nf")):
        text = module.read_text()
        starts = [m.start() for m in re.finditer(r"^process\s+\w+\s*\{", text, re.M)]
        for start, end in zip(starts, starts[1:] + [len(text)]):
            block = text[start:end]
            yield module.name, re.match(r"^process\s+(\w+)", block).group(1), block


class DigestCoverage(unittest.TestCase):
    def test_every_process_carries_the_digests_it_depends_on(self):
        missing = []
        for module, process, block in process_blocks():
            for param, key in TRACKED.items():
                if not re.search(rf"params\.{param}\b", block):
                    continue
                if f"params.dep_digest?.{key}" not in block:
                    missing.append(f"{module}:{process} uses params.{param} without its digest")
        self.assertEqual(missing, [], "\n".join(missing))

    def test_no_process_carries_a_digest_it_does_not_use(self):
        # Over-tagging is not harmless: an unrelated database change would invalidate
        # per-sample work across the cohort.
        extra = []
        for module, process, block in process_blocks():
            for param, key in TRACKED.items():
                if f"params.dep_digest?.{key}" not in block:
                    continue
                if not re.search(rf"params\.{param}\b", block):
                    extra.append(f"{module}:{process} carries the {key} digest but does not use it")
        self.assertEqual(extra, [], "\n".join(extra))

    def test_main_defines_a_digest_for_every_key_the_modules_reference(self):
        referenced = set()
        for _module, _process, block in process_blocks():
            referenced.update(re.findall(r"params\.dep_digest\?\.(\w+)", block))
        defined = set(re.findall(r"^\s{4}(\w+)\s*:\s*Provenance\.", MAIN, re.M))
        self.assertEqual(referenced - defined, set(),
                         f"referenced but not computed in main.nf: {sorted(referenced - defined)}")

    def test_code_and_data_modes_are_used_for_the_right_kinds_of_dependency(self):
        # Content for the small things, metadata for indexes and databases that cannot
        # be read on every launch.
        for key in ("scripts", "hmm_model", "clb_protein_hmm"):
            self.assertRegex(MAIN, rf"{key}\s*:\s*Provenance\.code\(")
        for key in ("hg38_db", "kraken_db", "gtdbtk_db", "pks_genome"):
            self.assertRegex(MAIN, rf"{key}\s*:\s*Provenance\.data\(")

    def test_the_helper_exists_and_offers_both_modes(self):
        self.assertTrue(PROVENANCE.exists(), "lib/Provenance.groovy is missing")
        text = PROVENANCE.read_text()
        self.assertIn("static String code(", text)
        self.assertIn("static String data(", text)


if __name__ == "__main__":
    unittest.main()
