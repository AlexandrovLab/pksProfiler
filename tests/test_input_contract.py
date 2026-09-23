"""One input contract, enforced identically by preflight and by the runtime.

Ludmil, revised report finding 10. Several places where preflight.py's sample-sheet
contract and the actual runtime code (extractReads.nf, filterReads.nf, map_reads.nf)
disagreed about what is valid:

- preflight accepted a .sam alignment as a normal name; extraction's own htsfile
  check rejects anything that is not BAM or CRAM, unconditionally.
- preflight accepted a plain, uncompressed .fastq/.fq; filterReads.nf runs
  `gzip -t` on the raw input before fastp ever starts, so an uncompressed file
  fails there every time.
- preflight accepted a `cram` sample-sheet column as satisfying the alignment
  requirement; main.nf resolves the alignment path from `row.alignment ?:
  row.bam` only and never reads `cram` at all.
- preflight validated --pangenome_db as a file only; map_reads.nf documents and
  implements support for a directory of .mmi files.
- preflight required --cram_reference whenever a sample looked like CRAM input;
  extractReads.nf, on a real (htsfile-detected) CRAM, used to fall through to a
  referenceless decode when --cram_reference was unset instead of failing the
  same way preflight already would have.

His acceptance test: CI must test every documented valid form and every
explicitly rejected form. tests/test_preflight.py covers the preflight side of
each disagreement (now resolved by matching preflight to the runtime contract,
or the runtime to preflight's stricter one for CRAM references); this file
checks the runtime side stays in sync.
"""
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXTRACT = (ROOT / "Modules/extract_reads.nf").read_text()
FILTER = (ROOT / "Modules/filter_reads.nf").read_text()
MAP_READS = (ROOT / "Modules/map_reads.nf").read_text()
MAIN = (ROOT / "main.nf").read_text()
PREFLIGHT = (ROOT / "scripts/preflight.py").read_text()


class ExtractionOnlySupportsBamAndCram(unittest.TestCase):
    def test_htsfile_detection_rejects_anything_else(self):
        self.assertIn('*CRAM*) DETECTED=cram ;;', EXTRACT)
        self.assertIn('*BAM*)  DETECTED=bam ;;', EXTRACT)
        self.assertIn("is neither BAM nor CRAM according to htsfile", EXTRACT)

    def test_preflight_no_longer_treats_sam_as_a_normal_name(self):
        self.assertNotIn('".sam"', PREFLIGHT.split("ALIGNMENT_SUFFIXES")[1].split("\n")[0])
        self.assertIn("extraction only supports BAM/CRAM", PREFLIGHT)


class FilterReadsRequiresGzip(unittest.TestCase):
    def test_gzip_dash_t_gates_every_input_before_fastp_runs(self):
        pre_fastp = FILTER.split("fastp \\", 1)[0]
        self.assertIn("gzip -t", pre_fastp)
        self.assertIn("Corrupt gzip input", pre_fastp)

    def test_preflight_now_requires_the_gzip_suffix(self):
        self.assertIn("GZIP_FASTQ_SUFFIXES", PREFLIGHT)
        self.assertIn("gzip-compressed", PREFLIGHT)


class OnlyAlignmentAndBamColumnsAreEverRead(unittest.TestCase):
    def test_main_nf_never_reads_a_cram_column(self):
        self.assertIn(
            'def alignment = row.alignment?.toString()?.trim() ?: row.bam?.toString()?.trim()',
            MAIN)
        self.assertNotIn("row.cram", MAIN)

    def test_preflight_no_longer_accepts_cram_as_satisfying_the_requirement(self):
        sheet_check = PREFLIGHT[PREFLIGHT.index("def check_sample_sheet"):
                                PREFLIGHT.index("def check_references")]
        # "cram" still appears as a legitimate --input_data_type *value*; the
        # thing that must be gone is the old 3-element column-satisfying set.
        self.assertIn('{"alignment", "bam"} & columns', sheet_check)
        self.assertNotIn('{"alignment", "bam", "cram"}', sheet_check)
        self.assertNotIn('row.get("cram"', sheet_check)


class PangenomeDbAcceptsAFileOrADirectory(unittest.TestCase):
    def test_map_reads_branches_on_both_forms(self):
        self.assertIn('-d "${params.pangenome_db}"', MAP_READS)
        self.assertIn('-f "${params.pangenome_db}"', MAP_READS)
        self.assertIn("*.mmi", MAP_READS)

    def test_preflight_now_accepts_both_forms_too(self):
        self.assertIn('kind="file_or_dir"', PREFLIGHT)
        self.assertIn("--pangenome_db {pangenome} is a directory", PREFLIGHT)


class CramReferenceIsRequiredEverywhereRealCramIsDetected(unittest.TestCase):
    def test_extraction_no_longer_falls_through_to_a_referenceless_decode(self):
        cram_block = EXTRACT[EXTRACT.index('if [[ "\\$DETECTED" == "cram" ]]'):
                             EXTRACT.index('samtools quickcheck')]
        self.assertIn("CRAM input requires --cram_reference", cram_block)
        self.assertNotIn('-n "${params.cram_reference', cram_block)

    def test_preflight_already_required_it_and_still_does(self):
        self.assertIn("CRAM input needs --cram_reference", PREFLIGHT)


if __name__ == "__main__":
    unittest.main()
