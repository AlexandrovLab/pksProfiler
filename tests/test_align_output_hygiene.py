"""F13 and F18 (Ludmil, GATE 4): pksProfiler_align.nf's coverage track and depth
computation, guarded against silently regressing back to their pre-fix behaviour.

Both fixes landed in commit e6d91f7 with no dedicated regression test -- the tracker
records the reasoning (handoff/PKSPROFILER_FIX_TRACKER.md, GATE 4 / F13 and GATE 5 /
F18) but nothing here would fail if either reverted.

F13. `bamCoverage --normalizeUsing RPKM` ran over the island-aligned BAM, so its "per
million mapped reads" denominator was reads mapping to the island at MAPQ >= 40 -- not
the library -- while the unit name claimed library-normalised abundance. The fix is
`--normalizeUsing None` (raw depth); this fails if RPKM (or any other library-scale
normalisation keyword) reappears.

F18. `bedtools genomecov -ibam -d` wrote one line per base of the whole ~5.1 Mb
reference, of which the pks island is ~1%, and nothing downstream read it (the file was
recomputed from the BAM anyway). The fix is a `samtools depth` call restricted to the
island region (`-r <contig>:<start>-<end>`), the same tool and filters used everywhere
else in this codebase for this purpose. This also fails if the process resumes writing a
`.sam` file: bowtie2 used to write one to disk, published it, and nothing read it either
(main.nf destructured and dropped it) -- the same "detailed alignment retention" the
finding's Fix line asks to make optional, closed here by removing it outright rather than
gating it.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ALIGN_NF = ROOT / "Modules/pksProfiler_align.nf"
ALIGN = ALIGN_NF.read_text()


def code_only(text):
    return "\n".join(line for line in text.splitlines()
                      if not line.strip().startswith(("#", "//")))


class CoverageTrackIsRawDepthNotRpkm(unittest.TestCase):
    """F13."""

    def test_bamcoverage_normalizes_to_none(self):
        self.assertIn("--normalizeUsing None", code_only(ALIGN))

    def test_no_library_scale_normalization_keyword_reappears(self):
        code = code_only(ALIGN)
        for keyword in ("RPKM", "CPM", "BPM", "RPGC"):
            self.assertNotIn(keyword, code,
                              f"{keyword} normalisation reappeared in pksProfiler_align.nf "
                              "(F13: it is not a library-normalised abundance estimate "
                              "over an island-restricted BAM)")


class DepthIsIslandRestrictedAndNoSamIsWritten(unittest.TestCase):
    """F18."""

    def test_no_genome_wide_genomecov_call(self):
        self.assertNotIn("genomecov", code_only(ALIGN),
                          "bedtools genomecov reappeared -- F18 replaced the genome-wide "
                          "(~5.1 Mb) depth file with samtools depth restricted to the island")

    def test_depth_is_computed_with_samtools_depth_over_the_island_region(self):
        code = code_only(ALIGN)
        self.assertIn("samtools depth", code)
        self.assertIn('-r "${params.pks_contig}:${params.pks_start_1based}-${params.pks_end_1based}"',
                       code)

    def test_no_sam_file_is_declared_or_published(self):
        # A literal ".sam" filename ending (e.g. "foo.sam"), not the substring inside
        # unrelated identifiers like "params.sample_dir".
        self.assertIsNone(re.search(r'\.sam["\'\s]', ALIGN),
                           "pksProfiler_align.nf used to write, declare and publish a "
                           ".sam file that nothing read; F18 removed it by streaming "
                           "bowtie2 straight into the sorter")

    def test_bowtie2_streams_into_the_sorter_without_an_intermediate_dash_s_flag(self):
        code = code_only(ALIGN)
        # bowtie2's own -S flag writes a SAM file to disk; the fixed invocation pipes
        # stdout straight into samtools view | samtools sort instead, across the
        # backslash-continued lines of the same shell command.
        self.assertRegex(code, re.compile(r"bowtie2\s+-x.*?\|\s*samtools view", re.DOTALL))
        self.assertNotRegex(code, r"bowtie2\b(?:(?!\n\n)[^\n])*\s-S\s")


if __name__ == "__main__":
    unittest.main()
