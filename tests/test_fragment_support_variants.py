"""Real, executing coverage for scripts/compute_fragment_support_variants.py.

The PI's ask: quantify raw reads, duplicate-marked support, unique
fragments/templates, and duplicate-collapsed support for a pks-positive tumour
sample, and confirm the read-evidence tier (multi_gene / broad_island /
extensive_island) is stable across all four. This builds a tiny synthetic
sorted BAM (hand-placed alignment records against a tiny reference with three
fake clb genes, plus 16 placeholder clb genes clb_gene_bed.py otherwise
requires) and runs REAL samtools (markdup, depth), REAL bedtools (bamtobed,
intersect) and REAL featureCounts against it -- exactly the tools and
invocations pksProfiler_align.nf and Modules/pks_taxa.nf already use -- then
asserts on the script's real output, not hand-constructed intermediates
standing in for what those tools would have produced.

Two fixtures:

  * MECHANICS: proves (a) two mates of one fragment landing in the SAME gene
    collapse to one fragment, not two reads; (b) a fragment whose two mates
    land in DIFFERENT genes counts toward both genes' fragment totals (see
    compute_fragment_support_variants.py's module docstring for why); and
    (c) samtools markdup, run directly on this no-real-pairing BAM with no
    collate/fixmate step, actually flags an exact-duplicate pair 0x400.

  * BOUNDARY: a minimal case engineered so the tier the script recomputes
    for dup_marked/dup_collapsed genuinely differs from raw -- proving the
    stability check catches a real divergence and names the right variants,
    not just that it always prints PASS.
"""
import csv
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/compute_fragment_support_variants.py"
CLB_GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]

CONTIG = "ctg1"
CONTIG_LEN = 1000
ISLAND_START = 1
ISLAND_END = 200  # island_length = 200; clbA/clbB/clbC (the only genes any read touches) sit well inside it.


def _resolve_bin(tool, env_yml):
    """Prefer an explicit override, then the pinned pks_align_env conda build
    already cached for this repo (conda_envs/pks_align_env.yml -- samtools
    1.24, bedtools 2.31.1, subread/featureCounts 2.1.1), then PATH. Mirrors
    tests/check_e2e_execution.sh's resolve_minimap2()."""
    override = os.environ.get(f"{tool.upper()}_BIN")
    if override and os.path.isfile(override):
        return override
    cache_dir = os.environ.get(
        "CONDA_CACHE_DIR",
        "/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache",
    )
    for history in glob.glob(os.path.join(cache_dir, "env-*", "conda-meta", "history")):
        try:
            text = Path(history).read_text()
        except OSError:
            continue
        if env_yml in text:
            candidate = Path(history).parents[1] / "bin" / tool
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)
    return shutil.which(tool)


SAMTOOLS = _resolve_bin("samtools", "conda_envs/pks_align_env.yml")
BEDTOOLS = _resolve_bin("bedtools", "conda_envs/pks_align_env.yml")
FEATURECOUNTS = _resolve_bin("featureCounts", "conda_envs/pks_align_env.yml")
TOOLS_AVAILABLE = all((SAMTOOLS, BEDTOOLS, FEATURECOUNTS))


def gene_interval(index):
    """The i-th clb gene's 1-based [start, end] under clb_gene_bed.py's own
    convention: 19 genes, 30 bp each, 50 bp apart, starting at position 11."""
    start = 11 + index * 50
    return start, start + 29


def annotation_gff():
    """All 19 clb genes clb_gene_bed.py requires. Only clbA/clbB/clbC (indices
    0-2) sit inside the island any read or gene interval below touches; the
    rest are placeholders so clb_gene_bed.py's completeness check passes,
    exactly as a real --pks_genome_annotation always carries all 19."""
    lines = ["##gff-version 3"]
    for index, gene in enumerate(CLB_GENES):
        start, end = gene_interval(index)
        lines.append(f"{CONTIG}\tsrc\tgene\t{start}\t{end}\t.\t+\t.\tID=g{gene};Name={gene}")
    return "\n".join(lines) + "\n"


def sam_record(read_id, pos, mapq=60, length=20):
    seq = ("ACGT" * ((length // 4) + 1))[:length]
    qual = "I" * length
    return f"{read_id}\t0\t{CONTIG}\t{pos}\t{mapq}\t{length}M\t*\t0\t0\t{seq}\t{qual}"


def build_bam(directory, reads):
    """reads: [(read_id, pos), ...]. Returns the coordinate-sorted BAM path,
    built with real samtools (view+sort), matching pksProfiler_align.nf's own
    `samtools sort -o sample.sorted.bam` output shape."""
    sam = directory / "in.sam"
    header = [
        "@HD\tVN:1.6\tSO:unsorted",
        f"@SQ\tSN:{CONTIG}\tLN:{CONTIG_LEN}",
    ]
    body = [sam_record(read_id, pos) for read_id, pos in reads]
    sam.write_text("\n".join(header + body) + "\n")
    bam = directory / "sample.sorted.bam"
    unsorted = directory / "in.bam"
    subprocess.run([SAMTOOLS, "view", "-bS", "-o", str(unsorted), str(sam)],
                   check=True, capture_output=True)
    subprocess.run([SAMTOOLS, "sort", "-o", str(bam), str(unsorted)],
                   check=True, capture_output=True)
    subprocess.run([SAMTOOLS, "index", str(bam)], check=True, capture_output=True)
    return bam


def build_variant1_inputs(directory, bam, annotation_path):
    """Real featureCounts + real samtools depth, producing exactly the files
    pksProfiler_align.nf publishes -- counts.txt, alignment.qc.tsv,
    coverage.txt -- so variant 1 in the script under test is reading real
    tool output, not a hand-authored stand-in."""
    counts = directory / "sample.counts.txt"
    subprocess.run(
        [FEATURECOUNTS, "-a", str(annotation_path), "-o", str(counts),
         "-t", "gene", "-F", "GFF", "-g", "Name", "--largestOverlap", str(bam)],
        check=True, capture_output=True, text=True,
    )

    mapped = subprocess.run([SAMTOOLS, "view", "-c", "-F", "4", str(bam)],
                             check=True, capture_output=True, text=True).stdout.strip()
    clb_reads = 0
    genes_detected = 0
    with counts.open() as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] in CLB_GENES:
                value = float(fields[-1])
                clb_reads += value
                if value >= 1:
                    genes_detected += 1

    qc = directory / "sample.alignment.qc.tsv"
    with qc.open("w") as handle:
        handle.write("Sample\tMetric\tValue\n")
        handle.write(f"sample\treads_mapped_ihe3034\t{mapped}\n")
        handle.write(f"sample\treads_clb_genes_align\t{int(clb_reads)}\n")
        handle.write(f"sample\tnum_clb_genes_align\t{genes_detected}\n")

    coverage = directory / "sample.coverage.txt"
    depth = subprocess.run(
        [SAMTOOLS, "depth", "-aa", "-s", "-Q", "40",
         "-r", f"{CONTIG}:{ISLAND_START}-{ISLAND_END}", str(bam)],
        check=True, capture_output=True, text=True,
    ).stdout
    coverage.write_text(depth)

    return counts, qc, coverage


def run_script(directory, bam, counts, qc, coverage, annotation_path, extra_args=()):
    output = directory / "variants.tsv"
    cmd = [
        sys.executable, str(SCRIPT),
        "--sample", "sample",
        "--bam", str(bam),
        "--counts", str(counts),
        "--qc", str(qc),
        "--coverage", str(coverage),
        "--annotation", str(annotation_path),
        "--output", str(output),
        "--samtools", SAMTOOLS,
        "--bedtools", BEDTOOLS,
        "--pks-contig", CONTIG,
        "--pks-start", str(ISLAND_START),
        "--pks-end", str(ISLAND_END),
        *extra_args,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result, output


def parse_report(path):
    rows = {}
    verdict = None
    diverged = []
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if not fields:
                continue
            if fields[0] == "sample" and fields[1] == "variant":
                continue
            if fields[0].startswith("# verdict"):
                verdict = fields[1]
                continue
            if fields[0].startswith("# diverged"):
                diverged.append(fields[1])
                continue
            _sample, variant, count, genes_detected, breadth, tier = fields
            rows[variant] = {
                "count": int(count),
                "genes_detected": int(genes_detected),
                "breadth_1x": float(breadth),
                "tier": tier,
            }
    return rows, verdict, diverged


@unittest.skipUnless(TOOLS_AVAILABLE, "real samtools/bedtools/featureCounts not found on PATH or in CONDA_CACHE_DIR")
class MechanicsTests(unittest.TestCase):
    """One fixture; fragment collapsing, gene-disagreement, and a real
    markdup-flagged duplicate pair, all in one small BAM.

        F1/1@12, F1/2@16   -- one fragment, both mates in clbA (11-40)
        F2/1@20            -- clbA
        F2/2@65            -- clbB (61-90) -- F2's mates DISAGREE on gene
        dupA/1@115,dupB/1@115 -- clbC (111-140), identical position/orientation:
                                  an exact-duplicate PAIR for markdup to flag
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        directory = Path(cls.tmp.name)
        cls.annotation = directory / "annotation.gff"
        cls.annotation.write_text(annotation_gff())

        reads = [
            ("F1/1", 12), ("F1/2", 16),
            ("F2/1", 20), ("F2/2", 65),
            ("dupA/1", 115), ("dupB/1", 115),
        ]
        cls.bam = build_bam(directory, reads)
        cls.counts, cls.qc, cls.coverage = build_variant1_inputs(directory, cls.bam, cls.annotation)
        cls.result, cls.output = run_script(
            directory, cls.bam, cls.counts, cls.qc, cls.coverage, cls.annotation)

    def test_script_exits_cleanly(self):
        self.assertEqual(self.result.returncode, 0, self.result.stderr)

    def test_variant1_matches_the_published_files_unchanged(self):
        rows, _verdict, _diverged = parse_report(self.output)
        # 3 (clbA) + 1 (clbB) + 2 (clbC) = 6 raw reads, straight off real featureCounts/qc.tsv.
        self.assertEqual(rows["raw"]["count"], 6)
        self.assertEqual(rows["raw"]["genes_detected"], 3)

    def test_a_same_gene_fragment_collapses_to_one_not_two(self):
        # F1's two mates both land in clbA: 3 raw clbA reads (F1/1, F1/2, F2/1)
        # but only 2 distinct clbA fragments (F1, F2) once F1's mates collapse.
        rows, _verdict, _diverged = parse_report(self.output)
        self.assertEqual(rows["raw"]["count"], 6)
        self.assertLess(rows["fragment"]["count"], rows["raw"]["count"])
        # Total distinct fragments touching any gene: F1, F2, dupA, dupB = 4.
        self.assertEqual(rows["fragment"]["count"], 4)
        # Incidentally (this fixture is about mechanics, not tier stability),
        # 4 fragments falls below the multi_gene read/fragment floor of 5, so
        # this fixture's own fragment/dup_collapsed tiers legitimately diverge
        # from raw -- TierBoundaryTests below is the deliberately-engineered
        # case for asserting on that divergence; this one just checks counts.

    def test_markdup_flags_exactly_one_of_the_exact_duplicate_pair(self):
        dup_ids = set()
        with tempfile.TemporaryDirectory() as tmp:
            markdup_bam = Path(tmp) / "markdup.bam"
            subprocess.run([SAMTOOLS, "markdup", "-s", str(self.bam), str(markdup_bam)],
                           check=True, capture_output=True, text=True)
            flagged = subprocess.run([SAMTOOLS, "view", "-f", "0x400", str(markdup_bam)],
                                      check=True, capture_output=True, text=True).stdout
            for line in flagged.splitlines():
                if line:
                    dup_ids.add(line.split("\t", 1)[0])
        self.assertEqual(dup_ids, {"dupA/1"} if "dupA/1" in dup_ids else {"dupB/1"})
        self.assertEqual(len(dup_ids), 1)
        self.assertTrue(dup_ids <= {"dupA/1", "dupB/1"})

    def test_dup_marked_reads_drop_by_exactly_the_flagged_duplicate(self):
        rows, _verdict, _diverged = parse_report(self.output)
        self.assertEqual(rows["dup_marked"]["count"], rows["raw"]["count"] - 1)

    def test_a_fragment_split_across_two_genes_counts_toward_both(self):
        # F2/1 -> clbA, F2/2 -> clbB: gene-level fragment counts (clbA=2, clbB=1,
        # clbC=2 -> sum 5) exceed the 4 distinct total fragments, because F2 is
        # counted once under each gene it touched. This is the documented,
        # intentional handling of a fragment whose mates disagree.
        import importlib.util
        spec = importlib.util.spec_from_file_location("compute_fragment_support_variants", SCRIPT)
        cfsv = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cfsv)
        with tempfile.TemporaryDirectory() as tmp:
            clb_bed = Path(tmp) / "clb_genes.bed"
            cfsv.build_clb_bed(self.annotation, clb_bed)
            assignment = cfsv.bam_to_gene_assignment(self.bam, clb_bed, BEDTOOLS, 1)
            gene_counts, total = cfsv.gene_counts_from_assignment(assignment, by_fragment=True)
        self.assertEqual(total, 4)
        self.assertEqual(sum(gene_counts.values()), 5)
        self.assertIn("clbA", gene_counts)
        self.assertIn("clbB", gene_counts)


@unittest.skipUnless(TOOLS_AVAILABLE, "real samtools/bedtools/featureCounts not found on PATH or in CONDA_CACHE_DIR")
class TierBoundaryTests(unittest.TestCase):
    """A minimal case where excluding duplicates genuinely flips the tier:

        a1/1@15   -- clbA, 1 read
        b1/1@65   -- clbB, 1 read
        c1..c4/1@115 -- clbC, 4 reads, all identical position/orientation

    Raw: 6 reads, 3 genes, breadth well above 1% -> multi_gene (>=5 reads).
    markdup collapses the 4 identical clbC reads to 1 survivor + 3 flagged
    duplicates: dup-marked/dup-collapsed reads drop to 3 (< 5) ->
    localized_indeterminate. Fragment grouping alone does not collapse
    anything here (c1..c4 are 4 distinct fragment IDs, not one fragment
    sequenced four times), so the fragment variant alone stays multi_gene --
    proving dedup and fragment-grouping are genuinely independent axes and the
    stability check flags only the variants that actually diverged.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        directory = Path(cls.tmp.name)
        cls.annotation = directory / "annotation.gff"
        cls.annotation.write_text(annotation_gff())

        reads = [("a1/1", 15), ("b1/1", 65),
                 ("c1/1", 115), ("c2/1", 115), ("c3/1", 115), ("c4/1", 115)]
        cls.bam = build_bam(directory, reads)
        cls.counts, cls.qc, cls.coverage = build_variant1_inputs(directory, cls.bam, cls.annotation)
        cls.result, cls.output = run_script(
            directory, cls.bam, cls.counts, cls.qc, cls.coverage, cls.annotation)

    def test_script_exits_cleanly(self):
        self.assertEqual(self.result.returncode, 0, self.result.stderr)

    def test_raw_is_multi_gene(self):
        rows, _verdict, _diverged = parse_report(self.output)
        self.assertEqual(rows["raw"]["count"], 6)
        self.assertEqual(rows["raw"]["genes_detected"], 3)
        self.assertEqual(rows["raw"]["tier"], "multi_gene")

    def test_fragment_alone_stays_stable(self):
        rows, _verdict, _diverged = parse_report(self.output)
        self.assertEqual(rows["fragment"]["count"], 6)  # no shared fragment IDs to collapse
        self.assertEqual(rows["fragment"]["tier"], "multi_gene")

    def test_dedup_pushes_reads_below_the_multi_gene_floor(self):
        rows, _verdict, _diverged = parse_report(self.output)
        self.assertEqual(rows["dup_marked"]["count"], 3)
        self.assertEqual(rows["dup_marked"]["tier"], "localized_indeterminate")
        self.assertEqual(rows["dup_collapsed"]["count"], 3)
        self.assertEqual(rows["dup_collapsed"]["tier"], "localized_indeterminate")

    def test_verdict_is_flag_naming_the_diverging_variants_and_direction(self):
        _rows, verdict, diverged = parse_report(self.output)
        self.assertEqual(verdict, "FLAG")
        joined = " ".join(diverged)
        self.assertIn("dup_marked", joined)
        self.assertIn("dup_collapsed", joined)
        self.assertIn("multi_gene -> localized_indeterminate", joined)
        self.assertNotIn("fragment:", joined)


@unittest.skipUnless(TOOLS_AVAILABLE, "real samtools/bedtools/featureCounts not found on PATH or in CONDA_CACHE_DIR")
class StableCaseIsPass(unittest.TestCase):
    """A plain case with no duplicates and no shared fragment IDs: all four
    variants must agree, and the verdict must say so."""

    def test_no_duplicates_no_shared_fragments_is_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            annotation = directory / "annotation.gff"
            annotation.write_text(annotation_gff())
            reads = [("r1/1", 15), ("r2/1", 20), ("r3/1", 65), ("r4/1", 70), ("r5/1", 115)]
            bam = build_bam(directory, reads)
            counts, qc, coverage = build_variant1_inputs(directory, bam, annotation)
            result, output = run_script(directory, bam, counts, qc, coverage, annotation)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows, verdict, diverged = parse_report(output)
        self.assertEqual(verdict, "PASS")
        self.assertEqual(diverged, ["none -- tier is stable across all four variants"])
        tiers = {rows[v]["tier"] for v in ("raw", "dup_marked", "fragment", "dup_collapsed")}
        self.assertEqual(len(tiers), 1)


if __name__ == "__main__":
    unittest.main()
