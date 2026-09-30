"""Ludmil, 25 Sep A1/A2: targeted-recruitment specificity and recruitment classes.

A1. "A long read can pass with a short high-identity local match and extensive soft
clipping. Add an explicit query-coverage criterion for internal island reads ... Use a
separate boundary-evidence rule that allows lower query coverage when the alignment
reaches/crosses the canonical island boundary within a defined tolerance. Consider
ordinary best-alignment reporting and an explicit MAPQ criterion rather than relying on
-k 1." Acceptance: "A 150-bp read with only a 60-bp internal perfect local match must
fail; a true boundary-spanning fragment should pass the dedicated boundary rule."

A2. "Report three recruitment classes separately and drive the 'why did assembly fail?'
diagnostic from island-hitting plus boundary-spanning fragments only." Acceptance: "A
pks-negative flank-only fixture must show flank-only recruitment but zero
island-supporting fragments."

Every fixture here is a hand-written synthetic SAM record against the public IHE3034
island+flank recruitment reference's coordinates (indices/pks_island); no sample data.
"""
import csv
import gzip
import importlib.util
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/recover_recruited_mates.py"
MAIN = (ROOT / "main.nf").read_text()
TARGETED = (ROOT / "Modules/pks_targeted.nf").read_text()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


recover = load("recover_recruited_mates", SCRIPT)

# The shipped recruitment reference: 5 kb flank + 50,768 bp island + 5 kb flank.
FLANK, ISLAND_LEN = 5000, 50768
ISLAND_START, ISLAND_END, REF_LEN = FLANK + 1, FLANK + ISLAND_LEN, ISLAND_LEN + 2 * FLANK
REF = "NC_017628.1_pks_island_flank5000"


def rule(**overrides):
    kwargs = dict(island_start=ISLAND_START, island_end=ISLAND_END, reference_length=REF_LEN)
    kwargs.update(overrides)
    return recover.RecruitmentRule(**kwargs)


def classify(pos, cigar, mapq=42, nm=0, **overrides):
    return recover.classify_alignment(pos, cigar, mapq, nm, rule(**overrides))


def sam_record(name, pos, cigar, mapq=42, nm=0, flag=0):
    length = sum(int(n) for n, op in re.findall(r"(\d+)([MIS=X])", cigar))
    return (f"{name}\t{flag}\t{REF}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{'A' * length}\t"
            f"{'I' * length}\tAS:i:100\tNM:i:{nm}\n")


def write_sam(path, records, length=REF_LEN):
    path.write_text(f"@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:{REF}\tLN:{length}\n" + "".join(records))


class InternalReadsNeedQueryCoverage(unittest.TestCase):
    """A1: the acceptance test and its neighbours."""

    def test_150bp_read_with_60bp_internal_local_match_fails(self):
        # 60 perfect aligned bases, 90 soft-clipped, well inside the island: passes the
        # old >=60 bp / >=90% rule, fails the new query-coverage rule (0.40 < 0.80).
        self.assertEqual(classify(20000, "60M90S"), ("rejected", "low_query_coverage"))
        self.assertEqual(classify(20000, "90S60M"), ("rejected", "low_query_coverage"))
        self.assertEqual(classify(20000, "45S60M45S"), ("rejected", "low_query_coverage"))

    def test_the_old_rule_would_have_accepted_it(self):
        # Same record under the pre-A1 rule (no coverage floor) -- the gap A1 closes.
        self.assertEqual(classify(20000, "60M90S", min_query_coverage=0.01),
                         ("island_internal", None))

    def test_fully_aligned_internal_read_passes(self):
        self.assertEqual(classify(20000, "150M"), ("island_internal", None))
        self.assertEqual(classify(20000, "10S140M"), ("island_internal", None))   # 0.93

    def test_coverage_counts_insertions_as_query_in_the_alignment(self):
        self.assertEqual(classify(20000, "70M5I70M5S", nm=5), ("island_internal", None))

    def test_the_threshold_is_the_parameter(self):
        self.assertEqual(classify(20000, "120M30S"), ("island_internal", None))   # 0.80
        self.assertEqual(classify(20000, "119M31S"), ("rejected", "low_query_coverage"))
        self.assertEqual(classify(20000, "119M31S", min_query_coverage=0.75),
                         ("island_internal", None))

    def test_identity_and_length_rules_still_apply(self):
        self.assertEqual(classify(20000, "150M", nm=20), ("rejected", "low_identity"))
        self.assertEqual(classify(20000, "59M"), ("rejected", "short_alignment"))


class BoundaryRule(unittest.TestCase):
    """A1: boundary-spanning reads pass with lower query coverage."""

    def test_read_crossing_the_left_boundary(self):
        # 4926..5075: 75 bp of flank, 75 bp of island.
        self.assertEqual(classify(ISLAND_START - 75, "150M"), ("boundary_spanning", None))

    def test_read_crossing_the_right_boundary(self):
        self.assertEqual(classify(ISLAND_END - 68, "150M"), ("boundary_spanning", None))

    def test_crossing_read_needs_no_query_coverage(self):
        # 70 aligned across the boundary, 80 clipped: coverage 0.47, still boundary.
        self.assertEqual(classify(ISLAND_START - 35, "70M80S"), ("boundary_spanning", None))

    def test_island_end_reached_with_the_clip_pointing_out(self):
        # A true junction read where the host locus differs from IHE3034's flank: the
        # island-side alignment stops at the boundary, the rest of the read is clipped.
        self.assertEqual(classify(ISLAND_START, "90S60M"), ("boundary_spanning", None))
        self.assertEqual(classify(ISLAND_END - 59, "60M90S"), ("boundary_spanning", None))

    def test_within_tolerance_of_the_boundary(self):
        self.assertEqual(classify(ISLAND_START + 10, "90S60M"), ("boundary_spanning", None))
        self.assertEqual(classify(ISLAND_START - 5, "90S60M"), ("boundary_spanning", None))
        # 11 bp in from the boundary is internal: the clip is not at the boundary.
        self.assertEqual(classify(ISLAND_START + 11, "90S60M"), ("rejected", "low_query_coverage"))
        self.assertEqual(classify(ISLAND_START + 11, "90S60M", boundary_tolerance_bp=20),
                         ("boundary_spanning", None))

    def test_clip_pointing_into_the_island_is_not_boundary_evidence(self):
        # Right soft clip at the LEFT boundary: the island end is not what was reached.
        self.assertEqual(classify(ISLAND_START, "60M90S"), ("rejected", "low_query_coverage"))

    def test_a_short_clip_is_not_boundary_evidence(self):
        self.assertEqual(classify(ISLAND_START, "5S145M"), ("island_internal", None))

    def test_overhang_within_tolerance_is_not_a_crossing(self):
        # 6 bp of flank, 144 of island: an internal read, judged on coverage.
        self.assertEqual(classify(ISLAND_START - 6, "150M"), ("island_internal", None))


class FlankOnly(unittest.TestCase):
    def test_fully_aligned_flank_read(self):
        self.assertEqual(classify(100, "150M"), ("flank_only", None))
        self.assertEqual(classify(ISLAND_END + 200, "150M"), ("flank_only", None))

    def test_flank_read_clipped_into_the_island_is_never_island_support(self):
        # What a pks-negative relative sharing the flank looks like at the junction.
        self.assertEqual(classify(ISLAND_START - 60, "60M90S"), ("rejected", "low_query_coverage"))
        self.assertEqual(classify(ISLAND_START - 145, "145M5S"), ("flank_only", None))

    def test_touching_the_island_by_less_than_the_tolerance_is_still_flank(self):
        self.assertEqual(classify(ISLAND_START - 145, "150M"), ("flank_only", None))


class MapqFilter(unittest.TestCase):
    def test_low_mapq_is_rejected_before_anything_else(self):
        self.assertEqual(classify(20000, "150M", mapq=5), ("rejected", "low_mapq"))
        self.assertEqual(classify(ISLAND_START - 75, "150M", mapq=19), ("rejected", "low_mapq"))
        self.assertEqual(classify(20000, "150M", mapq=20), ("island_internal", None))

    def test_the_threshold_is_the_parameter(self):
        self.assertEqual(classify(20000, "150M", mapq=5, min_mapq=0), ("island_internal", None))

    def test_secondary_and_supplementary_records_are_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            sam = Path(tmp) / "r.sam"
            write_sam(sam, [sam_record("a/1", 20000, "150M", flag=256),
                            sam_record("b/1", 20000, "150M", flag=2048)])
            best, reads = recover.classify_sam(sam, rule())
        self.assertEqual(best, {})
        self.assertEqual(reads["aligned_reads"], 0)


class FragmentClasses(unittest.TestCase):
    def test_fragment_takes_the_best_class_of_its_mates(self):
        with tempfile.TemporaryDirectory() as tmp:
            sam = Path(tmp) / "r.sam"
            write_sam(sam, [
                sam_record("f1/1", 100, "150M"), sam_record("f1/2", 20000, "150M"),     # flank + island
                sam_record("f2/1", 100, "150M"), sam_record("f2/2", ISLAND_START, "90S60M"),  # flank + boundary
                sam_record("f3/1", ISLAND_START - 75, "150M"), sam_record("f3/2", 20000, "150M"),
                sam_record("f4/1", 20000, "60M90S"), sam_record("f4/2", 300, "150M"),  # rejected + flank
                sam_record("f5/1", 20000, "60M90S"),                                   # rejected only
            ])
            best, _reads = recover.classify_sam(sam, rule())
        self.assertEqual(best, {"f1": "island_internal", "f2": "boundary_spanning",
                                "f3": "island_internal", "f4": "flank_only", "f5": "rejected"})

    def test_reference_length_mismatch_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            sam = Path(tmp) / "r.sam"
            write_sam(sam, [sam_record("a/1", 20000, "150M")], length=REF_LEN + 2000)
            with self.assertRaises(SystemExit):
                recover.classify_sam(sam, rule())

    def test_rule_rejects_nonsense(self):
        for bad in (dict(island_start=0), dict(island_end=REF_LEN + 1),
                    dict(min_query_coverage=0), dict(min_query_coverage=1.5),
                    dict(boundary_tolerance_bp=-1), dict(min_mapq=-1)):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                rule(**bad)


def fastq(records):
    return "".join(f"@{name}\n{'A' * 150}\n+\n{'I' * 150}\n" for name in records)


def run_cli(tmp, sam_records, names, *extra):
    tmp = Path(tmp)
    sam, reads = tmp / "r.sam", tmp / "reads.fastq.gz"
    write_sam(sam, sam_records)
    with gzip.open(reads, "wt") as h:
        h.write(fastq(names))
    out = {k: tmp / f"{k}.fastq.gz" for k in ("r1", "r2", "single")}
    stats = tmp / "stats.tsv"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--sam", str(sam), "--fastq", str(reads),
         "--r1", str(out["r1"]), "--r2", str(out["r2"]), "--single", str(out["single"]),
         "--stats", str(stats), "--island-start", str(ISLAND_START), "--island-end", str(ISLAND_END),
         "--reference-length", str(REF_LEN), *extra],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    with stats.open() as h:
        metrics = {r["metric"]: r["value"] for r in csv.DictReader(h, delimiter="\t")}
    with gzip.open(out["r1"], "rt") as h:
        pairs = sum(1 for line in h if line.startswith("@"))
    return metrics, pairs


class PksNegativeFlankOnlyFixture(unittest.TestCase):
    """A2 acceptance: a relative sharing the IHE3034 flanks but carrying no island."""

    RECORDS = [
        sam_record("n1/1", 120, "150M"), sam_record("n1/2", 380, "150M"),
        sam_record("n2/1", 2500, "150M"), sam_record("n2/2", 2760, "150M"),
        sam_record("n3/1", ISLAND_END + 500, "150M"), sam_record("n3/2", ISLAND_END + 760, "150M"),
        # the relative's own junction: flank aligned, the rest clipped where the island isn't
        sam_record("n4/1", ISLAND_START - 60, "60M90S"), sam_record("n4/2", ISLAND_START - 400, "150M"),
        sam_record("n5/1", ISLAND_END + 1, "90S60M"), sam_record("n5/2", ISLAND_END + 300, "150M"),
    ]
    NAMES = [f"n{i}/{m}" for i in range(1, 6) for m in (1, 2)]

    def test_flank_only_recruitment_but_zero_island_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            metrics, pairs = run_cli(tmp, self.RECORDS, self.NAMES)
        self.assertEqual(metrics["fragments_flank_only"], "5")
        self.assertEqual(metrics["island_supporting_fragments"], "0")
        self.assertEqual(metrics["fragments_island_internal"], "0")
        self.assertEqual(metrics["fragments_boundary_spanning"], "0")
        self.assertEqual(metrics["reads_flank_only"], "8")
        self.assertEqual(metrics["reads_rejected"], "2")
        self.assertEqual(metrics["reads_rejected_low_query_coverage"], "2")
        # flank-only fragments still go to the assemblers, for context
        self.assertEqual(metrics["recruited_fragment_ids"], "5")
        self.assertEqual(metrics["flank_only_fragments_assembled"], "5")
        self.assertEqual(pairs, 5)

    def test_flank_only_can_be_left_out_of_assembly(self):
        with tempfile.TemporaryDirectory() as tmp:
            metrics, pairs = run_cli(tmp, self.RECORDS, self.NAMES, "--exclude-flank-only")
        self.assertEqual(metrics["fragments_flank_only"], "5")
        self.assertEqual(metrics["recruited_fragment_ids"], "0")
        self.assertEqual(metrics["flank_only_fragments_assembled"], "0")
        self.assertEqual(metrics["include_flank_only"], "false")
        self.assertEqual(pairs, 0)

    def test_the_diagnostic_reads_island_support_only(self):
        summarize = load("summarize_a2", ROOT / "scripts/summarize_tumor_pks_contigs.py")
        with tempfile.TemporaryDirectory() as tmp:
            metrics, _ = run_cli(tmp, self.RECORDS, self.NAMES)
        self.assertEqual(summarize.assembly_input_diagnostic(metrics, 0.0),
                         "no_island_supporting_fragments")


class PksPositiveFixture(unittest.TestCase):
    def test_classes_and_stats_row_names(self):
        records = [sam_record("p1/1", 20000, "150M"), sam_record("p1/2", 20300, "150M"),
                   sam_record("p2/1", ISLAND_START - 75, "150M"), sam_record("p2/2", ISLAND_START + 200, "150M"),
                   sam_record("p3/1", ISLAND_END - 59, "60M90S"), sam_record("p3/2", 30000, "60M90S"),
                   sam_record("p4/1", 100, "150M"), sam_record("p4/2", 400, "150M"),
                   sam_record("p5/1", 20000, "150M", mapq=3)]
        names = [f"p{i}/{m}" for i in range(1, 5) for m in (1, 2)] + ["p5/1"]
        with tempfile.TemporaryDirectory() as tmp:
            metrics, pairs = run_cli(tmp, records, names)
        self.assertEqual(metrics["fragments_island_internal"], "2")
        self.assertEqual(metrics["fragments_boundary_spanning"], "1")
        self.assertEqual(metrics["fragments_flank_only"], "1")
        self.assertEqual(metrics["fragments_rejected"], "1")
        self.assertEqual(metrics["island_supporting_fragments"], "3")
        self.assertEqual(metrics["reads_rejected_low_mapq"], "1")
        self.assertEqual(metrics["recruited_fragment_ids"], "4")
        self.assertEqual(metrics["island_start_in_reference"], str(ISLAND_START))
        self.assertEqual(metrics["island_end_in_reference"], str(ISLAND_END))
        self.assertEqual(metrics["min_query_coverage"], "0.8")
        self.assertEqual(metrics["min_mapq"], "20")
        self.assertEqual(metrics["boundary_tolerance_bp"], "10")
        self.assertEqual(pairs, 4)


class AssemblyDiagnostic(unittest.TestCase):
    summarize = load("summarize_a2_diag", ROOT / "scripts/summarize_tumor_pks_contigs.py")

    def test_values(self):
        d = self.summarize.assembly_input_diagnostic
        self.assertEqual(d({"island_supporting_fragments": "0", "fragments_flank_only": "40"}, 0),
                         "no_island_supporting_fragments")
        self.assertEqual(d({"island_supporting_fragments": "12"}, 0), "island_fragments_no_contig")
        self.assertEqual(d({"island_supporting_fragments": "12"}, 0.4), "island_contigs_recovered")
        self.assertEqual(d({"recruited_fragment_ids": "9"}, 0), "recruitment_classes_unavailable")

    def test_final_evidence_carries_the_classes(self):
        src = (ROOT / "scripts/summarize_tumor_pks_contigs.py").read_text()
        for name in ("island_supporting_fragments", "island_internal_fragments",
                     "boundary_spanning_fragments", "flank_only_fragments",
                     "assembly_input_diagnostic"):
            self.assertIn(f'"{name}"', src)

    def test_master_summary_and_report_carry_them(self):
        master = load("master_a2", ROOT / "scripts/build_master_summary.py")
        for name in ("island_supporting_fragments", "boundary_spanning_fragments",
                     "flank_only_fragments", "assembly_input_diagnostic"):
            self.assertIn(name, master.CONTIG_COLUMNS)
        report = load("report_a2", ROOT / "scripts/build_cohort_report.py")
        tip = report.recruitment_classes_tip({"island_supporting_fragments": "0",
                                              "boundary_spanning_fragments": "0",
                                              "flank_only_fragments": "7",
                                              "assembly_input_diagnostic": "no_island_supporting_fragments"})
        self.assertIn("0 island-supporting", tip)
        self.assertIn("7 flank-only", tip)
        self.assertIn("no_island_supporting_fragments", tip)
        self.assertEqual(report.recruitment_classes_tip({}), "")


class PipelineWiring(unittest.TestCase):
    def recruit_process(self):
        return TARGETED.split("process targetedPksRecruit", 1)[1].split("\nprocess ", 1)[0]

    def test_bowtie2_uses_best_alignment_reporting(self):
        command = next(l for l in self.recruit_process().splitlines() if l.strip().startswith("bowtie2"))
        self.assertNotIn("-k 1", command)
        self.assertNotRegex(command, r"\s-k\s")
        self.assertIn("--very-sensitive-local", command)

    def test_the_script_gets_every_rule_parameter(self):
        command = next(l for l in self.recruit_process().splitlines()
                       if l.strip().startswith("python") and "recover_recruited_mates.py" in l)
        for flag, param in (("--min-query-coverage", "recruit_min_query_coverage"),
                            ("--min-mapq", "recruit_min_mapq"),
                            ("--boundary-tolerance-bp", "recruit_boundary_tolerance_bp"),
                            ("--min-aligned-bases", "pks_recruit_min_aligned"),
                            ("--min-identity", "pks_recruit_min_identity")):
            self.assertIn(f'{flag} "${{params.{param}}}"', command)
        for flag in ("--island-start", "--island-end", "--reference-length"):
            self.assertIn(flag, command)
        self.assertIn("params.recruit_include_flank_only.toString().toBoolean()", self.recruit_process())

    def test_defaults(self):
        for line in ("params.recruit_min_query_coverage    = 0.80",
                     "params.recruit_boundary_tolerance_bp = 10",
                     "params.recruit_min_mapq              = 20",
                     "params.recruit_include_flank_only    = true",
                     "params.pks_recruit_flank_bp          = 5000",
                     "params.pks_recruit_min_aligned  = 60",
                     "params.pks_recruit_min_identity = 0.90"):
            self.assertIn(line, MAIN)

    def test_island_placement_matches_the_shipped_index(self):
        # The module computes flank+1 .. flank+island_len; the public index FASTA's own
        # header must say the same island and flank, and its length must be 2*flank+island.
        fasta = (ROOT / "indices/pks_island/pks_island.fna").read_text().splitlines()
        header = fasta[0]
        m = re.search(r"island=(\d+)-(\d+) flank=(\d+)", header)
        start = int(re.search(r"params\.pks_start_1based\s*=\s*(\d+)", MAIN).group(1))
        end = int(re.search(r"params\.pks_end_1based\s*=\s*(\d+)", MAIN).group(1))
        self.assertEqual((int(m.group(1)), int(m.group(2)), int(m.group(3))), (start, end, FLANK))
        self.assertEqual(sum(len(l) for l in fasta[1:]), REF_LEN)
        self.assertEqual(end - start + 1, ISLAND_LEN)
        for needle in ("def recruit_island_start = recruit_flank + 1",
                       "def recruit_island_end = recruit_flank + params.pks_island_len_1based.toString().toInteger()",
                       "def recruit_reference_length = recruit_island_end + recruit_flank"):
            self.assertIn(needle, TARGETED)

    def test_main_validates_the_parameters(self):
        for message in ("--recruit_min_query_coverage must be a fraction in (0,1]",
                        "--recruit_boundary_tolerance_bp must be less than half of --pks_recruit_min_aligned",
                        "--recruit_include_flank_only must be true or false",
                        "recruitment classes would be misplaced"):
            self.assertIn(message, MAIN)
        self.assertIn("recruit_boundary_tolerance_bp: params.recruit_boundary_tolerance_bp", MAIN)
        self.assertIn("recruit_min_mapq: params.recruit_min_mapq", MAIN)

    def test_run_report_records_them(self):
        prov = load("finalize_provenance_a1", ROOT / "scripts/finalize_provenance.py")
        names = {n for _group, members in prov.PARAMETER_GROUPS for n in members}
        for name in ("pks_recruit_min_aligned", "pks_recruit_min_identity",
                     "recruit_min_query_coverage", "recruit_boundary_tolerance_bp",
                     "recruit_min_mapq", "recruit_include_flank_only", "pks_recruit_flank_bp"):
            self.assertIn(name, names)

    def test_docs_describe_them(self):
        docs = (ROOT / "docs/running/assembly.md").read_text()
        for name in ("recruit_min_query_coverage", "recruit_boundary_tolerance_bp",
                     "recruit_min_mapq", "recruit_include_flank_only", "island_internal",
                     "boundary_spanning", "flank_only", "island_supporting_fragments",
                     "assembly_input_diagnostic"):
            self.assertIn(name, docs)
        output = (ROOT / "docs/output.md").read_text()
        self.assertIn("island_supporting_fragments", output)


if __name__ == "__main__":
    unittest.main()
