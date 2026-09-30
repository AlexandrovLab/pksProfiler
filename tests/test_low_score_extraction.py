"""E1 (Ludmil, 25 Sep): extraction also takes primary mapped records with AS:i < 30.

GDC aligned TCGA with `bwa mem -T 0` (read from the real @PG CL); Mutographs with -T 30,
bwa's default. bwa reports a read as unmapped only when its best local score is below
-T, so under -T 0 about a third of genuine pks reads come out "mapped" to a short
spurious hit, and `samtools fastq -f 4` never extracts them. AS < 30 among the mapped
primaries is exactly the set -T 30 would have left unmapped. Validated in
projects/pksprofiler/validation/E1 (synthetic: 99.9-100% of the hidden reads recovered;
20-sample real pilot: tumours gained clb evidence, matched normals gained none).

The extraction tests below run the on-branch shell pipeline exactly as
Modules/extract_reads.nf renders it (the Groovy escapes undone, the three interpolated
values substituted) through real samtools on hand-written SAM records. Every record
is synthetic: invented names, invented sequence, one invented contig.
"""
import csv
import gzip
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXTRACT = (ROOT / "Modules/extract_reads.nf").read_text()
MAIN = (ROOT / "main.nf").read_text()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


hp = load("header_provenance_e1", ROOT / "scripts/header_provenance.py")
bqs = load("build_qc_summary_e1", ROOT / "scripts/build_qc_summary.py")
fp = load("finalize_provenance_e1", ROOT / "scripts/finalize_provenance.py")


def find_samtools():
    """samtools with -e expressions and exists() (>= 1.16). Beside this python first."""
    for candidate in (Path(sys.executable).parent / "samtools", shutil.which("samtools")):
        if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            out = subprocess.run([str(candidate), "--version"], capture_output=True, text=True).stdout
            m = re.match(r"samtools (\d+)\.(\d+)", out)
            if m and (int(m.group(1)), int(m.group(2))) >= (1, 16):
                return str(candidate)
    return None


SAMTOOLS = find_samtools()

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def seq(i, n=60):
    """Deterministic, distinct invented sequence."""
    bases = "ACGT"
    return "".join(bases[(i * 7 + k * (i + 3) + k // 5) % 4] for k in range(n))


QUAL = "I" * 59 + "#"   # asymmetric, so a reversed QUAL is visible

# (name, flag, mapped?, AS or None, original read sequence, expected class)
# flag bits: 1 paired, 4 unmapped, 16 reverse, 64 first, 128 second, 256 secondary, 2048 supplementary
RECORDS = [
    ("unm_pair", 77, False, None, seq(1), "flag4"),        # both unmapped, /1
    ("unm_pair", 141, False, None, seq(2), "flag4"),       # /2
    ("low_fwd", 65 + 32, True, 29, seq(3), "low"),        # AS 29: kept, /1
    ("low_fwd", 145, True, 0, seq(4), "low"),             # AS 0, reverse strand: kept, /2
    ("at_threshold", 99, True, 30, seq(5), "excluded"),    # AS 30: not below 30
    ("high", 65, True, 58, seq(6), "excluded"),            # AS 58
    ("no_as", 65, True, None, seq(7), "no_as"),            # mapped, no AS: never added, counted
    ("secondary", 65 + 256, True, 5, seq(8), "excluded"),  # secondary: never
    ("supplementary", 65 + 2048, True, 5, seq(9), "excluded"),  # supplementary: never
    ("unm_mate_mapped", 133 + 32, False, None, seq(10), "flag4"),  # FLAG 4 with a mapped mate
    ("unpaired_low", 16, True, 12, seq(11), "low"),        # unpaired, reverse, AS 12: no suffix
    ("unpaired_unm", 4, False, None, seq(12), "flag4"),    # unpaired FLAG 4: READ_OTHER
]


def sam_text(records=RECORDS, pg=None):
    lines = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:synthetic_chr1\tLN:5000"]
    if pg:
        lines.append(pg)
    for name, flag, mapped, score, original, _ in records:
        stored = revcomp(original) if flag & 16 else original
        qual = QUAL[::-1] if flag & 16 else QUAL
        rname, pos, cigar, mapq = ("synthetic_chr1", 100, "20M40S", 0) if mapped else ("*", 0, "*", 0)
        if not mapped and flag & 1 and not flag & 8:
            rname, pos = "synthetic_chr1", 100   # placed at its mapped mate, as bwa does
        fields = [name, str(flag), rname, str(pos), str(mapq), cigar, "*", "0", "0", stored, qual]
        if score is not None:
            fields += ["NM:i:0", f"AS:i:{score}", "XS:i:0"]
        lines.append("\t".join(fields))
    return "\n".join(lines) + "\n"


def expected_fastq(records=RECORDS, classes=("flag4", "low")):
    """What samtools fastq -N writes: /1 /2 from the mate bits, reverse strand undone."""
    out = []
    for name, flag, _mapped, _score, original, cls in records:
        if cls not in classes:
            continue
        suffix = "/1" if flag & 64 and not flag & 128 else "/2" if flag & 128 and not flag & 64 else ""
        out.append((name + suffix, original, QUAL))
    return sorted(out)


def on_branch_script(min_as=30):
    """The fifo counter + the on-branch pipeline + the class check, as bash would see it."""
    start = EXTRACT.index("    mkfifo extract.count.fifo")
    on = EXTRACT.index('    if [[ "\\$LOW_SCORE" == "on" ]]; then', start)
    off = EXTRACT.index("\n    else\n", on)
    wait = EXTRACT.index('		wait "\\$EXTRACT_COUNTER"', off)
    classes_end = EXTRACT.index("\n\t# F09", wait)
    block = (EXTRACT[start:on] + EXTRACT[on:off] + "\n    fi\n" + EXTRACT[wait:classes_end])
    block = block.replace("\\\\", "\x00").replace("\\$", "$").replace("\x00", "\\")
    for groovy, value in (("${task.cpus}", "1"), ('${alignment}', '$ALN'),
                          ("${params.extract_min_as}", str(min_as)), ("${sampleID}", "S")):
        block = block.replace(groovy, value)
    assert "${" not in block.replace("${REFERENCE_ARGS[@]}", ""), block
    return ("set -euo pipefail\nLOW_SCORE=on\nREFERENCE_ARGS=()\nREADS=out.fastq.gz\n" + block +
            '\nprintf "%s\\t%s\\t%s\\t%s\\n" "$UNMAPPED_READS" "$LOW_SCORE_READS" "$NO_AS_RECORDS" "$EXTRACTED_READS" > result.tsv\n')


@unittest.skipUnless(SAMTOOLS, "samtools >= 1.16 not found beside python or on PATH")
class OnBranchOnSyntheticSam(unittest.TestCase):
    def run_branch(self, sam, min_as=30):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        (tmp / "in.sam").write_text(sam)
        env = dict(os.environ, ALN=str(tmp / "in.sam"),
                   PATH=f"{Path(SAMTOOLS).parent}:{os.environ.get('PATH', '')}")
        (tmp / "run.sh").write_text(on_branch_script(min_as))
        result = subprocess.run(["bash", "run.sh"], cwd=tmp, env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        unmapped, low, no_as, total = (tmp / "result.tsv").read_text().split()
        with gzip.open(tmp / "out.fastq.gz", "rt") as fh:
            lines = fh.read().splitlines()
        reads = sorted((lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4))
        return int(unmapped), int(low), int(no_as), int(total), reads, result.stderr

    def test_the_rule_keeps_exactly_flag4_and_primary_mapped_as_below_30(self):
        unmapped, low, no_as, total, reads, stderr = self.run_branch(sam_text())
        self.assertEqual(reads, expected_fastq())
        self.assertEqual((unmapped, low, no_as, total), (4, 3, 1, 7))
        self.assertIn("1 mapped primary records had no AS tag and were not added", stderr)

    def test_as_at_or_above_the_threshold_is_excluded(self):
        names = {r[0].split("/")[0] for r in self.run_branch(sam_text())[4]}
        self.assertNotIn("at_threshold", names)
        self.assertNotIn("high", names)

    def test_secondary_and_supplementary_never_extracted(self):
        names = {r[0].split("/")[0] for r in self.run_branch(sam_text())[4]}
        self.assertNotIn("secondary", names)
        self.assertNotIn("supplementary", names)

    def test_missing_as_is_counted_not_added(self):
        unmapped, low, no_as, total, reads, _ = self.run_branch(sam_text())
        self.assertEqual(no_as, 1)
        self.assertNotIn("no_as", {r[0].split("/")[0] for r in reads})

    def test_orientation_and_names_match_samtools_fastq_minus_n(self):
        reads = dict((n, (s, q)) for n, s, q in self.run_branch(sam_text())[4])
        self.assertEqual(reads["low_fwd/2"], (seq(4), QUAL))        # reverse strand undone
        self.assertEqual(reads["unpaired_low"], (seq(11), QUAL))    # no mate bit, no suffix
        self.assertEqual(reads["low_fwd/1"], (seq(3), QUAL))

    def test_flag4_set_is_what_the_off_branch_extracts(self):
        """extracted_unmapped_reads keeps its meaning: the FLAG 4 reads, same FASTQ records."""
        unmapped, _, _, _, reads, _ = self.run_branch(sam_text())
        flag4 = expected_fastq(classes=("flag4",))
        self.assertEqual(unmapped, len(flag4))
        self.assertTrue(set(flag4) <= set(reads))
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        (tmp / "in.sam").write_text(sam_text())
        off = subprocess.run([SAMTOOLS, "fastq", "-N", "-f", "4", "-F", "2304", str(tmp / "in.sam")],
                             capture_output=True, text=True, check=True).stdout.splitlines()
        off_reads = sorted((off[i][1:], off[i + 1], off[i + 3]) for i in range(0, len(off), 4))
        self.assertEqual(off_reads, flag4)

    def test_the_threshold_is_the_parameter(self):
        unmapped, low, no_as, total, reads, _ = self.run_branch(sam_text(), min_as=1)
        self.assertEqual(low, 1)   # only AS 0
        self.assertEqual(total, unmapped + low)

    def test_no_records_at_all(self):
        unmapped, low, no_as, total, reads, _ = self.run_branch(sam_text(records=[]))
        self.assertEqual((unmapped, low, no_as, total, reads), (0, 0, 0, 0, []))


PG_GDC = "@PG\tID:bwa\tPN:bwa\tVN:0.7.15-r1140\tCL:bwa mem -t 36 -T 0 -R @RG\\tID:SECRET_RG\\tSM:SECRET_SAMPLE /ref/GRCh38.d1.vd1.fa SECRET_1.fq SECRET_2.fq"
PG_SANGER = "@PG\tID:bwa_mem\tPN:bwa\tVN:0.7.17\tCL:/opt/bin/bwa mem -Y -K 100000000 -t 8 -T30 ref.fa in.fq"
PG_DEFAULT = "@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:2.2.1\tCL:bwa-mem2 mem -t 8 ref.fa a.fq b.fq"
PG_STAR = "@PG\tID:STAR\tPN:STAR\tVN:2.7.5c\tCL:STAR --runThreadN 8 --genomeDir /ref"
PG_ALN = "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa sampe ref.fa a.sai b.sai a.fq b.fq"
PG_NO_CL = "@PG\tID:bwa\tPN:bwa\tVN:0.7.17"
PG_MM2 = "@PG\tID:minimap2\tPN:minimap2\tVN:2.28\tCL:minimap2 -ax sr ref.mmi a.fq"


def header(*pg):
    return "\n".join(["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:248956422"] + list(pg)) + "\n"


class AutoDetection(unittest.TestCase):
    def test_bwa_mem_turns_auto_on_and_reads_minus_t(self):
        row = hp.parse_header(header(PG_GDC, "@PG\tID:samtools\tPN:samtools\tPP:bwa\tVN:1.3"), "A")
        self.assertEqual((row["aligner_bwa_mem"], row["bwa_min_score"]), ("yes", "0"))
        self.assertEqual(row["extract_low_score_mapped"], "on_auto_bwa_mem")
        self.assertEqual(row["extract_min_as"], 30)
        for value in row.values():   # the CL is read, never recorded
            self.assertNotIn("SECRET", str(value))

    def test_minus_t_forms_and_bwa_mem2(self):
        self.assertEqual(hp.bwa_mem_evidence(header(PG_SANGER)), ("yes", "30"))
        self.assertEqual(hp.bwa_mem_evidence(header(PG_DEFAULT)), ("yes", "default"))
        self.assertEqual(hp.bwa_mem_evidence(header("@PG\tID:x\tPN:sentieon\tCL:sentieon bwa mem -T 0 r.fa")),
                         ("yes", "0"))

    def test_star_minimap2_and_no_pg_leave_auto_off(self):
        for pg in ((PG_STAR,), (PG_MM2,), ()):
            with self.subTest(pg=pg):
                row = hp.parse_header(header(*pg), "A")
                self.assertEqual(row["aligner_bwa_mem"], "no")
                self.assertEqual(row["extract_low_score_mapped"], "off_auto_not_bwa_mem")
                self.assertEqual(row["extract_min_as"], "NA")

    def test_bwa_without_a_mem_command_line_is_unconfirmed_and_off(self):
        for pg in (PG_ALN, PG_NO_CL):
            with self.subTest(pg=pg):
                row = hp.parse_header(header(pg), "A")
                self.assertEqual(row["aligner_bwa_mem"], "unconfirmed")
                self.assertEqual(row["extract_low_score_mapped"], "off_auto_not_bwa_mem")

    def test_explicit_on_and_off_override_detection(self):
        self.assertEqual(hp.parse_header(header(PG_STAR), "A", "on")["extract_low_score_mapped"], "on_requested")
        self.assertEqual(hp.parse_header(header(PG_GDC), "A", "off")["extract_low_score_mapped"], "off_requested")
        with self.assertRaises(ValueError):
            hp.parse_header(header(PG_GDC), "A", "maybe")

    def test_every_state_is_a_qc_text_metric_value(self):
        self.assertEqual(set(hp.LOW_SCORE_STATES), bqs.TEXT_METRICS["extract_low_score_mode"])

    def test_the_decision_file_the_shell_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            for pg, first, state in ((PG_GDC, "on", "on_auto_bwa_mem"), (PG_STAR, "off", "off_auto_not_bwa_mem"),
                                     (PG_SANGER, "on", "on_auto_bwa_mem")):
                (d / "h.sam").write_text(header(pg))
                r = subprocess.run([sys.executable, str(ROOT / "scripts/header_provenance.py"), "sample",
                                    "--header", str(d / "h.sam"), "--sample", "A", "--output", str(d / "o.tsv"),
                                    "--low-score-mode", "auto", "--min-as", "30", "--decision", str(d / "dec")],
                                   capture_output=True, text=True)
                self.assertEqual(r.returncode, 0, r.stderr)
                fields = (d / "dec").read_text().rstrip("\n").split("\t")
                self.assertEqual(fields[:2], [first, state])
                self.assertNotIn("SECRET", fields[2])
                if pg is PG_SANGER:
                    self.assertIn("none should be added", fields[2])
            r = subprocess.run([sys.executable, str(ROOT / "scripts/header_provenance.py"), "sample",
                                "--header", str(d / "h.sam"), "--sample", "A", "--output", str(d / "o.tsv"),
                                "--min-as", "0"], capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)

    def test_summary_counts_states_and_run_report_renders_them(self):
        rows = [hp.parse_header(header(PG_GDC), "A"), hp.parse_header(header(PG_GDC), "B"),
                hp.parse_header(header(PG_STAR), "C")]
        summary = hp.summarize(rows)
        states = {r["key"]: r["n_files"] for r in summary if r["section"] == "extract_low_score"}
        self.assertEqual(states, {"on_auto_bwa_mem": 2, "off_auto_not_bwa_mem": 1})
        report = fp.render({"status": "completed", "code": {}, "invocation": {"params": {}}, "run": {},
                            "input_headers": [{k: str(v) for k, v in r.items()} for r in summary]})
        self.assertIn("low-score mapped extraction", report)
        self.assertIn("on_auto_bwa_mem", report)


class QcSummary(unittest.TestCase):
    def row(self, **values):
        base = {c: None for c in bqs.OUTPUT_COLUMNS}
        base.update(values)
        return base

    def test_columns_exist(self):
        for column in ("extracted_unmapped_reads", "extracted_low_score_mapped_reads", "extracted_reads",
                       "extract_low_score_no_as_records", "extract_low_score_mode"):
            self.assertIn(column, bqs.OUTPUT_COLUMNS)

    def test_extracted_reads_is_the_chain_step_so_filter_input_may_exceed_flag4(self):
        out = bqs.output_row("S", self.row(input_alignment_records=1000, extracted_unmapped_reads=10,
                                           extracted_low_score_mapped_reads=5, extracted_reads=15,
                                           filter_input_reads=15, reads_after_fastp=14))
        self.assertEqual((out["extracted_reads"], out["filter_input_reads"]), (15, 15))

    def test_parts_must_add_up(self):
        with self.assertRaisesRegex(ValueError, "extracted_unmapped_reads \\+ extracted_low_score_mapped_reads"):
            bqs.output_row("S", self.row(extracted_unmapped_reads=10, extracted_low_score_mapped_reads=5,
                                         extracted_reads=16))

    def test_the_total_is_still_bounded_by_the_input_and_bounds_fastp(self):
        with self.assertRaises(ValueError):
            bqs.output_row("S", self.row(input_alignment_records=10, extracted_unmapped_reads=8,
                                         extracted_low_score_mapped_reads=5, extracted_reads=13))
        with self.assertRaises(ValueError):
            bqs.output_row("S", self.row(extracted_unmapped_reads=8, extracted_low_score_mapped_reads=0,
                                         extracted_reads=8, filter_input_reads=9))

    def test_a_pre_e1_fragment_still_chains_on_extracted_unmapped_reads(self):
        with self.assertRaises(ValueError):
            bqs.output_row("S", self.row(extracted_unmapped_reads=8, filter_input_reads=9))

    def test_mode_is_a_checked_text_metric(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = Path(tmp) / "good.tsv"
            good.write_text("Sample\tMetric\tValue\nS\textract_low_score_mode\ton_auto_bwa_mem\n")
            self.assertEqual(bqs.load_fragments([good])["S"]["extract_low_score_mode"], "on_auto_bwa_mem")
            bad = Path(tmp) / "bad.tsv"
            bad.write_text("Sample\tMetric\tValue\nS\textract_low_score_mode\tyes\n")
            with self.assertRaises(ValueError):
                bqs.load_fragments([bad])

    def test_written_by_extraction(self):
        for metric in ("extracted_unmapped_reads", "extracted_low_score_mapped_reads", "extracted_reads",
                       "extract_low_score_no_as_records", "extract_low_score_mode"):
            self.assertIn(f"\\t{metric}\\t", EXTRACT)


class Wiring(unittest.TestCase):
    def test_params_and_validation(self):
        self.assertIn('params.extract_low_score_mapped = "auto"', MAIN)
        self.assertRegex(MAIN, r"params\.extract_min_as\s+= 30")
        self.assertIn("--extract_low_score_mapped must be auto, on or off", MAIN)
        self.assertIn("--extract_min_as must be an integer >= 1", MAIN)

    def test_recorded_in_run_report_parameters(self):
        names = [n for _t, group in fp.PARAMETER_GROUPS for n in group]
        self.assertIn("extract_low_score_mapped", names)
        self.assertIn("extract_min_as", names)

    def test_decided_from_the_header_before_the_decode(self):
        self.assertIn('--low-score-mode "${params.extract_low_score_mapped}"', EXTRACT)
        self.assertLess(EXTRACT.index("--decision extract_low_score.decision"), EXTRACT.index("mkfifo"))

    def test_one_decode_per_branch(self):
        on = EXTRACT[EXTRACT.index('if [[ "\\$LOW_SCORE" == "on" ]]; then'):EXTRACT.index("\n    else\n")]
        self.assertEqual(on.count('"${alignment}"'), 1)
        self.assertIn("-F 2304", on)
        self.assertIn('"\\${REFERENCE_ARGS[@]}"', on)   # CRAM decodes with --reference

    def test_the_e2e_check_runs_it_for_real(self):
        """The live proof: tests/check_e2e_execution.sh's seventh run (bwa -T 0 style BAM/CRAM)."""
        script = (ROOT / "tests/check_e2e_execution.sh").read_text()
        code = "\n".join(l for l in script.splitlines() if not l.strip().startswith("#"))
        for needle in ("generate_low_score_bam_fixtures.py", "assert_low_score_outputs.py",
                       '"$low_fx/bwa_T0.cram"', '"$low_fx/star.bam"', "--input_data_type auto",
                       "--save_intermediates true"):
            self.assertIn(needle, code)
        self.assertNotIn("--extract_low_score_mapped", code)   # the default, auto, is what runs
        for path in ("tests/fixtures/generate_low_score_bam_fixtures.py",
                     "tests/fixtures/assert_low_score_outputs.py"):
            self.assertTrue((ROOT / path).stat().st_mode & 0o100, path)
        generator = (ROOT / "tests/fixtures/generate_low_score_bam_fixtures.py").read_text()
        self.assertIn("CL:bwa mem -t 8 -T 0", generator)
        self.assertIn("PN:STAR", generator)

    def test_documented(self):
        for path in ("README.md", "docs/running/bam.md", "CHANGELOG.md", "docs/output.md"):
            with self.subTest(path=path):
                text = (ROOT / path).read_text()
                self.assertIn("extract_low_score_mapped", text)
        self.assertIn("-T 0", (ROOT / "docs/running/bam.md").read_text())


if __name__ == "__main__":
    unittest.main()
