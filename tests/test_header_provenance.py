"""E3 (Ludmil, 25 Sep): input reference and upstream-aligner provenance, in the pipeline.

"Hash a canonicalized @SQ SN/LN/M5 list for every BAM/CRAM and capture the @PG chain
(program/version) when present. Store these in run/QC provenance. Acceptance test:
Equivalent BAM/CRAMs using the same reference produce the same SQ digest; different
reference builds produce different digests."

Every header here is synthetic: invented contig lists with the public GRCh37/GRCh38 chr1
lengths, invented M5 values, invented @PG records.
"""
import csv
import datetime
import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/header_provenance.py"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


hp = load("header_provenance", SCRIPT)
fp = load("finalize_provenance", ROOT / "scripts/finalize_provenance.py")

GRCH38 = [("chr1", 248956422, "6aef897c3d6ff0c78aff06ac189178dd"),
          ("chr2", 242193529, "f98db672eb0993dcfdabafe2a882905c"),
          ("chrM", 16569, "c68f52674c9fb33aef52dcf399755519")]
GRCH37 = [("1", 249250621, "1b22b98cdeb4a9304cb5d48026a85128"),
          ("2", 243199373, "a0d9851da00400dec1098a9255ac712e"),
          ("MT", 16569, "c68f52674c9fb33aef52dcf399755519")]
PG_BWA = ["@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem /data/SECRET_SAMPLE_42_R1.fq",
          "@PG\tID:MarkDuplicates\tPN:MarkDuplicates\tPP:bwa\tVN:2.23\tCL:java -jar picard.jar I=SECRET_SAMPLE_42.bam",
          "@PG\tID:samtools\tPN:samtools\tPP:MarkDuplicates\tVN:1.21\tCL:samtools view -C"]


def header(contigs, pg=(), m5=True, rename=None):
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    for sn, ln, md5 in contigs:
        sn = rename(sn) if rename else sn
        lines.append(f"@SQ\tSN:{sn}\tLN:{ln}" + (f"\tM5:{md5}" if m5 else ""))
    lines.append("@RG\tID:rg1\tSM:SECRET_SAMPLE_42\tPL:ILLUMINA")
    return "\n".join(lines + list(pg)) + "\n"


class Fingerprints(unittest.TestCase):
    def test_same_reference_same_digest_whatever_else_differs(self):
        bam = hp.parse_header(header(GRCH38, PG_BWA), "A")
        # a CRAM of the same alignment: extra @PG, same @SQ
        cram = hp.parse_header(header(GRCH38, PG_BWA + ["@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.22"]), "B")
        self.assertEqual(bam["sq_fingerprint"], cram["sq_fingerprint"])
        self.assertEqual(bam["m5_set_fingerprint"], cram["m5_set_fingerprint"])

    def test_different_builds_differ_and_are_named(self):
        a = hp.parse_header(header(GRCH38), "A")
        b = hp.parse_header(header(GRCH37), "B")
        self.assertNotEqual(a["sq_fingerprint"], b["sq_fingerprint"])
        self.assertNotEqual(a["m5_set_fingerprint"], b["m5_set_fingerprint"])
        self.assertEqual((a["build_guess"], a["chr1_name"], a["chr1_length"]), ("GRCh38", "chr1", 248956422))
        self.assertEqual((b["build_guess"], b["chr1_name"]), ("GRCh37", "1"))
        self.assertEqual(a["n_contigs"], 3)

    def test_the_digest_is_md5_over_ordered_sn_ln(self):
        import hashlib
        row = hp.parse_header(header(GRCH38), "A")
        expected = hashlib.md5("\n".join(f"{sn}:{ln}" for sn, ln, _ in GRCH38).encode()).hexdigest()
        self.assertEqual(row["sq_fingerprint"], expected)
        reordered = hp.parse_header(header(list(reversed(GRCH38))), "A")
        self.assertNotEqual(reordered["sq_fingerprint"], row["sq_fingerprint"])
        self.assertEqual(reordered["m5_set_fingerprint"], row["m5_set_fingerprint"])

    def test_m5_set_survives_contig_renaming_and_is_na_without_m5(self):
        named = hp.parse_header(header(GRCH38), "A")
        renamed = hp.parse_header(header(GRCH38, rename=lambda s: s[3:]), "A")
        self.assertNotEqual(named["sq_fingerprint"], renamed["sq_fingerprint"])
        self.assertEqual(named["m5_set_fingerprint"], renamed["m5_set_fingerprint"])
        bare = hp.parse_header(header(GRCH38, m5=False), "A")
        self.assertEqual((bare["m5_set_fingerprint"], bare["n_contigs_with_m5"]), ("NA", 0))

    def test_no_sq_is_refused(self):
        with self.assertRaises(ValueError):
            hp.parse_header("@HD\tVN:1.6\n@PG\tID:x\n", "A")


class AlignerChain(unittest.TestCase):
    def test_chain_follows_pp_not_header_order_and_never_carries_cl(self):
        shuffled = [PG_BWA[2], PG_BWA[0], PG_BWA[1]]
        row = hp.parse_header(header(GRCH38, shuffled), "A")
        self.assertEqual(row["pg_chain"], "bwa|bwa|0.7.17 > MarkDuplicates|MarkDuplicates|2.23 > samtools|samtools|1.21")
        self.assertEqual(row["pg_programs"], "bwa:0.7.17 > MarkDuplicates:2.23 > samtools:1.21")
        self.assertEqual(row["n_pg"], 3)
        for value in row.values():
            self.assertNotIn("SECRET_SAMPLE_42", str(value))
            self.assertNotIn("CL:", str(value))

    def test_no_pg(self):
        row = hp.parse_header(header(GRCH38), "A")
        self.assertEqual((row["n_pg"], row["pg_chain"], row["pg_programs"]), (0, "NA", "NA"))


class CommandLineAndSummary(unittest.TestCase):
    def test_sample_then_summarize_counts_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            tables = []
            for name, contigs, pg in (("SAMPLE_ALPHA", GRCH38, PG_BWA), ("SAMPLE_BETA", GRCH38, PG_BWA),
                                      ("SAMPLE_GAMMA", GRCH37, [])):
                (d / f"{name}.sam").write_text(header(contigs, pg))
                out = d / f"{name}.tsv"
                result = subprocess.run([sys.executable, str(SCRIPT), "sample", "--header", str(d / f"{name}.sam"),
                                         "--sample", name, "--output", str(out)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                tables.append(out)
            (d / "list.txt").write_text("".join(f"{t}\n" for t in tables))
            result = subprocess.run([sys.executable, str(SCRIPT), "summarize", "--inputs-file", str(d / "list.txt"),
                                     "--table", str(d / "all.tsv"), "--summary", str(d / "summary.tsv")],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (d / "all.tsv").open() as fh:
                self.assertEqual([r["sample"] for r in csv.DictReader(fh, delimiter="\t")],
                                 ["SAMPLE_ALPHA", "SAMPLE_BETA", "SAMPLE_GAMMA"])
            summary_text = (d / "summary.tsv").read_text()
            for secret in ("SAMPLE_ALPHA", "SAMPLE_BETA", "SAMPLE_GAMMA", "SECRET_SAMPLE_42", "|"):
                self.assertNotIn(secret, summary_text)
            with (d / "summary.tsv").open() as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            refs = [r for r in rows if r["section"] == "reference"]
            self.assertEqual(sorted(int(r["n_files"]) for r in refs), [1, 2])
            chains = {r["key"]: int(r["n_files"]) for r in rows if r["section"] == "aligner_chain"}
            self.assertEqual(chains, {"bwa:0.7.17 > MarkDuplicates:2.23 > samtools:1.21": 2, "(no @PG)": 1})

            # RUN_REPORT's INPUT ALIGNMENTS section renders it -- if this session wrote it.
            cohort = d / "cohort"
            (cohort / "provenance").mkdir(parents=True)
            (cohort / "provenance" / "input_header_summary.tsv").write_text(summary_text)
            now = datetime.datetime.now().astimezone()
            fresh = fp.input_header_summary(cohort, (now - datetime.timedelta(hours=1)).isoformat())
            report = fp.render({"status": "completed", "code": {}, "invocation": {"params": {}},
                                "run": {}, "input_headers": fresh})
            self.assertIn("INPUT ALIGNMENTS", report)
            self.assertIn("distinct @SQ fingerprints  2", report)
            self.assertIn("distinct @PG chains     2", report)
            self.assertIn("build=GRCh37", report)
            self.assertNotIn("SAMPLE_ALPHA", report)
            stale = fp.input_header_summary(cohort, (now + datetime.timedelta(hours=1)).isoformat())
            self.assertIsNone(stale)

    def test_an_unaligned_header_fails_the_sample_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "h.sam"
            path.write_text("@HD\tVN:1.6\n")
            result = subprocess.run([sys.executable, str(SCRIPT), "sample", "--header", str(path),
                                     "--sample", "S", "--output", str(Path(tmp) / "o.tsv")],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("no @SQ lines", result.stderr)


class Wiring(unittest.TestCase):
    EXTRACT = (ROOT / "Modules/extract_reads.nf").read_text()
    MAIN = (ROOT / "main.nf").read_text()

    def test_extract_reads_writes_and_publishes_it(self):
        self.assertIn('header_provenance.py" sample', self.EXTRACT)
        self.assertIn("--header input.header.sam", self.EXTRACT)
        self.assertIn('path("${sampleID}.header_provenance.tsv"), emit: provenance', self.EXTRACT)
        self.assertIn('{ "${params.sample_dir}/${sampleID}/provenance" }', self.EXTRACT)
        self.assertIn('saveAs: { "input_header.tsv" }', self.EXTRACT)
        # after the @SQ check, before any record is decoded
        self.assertLess(self.EXTRACT.index("grep -q '^@SQ' input.header.sam"),
                        self.EXTRACT.index('header_provenance.py" sample'))
        self.assertLess(self.EXTRACT.index('header_provenance.py" sample'),
                        self.EXTRACT.index("samtools fastq"))

    def test_the_cohort_summary_is_built_and_rendered(self):
        self.assertIn("include { extractReads; inputHeaderSummary; lowScoreRescueProfile }", self.MAIN)
        self.assertIn("inputHeaderSummary(", self.MAIN)
        self.assertIn('publishDir "${params.cohort_dir}/provenance"', self.EXTRACT)
        self.assertIn('header_provenance.py" summarize', self.EXTRACT)
        finalize = (ROOT / "scripts/finalize_provenance.py").read_text()
        self.assertIn('record["input_headers"] = input_header_summary(', finalize)

    def test_documented(self):
        docs = (ROOT / "docs/output.md").read_text()
        self.assertIn("provenance/input_header.tsv", docs)
        self.assertIn("input_header_summary.tsv", docs)



class UnconfirmedIsTwoDifferentThings(unittest.TestCase):
    """Ludmil's v3 review, Tier A item 1.

    low_score_decision used to send every non-`yes` bwa_mem to off_auto_not_bwa_mem. That
    lumped together an input we know cannot be rescued with one we simply cannot read, and
    the second silently lowers sensitivity. C1 showed the second never occurred in the TCGA
    cohorts -- every file carries a CL -- but the code path was live.
    """

    PG_SAMPE = ('@PG\tID:bwa\tPN:bwa\tVN:0.7.15\t'
                'CL:bwa sampe -r @RG\\tID:x ref.fa a.sai b.sai a.fq b.fq\n')
    PG_NO_CL = '@PG\tID:bwa\tPN:bwa\tVN:0.7.15\n'
    PG_MEM = '@PG\tID:bwa\tPN:bwa\tVN:0.7.15\tCL:bwa mem -T 0 ref.fa a.fq b.fq\n'
    PG_STAR = '@PG\tID:STAR\tPN:STAR\tVN:2.7.5c\tCL:STAR --runMode alignReads\n'

    def test_sampe_is_not_mem_and_is_not_the_same_as_an_unreadable_header(self):
        self.assertEqual(hp.bwa_mem_evidence(self.PG_SAMPE)[0], "not_mem")
        self.assertEqual(hp.bwa_mem_evidence(self.PG_NO_CL)[0], "no_cl")

    def test_mem_and_a_non_bwa_aligner_are_unchanged(self):
        self.assertEqual(hp.bwa_mem_evidence(self.PG_MEM)[0], "yes")
        self.assertEqual(hp.bwa_mem_evidence(self.PG_STAR)[0], "no")

    def test_sampe_resolves_to_its_own_off_state_and_the_run_continues(self):
        # 17 of 1,304 TCGA WES samples are this. They are correctly declined -- sampe
        # emits no AS tag (tests/check_bwa_sampe_tags.sh) -- so the run must not fail.
        self.assertEqual(hp.low_score_decision("auto", "not_mem"), "off_auto_no_as_tag")

    def test_a_missing_command_line_refuses_rather_than_guessing(self):
        with self.assertRaises(ValueError) as caught:
            hp.low_score_decision("auto", "no_cl")
        msg = str(caught.exception)
        self.assertIn("--extract_low_score_mapped", msg,
                      "the error has to say how to resolve it")
        self.assertIn("sensitivity", msg,
                      "and why guessing off is not safe")

    def test_an_explicit_request_still_wins_over_an_unreadable_header(self):
        # The escape hatch the error points at must actually work.
        self.assertEqual(hp.low_score_decision("on", "no_cl"), "on_requested")
        self.assertEqual(hp.low_score_decision("off", "no_cl"), "off_requested")

    def test_a_genuinely_different_aligner_is_still_off_auto_not_bwa_mem(self):
        # TCGA RNA: 701/701 STAR. Unchanged behaviour, and distinct from the sampe state.
        self.assertEqual(hp.low_score_decision("auto", "no"), "off_auto_not_bwa_mem")

    def test_every_state_the_decision_can_return_is_declared(self):
        for bwa_mem in ("yes", "not_mem", "no"):
            self.assertIn(hp.low_score_decision("auto", bwa_mem), hp.LOW_SCORE_STATES)
        for request in ("on", "off"):
            self.assertIn(hp.low_score_decision(request, "no"), hp.LOW_SCORE_STATES)

    def test_each_state_has_a_reason_a_person_can_read(self):
        for bwa_mem, state in (("not_mem", "off_auto_no_as_tag"),
                               ("no", "off_auto_not_bwa_mem"),
                               ("yes", "on_auto_bwa_mem")):
            with self.subTest(bwa_mem=bwa_mem):
                reason = hp.low_score_reason(state, bwa_mem, "NA", 30)
                self.assertTrue(reason.strip())
                self.assertIn(state, reason)
        self.assertIn("no AS tag", hp.low_score_reason("off_auto_no_as_tag", "not_mem", "NA", 30))


if __name__ == "__main__":
    unittest.main()
