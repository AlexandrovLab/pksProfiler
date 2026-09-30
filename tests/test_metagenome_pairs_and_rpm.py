"""Metagenome paired assembly and pks reads per million (Ludmil, 26 Sep mode memo section 5).

"Restore pairs before metagenome assembly ... reconstruct R1/R2 plus singletons, then feed
paired inputs to MEGAHIT" and "Make normalized abundance the primary metagenome output.
Lead with pks reads per million host-depleted reads; retain raw pks reads."

Plus A7 (25 Sep): mate-tag self-diagnostics -- "current valid paired runs report zero; a
deliberately broken QNAME fixture triggers the diagnostic."

And the invariant the pairing must not break: profiling narrows, assembly never does.
Every host-depleted read reaches the assembler, as a mate or as a singleton.

All reads here are synthetic.
"""
import csv
import gzip
import importlib.util
import random
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/recover_recruited_mates.py"
MAIN = (ROOT / "main.nf").read_text()
MAG = (ROOT / "Modules/pks_mag.nf").read_text()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


classify = load("classify_tumor_pks_evidence", ROOT / "scripts/classify_tumor_pks_evidence.py")
recover = load("recover_recruited_mates", SCRIPT)


def record(name, seq):
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n"


def fragments(n, seed=7):
    rng = random.Random(seed)
    return [(f"frag{i:04d}", "".join(rng.choices("ACGT", k=60)), "".join(rng.choices("ACGT", k=60)))
            for i in range(n)]


def write_fastq(path, records):
    with gzip.open(path, "wt") as fh:
        fh.write("".join(records))


def read_fastq(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    return ["\n".join(lines[i:i + 4]) + "\n" for i in range(0, len(lines), 4)]


def run_pairing(d, records, extra=()):
    fastq = d / "in.fastq.gz"
    write_fastq(fastq, records)
    out = {k: d / f"{k}.fastq.gz" for k in ("r1", "r2", "single")}
    result = subprocess.run([sys.executable, str(SCRIPT), "--all-reads", "--fastq", str(fastq),
                             "--r1", str(out["r1"]), "--r2", str(out["r2"]), "--single", str(out["single"]),
                             "--stats", str(d / "stats.tsv"), "--tmpdir", str(d), *extra],
                            capture_output=True, text=True)
    if result.returncode != 0:
        return result, None, None
    with (d / "stats.tsv").open() as fh:
        stats = {row["metric"]: row["value"] for row in csv.DictReader(fh, delimiter="\t")}
    return result, stats, {k: read_fastq(v) for k, v in out.items()}


class Pairing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.d = Path(self.tmp.name)

    def test_block_ordered_and_interleaved_inputs_give_the_same_pairs(self):
        frags = fragments(50)
        blocks = [record(f"{n}/1", a) for n, a, _ in frags] + [record(f"{n}/2", b) for n, _, b in frags]
        interleaved = [r for n, a, b in frags for r in (record(f"{n}/1", a), record(f"{n}/2", b))]
        (self.d / "a").mkdir(); (self.d / "b").mkdir()
        _, stats_a, out_a = run_pairing(self.d / "a", blocks)
        _, stats_b, out_b = run_pairing(self.d / "b", interleaved)
        self.assertEqual(out_a, out_b)
        self.assertEqual(stats_a["paired_fragments"], "50")
        self.assertEqual(out_a["single"], [])
        # R1[i] and R2[i] are mates
        for one, two in zip(out_a["r1"], out_a["r2"]):
            self.assertEqual(one.split("\n")[0][:-2], two.split("\n")[0][:-2])
        # A7 acceptance: a valid paired run reports zero anomalies
        for metric in ("fragments_without_mate_tag", "orphan_reads", "duplicate_mate_slot_overwrites"):
            self.assertEqual(stats_a[metric], "0", metric)
        self.assertEqual((stats_a["mate_tag_1"], stats_a["mate_tag_2"]), ("50", "50"))

    def test_assembly_never_narrows_every_read_is_written_once(self):
        frags = fragments(30)
        reads = [record(f"{n}/1", a) for n, a, _ in frags[:20]]            # 20 with mates below
        reads += [record(f"{n}/2", b) for n, _, b in frags[:15]]           # 5 orphans /1
        reads += [record(f"{n}/2", b) for n, _, b in frags[20:25]]         # 5 orphans /2
        reads += [record(f"single{i}", a) for i, (_, a, _) in enumerate(frags[25:])]  # 5 unsuffixed
        reads += [record(f"{frags[0][0]}/1", "ACGT" * 15)]                 # one duplicate slot
        result, stats, out = run_pairing(self.d, reads, ("--max-duplicate-fraction", "0.5"))
        self.assertEqual(result.returncode, 0, result.stderr)
        written = out["r1"] + out["r2"] + out["single"]
        self.assertEqual(Counter(written), Counter(reads), "pairs + singletons must be the input exactly")
        self.assertEqual(stats["input_reads"], str(len(reads)))
        self.assertEqual(stats["paired_fragments"], "15")
        self.assertEqual(stats["orphan_reads"], "10")
        self.assertEqual(stats["fragments_without_mate_tag"], "5")
        self.assertEqual(stats["duplicate_mate_slot_overwrites"], "1")
        self.assertIn("duplicate read names", result.stderr)

    def test_a_broken_qname_fixture_trips_the_diagnostic_and_the_tolerance(self):
        # every read name collapsed to the same token: the malformed case A7 names
        reads = [record("read", a) for _, a, _ in fragments(20)]
        result, _stats, _out = run_pairing(self.d, reads)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("have duplicate names", result.stderr)

    def test_single_end_input_is_all_singletons(self):
        reads = [record(n, a) for n, a, _ in fragments(10)]
        _, stats, out = run_pairing(self.d, reads)
        self.assertEqual((stats["paired_fragments"], stats["singleton_reads"]), ("0", "10"))
        self.assertEqual(Counter(out["single"]), Counter(reads))

    def test_casava_second_field_is_a_mate_tag_but_a_number_is_not(self):
        self.assertEqual(recover.base_and_mate("@r1 1:N:0:ACGT"), ("r1", "1"))
        self.assertEqual(recover.base_and_mate("@r1 2:N:0:ACGT"), ("r1", "2"))
        self.assertEqual(recover.base_and_mate("@r1/2 extra"), ("r1", "2"))
        self.assertEqual(recover.base_and_mate("@r1 150"), ("r1", "0"))

    def test_recruitment_keeps_its_behaviour_and_gains_the_census(self):
        frags = fragments(4)
        reads = [r for n, a, b in frags for r in (record(f"{n}/1", a), record(f"{n}/2", b))]
        fastq = self.d / "in.fastq.gz"
        write_fastq(fastq, reads)
        sam = self.d / "rec.sam"
        # Ludmil, 25 Sep A1: recruitment now filters on MAPQ and needs the island's span
        # in the recruitment reference; these alignments are unique (MAPQ 42) island hits.
        sam.write_text("".join(f"{n}/1\t0\tref\t1\t42\t60M\t*\t0\t0\t{a}\t{'I' * 60}\tNM:i:0\n"
                               for n, a, _ in frags[:2]))
        stats_path = self.d / "stats.tsv"
        subprocess.run([sys.executable, str(SCRIPT), "--sam", str(sam), "--fastq", str(fastq),
                        "--r1", str(self.d / "r1.fq.gz"), "--r2", str(self.d / "r2.fq.gz"),
                        "--single", str(self.d / "s.fq.gz"), "--stats", str(stats_path),
                        "--island-start", "1", "--island-end", "1000"], check=True)
        with stats_path.open() as fh:
            stats = {row["metric"]: row["value"] for row in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual((stats["recruited_fragment_ids"], stats["paired_fragments"],
                          stats["singleton_reads"]), ("2", "2", "0"))
        self.assertEqual(stats["duplicate_mate_slot_overwrites"], "0")
        self.assertEqual(len(read_fastq(self.d / "r1.fq.gz")), 2)


class Wiring(unittest.TestCase):
    def test_the_mag_lane_pairs_the_complete_host_depleted_stream(self):
        # main.nf still hands pksMAG the complete stream; the pairing is inside the lane
        self.assertIn("pksMAG(MAG_ASSEMBLY_READS)", MAIN)
        workflow = MAG[MAG.index("workflow pksMAG"):]
        self.assertIn("pairDepletedReads(reads)", workflow)
        block = MAG[MAG.index("process pairDepletedReads"):MAG.index("process megahitAssemble")]
        self.assertIn("--all-reads", block)
        self.assertIn('saveAs: { "mate_pairing.tsv" }', block)
        self.assertIn("mate_pairing = pairDepletedReads.out.stats", MAG)

    def test_profiling_is_unchanged(self):
        self.assertIn("pksProfilerAlign(PROFILING_READS)", MAIN)
        self.assertNotIn("pairDepletedReads", MAIN)


class ReadsPerMillion(unittest.TestCase):
    def test_the_rate_and_its_denominator(self):
        self.assertEqual(classify.reads_per_million(500, 2_500_000), "200.0000")
        self.assertEqual(classify.reads_per_million(0, 10), "0.0000")
        self.assertEqual(classify.reads_per_million(0, 0), "NA")
        with tempfile.TemporaryDirectory() as tmp:
            qc = Path(tmp) / "dep.tsv"
            qc.write_text("Sample\tMetric\tValue\nS\treads_after_t2t_phix\t100\nS\treads_after_pangenome\t80\n")
            self.assertEqual(classify.depleted_reads(qc), 80)
            qc.write_text("Sample\tMetric\tValue\nS\treads_after_t2t_phix\t100\nS\treads_after_pangenome\tNA\n")
            self.assertEqual(classify.depleted_reads(qc), 100)

    def test_only_metagenome_rows_carry_it_and_the_raw_count_stays(self):
        self.assertEqual(classify.SCHEMAS["metagenome"][-2:], ["host_depleted_reads", "pks_reads_per_million"])
        self.assertIn("pks_reads", classify.SCHEMAS["metagenome"])
        for mode in ("wgs", "wes", "rna"):
            self.assertNotIn("pks_reads_per_million", classify.SCHEMAS[mode])

    def test_the_denominator_is_host_depletion_not_the_prefiltered_input(self):
        block = MAIN[MAIN.index("read_evidence_input_ch = ALIGN_OUT.profile"):]
        block = block[:block.index("classifyPksReadEvidence(read_evidence_input_ch)")]
        self.assertIn(".join(MAP_OUT.qc, by: 0)", block)
        targeted = (ROOT / "Modules/pks_targeted.nf").read_text()
        self.assertIn('--depletion-qc "${depletion_qc}"', targeted)

    def test_cohort_report_and_master_summary_show_it(self):
        report = (ROOT / "scripts/build_cohort_report.py").read_text()
        self.assertIn('"pks_reads_per_million", "pks_reads"', report)
        source = (ROOT / "scripts/build_master_summary.py").read_text()
        self.assertIn('"pks_reads_per_million", "pks_reads", "host_depleted_reads"', source)
        docs = (ROOT / "docs/output.md").read_text()
        self.assertIn("pks_reads_per_million", docs)
        self.assertIn("mate_pairing.tsv", docs)

    def test_master_summary_end_to_end(self):
        sys.path.insert(0, str(ROOT / "scripts"))
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            ev = d / "inputs/1/read_evidence.tsv"
            ev.parent.mkdir(parents=True)
            fields = classify.SCHEMAS["metagenome"]
            row = classify.evidence_row("metagenome", "M1", 12, 3, (0.1, 0.05, 0.01),
                                        {"multi_gene": (5, 3, .01), "broad_island": (30, 8, .075),
                                         "extensive_island": (100, 10, .15)}, depleted=3_000_000)
            with ev.open("w", newline="") as fh:
                w = csv.DictWriter(fh, fields, delimiter="\t"); w.writeheader(); w.writerow(row)
            (d / "manifest.tsv").write_text("sample\tkind\tpath\nM1\tread_evidence\tinputs/1/read_evidence.tsv\n")
            out = d / "master.tsv"
            subprocess.run([sys.executable, str(ROOT / "scripts/build_master_summary.py"),
                            "--inputs-manifest", str(d / "manifest.tsv"), "--output", str(out),
                            "--sample-type", "metagenome"], check=True, capture_output=True)
            with out.open() as fh:
                got = list(csv.DictReader(fh, delimiter="\t"))[0]
        self.assertEqual((got["pks_reads_per_million"], got["pks_reads"], got["host_depleted_reads"]),
                         ("4.0000", "12", "3000000"))


if __name__ == "__main__":
    unittest.main()
