"""A8 (Ludmil, 25 Sep): targetedPksRecruit streams bowtie2 into recover_recruited_mates.py.

It used to write the recruitment SAM uncompressed to scratch (-S) and read it back. The
classifier reads the SAM exactly once, in order, so it can take it on stdin. Acceptance:
recruitment IDs and stats byte-identical before and after. Here: the same synthetic SAM
(hand-written records on the public island+flank reference's coordinates) given as a
file and on stdin must give byte-identical stats and identical recovered reads, and the
module must pipe rather than write a SAM. The live comparison against a pre-A8 run is in
the A8 commit message.
"""
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_recruitment_classes as rc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TARGETED = (ROOT / "Modules/pks_targeted.nf").read_text()

RECORDS = [rc.sam_record("p1/1", 20000, "150M"), rc.sam_record("p1/2", 20300, "150M"),
           rc.sam_record("p2/1", rc.ISLAND_START - 75, "150M"),
           rc.sam_record("p2/2", rc.ISLAND_START + 200, "150M"),
           rc.sam_record("p3/1", rc.ISLAND_END - 59, "60M90S"), rc.sam_record("p3/2", 30000, "60M90S"),
           rc.sam_record("p4/1", 100, "150M"), rc.sam_record("p4/2", 400, "150M"),
           rc.sam_record("p5/1", 20000, "150M", mapq=3),
           rc.sam_record("p6/1", 21000, "60M90S"),                 # rejected: low query coverage
           rc.sam_record("p7/1", 22000, "150M", flag=256)]          # secondary: ignored
NAMES = [f"p{i}/{m}" for i in range(1, 5) for m in (1, 2)] + ["p5/1", "p6/1", "p6/2", "p7/1", "p8/1"]


def recruit(tmp, sam, use_stdin):
    tmp = Path(tmp)
    out = {k: tmp / f"{k}.fastq.gz" for k in ("r1", "r2", "single")}
    args = [sys.executable, str(rc.SCRIPT), "--sam", "-" if use_stdin else str(sam),
            "--fastq", str(tmp.parent / "reads.fastq.gz"),
            "--r1", str(out["r1"]), "--r2", str(out["r2"]), "--single", str(out["single"]),
            "--stats", str(tmp / "stats.tsv"), "--island-start", str(rc.ISLAND_START),
            "--island-end", str(rc.ISLAND_END), "--reference-length", str(rc.REF_LEN)]
    with open(sam) as handle:
        result = subprocess.run(args, stdin=handle if use_stdin else subprocess.DEVNULL,
                                capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    reads = {k: gzip.decompress(p.read_bytes()) for k, p in out.items()}
    return (tmp / "stats.tsv").read_bytes(), reads


class StreamedSamIsTheSameSam(unittest.TestCase):
    def test_stdin_and_file_give_byte_identical_stats_and_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sam = tmp / "r.sam"
            rc.write_sam(sam, RECORDS)
            with gzip.open(tmp / "reads.fastq.gz", "wt") as handle:
                handle.write(rc.fastq(NAMES))
            (tmp / "file").mkdir()
            (tmp / "stdin").mkdir()
            file_stats, file_reads = recruit(tmp / "file", sam, use_stdin=False)
            pipe_stats, pipe_reads = recruit(tmp / "stdin", sam, use_stdin=True)
        self.assertEqual(pipe_stats, file_stats)
        self.assertEqual(pipe_reads, file_reads)
        # Not vacuous: fragments were recruited and classes counted.
        self.assertIn(b"recruited_fragment_ids\t4\n", file_stats)
        self.assertIn(b"reads_rejected_low_query_coverage\t2\n", file_stats)
        self.assertEqual(file_reads["r1"].count(b"@p"), 4)

    def test_the_reference_length_check_still_runs_on_a_stream(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sam = tmp / "r.sam"
            rc.write_sam(sam, RECORDS, length=rc.REF_LEN + 1)
            with gzip.open(tmp / "reads.fastq.gz", "wt") as handle:
                handle.write(rc.fastq(NAMES))
            (tmp / "stdin").mkdir()
            with self.assertRaises(AssertionError) as caught:
                recruit(tmp / "stdin", sam, use_stdin=True)
        self.assertIn("do not match the expected island+flank reference length", str(caught.exception))


class TheModulePipes(unittest.TestCase):
    def block(self):
        return TARGETED.split("process targetedPksRecruit", 1)[1].split("\nprocess ", 1)[0]

    def test_bowtie2_streams_into_the_classifier(self):
        block = self.block()
        self.assertRegex(block, r'-U "\$\{reads\}" \|\n\s*python "\$\{params\.scripts\}/recover_recruited_mates\.py" --sam - ')
        self.assertNotIn("-S ", block)
        self.assertNotIn("pks_recruitment.sam", block)
        self.assertIn("set -euo pipefail", block)


if __name__ == "__main__":
    unittest.main()
