"""L5: the AS boundary of E1's low-score rescue, named on its own.

Ludmil, 1 Oct: an explicitly named regression test at AS 29 vs AS 30. E1 extracts primary
mapped records with AS:i strictly BELOW --extract_min_as (default 30, bwa mem's default -T):
AS 29 is extracted, AS 30 is not. Secondary/supplementary and missing-AS records are in
tests/test_low_score_extraction.py; this file is only the boundary, on both mates and an
unpaired record, through the module's own on-branch pipeline (its samtools -e filter, its
awk, samtools fastq -N), and through the L4 rescued-name list the same awk writes.
"""
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_low_score_extraction as e1  # noqa: E402
from test_fragment_support_variants import _resolve_bin  # noqa: E402

SAMTOOLS = e1.SAMTOOLS or _resolve_bin("samtools", "conda_envs/samtools_env.yml")

# (name, flag, mapped?, AS, sequence, expected class) -- e1.sam_text()'s record layout.
BOUNDARY = [
    ("as29_pair", 65 + 32, True, 29, e1.seq(21), "low"),        # /1, AS 29: extracted
    ("as29_pair", 129 + 16, True, 30, e1.seq(22), "excluded"),  # /2, AS 30: not
    ("as30_pair", 65, True, 30, e1.seq(23), "excluded"),        # /1, AS 30: not
    ("as30_pair", 129, True, 29, e1.seq(24), "low"),            # /2, AS 29: extracted
    ("as29_single", 0, True, 29, e1.seq(25), "low"),            # unpaired AS 29: extracted
    ("as30_single", 16, True, 30, e1.seq(26), "excluded"),      # unpaired AS 30: not
]


def run_on_branch(records, min_as=30):
    tmp = Path(tempfile.mkdtemp())
    try:
        (tmp / "in.sam").write_text(e1.sam_text(records))
        (tmp / "run.sh").write_text(e1.on_branch_script(min_as))
        env = dict(os.environ, ALN=str(tmp / "in.sam"),
                   PATH=f"{Path(SAMTOOLS).parent}:{os.environ.get('PATH', '')}")
        result = subprocess.run(["bash", "run.sh"], cwd=tmp, env=env, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        unmapped, low, no_as, total = (tmp / "result.tsv").read_text().split()
        with gzip.open(tmp / "out.fastq.gz", "rt") as fh:
            lines = fh.read().splitlines()
        names_file = tmp / "S.low_score_rescued.names"
        names = names_file.read_text().split() if names_file.exists() else []
        return int(low), sorted(lines[i][1:] for i in range(0, len(lines), 4)), sorted(names)
    finally:
        shutil.rmtree(tmp)


@unittest.skipUnless(SAMTOOLS, "samtools >= 1.16 not found (set SAMTOOLS_BIN or CONDA_CACHE_DIR)")
class AsBoundary29Vs30(unittest.TestCase):
    def test_as29_is_extracted_and_as30_is_not(self):
        low, reads, names = run_on_branch(BOUNDARY)
        self.assertEqual(reads, ["as29_pair/1", "as29_single", "as30_pair/2"])
        self.assertEqual(low, 3)

    def test_the_rescued_name_list_is_exactly_the_as29_reads(self):
        _low, reads, names = run_on_branch(BOUNDARY)
        self.assertEqual(names, reads)

    def test_the_boundary_moves_with_the_parameter(self):
        # --extract_min_as 29: AS 29 is no longer below it.
        low, reads, names = run_on_branch(BOUNDARY, min_as=29)
        self.assertEqual((low, reads, names), (0, [], []))
        low, reads, _ = run_on_branch(BOUNDARY, min_as=31)
        self.assertEqual(low, 6)

    def test_nothing_rescued_writes_no_name_list(self):
        # L4: the optional rescue_names output is absent, so nothing downstream changes.
        low, _reads, names = run_on_branch([r for r in BOUNDARY if r[3] == 30])
        self.assertEqual((low, names), (0, []))


if __name__ == "__main__":
    unittest.main()
