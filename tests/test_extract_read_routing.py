"""Extraction must route every extracted primary record into the one output stream.

F01, from Ludmil's 2026-09-19 report. `samtools fastq` splits its output by mate bit:
-o takes READ1/READ2, -0 takes READ_OTHER, -s takes singletons. Naming any of them
sends the categories you did not name somewhere else -- `-0 /dev/null` discarded every
READ_OTHER record, and the earlier `-o FILE` form leaked them to the task's stdout, so
read sequences reached the task logs. Naming none of them puts all four categories on
stdout, which the pipeline captures.

These are static checks on the invocation. The empirical proof that each category
behaves as described is handoff/f01_check/01_synthetic_routing_test.sh, which builds
its own SAM records and runs real samtools.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXTRACT = (ROOT / "Modules/extract_reads.nf").read_text()

# The invocations, with their backslash continuations joined: from "samtools fastq" up
# to and including the line that ends the pipeline stage. Comment lines mention the
# command too, so only non-comment occurrences count. E1 (Ludmil, 25 Sep) added a second
# one: with low-score mapped extraction on, samtools view selects the records and
# samtools fastq reads them from stdin; off, samtools fastq decodes the input itself.
def fastq_invocations(text):
    lines = text.splitlines()
    found = []
    for start, line in enumerate(lines):
        if "samtools fastq" not in line or line.strip().startswith("#"):
            continue
        out = []
        for l in lines[start:]:
            out.append(l)
            if not l.rstrip().endswith("\\"):
                break
        found.append(" ".join(l.strip().rstrip("\\").strip() for l in out))
    return found


INVOCATIONS = fastq_invocations(EXTRACT)
# The off branch: samtools fastq decodes the alignment itself, FLAG 4 only (pre-E1).
INVOCATION = next(i for i in INVOCATIONS if "-f 4" in i)
# The on branch: samtools fastq reads what samtools view selected (FLAG 4 + AS < min).
ON_INVOCATION = next(i for i in INVOCATIONS if i.rstrip(" |").endswith(" -"))
ON_BRANCH = EXTRACT[EXTRACT.index('if [[ "\\$LOW_SCORE" == "on" ]]; then'):EXTRACT.index("\n    else\n")]

# -1/-2 are the paired equivalents of -o; -s is singletons. Each one, named alone,
# silently drops the categories it does not cover.
ROUTING_FLAGS = ("-o", "-0", "-1", "-2", "-s")


class ReadRoutingTests(unittest.TestCase):
    def test_there_are_exactly_the_two_invocations(self):
        self.assertEqual(len(INVOCATIONS), 2, INVOCATIONS)

    def test_no_output_routing_flag_is_named(self):
        for invocation in INVOCATIONS:
            named = [f for f in ROUTING_FLAGS if re.search(rf"(?<!\S){re.escape(f)}(?!\S)", invocation)]
            self.assertEqual(named, [], f"output-routing flags drop read categories: {named} in {invocation!r}")

    def test_nothing_is_routed_to_dev_null(self):
        for invocation in INVOCATIONS:
            self.assertNotIn("/dev/null", invocation)

    def test_the_stream_is_captured_not_left_on_stdout(self):
        # Everything samtools writes must reach the compressed output; a bare invocation
        # with no pipe would put read sequences in the task log instead.
        for invocation in INVOCATIONS:
            self.assertTrue(invocation.rstrip().endswith("|"), f"not piped: {invocation!r}")
        for marker in ("samtools fastq \\\n        -@", "samtools fastq \\\n            -@"):
            tail = EXTRACT[EXTRACT.index(marker):]
            self.assertRegex(tail, r"tee extract\.count\.fifo \|\s*bgzip[^\n]*-c\s*>\s*\"\\?\$READS\"")

    def test_only_primary_unmapped_records_are_taken(self):
        for flag in ("-f 4", "-F 2304", "-N"):
            self.assertIn(flag, INVOCATION, f"missing {flag}: {INVOCATION!r}")

    def test_the_on_branch_takes_primary_flag4_or_low_score_records_only(self):
        for flag in ("-F 2304", "-N"):
            self.assertIn(flag, ON_INVOCATION)
        self.assertIn("samtools view -h", ON_BRANCH)
        self.assertIn("-F 2304", ON_BRANCH[ON_BRANCH.index("samtools view -h"):ON_BRANCH.index("awk -F")])
        self.assertIn("-e '(flag & 4) || !exists([AS]) || [AS] < ${params.extract_min_as}'", ON_BRANCH)

    def test_cram_reference_arguments_survive_the_routing_fix(self):
        self.assertIn("REFERENCE_ARGS", INVOCATION)
        # on: the decode is samtools view's, so it carries the reference
        self.assertIn("REFERENCE_ARGS", ON_BRANCH[ON_BRANCH.index("samtools view -h"):ON_BRANCH.index("awk -F")])


if __name__ == "__main__":
    unittest.main()
