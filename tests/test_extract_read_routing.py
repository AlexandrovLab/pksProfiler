"""Extraction must route every primary unmapped record into the one output stream.

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

# The invocation, with its backslash continuations joined: from "samtools fastq" up to
# and including the line that ends the pipeline stage. Comment lines mention the command
# too, so the first non-comment occurrence is the real one.
def fastq_invocation(text):
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines)
                 if "samtools fastq" in l and not l.strip().startswith("#"))
    out = []
    for line in lines[start:]:
        out.append(line)
        if not line.rstrip().endswith("\\"):
            break
    return " ".join(l.strip().rstrip("\\").strip() for l in out)


INVOCATION = fastq_invocation(EXTRACT)

# -1/-2 are the paired equivalents of -o; -s is singletons. Each one, named alone,
# silently drops the categories it does not cover.
ROUTING_FLAGS = ("-o", "-0", "-1", "-2", "-s")


class ReadRoutingTests(unittest.TestCase):
    def test_no_output_routing_flag_is_named(self):
        named = [f for f in ROUTING_FLAGS if re.search(rf"(?<!\S){re.escape(f)}(?!\S)", INVOCATION)]
        self.assertEqual(named, [], f"output-routing flags drop read categories: {named} in {INVOCATION!r}")

    def test_nothing_is_routed_to_dev_null(self):
        self.assertNotIn("/dev/null", INVOCATION)

    def test_the_stream_is_captured_not_left_on_stdout(self):
        # Everything samtools writes must reach the compressed output; a bare invocation
        # with no pipe would put read sequences in the task log instead.
        self.assertTrue(INVOCATION.rstrip().endswith("|"), f"not piped: {INVOCATION!r}")
        tail = EXTRACT[EXTRACT.index(INVOCATION.split("samtools fastq")[0] + "samtools fastq"):]
        self.assertRegex(tail, r"bgzip[^\n]*-c\s*>\s*\"\\?\$READS\"")

    def test_only_primary_unmapped_records_are_taken(self):
        for flag in ("-f 4", "-F 2304", "-N"):
            self.assertIn(flag, INVOCATION, f"missing {flag}: {INVOCATION!r}")

    def test_cram_reference_arguments_survive_the_routing_fix(self):
        self.assertIn("REFERENCE_ARGS", INVOCATION)


if __name__ == "__main__":
    unittest.main()
