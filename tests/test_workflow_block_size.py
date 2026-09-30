"""Nextflow 24.10 (the manifest's floor) keeps the main workflow block's source text as one string,
and Groovy caps a string constant at 65,535 characters. The combined rc1 branch hit that limit
(70,798) and failed to compile ("String too long") while every branch alone was under it.
This fails early, with headroom, so the next addition moves code into a function or module
instead of breaking compilation."""
import re
import unittest
from pathlib import Path

MAIN = (Path(__file__).resolve().parents[1] / "main.nf").read_text()
LIMIT = 65535
HEADROOM = 3000


def workflow_block_chars(text):
    i = text.index("\nworkflow {")
    j = text.index("{", i)
    depth = 0
    for k in range(j, len(text)):
        if text[k] == "{":
            depth += 1
        elif text[k] == "}":
            depth -= 1
            if depth == 0:
                return k - j + 1
    raise AssertionError("unbalanced workflow block")


class WorkflowBlockFitsInOneGroovyString(unittest.TestCase):
    def test_main_workflow_block_leaves_headroom_under_65535(self):
        n = workflow_block_chars(MAIN)
        self.assertLess(n, LIMIT - HEADROOM,
                        f"main workflow block is {n} chars; Nextflow 24.10 fails at {LIMIT}. "
                        "Move self-contained logic (validation, helpers) into a top-level def or a module.")

    def test_moved_validation_is_still_called_once_from_the_workflow(self):
        # Same position as before the move. Process calls in the workflow block only build the
        # dataflow graph; no task runs until the block has been evaluated, so these checks
        # still fail the run before anything executes.
        body = MAIN.split("\nworkflow {", 1)[1]
        for fn in ("validatePrefilterParams(", "validateStageParams("):
            self.assertEqual(body.count(fn), 1, fn)
            self.assertIn(f"def {fn}Map flags)", MAIN)


if __name__ == "__main__":
    unittest.main()
