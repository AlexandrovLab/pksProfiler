"""No unescaped `$` before a digit inside a Nextflow script block.

On 21 September the F07 parser never ran. Three of five pksProfilerHMM tasks died with

    python3: can't open file '<sample>.nhmmscan.tblout/hmm_best_hit.py'

because a comment in the script block -- written to explain the very bug F07 fixes --
contained a bare `$13`:

    # F07: a typed parser, not awk coercion. `$13+0` turned a non-numeric E-value

Groovy interpolates inside a triple-quoted string, so `$1` was consumed and every later
`${...}` in that block shifted one position: `${params.scripts}` rendered as the tblout
filename, `${tblout}` as the counts file, and so on. Reproduced in four lines and fixed by
escaping it to `\\$13`.

It was invisible everywhere it mattered. The script compiles. `nextflow run --help`
passes. `-preview` builds the graph without rendering task scripts. The unit suite reads
the module as text. And `params.sample_failure_strategy` defaults to `ignore` for the
per-sample lane, so the run logged "Error is ignored" and carried on -- a validation run
can report success with the HMM lane having produced nothing.

Awk field references (`$1`, `$13`, `$NF`) are the common case in this codebase, so the
hazard recurs every time one is written into a script block or a comment inside one.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NF_FILES = sorted(ROOT.glob("Modules/*.nf")) + [ROOT / "main.nf"]
# `script:`/`shell:` up to the closing triple quote.
BLOCK = re.compile(r'(?:script|shell):.*?"""(.*?)"""', re.S)
# An unescaped $ followed by a digit: Groovy eats it and shifts everything after.
HAZARD = re.compile(r'(?<!\\)\$\d')


class ScriptBlocksDoNotShiftTheirInterpolation(unittest.TestCase):
    def test_no_unescaped_dollar_before_a_digit(self):
        offenders = []
        for nf in NF_FILES:
            text = nf.read_text()
            for m in BLOCK.finditer(text):
                block, start = m.group(1), m.start(1)
                for hit in HAZARD.finditer(block):
                    line = text[: start + hit.start()].count("\n") + 1
                    frag = block[max(0, hit.start() - 50):hit.start() + 20]
                    offenders.append(f"{nf.name}:{line}  ...{frag.strip()}...")
        self.assertEqual(
            offenders, [],
            "\nUnescaped $<digit> in a script block shifts every later ${...}:\n"
            + "\n".join(offenders)
            + "\nWrite \\$1 for an awk field, or move the note outside the block.")

    def test_the_f07_comment_is_escaped(self):
        text = (ROOT / "Modules/pksProfiler_hmm.nf").read_text()
        self.assertIn(r"`\$13+0`", text)
        self.assertNotIn("`$13+0`", text)

    def test_the_guard_catches_the_original(self):
        sample = 'script:\n    """\n    # `$13+0` coerced it\n    echo "${x}"\n    """'
        self.assertTrue(any(HAZARD.search(m.group(1)) for m in BLOCK.finditer(sample)))


if __name__ == "__main__":
    unittest.main()
