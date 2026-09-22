"""A process must write its outputs where its output block says they are.

`cohortReport` did `mkdir -p cohort` and wrote `cohort/pks_cohort_report.html`, while
declaring `path "pks_cohort_report.html"`. The script exited 0 and the file existed, one
directory below where Nextflow looked, so the task failed with "Missing output file(s)"
after everything upstream of it had already run -- twenty minutes into a two-sample run.

No unit test runs Nextflow, so this compares the literal filenames in each output block
against the paths the script writes, and flags a declaration whose name appears in the
script only with a directory in front of it.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NF_FILES = sorted(ROOT.glob("Modules/*.nf")) + [ROOT / "main.nf"]


def processes():
    for nf in NF_FILES:
        text = nf.read_text()
        for m in re.finditer(r"^process\s+(\w+)\s*\{", text, re.M):
            body = text[m.end():]
            nxt = body.find("\nprocess ")
            body = body[:nxt] if nxt > 0 else body
            om = re.search(r"^\s*output:\s*$(.*?)^\s*(script|shell|exec):",
                           body, re.S | re.M)
            sm = re.search(r"^\s*(?:script|shell):.*?\"\"\"(.*?)\"\"\"", body, re.S | re.M)
            if om and sm:
                yield nf.name, m.group(1), om.group(1), sm.group(1)


class DeclaredOutputsMatchTheScript(unittest.TestCase):
    def test_no_output_is_written_one_directory_below_its_declaration(self):
        offenders = []
        for nf_name, process, outputs, script in processes():
            for decl in re.findall(r'path\s*\(?\s*"([^"${}]+\.[a-z]{2,5})"', outputs):
                if "/" in decl:
                    continue                       # declares a path; nothing to check
                # written with a directory in front of it, and never at the task root
                nested = re.search(rf"[\w.]+/{re.escape(decl)}\b", script)
                bare = re.search(rf"(?<![\w/]){re.escape(decl)}\b", script)
                if nested and not bare:
                    offenders.append(
                        f"{nf_name}:{process} declares '{decl}' but the script writes "
                        f"'{nested.group(0)}'")
        self.assertEqual(offenders, [], "\n" + "\n".join(offenders))

    def test_cohort_report_writes_to_the_task_root(self):
        text = (ROOT / "Modules/plotting.nf").read_text()
        block = text[text.index("process cohortReport"):]
        block = block[:block.index("\n}")]
        self.assertNotIn("mkdir -p cohort", block)
        self.assertIn("--output pks_cohort_report.html", block)
        self.assertNotIn("cohort/pks_cohort_report.html", block)


if __name__ == "__main__":
    unittest.main()
