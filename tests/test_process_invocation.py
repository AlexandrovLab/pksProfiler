"""No process or workflow may be invoked twice in main.nf.

Nextflow rejects a second invocation of the same component in one workflow scope with
"Process 'X' has been already used". It is a launch-time failure, invisible to unit tests
that only read the script, and it only surfaces on the specific flag combination that
reaches both call sites -- which is how it escaped into a real run.
"""
import collections
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()


def strip_comments(text):
    """Comments are not invocations.

    A comment explaining a call -- "`Bracken(...).set { }` captured the multi-channel
    object" -- counted as a second invocation and failed this test. The rule is about
    executable code, so the text it reads has to be executable code.
    """
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return "\n".join(re.sub(r"//.*$", "", line) for line in text.splitlines())


BODY = strip_comments(MAIN.split("workflow {", 1)[1])


def included_names():
    names = set()
    for block in re.finditer(r"include\s*\{([^}]*)\}", MAIN):
        for part in block.group(1).split(";"):
            part = part.strip()
            if part:
                names.add(part.split(" as ")[-1].strip())
    return names


def invocation_counts():
    return {n: len(re.findall(rf"(?<![\w.]){re.escape(n)}\s*\(", BODY)) for n in included_names()}


class InvocationTests(unittest.TestCase):
    def test_no_component_is_invoked_more_than_once(self):
        repeated = {n: c for n, c in invocation_counts().items() if c > 1}
        self.assertEqual(repeated, {}, f"invoked more than once: {repeated}")

    def test_the_two_bracken_plots_use_distinct_aliases(self):
        self.assertIn("plotBrackenTaxa as plotPKSTaxa", MAIN)
        self.assertIn("plotBrackenTaxa as plotCommunityTaxa", MAIN)


class PublishDirTests(unittest.TestCase):
    def test_every_publishdir_param_is_defined(self):
        defined = set(re.findall(r"^params\.(\w+)\s*=", MAIN, re.M))
        missing = collections.defaultdict(list)
        for module in (ROOT / "Modules").glob("*.nf"):
            for line in module.read_text().splitlines():
                if "publishDir" not in line:
                    continue
                for name in re.findall(r"params\.(\w+)", line):
                    if name not in defined:
                        missing[module.name].append(name)
        self.assertEqual(dict(missing), {}, f"publishDir references undefined params: {dict(missing)}")


if __name__ == "__main__":
    unittest.main()
