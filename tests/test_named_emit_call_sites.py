"""A process with named outputs must be addressed by name where it is called.

The bug this exists for. F10 added a second output to Bracken:

    tuple val(sampleID), path(...)                        # the reports
    tuple val(sampleID), path("...taxonomy.qc.tsv"), emit: qc

which turned `Bracken(...)` from a channel into a multi-channel output object. The call
site still did

    Bracken(PKS_ISLAND_FASTQ).set { BRACKEN_PER_SAMPLE }

and the next `.map` on that object failed with "Multi-channel output cannot be applied to
operator map for which argument is already provided". That is every --pks_taxa run, and
`nextflow run --help` does not catch it: the script compiles, and the failure happens when
the workflow is constructed.

Nothing in the suite runs Nextflow, so the check has to be static. Two rules:

  1. a process whose outputs are named must not be consumed as a whole object -- no
     `Name(...).set`, `.map`, `.join` and so on directly off the call
  2. the set of processes that name only some of their outputs is pinned. Mixing is legal
     Nextflow and the existing ones work, but a half-named output block is what makes
     rule 1 easy to get wrong, so a new one should be a deliberate decision rather than
     something that appears unnoticed.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NF_FILES = sorted(ROOT.glob("Modules/*.nf")) + [ROOT / "main.nf"]
OPERATORS = ("set", "map", "filter", "join", "mix", "multiMap", "view", "flatMap",
             "collect", "groupTuple", "branch", "combine", "transpose")


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return "\n".join(re.sub(r"//.*$", "", line) for line in text.splitlines())


def output_blocks():
    """{process name: (named count, total declarations)}"""
    found = {}
    for nf in NF_FILES:
        text = strip_comments(nf.read_text())
        for m in re.finditer(r"^process\s+(\w+)\s*\{", text, re.M):
            body = text[m.end():]
            nxt = body.find("\nprocess ")
            body = body[:nxt] if nxt > 0 else body
            om = re.search(r"^\s*output:\s*$(.*?)^\s*(script|shell|exec):",
                           body, re.S | re.M)
            if not om:
                continue
            lines = om.group(1).splitlines()
            # A tuple output spans several lines and its `emit:` sits on the last one,
            # so declarations are grouped by indent and searched as blocks. Counting
            # per line made an 11-line tuple look like eleven unnamed outputs.
            starts = [i for i, l in enumerate(lines)
                      if re.match(r"\s*(tuple|path|val|env|stdout)\b", l)]
            if not starts:
                continue
            base = min(len(lines[i]) - len(lines[i].lstrip()) for i in starts)
            heads = [i for i in starts
                     if len(lines[i]) - len(lines[i].lstrip()) == base]
            blocks = [lines[a:b] for a, b in zip(heads, heads[1:] + [len(lines)])]
            named = sum(any("emit:" in l for l in blk) for blk in blocks)
            found[m.group(1)] = (named, len(blocks))
    return found


def call_then_operator(text, process):
    """`process ( ... ) . operator` with balanced parentheses, comments removed."""
    hits = []
    for m in re.finditer(rf"(?<![\w.]){re.escape(process)}\s*\(", text):
        i, depth = m.end() - 1, 0
        while i < len(text):
            if text[i] == "(":
                depth += 1
            elif text[i] == ")":
                depth -= 1
                if depth == 0:
                    break
            i += 1
        tail = text[i + 1:i + 40]
        op = re.match(r"\s*\.\s*(\w+)", tail)
        if op and op.group(1) in OPERATORS:
            hits.append(op.group(1))
    return hits


class NamedOutputsAreAddressedByName(unittest.TestCase):
    def setUp(self):
        self.procs = output_blocks()
        self.assertIn("Bracken", self.procs, "process scan found nothing")

    def test_no_process_with_named_outputs_is_consumed_whole(self):
        offenders = []
        for nf in NF_FILES:
            text = strip_comments(nf.read_text())
            for process, (named, _total) in self.procs.items():
                if not named:
                    continue
                for op in call_then_operator(text, process):
                    offenders.append(f"{nf.name}: {process}(...).{op} "
                                     f"-- use {process}.out.<name>")
        self.assertEqual(offenders, [], "\n" + "\n".join(offenders))

    # Known and accepted: each is addressed by name at every call site, which rule 1
    # enforces. Add to this only with a reason.
    KNOWN_MIXED = {"alignToContigs", "pksProfiler_align", "pksProfiler_hmm"}

    def test_no_new_process_mixes_named_and_unnamed_outputs(self):
        mixed = {p for p, (named, total) in self.procs.items()
                 if named and total > named}
        new = sorted(mixed - self.KNOWN_MIXED)
        self.assertEqual(new, [], f"\nnew half-named output block(s): {new}. "
                                  "Name every output, or add it to KNOWN_MIXED with a "
                                  "reason and check every call site uses .out.<name>.")

    def test_bracken_no_longer_mixes(self):
        named, total = self.procs["Bracken"]
        self.assertEqual((named, total), (2, 2),
                         "Bracken should declare exactly two fully named outputs")

    def test_bracken_reports_channel_is_named(self):
        taxa = (ROOT / "Modules/pks_taxa.nf").read_text()
        self.assertIn("emit: reports", taxa)
        self.assertIn("emit: qc", taxa)
        main = strip_comments((ROOT / "main.nf").read_text())
        self.assertIn("Bracken.out.reports", main)


if __name__ == "__main__":
    unittest.main()


class ChannelsUsedUnconditionallyAreDeclaredUnconditionally(unittest.TestCase):
    """A variable the workflow always reads must not be assigned only inside a branch.

    `cohort_report_gate` was assigned inside `if (do_align)` and read unconditionally by
    `cohortReport`. Under `--profiling_method hmm` it was never assigned, and the run died
    while the workflow was still being built: "No such variable: cohort_report_gate".
    Nothing caught it -- the script compiles, and no combination in the suite ran the
    HMM-only path. tests/check_workflow_construction.sh is the real check; this pins the
    specific variable so a later edit cannot quietly push it back inside the branch.
    """

    def test_the_cohort_report_gate_is_declared_outside_the_alignment_branch(self):
        body = strip_comments((ROOT / "main.nf").read_text()).split("workflow {", 1)[1]
        declare = body.index("cohort_report_gate = classifyPksReadEvidence_gate()")
        branch = body.index("if (do_align)")
        self.assertLess(declare, branch,
                        "cohort_report_gate must be declared before `if (do_align)`; "
                        "assigned only inside it, --profiling_method hmm crashes at "
                        "workflow construction.")
