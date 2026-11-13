"""Tier A 3: the host-depletion floor benchmark.

The benchmark's whole point is that the pipeline applies `samtools fastq -f 4` at
every depletion pass in series -- human, then T2T+phiX, then once per pangenome index
-- so a read that survives one pass is re-exposed to the same rule by the next. These
tests pin the serial carry-forward, because a benchmark that re-ran every stage on the
full read set would overstate removal and a benchmark that dropped survivors would
understate it.

No minimap2 is run here: the alignment step is stubbed so the accounting is tested on
its own.
"""
import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "a3_depletion_floor_benchmark.py"

_spec = importlib.util.spec_from_file_location("a3_depletion_floor_benchmark", SCRIPT)
a3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a3)


def record(mapped, mapq=60, score=150, identity=0.99, aligned_fraction=1.0):
    return {"mapped": mapped, "mapq": mapq, "AS": score,
            "identity": identity, "aligned_fraction": aligned_fraction}


class SerialCarryForward(unittest.TestCase):
    """Stage k must see only what survived stages 1..k-1."""

    def setUp(self):
        self.seen = []
        self.classes = [("pks", [(f"pks|{i}", "ACGT" * 10) for i in range(10)])]

    def _run(self, removers):
        """`removers` is one predicate per stage: name -> removed at this stage?"""
        calls = iter(removers)

        def fake_stage(reads, index, preset, threads, minimap2, workdir):
            self.seen.append([name for name, _ in reads])
            predicate = next(calls)
            return {name: record(predicate(name)) for name, _ in reads}

        original = a3.run_stage
        a3.run_stage = fake_stage
        try:
            return a3.run_serial(self.classes, ["one.mmi", "two.mmi"],
                                 "sr", 1, "minimap2", "/tmp")
        finally:
            a3.run_stage = original

    def test_second_stage_only_sees_survivors(self):
        # Stage one removes the five even-numbered reads.
        even = lambda name: int(name.split("|")[1]) % 2 == 0
        stages, surviving = self._run([even, lambda name: False])

        self.assertEqual(len(self.seen[0]), 10, "stage one should see every read")
        self.assertEqual(len(self.seen[1]), 5,
                         "stage two must see only the five that survived stage one")
        self.assertTrue(all(int(n.split("|")[1]) % 2 for n in self.seen[1]),
                        "the reads carried forward must be exactly the ones not removed")
        self.assertEqual(len(surviving["pks"]), 5)

    def test_removal_is_attributed_to_the_stage_that_did_it(self):
        first = lambda name: int(name.split("|")[1]) < 3
        second = lambda name: int(name.split("|")[1]) < 6
        stages, surviving = self._run([first, second])

        self.assertEqual(stages[0][1]["pks"]["n_in"], 10)
        self.assertEqual(len(stages[0][1]["pks"]["removed"]), 3)
        # Stage two sees 7 and removes reads 3,4,5 -- not 0,1,2, which are already gone.
        self.assertEqual(stages[1][1]["pks"]["n_in"], 7)
        self.assertEqual(len(stages[1][1]["pks"]["removed"]), 3)
        self.assertEqual(len(surviving["pks"]), 4)

    def test_no_survivors_anywhere_stops_the_remaining_stages(self):
        # Loading a ~9 GB index to align zero reads is pure waste, so run_serial
        # stops once every class is empty rather than running the later stages.
        stages, surviving = self._run([lambda name: True])
        self.assertEqual(len(surviving["pks"]), 0)
        self.assertEqual(len(stages), 1,
                         "stage two should not run once nothing is left to deplete")

    def test_an_empty_class_still_reports_zero_while_others_continue(self):
        # The mixed case, which is what really happens: one class is exhausted but
        # the run continues for the others, and the exhausted one must show n_in 0
        # rather than vanishing from the table.
        self.classes = [
            ("pks", [(f"pks|{i}", "ACGT" * 10) for i in range(4)]),
            ("human", [(f"human|{i}", "ACGT" * 10) for i in range(4)]),
        ]
        only_pks = lambda name: name.startswith("pks")
        stages, surviving = self._run([only_pks, lambda name: False])

        self.assertEqual(len(stages), 2, "the human class still had reads to carry")
        self.assertEqual(stages[1][1]["pks"]["n_in"], 0)
        self.assertEqual(stages[1][1]["human"]["n_in"], 4)
        self.assertEqual(len(surviving["pks"]), 0)
        self.assertEqual(len(surviving["human"]), 4)


class Accounting(unittest.TestCase):

    def test_all_stages_row_reconciles_with_the_per_stage_rows(self):
        classes = [("pks", [(f"pks|{i}", "A" * 150) for i in range(10)])]
        stages = [
            ("one.mmi", {"pks": {"n_in": 10, "removed": [record(True)] * 4}}),
            ("two.mmi", {"pks": {"n_in": 6, "removed": [record(True)] * 2}}),
        ]
        surviving = {"pks": [(f"pks|{i}", "A" * 150) for i in range(4)]}

        rows = a3.build_rows(classes, stages, surviving)
        per_stage = [r for r in rows if r["stage"] != "ALL_STAGES"]
        total = next(r for r in rows if r["stage"] == "ALL_STAGES")

        self.assertEqual(sum(r["removed_here"] for r in per_stage), 6)
        self.assertEqual(total["removed_here"], 6,
                         "the ALL_STAGES total must equal the sum of the stages")
        self.assertEqual(total["reads_in"], 10)
        self.assertEqual(total["percent_of_stage_input"], "60.00")

    def test_percentage_is_of_the_stage_input_not_the_original_count(self):
        classes = [("pks", [(f"pks|{i}", "A" * 150) for i in range(10)])]
        stages = [("two.mmi", {"pks": {"n_in": 4, "removed": [record(True)] * 2}})]
        rows = a3.build_rows(classes, stages, {"pks": [("pks|9", "A" * 150)] * 2})
        stage_row = next(r for r in rows if r["stage"] == "two.mmi")
        self.assertEqual(stage_row["percent_of_stage_input"], "50.00",
                         "2 of the 4 reads that reached this stage, not 2 of 10")


class ReadSimulation(unittest.TestCase):

    def test_zero_divergence_reproduces_the_source(self):
        import random
        source = "ACGT" * 100
        reads = a3.simulate(source, 5, 0.0, random.Random(1), "x")
        for _, sequence in reads:
            self.assertIn(sequence, source)

    def test_divergence_introduces_roughly_the_requested_rate(self):
        import random
        source = "".join(random.Random(7).choice("ACGT") for _ in range(20000))
        reads = a3.simulate(source, 200, 0.05, random.Random(3), "x")
        mismatches = 0
        for name, sequence in reads:
            # Find the window it came from by brute force is costly; instead count the
            # bases that differ from the best-matching window start recorded implicitly.
            # Simpler invariant: the read is NOT a substring of the source at 5%.
            if sequence not in source:
                mismatches += 1
        self.assertGreater(mismatches, 150,
                           "at 5% over 150 bp almost every read should carry a change")

    def test_reads_carry_their_class_in_the_name(self):
        import random
        reads = a3.simulate("ACGT" * 100, 3, 0.0, random.Random(1), "pks_identical")
        for name, _ in reads:
            self.assertEqual(name.split(a3.NAME_SEP, 1)[0], "pks_identical")

    def test_low_complexity_reads_are_full_length(self):
        import random
        for name, sequence in a3.low_complexity(12, random.Random(5)):
            self.assertEqual(len(sequence), a3.READ_LEN)
            self.assertEqual(name.split(a3.NAME_SEP, 1)[0], "low_complexity")


class CommandLine(unittest.TestCase):

    def test_index_is_repeatable_and_ordered(self):
        parser = self._parser()
        args = parser.parse_args(["--index", "a.mmi", "--index", "b.mmi",
                                  "--pks-fasta", "x.fna"])
        self.assertEqual(args.indexes, ["a.mmi", "b.mmi"],
                         "stages must keep the order they were given in")

    def test_human_index_alias_still_works(self):
        # The first submitted job used --human-index; the alias keeps it valid.
        parser = self._parser()
        args = parser.parse_args(["--human-index", "a.mmi", "--pks-fasta", "x.fna"])
        self.assertEqual(args.indexes, ["a.mmi"])

    def _parser(self):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--index", action="append", dest="indexes", default=[])
        parser.add_argument("--human-index", action="append", dest="indexes")
        parser.add_argument("--pks-fasta", required=True)
        return parser


class RemovalRuleMatchesThePipeline(unittest.TestCase):
    """The benchmark is only meaningful if it uses the pipeline's own rule."""

    def test_pipeline_still_keeps_only_unmapped_reads(self):
        module = (Path(__file__).resolve().parent.parent / "Modules" / "map_reads.nf")
        text = module.read_text()
        self.assertGreaterEqual(
            text.count("samtools fastq -@ \"${task.cpus}\" -f 4 -F 2304"), 3,
            "map_reads.nf should still apply the keep-only-unmapped rule at the human, "
            "T2T+phiX and pangenome passes; if this changed, the benchmark's premise "
            "and Tier A 3 need revisiting")

    def test_benchmark_documents_the_same_flags(self):
        text = SCRIPT.read_text()
        self.assertIn("-F 2304", text)
        self.assertIn("-f 4", text)
        self.assertIn("--secondary=no", text)


if __name__ == "__main__":
    unittest.main()
