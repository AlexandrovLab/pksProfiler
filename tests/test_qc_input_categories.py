"""E2 and the prefilter QC columns (Ludmil, 25 Sep).

E2: extraction kept only the idxstats grand total. It now keeps compact counts per
input-reference category, which must still add up to that total.
P1: the metagenome prefilter lane's Kraken routing, DIAMOND rescue and candidate
merge counts were written per sample but never reached the cohort QC table.
"""
import csv
import importlib.util
import re
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


idx = load("summarize_idxstats")
qc = load("build_qc_summary")

IDXSTATS = "\n".join([
    "chr1\t248956422\t1000\t10",
    "chrX\t156040895\t200\t2",
    "chrM\t16569\t50\t0",
    "chr6_GL000251v2_alt\t4795265\t30\t1",
    "chr1_KI270706v1_random\t175055\t7\t0",
    "chrUn_GL000220v1\t161802\t5\t0",
    "HLA-A*01:01:01:01\t3503\t4\t0",
    "chrEBV\t171823\t3\t0",
    "phiX174\t5386\t2\t0",
    "*\t0\t0\t99",
]) + "\n"


class Categories(unittest.TestCase):
    def test_names_map_to_their_category(self):
        cases = {"chr1": "canonical", "22": "canonical", "chrY": "canonical", "X": "canonical",
                 "chrM": "mito", "MT": "mito", "chr6_GL000251v2_alt": "alt",
                 "chr1_KI270706v1_random": "unplaced", "chrUn_GL000220v1": "unplaced",
                 "GL000192.1": "unplaced", "chrUn_JTFH01000001v1_decoy": "decoy",
                 "hs37d5": "decoy", "HLA-A*01:01:01:01": "decoy", "chrEBV": "viral",
                 "NC_007605": "viral", "chr23": "other", "phiX174": "other"}
        for name, expected in cases.items():
            with self.subTest(name=name):
                self.assertEqual(idx.category(name), expected)

    def test_categories_add_up_to_the_total(self):
        m = idx.summarize(IDXSTATS.splitlines())
        self.assertEqual(m["input_mapped_canonical"], 1200)
        self.assertEqual(m["input_mapped_mito"], 50)
        self.assertEqual(m["input_mapped_alt"], 30)
        self.assertEqual(m["input_mapped_unplaced"], 12)
        self.assertEqual(m["input_mapped_decoy"], 4)
        self.assertEqual(m["input_mapped_viral"], 3)
        self.assertEqual(m["input_mapped_other"], 2)
        self.assertEqual(m["input_unmapped_placed"], 13)
        self.assertEqual(m["input_unmapped_unplaced"], 99)
        parts = sum(v for k, v in m.items() if k != "input_alignment_records")
        self.assertEqual(parts, m["input_alignment_records"])
        # Same total the old single-awk version computed: columns 3 + 4 over every row.
        old = sum(int(l.split("\t")[2]) + int(l.split("\t")[3]) for l in IDXSTATS.splitlines())
        self.assertEqual(m["input_alignment_records"], old)

    def test_malformed_idxstats_is_refused(self):
        with self.assertRaises(ValueError):
            idx.summarize(["chr1\t100\tnot_a_number\t0"])


class CohortTable(unittest.TestCase):
    def table(self, metrics, sample="S1"):
        with tempfile.TemporaryDirectory() as directory:
            fragment = Path(directory) / "f.qc.tsv"
            fragment.write_text("Sample\tMetric\tValue\n" + "".join(
                f"{sample}\t{k}\t{v}\n" for k, v in metrics.items()))
            out = Path(directory) / "out.tsv"
            qc.write_summary(qc.load_fragments([fragment]), out)
            return list(csv.DictReader(out.open(), delimiter="\t"))[0]

    def test_input_categories_and_prefilter_counts_reach_the_table(self):
        m = idx.summarize(IDXSTATS.splitlines())
        m.update({"prefilter_input_reads": 100, "prefilter_primary_reads": 30,
                  "prefilter_non_target_reads": 70, "diamond_input_reads": 70,
                  "diamond_rescue_fragments": 4, "diamond_rescued_reads": 5,
                  "candidate_primary_reads": 30, "candidate_rescued_reads": 5,
                  "candidate_total_reads": 35})
        row = self.table(m)
        for key, value in m.items():
            with self.subTest(column=key):
                self.assertEqual(row[key], str(value))

    def test_categories_that_do_not_add_up_are_rejected(self):
        m = idx.summarize(IDXSTATS.splitlines())
        m["input_mapped_canonical"] += 1
        with self.assertRaisesRegex(ValueError, "input_alignment_records"):
            self.table(m)

    def test_candidate_merge_that_does_not_add_up_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "candidate_total_reads"):
            self.table({"candidate_primary_reads": 3, "candidate_rescued_reads": 1,
                        "candidate_total_reads": 5})

    def test_rescue_larger_than_its_pool_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "diamond_rescued_reads"):
            self.table({"diamond_input_reads": 3, "diamond_rescued_reads": 4})

    def test_a_partial_prefilter_group_is_reported_not_rejected(self):
        row = self.table({"prefilter_input_reads": 10})
        self.assertEqual(row["prefilter_primary_reads"], "NA")

    def test_every_new_column_is_written_by_some_producer(self):
        producers = "\n".join(p.read_text() for p in
                              list((ROOT / "Modules").glob("*.nf")) + list((ROOT / "scripts").glob("*.py")))
        for column in qc.OUTPUT_COLUMNS:
            if not column.startswith(("prefilter_", "diamond_", "candidate_")):
                continue
            with self.subTest(column=column):
                stem = column[len("prefilter_"):-len("_reads")] if column.startswith("prefilter_") else None
                self.assertTrue(
                    column in producers
                    or (stem and re.search(rf'"{stem}"', (ROOT / "scripts/route_kraken_reads.py").read_text())),
                    f"{column} is in the table but nothing writes it")


if __name__ == "__main__":
    unittest.main()


class BowtieRuleIsExplicit(unittest.TestCase):
    """A4 (25 Sep): the primary alignment's scoring rule is named and recorded."""

    MAIN = (ROOT / "main.nf").read_text()
    ALIGN = (ROOT / "Modules/pksProfiler_align.nf").read_text()

    def test_defaults_are_the_values_the_lane_always_used(self):
        for line in ('params.pks_align_bowtie2_mode      = "end-to-end"',
                     'params.pks_align_bowtie2_preset    = "very-sensitive"',
                     'params.pks_align_bowtie2_score_min = "L,-0.6,-0.6"',
                     'params.pks_align_min_mapq          = 40'):
            self.assertIn(line, self.MAIN)

    def test_the_command_uses_them(self):
        cmd = self.ALIGN.split("bowtie2 -x", 1)[1].split("samtools sort", 1)[0]
        for name in ("pks_align_bowtie2_mode", "pks_align_bowtie2_preset",
                     "pks_align_bowtie2_score_min", "pks_align_min_mapq"):
            self.assertIn(f"params.{name}", cmd)
        self.assertNotIn("--very-sensitive ", cmd)
        self.assertNotIn("-q 40", cmd)

    def test_run_report_records_them(self):
        prov = load("finalize_provenance")
        names = {n for _group, members in prov.PARAMETER_GROUPS for n in members}
        for name in ("pks_align_bowtie2_mode", "pks_align_bowtie2_preset",
                     "pks_align_bowtie2_score_min", "pks_align_min_mapq"):
            self.assertIn(name, names)

    def test_a_different_mapq_is_refused(self):
        self.assertIn('params.pks_align_min_mapq.toString() != "40"', self.MAIN)


class FragModeIsASwitch(unittest.TestCase):
    """E0 (25 Sep): --frag=no is available for the benchmark; the default is unchanged.

    Since the 26 Sep mode memo each pass's preset and flags come from
    lib/HostDepletion.groovy; the switch must still reach every call.
    """

    MAP = (ROOT / "Modules/map_reads.nf").read_text()
    MAIN = (ROOT / "main.nf").read_text()
    PLAN = (ROOT / "lib/HostDepletion.groovy").read_text()

    def test_every_depletion_call_takes_the_switch(self):
        calls = re.findall(r"minimap2 -2 -ax \$\{(\w+)\.preset\} \$\{\1\.args\}[^\n]*", self.MAP)
        self.assertEqual(len(calls), 3)
        self.assertEqual(set(calls), {"human", "t2t", "pangenome"})
        for plan in ("humanPass", "t2tPass", "pangenomePass"):
            with self.subTest(plan=plan):
                self.assertRegex(self.MAP, rf"HostDepletion\.{plan}\((params\.sample_type, )?"
                                           rf"params\.read_type, params\.minimap2_frag_mode")
                body = self.PLAN[self.PLAN.index(f"static Map {plan}("):]
                body = body[:body.index("\n    }\n")]
                self.assertIn("fragArg(rt, fragMode)", body)

    def test_the_default_is_frag_no_since_rc1(self):
        # E0 decision (29 Sep): route-invariant depletion by default; "default" stays
        # selectable for reproducing v0.0.1.
        self.assertIn('params.minimap2_frag_mode = "no"', self.MAIN)
        # --frag only concerns short paired reads, so long reads never get it.
        self.assertIn("(!isLong(readType) && fragMode?.toString() == 'no') ? '--frag=no' : ''", self.PLAN)

    def test_run_report_records_it(self):
        prov = load("finalize_provenance")
        names = {n for _group, members in prov.PARAMETER_GROUPS for n in members}
        self.assertIn("minimap2_frag_mode", names)
