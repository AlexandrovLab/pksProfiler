"""One sample failing must not stop the cohort reducers from ever launching.

Ludmil, revised report finding 8: only the original core per-sample processes
(extractReads, filterReads, mapReads, pksProfiler_align, pksProfiler_hmm) got
params.sample_failure_strategy through a withName selector in conf/base.config.
Every v0.0.2 sample stage added since -- targeted assembly, MAG binning and
classification, community taxonomy, strain typing -- inherited the process-wide
'finish' default, under which Nextflow submits no new tasks after any error. One
sample hitting a wall-time cap in, say, GTDB-Tk therefore took the whole cohort's
QC summary and gene-count tables with it, the exact failure mode the withName
selector already existed to prevent for the older lanes.

The fix is one shared label, `sample_stage`, carrying the same cohort-safe
errorStrategy, added to every independent per-sample process in those lanes --
and deliberately not to the cohort reducers that summarize across samples, or to
the one shared per-run resource (the DIAMOND clb database) every sample's rescue
step depends on.

His acceptance test: force one sample to fail in targeted assembly and another in
MAG/taxonomy -- other samples must complete and cohort reducers must still run.
These tests check the structural half of that (every relevant process carries the
label, the label carries the policy, the reducers/singleton do not); actually
forcing a task to fail needs a live Nextflow run, out of scope for a unit test.

----------------------------------------------------------------------------
rc2 review section 3, item 3.2 -- why this file was rewritten
----------------------------------------------------------------------------
The first version of this test carried a hand-maintained map of module file ->
process names. Three things follow from that, and all three bit us:

  * a module file absent from the map was not checked at all. Modules/
    extract_reads.nf, Modules/map_reads.nf, Modules/pks_carrier.nf,
    Modules/pks_specificity.nf and Modules/plot_bracken_taxa.nf were all absent;
  * a process absent from the map inside a file that WAS in the map was likewise
    not checked (pks_mag.nf's pairDepletedReads, for one);
  * a process added later was covered by nothing until somebody remembered to
    edit this list.

That is how `lowScoreRescueProfile` (Modules/extract_reads.nf -- a per-sample
process, `tag "$sampleID"`, one task per sample) shipped in rc1 and rc2 with no
cohort-safe failure policy at all: under the strict default, one sample failing
there stopped new task submission and took the cohort QC summary with it. The
allowlist could not catch it because extract_reads.nf was never in the list.

So the per-sample processes are now DISCOVERED from the source instead. The rule
is the one the pipeline itself uses: a process is per-sample when its input block
declares a per-sample identifier (`val(sampleID)` / `val sampleID`), because that
is exactly what makes Nextflow instantiate one task per sample. A cohort reducer
takes collected `path(...)` lists and no sample id, so it is never discovered.

Every discovered process must then be covered, by either:

  * `label 'sample_stage'`, or
  * the older `withName:` selector in conf/base.config, which carries the same
    errorStrategy for the five original core stages,

or appear in EXCLUSIONS below with a written reason. EXCLUSIONS is deliberately
empty today -- every per-sample process in the tree is covered -- and anything
added to it is a decision someone has to justify in the diff.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_DIR = ROOT / "Modules"
BASE_CONFIG = (ROOT / "conf/base.config").read_text()

PROCESS_RE = re.compile(r"^process\s+(\w+)\s*\{", re.M)
# `tuple val(sampleID), ...`, `val sampleID`, `val(sampleID)`. The pipeline spells the
# per-sample key `sampleID` everywhere; a new spelling is a deliberate change and should
# break this test rather than quietly drop a process out of the discovered set.
SAMPLE_KEY_RE = re.compile(r"\bval\s*\(?\s*sampleID\b")
LABEL = "label 'sample_stage'"

# Per-sample processes that are deliberately NOT covered by the cohort-safe policy.
# {process name: the reason, which must be a real one}. Empty is the correct state:
# every per-sample process in the tree today is covered. A process goes in here only
# when carrying on past its failure would produce wrong results rather than a smaller
# cohort -- and then the reason has to say why.
EXCLUSIONS = {}

# Processes that are not per-sample and must therefore keep the strict default: the
# cohort reducers that summarize across every sample, and the one shared per-run
# resource build. Listed by name so that a reducer silently acquiring the label (which
# would let a bad cohort table through) is caught too, not only the other direction.
STRICT_BY_DESIGN = {
    # Modules/pks_prefilter.nf -- one shared DIAMOND clb database every sample's rescue
    # depends on. Its failure genuinely should stop the run rather than let every sample
    # silently skip clb rescue.
    "buildClbDiamondDb",
    # Modules/pks_taxa.nf -- cohort reducers.
    "process_bracken", "combineClbTaxonomySupport",
    # Modules/plotting.nf -- cohort reducers.
    "masterTableAlign", "masterTableHMM", "masterQCSummary", "cohortReport",
    "masterSummary",
    # Modules/extract_reads.nf -- one cohort table of every input's header provenance.
    "inputHeaderSummary",
}


def process_blocks(text):
    """{process name: its source block} for one module file."""
    starts = [(m.start(), m.group(1)) for m in PROCESS_RE.finditer(text)]
    blocks = {}
    for index, (start, name) in enumerate(starts):
        end = starts[index + 1][0] if index + 1 < len(starts) else len(text)
        blocks[name] = text[start:end]
    return blocks


def input_block(block):
    """The process's input: declaration, or "" when it has none."""
    match = re.search(r"\n\s*input:\s*\n(.*?)\n\s*(output:|script:|exec:|shell:|stub:)",
                      block, re.S)
    return match.group(1) if match else ""


def all_processes():
    """{process name: (module file name, source block)} across every Modules/*.nf."""
    found = {}
    for module in sorted(MODULE_DIR.glob("*.nf")):
        for name, block in process_blocks(module.read_text()).items():
            assert name not in found, f"two processes named {name}"
            found[name] = (module.name, block)
    return found


PROCESSES = all_processes()
PER_SAMPLE = {name: (fname, block) for name, (fname, block) in PROCESSES.items()
              if SAMPLE_KEY_RE.search(input_block(block))}


def withname_selector_names():
    """The process names conf/base.config's withName: errorStrategy selector covers."""
    names = set()
    for match in re.finditer(r"withName:\s*'([^']+)'\s*\{(.*?)\n\s*\}", BASE_CONFIG, re.S):
        if "params.sample_failure_strategy" in match.group(2):
            names |= {part.strip() for part in match.group(1).split("|")}
    return names


SELECTOR_NAMES = withname_selector_names()


class DiscoveryItselfStillWorks(unittest.TestCase):
    """If the regexes stop matching, every assertion below passes vacuously."""

    def test_every_module_file_is_parsed(self):
        files = {fname for fname, _block in PROCESSES.values()}
        on_disk = {p.name for p in MODULE_DIR.glob("*.nf")}
        self.assertEqual(files, on_disk,
                         "a Modules/*.nf file contributed no process: the parser missed it")

    def test_a_plausible_number_of_processes_was_found(self):
        self.assertGreaterEqual(len(PROCESSES), 50, "process discovery found too few")
        self.assertGreaterEqual(len(PER_SAMPLE), 40, "per-sample discovery found too few")

    def test_known_per_sample_processes_are_discovered(self):
        # One from each lane, including the three the old hardcoded list never saw.
        for name in ("extractReads", "mapReads", "lowScoreRescueProfile", "plotPKS",
                     "megahitAssemble", "pairDepletedReads", "krakenPrefilter",
                     "classifyPksReadEvidence", "mlstTypeAssembly", "wesInputQC",
                     "alignCarrierPanel", "alignConcatenatedIndex", "plotBrackenTaxa"):
            with self.subTest(process=name):
                self.assertIn(name, PER_SAMPLE)

    def test_known_cohort_reducers_are_not_discovered_as_per_sample(self):
        for name in sorted(STRICT_BY_DESIGN):
            with self.subTest(process=name):
                self.assertIn(name, PROCESSES, f"{name} no longer exists; update the list")
                self.assertNotIn(name, PER_SAMPLE,
                                 f"{name} now takes a sample id: it is no longer a reducer")

    def test_the_withname_selector_is_found_and_names_real_processes(self):
        self.assertTrue(SELECTOR_NAMES, "no withName: sample_failure_strategy selector found")
        for name in sorted(SELECTOR_NAMES):
            with self.subTest(process=name):
                self.assertIn(name, PROCESSES)


class EveryPerSampleProcessCarriesTheCohortSafePolicy(unittest.TestCase):
    def test_every_discovered_per_sample_process_is_covered_or_excluded(self):
        uncovered = []
        for name, (fname, block) in sorted(PER_SAMPLE.items()):
            if name in EXCLUSIONS or name in SELECTOR_NAMES:
                continue
            if LABEL not in block:
                uncovered.append(f"{fname}:{name}")
        self.assertEqual(uncovered, [],
                         "per-sample process(es) with no cohort-safe failure policy: "
                         "add label 'sample_stage', or an EXCLUSIONS entry saying why not")

    def test_the_label_is_declared_at_most_once_per_process(self):
        for name, (fname, block) in sorted(PER_SAMPLE.items()):
            with self.subTest(file=fname, process=name):
                self.assertLessEqual(block.count(LABEL), 1,
                                     f"{name} in {fname} declares the label twice")

    def test_every_exclusion_is_a_real_per_sample_process_with_a_reason(self):
        for name, reason in sorted(EXCLUSIONS.items()):
            with self.subTest(process=name):
                self.assertIn(name, PER_SAMPLE, f"{name} is not a per-sample process")
                self.assertTrue(str(reason).strip(), f"{name} is excluded with no reason")

    def test_a_process_is_not_covered_twice(self):
        for name in sorted(SELECTOR_NAMES):
            with self.subTest(process=name):
                self.assertNotIn(LABEL, PROCESSES[name][1],
                                 f"{name} is in the withName selector AND labelled")


class ReducersAndTheSharedSingletonStayStrict(unittest.TestCase):
    def test_cohort_reducers_and_the_shared_db_build_do_not_carry_the_label(self):
        for name in sorted(STRICT_BY_DESIGN):
            fname, block = PROCESSES[name]
            with self.subTest(file=fname, process=name):
                self.assertNotIn(LABEL, block,
                                 f"{name} in {fname} should stay on the strict default")

    def test_no_process_without_a_sample_id_carries_the_label(self):
        """The label means "one sample's failure"; on a reducer it would hide a bad table."""
        wrong = sorted(f"{fname}:{name}" for name, (fname, block) in PROCESSES.items()
                       if name not in PER_SAMPLE and LABEL in block)
        self.assertEqual(wrong, [])


class TheLabelCarriesTheCohortSafePolicy(unittest.TestCase):
    def test_withLabel_sample_stage_block_exists(self):
        self.assertIn("withLabel:sample_stage", BASE_CONFIG)

    def test_it_uses_the_same_retry_then_configurable_pattern_as_the_core_lane(self):
        core = BASE_CONFIG.split("withName:", 1)[1]
        core_strategy = re.search(
            r"errorStrategy\s*=\s*\{[^}]*\}", core).group(0)
        stage = BASE_CONFIG.split("withLabel:sample_stage", 1)[1]
        stage_strategy = re.search(r"errorStrategy\s*=\s*\{[^}]*\}", stage).group(0)
        self.assertEqual(core_strategy, stage_strategy)
        self.assertIn("params.sample_failure_strategy", stage_strategy)

    def test_retryable_exit_codes_still_retry_before_falling_back(self):
        stage = BASE_CONFIG.split("withLabel:sample_stage", 1)[1]
        strategy = re.search(r"errorStrategy\s*=\s*\{([^}]*)\}", stage).group(1)
        self.assertIn("'retry'", strategy)
        self.assertIn("130..145", strategy)


if __name__ == "__main__":
    unittest.main()
