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
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_CONFIG = (ROOT / "conf/base.config").read_text()

MODULE_NAMES = ["pks_targeted.nf", "pks_mag.nf", "pks_prefilter.nf", "pks_taxa.nf",
                "pks_typing.nf", "plotting.nf"]
MODULES = {name: (ROOT / "Modules" / name).read_text() for name in MODULE_NAMES}

SAMPLE_STAGE_PROCESSES = {
    # plotPKS was missed by the original p1-2 fix: it's a per-sample process (one
    # circos plot per sample, no cross-sample input) sitting in plotting.nf
    # alongside the cohort reducers, and it inherited the strict default like they
    # did -- so a single node failure in one sample's plot took the whole cohort's
    # run down with it (observed live 2026-09-23, node failure on job 12331589).
    "plotting.nf": ["plotPKS"],
    "pks_targeted.nf": [
        "classifyPksReadEvidence", "targetedPksRecruit", "targetedPksMegahit",
        "targetedPksSpades", "alignTargetedContigsToCanonicalReference",
        "summarizeTargetedPksEvidence",
    ],
    "pks_mag.nf": [
        "megahitAssemble", "alignToContigs", "jgiContigDepths", "metabat2Bin",
        "checkm2Predict", "gtdbtkClassify", "prokkaAnnotate", "hmmsearchClb",
        "alignMagBinToCanonicalReference", "magBinLocusEvidence", "genomadProphages",
        "communityProphageSummary", "extractGenomicContext", "magSampleStatus",
        "magSummaryTable",
    ],
    "pks_prefilter.nf": ["krakenPrefilter", "diamondRescue", "mergePksCandidates",
                         "sampleBracken"],
    "pks_taxa.nf": ["extractPksIslandReads", "Bracken"],
    "pks_typing.nf": ["mlstTypeAssembly", "assignStrainType", "strainTypeSummary"],
}

# Cohort reducers (summarize across every sample) and the one shared per-run
# resource: failure there is not something to carry on past, so the default
# strict 'finish' behaviour is deliberate, not an oversight.
EXCLUDED = {
    # tumorEligibilityStatus used to be here (dead code, finding d-1); it has
    # since been removed from pks_mag.nf entirely, so there is nothing left to
    # exclude it from.
    "pks_prefilter.nf": ["buildClbDiamondDb"],  # shared DB every sample's rescue depends on
    "pks_taxa.nf": ["process_bracken", "combineClbTaxonomySupport"],
}


def process_block(text, name):
    start = text.index(f"process {name} {{")
    end = text.index("\nprocess ", start + 1) if "\nprocess " in text[start + 1:] else len(text)
    # search only for the next process boundary after start
    rest = text[start:]
    m = re.search(r"\nprocess \w+ \{", rest[1:])
    return rest if m is None else rest[:m.start() + 1]


class EveryIndependentPerSampleProcessCarriesTheLabel(unittest.TestCase):
    def test_every_named_process_declares_label_sample_stage(self):
        for fname, names in SAMPLE_STAGE_PROCESSES.items():
            text = MODULES[fname]
            for name in names:
                with self.subTest(file=fname, process=name):
                    block = process_block(text, name)
                    self.assertIn("label 'sample_stage'", block,
                                  f"{name} in {fname} is missing label 'sample_stage'")

    def test_the_label_is_declared_exactly_once_per_process(self):
        for fname, names in SAMPLE_STAGE_PROCESSES.items():
            text = MODULES[fname]
            for name in names:
                with self.subTest(file=fname, process=name):
                    block = process_block(text, name)
                    self.assertEqual(block.count("label 'sample_stage'"), 1)


class ReducersAndTheSharedSingletonStayStrict(unittest.TestCase):
    def test_cohort_reducers_and_the_shared_db_build_do_not_carry_the_label(self):
        for fname, names in EXCLUDED.items():
            text = MODULES[fname]
            for name in names:
                with self.subTest(file=fname, process=name):
                    block = process_block(text, name)
                    self.assertNotIn("label 'sample_stage'", block,
                                      f"{name} in {fname} should stay on the strict default")

    def test_plotting_reducers_are_also_unlabeled(self):
        plotting = (ROOT / "Modules/plotting.nf").read_text()
        for name in ("masterTableAlign", "masterTableHMM", "masterQCSummary", "cohortReport"):
            with self.subTest(process=name):
                self.assertNotIn("label 'sample_stage'", process_block(plotting, name))


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
