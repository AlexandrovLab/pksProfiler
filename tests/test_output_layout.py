"""Output layout: three trees, and a sample's artefacts all in one of them.

v0.0.2 publishes into exactly `by_sample/<sample>/`, `cohort/` and `runs/`. The layout
this replaced split each sample across six top-level trees with inconsistent nesting --
`community_context/` carried an extra `tumor_wgs/` level on the tumour side but not the
metagenome side, `strain_typing/` carried none, `pks_per_sample/` was flat -- and put
figures in three separate places. Finding one sample's results meant knowing which stage
produced each file.

These tests lock the shape so it cannot drift back.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "main.nf").read_text()
MODULES = {p.name: p.read_text() for p in sorted((ROOT / "Modules").glob("*.nf"))
           if ".bak" not in p.name}

# Trees that existed before and must not come back as publish targets.
RETIRED_PARAMS = ("pks_dir", "pks_summary_dir", "pks_mag_dir", "pks_community_dir",
                  "pks_typing_dir", "prefilter_dir", "pks_coverage_plots_dir",
                  "pks_taxonomy_plots_dir", "unmapped_bam_dir", "mapped_reads_dir")

PROCESS = re.compile(r"^process\s+(\w+)", re.M)


def publish_targets():
    """{process: [publishDir directive, ...]} across every module.

    A directive is its own line plus any continuation lines, which is how the multi-line
    `publishDir(...)` form and trailing `saveAs:` closures are written. Brace matching is
    deliberately avoided: `${params.sample_dir}` contains a closing brace, so a naive
    `\{[^}]*\}` truncates the expression and silently passes assertions.
    """
    out = {}
    for text in MODULES.values():
        bounds = [(m.group(1), m.start()) for m in PROCESS.finditer(text)]
        for i, (name, start) in enumerate(bounds):
            end = bounds[i + 1][1] if i + 1 < len(bounds) else len(text)
            lines = text[start:end].splitlines()
            directives = []
            for j, line in enumerate(lines):
                if "publishDir" not in line:
                    continue
                block = [line]
                for cont in lines[j + 1:]:
                    if re.match(r"\s*(saveAs|pattern|mode|enabled|overwrite|//|\)|\"|\{)", cont):
                        block.append(cont)
                    else:
                        break
                directives.append("\n".join(block))
            out[name] = directives
    return out


TARGETS = publish_targets()

# Processes whose output spans samples, so they belong in cohort/, not by_sample/.
# tumorEligibilityStatus used to be here too; removed as dead code (finding d-1),
# it was never invoked and its publishDir target no longer exists.
COHORT_PROCESSES = ("masterTableAlign", "masterTableHMM", "masterQCSummary",
                    "cohortReport", "masterSummary",
                    "process_bracken", "combineClbTaxonomySupport")


class LayoutParamsTests(unittest.TestCase):
    def test_the_three_trees_are_declared(self):
        for decl in ('params.sample_dir = "${params.outdir}/by_sample"',
                     'params.cohort_dir = "${params.outdir}/cohort"',
                     'params.runs_dir   = "${params.outdir}/runs"'):
            self.assertIn(decl, MAIN)

    def test_retired_trees_are_gone(self):
        for name in RETIRED_PARAMS:
            self.assertNotIn(f"params.{name}", MAIN,
                             f"params.{name} was retired when the layout was flattened")

    def test_cohort_subdirs_hang_off_cohort(self):
        for decl in ('params.pks_counts_dir = "${params.cohort_dir}/gene_counts"',
                     'params.pks_taxonomy_dir = "${params.cohort_dir}/taxonomy"',
                     'params.pks_qc_dir = "${params.cohort_dir}/qc"'):
            self.assertIn(decl, MAIN)


class PublishTargetTests(unittest.TestCase):
    def test_every_publish_target_is_one_of_the_three_trees(self):
        allowed = ("params.sample_dir", "params.cohort_dir", "params.runs_dir",
                   "params.pks_counts_dir", "params.pks_taxonomy_dir", "params.pks_qc_dir")
        for process, targets in TARGETS.items():
            for target in targets:
                self.assertTrue(any(a in target for a in allowed),
                                f"{process} publishes outside the three trees: {target}")

    def test_no_process_references_a_retired_tree(self):
        for process, targets in TARGETS.items():
            for target in targets:
                for name in RETIRED_PARAMS:
                    self.assertNotIn(f"params.{name}", target,
                                     f"{process} still publishes to params.{name}")

    def test_cohort_processes_publish_to_cohort(self):
        for process in COHORT_PROCESSES:
            self.assertIn(process, TARGETS, f"{process} not found")
            for target in TARGETS[process]:
                self.assertTrue("cohort_dir" in target or "pks_counts_dir" in target
                                or "pks_qc_dir" in target or "pks_taxonomy_dir" in target,
                                f"{process} should publish under cohort/: {target}")

    def test_per_sample_processes_publish_under_by_sample(self):
        for process, targets in TARGETS.items():
            if process in COHORT_PROCESSES or not targets:
                continue
            for target in targets:
                self.assertIn("params.sample_dir", target,
                              f"{process} should publish under by_sample/: {target}")


class ArmLevelTests(unittest.TestCase):
    """The tumour lane had an extra tumor_wgs/ level the metagenome lane did not."""

    TUMOUR = ("prokkaTumorContigs", "hmmsearchTumorContigs", "genomadTumorContigs",
              "tumorContigContext")
    METAGENOME = ("genomadProphages", "communityProphageSummary", "extractGenomicContext")

    def test_no_process_inserts_an_arm_level(self):
        for process in self.TUMOUR + self.METAGENOME:
            for target in TARGETS[process]:
                self.assertNotIn("tumor_wgs", target,
                                 f"{process} still nests under an arm directory")

    def test_both_lanes_write_community_results_to_the_same_place(self):
        for process in self.TUMOUR + self.METAGENOME:
            for target in TARGETS[process]:
                self.assertIn("/community", target,
                              f"{process} does not publish under by_sample/<sample>/community/")

    def test_nothing_from_the_tumour_lane_lands_in_the_genomes_tree(self):
        # mags/tumor_wgs was a misnomer: nothing in the tumour lane is a draft genome.
        for process in self.TUMOUR:
            for target in TARGETS[process]:
                self.assertNotIn("/genomes", target,
                                 f"{process} is not a draft genome")


class FigureTests(unittest.TestCase):
    """Plots used to sit in coverage_plots/, taxonomy/plots/ and final_evidence/."""

    def test_every_figure_goes_to_the_sample_figures_dir(self):
        for process in ("plotPKS", "plotBrackenTaxa"):
            self.assertTrue(any("figures" in t for t in TARGETS[process]),
                            f"{process} does not publish into figures/")

    def test_the_contig_validation_plot_is_a_figure_not_a_result(self):
        targets = TARGETS["summarizeTargetedPksEvidence"]
        self.assertTrue(any("figures" in t for t in targets),
                        "the contig validation SVG should publish alongside the other figures")
        self.assertTrue(any("final_evidence" in t for t in targets),
                        "the evidence table should still publish under contigs/final_evidence")


if __name__ == "__main__":
    unittest.main()
