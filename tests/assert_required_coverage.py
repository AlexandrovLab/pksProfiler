#!/usr/bin/env python3
"""Run the test suite and refuse to call it green if the consequential tests were skipped.

rc2 review item 2.3:

    Fail the release if critical tests are skipped. The fast test environment can skip
    rescue, fragment-support and competitive-realignment tests when external tools are
    absent. Add an explicit skip-count/required-test assertion to the release job so
    'tests pass' cannot exclude the most consequential coverage.

Our own rc2 CI run reported `OK (skipped=13)` and nothing looked at the 13. Those skips
are not incidental: the tests that need a real samtools, minimap2, bedtools or
featureCounts are exactly the ones that execute the AS<30 rescue, the fragment/duplicate
collapse behind `tier_robustness`, and the competitive-realignment score. A release whose
evidence of correctness is a suite that quietly did not run them is not evidence.

So this is `python -m unittest discover`, plus three assertions:

  1. every module in REQUIRED ran, and nothing inside it was skipped;
  2. every other skip is in ALLOWED_SKIPS, with a reason written down here;
  3. the total skip count is at or below --max-skips (default: the size of ALLOWED_SKIPS).

Usage
    python tests/assert_required_coverage.py [--max-skips N] [--json FILE] [-q]

It exits non-zero on a failure, an error, a missing required module, an unexpected skip,
or a skip count over budget, and it prints the exact counts either way. --json writes a
machine-readable summary, which is what the release job records as the run's evidence.

Making a required test pass takes real tools on PATH or a built conda cache; see
tests/README.md. On TSCC:

    export CONDA_CACHE_DIR=/.../pksProfiler/conda_cache
    python tests/assert_required_coverage.py
"""
import argparse
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"

# Modules that must RUN, with no skipped test inside them. These are the ones Ludmil
# named, plus the two that the same evidence flows through -- and, deliberately, the ones
# whose skip conditions make them the easiest to lose:
REQUIRED = {
    # ---- rescue (Ludmil, 25 Sep E1 / 1 Oct L4) ----
    # The AS<30 low-score rescue is the single change that moved TCGA prevalence. These
    # execute it against a real samtools rather than asserting on the shell text.
    "test_low_score_extraction": "executes the rescue branch with a real samtools",
    "test_low_score_as_boundary": "the AS 29 / AS 30 boundary, executed",
    "test_low_score_rescue_trace": "the rescued reads' survival to a clb gene",
    # ---- fragment support (V2 robustness, decision A, 1 Oct) ----
    # tier_robustness is computed here; rc2 section A3 turns on what this measures.
    "test_fragment_support_variants": "real samtools markdup / bedtools / featureCounts",
    "test_tier_robustness": "the tier each collapsed count implies",
    # ---- competitive realignment (V8a) ----
    "test_competitive_realign_score": "real minimap2 alignments, both arms",
    # ---- the declared Nextflow floor ----
    # Without a nextflow binary this silently stops checking that the version the
    # manifest declares can compile the workflow at all -- the m-1 failure.
    "test_nextflow_2410_compat": "the workflow compiles on the declared Nextflow floor",
}

# Skips that are tolerated, each with the reason it is tolerated. A skip anywhere else
# fails the run. Keyed by module name.
ALLOWED_SKIPS = {
    "test_coverage_plot_boundaries":
        "needs Rscript. It checks the per-sample coverage FIGURE's axis boundaries; no "
        "number in any cohort table, tier or QC column comes from it, and the figure is "
        "regenerated from published data. Install R to run it.",
}


def module_of(test):
    """`test_foo` from a TestCase id like `tests.test_foo.Klass.test_bar`."""
    parts = test.id().split(".")
    for part in parts:
        if part.startswith("test_") and part != parts[-1]:
            return part
    return parts[0]


class Collector(unittest.TextTestResult):
    """TextTestResult, plus which module each run/skipped test belonged to."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ran_modules = set()
        self.skipped_by_module = {}

    def startTest(self, test):
        self.ran_modules.add(module_of(test))
        super().startTest(test)

    def addSkip(self, test, reason):
        self.skipped_by_module.setdefault(module_of(test), []).append((test.id(), reason))
        super().addSkip(test, reason)


def run_suite(verbosity):
    loader = unittest.TestLoader()
    # Same discovery as `python -m unittest discover -s tests -p 'test_*.py'`: tests/ is
    # not a package (no __init__.py), so it is its own top level and ids read
    # `test_module.Class.test_name`.
    suite = loader.discover(start_dir=str(TESTS), pattern="test_*.py")
    if loader.errors:
        for error in loader.errors:
            print(error, file=sys.stderr)
        raise SystemExit("tests could not be loaded")
    runner = unittest.TextTestRunner(verbosity=verbosity, resultclass=Collector)
    return runner.run(suite)


def check(result, max_skips):
    """[] when the run may be called green, else the reasons it may not."""
    problems = []

    if result.failures:
        problems.append(f"{len(result.failures)} test failure(s)")
    if result.errors:
        problems.append(f"{len(result.errors)} test error(s)")

    for module, why in sorted(REQUIRED.items()):
        if module not in result.ran_modules:
            problems.append(f"required module {module} did not run at all ({why})")
        elif module in result.skipped_by_module:
            for test_id, reason in result.skipped_by_module[module]:
                problems.append(f"required test skipped: {test_id} -- {reason} ({why})")

    for module, skips in sorted(result.skipped_by_module.items()):
        if module in REQUIRED or module in ALLOWED_SKIPS:
            continue
        for test_id, reason in skips:
            problems.append(f"unexpected skip, not in ALLOWED_SKIPS: {test_id} -- {reason}")

    total = len(result.skipped)
    if total > max_skips:
        problems.append(f"{total} skipped test(s), budget is {max_skips}")

    return problems


def summarize(result, max_skips, problems):
    return {
        "tests_run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "skipped": len(result.skipped),
        "max_skips": max_skips,
        "required_modules": sorted(REQUIRED),
        "required_modules_ran": sorted(m for m in REQUIRED if m in result.ran_modules),
        "skips": {module: [{"test": t, "reason": r} for t, r in skips]
                  for module, skips in sorted(result.skipped_by_module.items())},
        "problems": problems,
        "ok": not problems,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--max-skips", type=int, default=len(ALLOWED_SKIPS),
                        help="total skipped tests tolerated (default: %(default)s)")
    parser.add_argument("--json", help="write a machine-readable summary here")
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args(argv)

    result = run_suite(verbosity=1 if args.quiet else 2)
    problems = check(result, args.max_skips)
    summary = summarize(result, args.max_skips, problems)

    print()
    print("=" * 72)
    print(f"tests run      {summary['tests_run']}")
    print(f"failures       {summary['failures']}")
    print(f"errors         {summary['errors']}")
    print(f"skipped        {summary['skipped']} (budget {summary['max_skips']})")
    for module, skips in summary["skips"].items():
        for entry in skips:
            allowed = "allowed" if module in ALLOWED_SKIPS else "NOT ALLOWED"
            print(f"  [{allowed}] {entry['test']}: {entry['reason']}")
    print(f"required modules ran  {len(summary['required_modules_ran'])}"
          f"/{len(summary['required_modules'])}")
    print("=" * 72)

    if args.json:
        Path(args.json).write_text(json.dumps(summary, indent=2) + "\n")

    if problems:
        print("\nREQUIRED COVERAGE NOT MET:")
        for problem in problems:
            print(f"  - {problem}")
        print("\nThis run must not be released. Install the missing tools (or set "
              "CONDA_CACHE_DIR to a built pipeline conda cache) and run it again.")
        return 1
    print("\nRequired coverage met.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
