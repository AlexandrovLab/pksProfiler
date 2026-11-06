# pksProfiler v1.0.0

The first release. It is the tree the rc2 review's reproducibility and operational items
were closed on, and it is what every number we report was produced with.

**No analytical parameter, threshold or positivity rule changed anywhere in the rc line.**
A sample's tier under v1.0.0 is the tier v1.0.0-rc1 computed for it. The work since rc1 is
reproducibility, failure handling and documentation.

## What it does

Detects and characterises the 19 *clb* genes of the *pks* island in gut metagenomes, WGS,
WES and RNA-seq, with a separate evidence report per assay. Read profiling is the core;
metagenomes can add community assembly and genome binning, and short-read WGS can add
targeted reassembly of the island from its own reads.

## Reproducibility

- **Release CI runs the real execution tests on the tag itself**, not only unit tests and
  graph checks. The two checks a GitHub runner cannot perform — a zero-skip unit run, which
  needs R, and the metagenome MAG check, which needs ~100 GB of real GTDB-Tk, CheckM2 and
  geNomad data — are recorded on the cluster against a content digest of everything that
  changes a run. The gate recomputes that digest from the tagged tree and refuses a
  mismatch, an unsigned record, a failure, or a check that asserted nothing.
- **A skipped check must be a signed, written-down waiver.** The gate prints it on every tag
  build and still refuses a skip of the unit gate, a skip with no stated reason, and a
  failure dressed as a waiver.
- **Databases are frozen** in `conf/databases.tsv` with a content checksum of each one's
  release-defining file; every run recomputes a fingerprint and prints frozen or CHANGED.
- **Environments are pinned and solved.** All 18 conda environments pin exact versions, and
  `conda_envs/solved/` records what the solver actually installed — every package, build,
  channel URL and md5 — for the 11 the recorded release checks built.
- **Host-depletion references are documented**, including that the index named `phix`
  contains no phiX.
- **Resume behaviour is documented**: a changed scripts digest means `-resume` will not
  reuse most tasks from an earlier release.

## Robustness

- One sample failing in `lowScoreRescueProfile` no longer stops the cohort; the
  failure-policy test now discovers per-sample processes from the source rather than a
  hand-maintained list.
- MAG analysis survives a failed CheckM2 or GTDB-Tk task, and a run in which *every* sample
  fails upstream now ends with its WHAT FAILED report instead of crashing.

## Interpretation limits

These are operational evidence tiers, not calibrated probabilities and not a clinical
result. Read them with the cohort's QC: a missing result is not a negative.

- **Depth decides what a negative means.** A sample with no *clb* reads is only a negative
  if the run could have detected the island had it been there. Minimum analysable depths at
  a 5-read gate and 80% power are roughly 0.56M (WES), 1.6M (WGS) and 6.6M (metagenome)
  non-host reads; a metagenome needs about four times a tumour genome's depth to earn the
  word. v1.0.0 does **not** compute this per sample — it is an analysis you run on the tag.
- **The one-read gene floor is not inflating cohorts.** A *clb* gene counts as present at
  one assigned read. On an 840-sample WGS cohort, requiring two reads moves positivity from
  14.3% to 13.9%, and the strongest tier is unchanged; the sensitivity sits in the marginal
  tiers, where it should.
- **Duplicate collapse changes nothing measurable; fragment collapse does.** On the same
  cohort, no positive is lost to deduplication, while fragment collapse costs 9 of 120.
- **Cohorts are not pooled.** Aligner, duplicate marking and upstream filtering differ
  between cohorts in ways that can move prevalence on their own, so raw prevalence is not
  combined across them.

## Not in this release

Per-sample detection power and a depth-adjusted tier were implemented and then removed
before the tag: the depth-adjusted ladder needs a per-level assessability floor that does
not exist yet, and shipping half of it would have been worse than shipping none. The
analysis stands on its own and can be merged back when that floor is defined.

V8a locus specificity and V8b carrier attribution remain off by default and labelled
experimental; they are not calibrated on real controls.
