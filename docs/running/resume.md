# `-resume`, and why an rc1 results tree will not resume into v1.0.0

rc2 review item 2.8: *"Document rc1-to-rc2 resume behaviour. Because the scripts digest
changed, `-resume` will not transparently reuse most rc1 tasks. State this explicitly so
production reruns do not silently span mixed code states."*

**The short version.** Point `-resume` at a results tree produced by `v1.0.0-rc1` and
every task re-runs. That is the correct behaviour and it is not a bug — but the warning
matters, because the failure it protects against is silent: a tree half-written by rc1
and half by v1.0.0, with one `RUN_REPORT.txt` describing only the later session, is a
cohort nobody can attribute to a version string.

---

## What Nextflow actually hashes

A task is reused when its hash is unchanged. The hash covers, among other things:

| Ingredient | Changed between rc1 and v1.0.0? |
|---|---|
| the rendered process script | for some processes, yes — and **every** script that interpolates `${params.dep_digest?.scripts}` changed, because `scripts/` changed |
| the task's conda environment (hashed from the `.yml`'s **content**) | yes, for 14 of the 18 environments — rc3 pinned every package to an exact version |
| each input file's identity (path, size, mtime) | yes, for anything produced by a task that re-ran |
| process name, `ext`, container, module | no |

It does **not** cover `conf/*.config` resources or `errorStrategy`, `publishDir`, or
`beforeScript`. Changing a time limit or a retry policy does not invalidate anything.

Two of those ingredients are what make rc1 → v1.0.0 a clean break.

**1. The scripts digest (F15).** `params.dep_digest.scripts` is a content digest of the
whole `scripts/` directory, interpolated as a comment into every task script that calls
a helper from it. It exists so that editing a helper invalidates exactly the tasks that
use it — `scripts/` is not a declared input, so without it a changed helper would be
invisible to Nextflow. Five files under `scripts/` changed since rc1
(`build_mag_summary.py`, `build_qc_summary.py`, `finalize_provenance.py`,
`record_database_versions.sh`, and the new `describe_minimap2_index.py`), so the digest
changed, so **every** task carrying it re-runs — including tasks that never touch the
changed helpers.

**2. The pinned conda environments (rc3, T6).** Pinning every package to an exact version
rewrote 14 of 18 `conda_envs/*.yml`. Nextflow hashes a conda environment from the file's
content, so every task in those 14 environments re-runs even where its script is
byte-identical.

## Which processes re-run

Counted from the source, `v1.0.0-rc1..HEAD`: **53 of the 60 processes** have a changed
hash ingredient from the code alone.

The seven that do not are `filterReads`, `mergePksCandidates`, `megahitAssemble`,
`targetedPksMegahit`, `metabat2Bin`, `jgiContigDepths` and `magSampleStatus` — the
`fastp`, `fastplong`, `megahit` and `metabat2` environments were already pinned, and
their scripts did not change. **Six of the seven re-run anyway**, because their inputs
are produced by tasks that did change, and a re-run upstream task writes its output to a
new work directory. The only one that can genuinely be reused is `filterReads` on a
**FASTQ** run, where its input is the sample sheet's own file rather than a pipeline
product — and everything downstream of it still re-runs.

In practice: **plan for a full re-run.** There is no partial-reuse path worth engineering
around.

| Class | Re-runs? | Why |
|---|---|---|
| extraction, host depletion, profiling (align + HMM) | yes | changed env and/or scripts digest, and `pksProfiler_align`/`pksProfiler_hmm` changed text in rc4 |
| prefilter lane (`krakenPrefilter`, `diamondRescue`, …) | yes | `prefilter_env` pinned; `diamondRescue` also changed |
| MAG lane | yes | `checkm2_env`, `gtdbtk_env`, `genomad_env`, `prokka_env` pinned; `stubFailedMagQuality` is new |
| targeted assembly, strain typing, taxonomy, specificity, carrier | yes | pinned environments |
| every cohort reducer | yes | its inputs are new files |
| `filterReads` on FASTQ input | **no** | input is the sheet's FASTQ; `fastp_env` unchanged |

## How to tell, on a run you already have

* **`runs/provenance/*.json`** — one record per session. More than one file, with
  different `code.commit` values, is a tree spanning code states. This is the definitive
  check, and it works after the fact.
* **`RUN_REPORT.txt`** — written per session and describes only that session. A resumed
  tree's report is *not* a description of the whole tree; `runs/reports/` keeps the
  earlier ones for exactly this reason.
* **`runs/trace.txt`** — `CACHED` versus `COMPLETED` per task, for the session that just
  ran. A resume across versions shows almost no `CACHED` rows.
* **`code.dirty` / `code.modified_files`** in the provenance record — a run launched from
  a working tree with uncommitted edits is not "commit X", and says so.

## What to do instead

**For a production cohort: run it into a new `--outdir`, with a new work directory, on
the tagged version, and do not resume across the tag.** The cost is a full re-run; the
alternative is a results tree whose provenance cannot be stated in one line.

`-resume` remains exactly right for what it is for: picking a run back up after a
timeout, a node failure or a scheduler problem, **on the same code**.

If a resume across versions has already happened, the tree is not ruined — but it must be
described as what it is. `runs/provenance/` lists every session's commit, and
`runs/reports/` every session's report; a cohort covered by more than one commit should
be re-run on one before anything is published from it.

## Within v1.0.0

After the tag, the same rules apply to any bug-fix release: a change under `scripts/`
re-runs every task carrying the scripts digest, and a change to a `conda_envs/*.yml`
re-runs every task in that environment. `tests/check_resume_reference_change.sh` is the
execution proof that the mechanism works — it changes a reference's content, resumes, and
asserts that exactly the dependent tasks re-ran.

---

Back to the [README](../../README.md).
