# Changelog

## v1.0.0-rc2 — 2026-10-14

Fixes on top of `v1.0.0-rc1` (4c73af5, including the tool bump). No analytical parameter or
positivity rule changed.

- **MAG lane no longer crashes when GTDB-Tk (or CheckM2) fails for a sample.** On the rc1 Mayo
  metagenome run (`--enable_mags --pks_taxa --pks_community_taxa`, prefilter balanced) several
  `gtdbtkClassify` tasks exited 1 under the per-sample `ignore` strategy and emitted no summary; a
  `remainder: true` join then produced `[sample, null, [bin.*.locus_evidence.tsv, ...]]` and
  Nextflow stopped the whole run after ~37 h (`Invalid method invocation 'call'`). Now every
  sample with bins but no CheckM2 report / GTDB-Tk summary gets a stand-in
  (`stubFailedMagQuality`) and is still reported in `genomes/pks_mag_summary.tsv` and the
  community tables: taxonomy reads `gtdbtk_failed` per bin, CheckM2 fields read `NA`. The failed
  task stays under RUN_REPORT "WHAT FAILED" and a warning names the sample. Join rows that still
  arrive left-orphaned are skipped with a warning, never destructured. Regression:
  `tests/check_mag_channel_regression.sh` (in `run_checks.sh`).
- **`pks_mag_summary.tsv` gains a last column, `mag_qc_status`**: `ok`, `checkm2_failed`,
  `gtdbtk_failed` or both, per bin (`NA` for the unbinned row); `gtdbtk_failed` is not flagged
  as an unexpected taxon. `tests/test_mag_tool_failure.py`. This is the only change under
  `scripts/`; because most tasks hash the scripts digest, a `-resume` across it re-runs them.
- **`diamondRescue` gets its own resource label** (`diamond_rescue`: 4 cpus / 24 GB / 8 h, each
  x attempt). It was `process_low` (8 GB) and hit exit 137 on several metagenomes; 137 retries
  (maxRetries 5) and each retry now asks for more memory.
- **F17, completed: Biowulf scheduler settings.** `conf/biowulf.config` now uses NIH's
  `queueStatInterval = '5 min'` and `submitRateLimit = '6/1min'` (were 2 min and 1 sec), as
  Ludmil asked; `pollInterval` stays 2 min.
- **CI runs on `integration/**` and `fix/**` branches and on `v*` tags.** The tagged rc1 commit
  had never run CI.
- **Legacy 50,767-bp island params removed.** `params.pks_shift` / `params.pks_island_len` were
  unused; the canonical 50,768-bp `pks_start_1based` / `pks_end_1based` /
  `pks_island_len_1based` are the only definition.
- **Docs.** README has a Software table of the exact pinned versions (tested against
  `conda_envs/*.yml`). README and `docs/prefilter.md` state that the default
  `--prefilter_mode auto` turns `balanced` on for every metagenome run given `--kraken_db`.
- Version 1.0.0-rc2 in the manifest, CITATION.cff, README badge and workflow figure.

## v1.0.0-rc1 — 2026-10-05

Extends the immutable v0.0.1 release with island reassembly, genome-resolved analysis and
strain typing. Full detail in `RELEASE_NOTES_v0.0.2-dev.md`.

**Tool bump** (`feature/tool-bump`, 2026-10-06; on `release/v1.0.0-rc1`; the manifest version is
unchanged). Detail, source references and A/B tables: `TOOL_BUMP_NOTES.md`.
- **Pins moved to the latest bioconda releases:** samtools 1.21/1.22 -> 1.24, bowtie2 2.5.4 -> 2.5.5, fastp
  0.24.0 -> 1.3.7, deeptools 3.5.6 -> 4.0.0, hmmer 3.3.2 -> 3.4, seqtk 1.4 -> 1.5, seqkit 2.6.1 -> 2.14.0,
  bracken 2.9 -> 3.1p1, spades 4.2.0 -> 4.3.0, metabat2 2.17 -> 2.18, mlst 2.23.0 -> 2.35.0, prokka 1.14.6 ->
  1.15.6, diamond 2.1.10/2.1.11 -> 2.2.8. minimap2 2.31, subread, bedtools, fastplong, krakenuniq,
  krakentools, megahit, checkm2, gtdbtk and genomad were already latest.
- **`pks_align_env` python 3.10 -> 3.12**: deeptools 4.0.0 ships only py312/py313 builds.
- **metabat2 is pinned `==2.18`**: `=2.18` resolved to the untagged `2.18_23_gc869c52` snapshot. Its version
  is now read from `metabat2 -h` / the depth tool's usage (it never answered `--version`, so every task
  recorded `version_unread`).
- **bracken 3.1p1 prints `v3.0.1`** (upstream v3.1 never updated VERSION). Strict `--tool_version_check`
  would have stopped every Bracken task. `scripts/tool_version_aliases.tsv` lists such known self-reports,
  one conda version each. `record_tool_versions.sh` and `finalize_provenance.py` both read it.
  `tests/test_tool_versions.py`.
- **No CLI change needed**: every flag the pipeline passes is unchanged in the new releases, and the fastp
  JSON layout is the same. Behaviour that can move counts: fastp 1.x turns poly-G trimming on for more
  first-read-name prefixes (any `@A...`, `@NA`, `@VL`, `@VH`, `@LH`, `@FS`, `@MN`, `@SH`; file input only), and
  its adapter matching tolerates one indel (with `--adapters`, or a detected SE adapter). The bracken G-level
  `-w` report's per-node column changed; the abundance and mpa tables are identical. On every synthetic e2e
  fixture (7 runs, 18 sample-runs) the QC, master, evidence and counts tables are cell-identical to the
  pinned tools.

**Fixed in rc1.1** (on top of `integration/v1.0.0-rc1`; the manifest version is unchanged)
- **Tool-version capture read `?` for a few tasks** (`version_unread`). On the rc1 TCGA WES
  run (1,304 samples) bamCoverage read `?` in 11 tasks, fastp 8, samtools 4 (one env),
  bowtie2 3, bgzip 1, and RNA python3 in 1 of 3, while every other task of the same tool
  and env read fine. The parse is deterministic, so in those tasks the one version call
  printed no version line: `scripts/record_tool_versions.sh` made a single attempt under a
  fixed 60 s `timeout`, with no retry, and threw away both the exit status and every line
  without a version number, so the record could not say why. The most likely cause is
  a cold start (exec, shared libraries and, for bamCoverage, ~30 Python modules from the
  conda env on the shared filesystem) that ran past the timeout while the whole cohort
  started at once, or a transient exec failure. Each read now captures stdout and stderr,
  runs under `timeout -k 5` (60 s, then 120 s, 180 s), and is retried up to 3 times with
  a 2 s / 4 s backoff when it exited non-zero or printed nothing. Output with no version
  and exit 0 is not retried, so a fast tool still costs one call. The record gains
  `attempts` and `read_exit`, and a `version_unread` row keeps the exit status and the
  first line of output in `raw`. Settings: `PKS_TOOL_VERSION_ATTEMPTS`, `_TIMEOUT`,
  `_BACKOFF`. Reading a version still never fails a task. `tests/test_tool_versions.py`
  (`AVersionReadIsRetried`: a tool that fails its first 1 or 2 calls, a hung first call,
  a transient exec error, a version on stderr, a tool that always fails).
- **Two processes published onto one path.** The `taxa_barplots.pdf` fix (6a03115) is now
  checked for every process: `tests/test_unique_publish_paths.py` works out each publishDir's
  published path templates (directory, output names, `pattern`, `saveAs`) and fails on any
  two that coincide, or on a process included under several aliases whose call sites
  do not each set its published name. It found the same overwrite one stage earlier:
  `Bracken` (`--pks_taxa`, island reads) and `sampleBracken` (`--pks_community_taxa`, whole
  community) both published `by_sample/<s>/taxonomy/bracken.{G,S}.report.txt` and
  `.krakenreport.txt`. sampleBracken's are now `community_bracken.*`. Nothing read those
  per-sample files by path; `docs/output.md` and `docs/running/taxonomy.md` list both.
- `tests/fixtures/assert_normal_parity.py` compares tool-version records without the new
  `attempts` / `read_exit` columns: a retried read is not a different binary.

**Added**
- **Prebuilt RNA references.** `scripts/download_rna_references.sh <dir>` fetches the GRCh38.d1.vd1
  `-x splice:sr` minimap2 index and the GENCODE v36 junction BED (~12 GB, public) from the shared Google Drive
  folder and checks fixed MD5s; `docs/running/wes_rna.md` documents it alongside the build commands.
- **Tier robustness columns** (V2; Ammal's decision A, 1 Oct: the tier stays on raw reads).
  `classifyPksReadEvidence` now runs `scripts/compute_fragment_support_variants.py` for every
  tiered (short-read wgs) sample with reads, with the tier's own thresholds, island depth and
  gene-overlap rule, and `read_evidence.tsv` gains `pks_fragments` / `tier_if_fragments`,
  `pks_reads_dedup` / `tier_if_dedup`, `pks_fragments_dedup` / `tier_if_fragments_dedup` and
  `tier_robustness` (stable | unstable). `read_evidence` and `pks_reads` are unchanged; the
  variant script's own raw tier must equal `read_evidence` or the task fails. Carried into
  `pks.qc.summary.tsv`, `pks.master_summary.tsv`, the cohort report (table and TSV) and
  RUN_REPORT (counts). `NA` when untiered. `tests/test_tier_robustness.py`.
- **Low-score rescue traced through every stage** (L4, Ludmil 1 Oct). `extractReads` lists the
  names of the reads E1 added (optional output, absent when none); `mapReads` counts them on
  the streams it already reads, after fastp, hg38 and T2T+phiX; the new `lowScoreRescueProfile`
  counts them in the profiling BAM and featureCounts-assigned to *clb* with the pipeline's own
  annotation filter and arguments. `low_score_rescued_*` columns in `pks.qc.summary.tsv` (only
  when some sample rescued a read, so rule-off and FASTQ tables are unchanged), checked to
  only decrease; RUN_REPORT **LOW-SCORE RESCUE (L4)**. Asserted on the bwa -T 0 e2e run;
  `tests/test_low_score_rescue_trace.py`.
- `tests/test_low_score_as_boundary.py` (L5): AS 29 is extracted and AS 30 is not, on both
  mates and an unpaired record, through the module's own on-branch pipeline and the L4 name
  list; the boundary moves with `--extract_min_as`.
- **Runtime tool versions, checked against their pins** (F16, "tool versions"). RUN_REPORT's
  SOFTWARE section listed what conda installed; on TSCC a user minimap2 2.24 earlier on PATH
  ran every host-depletion pass while it said 2.31. Every process that runs an external tool
  now starts with `scripts/record_tool_versions.sh` (rendered by `lib/ToolVersions.groovy`),
  which asks each binary on the task's PATH, after conda activation, for its version and
  location and writes `<work dir>/.tool_versions.tsv`. `--tool_version_check strict`
  (default) stops the run (exit 78, `terminate` in `conf/base.config`) when a tool ran from
  outside the task's conda env, is not the version its `conda_envs/*.yml` pins, or is not on
  PATH; `warn` records and carries on. `finalize_provenance.py` reads every task's record
  through the trace's new `workdir` column (cached tasks included), compares it again with
  the yml pins and the env's conda-meta, and writes **TOOL VERSIONS (as run)** and, when
  anything disagrees or a tool-running task left no record, **TOOL VERSION MISMATCH** to
  `RUN_REPORT.txt` (and `tool_versions` to the JSON record); it exits non-zero on a mismatch
  in strict mode. Checked on every mode and e2e execution run
  (`tests/fixtures/assert_tool_versions.py`); unit tests in `tests/test_tool_versions.py`.
- Matched normals are processed identically to tumours, now tested (Ludmil, 26 Sep mode
  memo). The pipeline carries no tumour/normal role: a row is `patient` plus its input, and
  the identifier is only a name. `tests/test_matched_normal_parity.py` proves that by
  construction (the sheet columns read, the per-sample tuple, no branch on an identifier, the
  run-wide `tumor_*` parameters, and an enumerated list of every per-sample divergence, each
  decided by the sample's own data under the same rule); `tests/check_mode_execution.sh`
  mode `normal_parity` runs one synthetic BAM as a tumour and as its normal in one run and
  requires identical processes, task scripts, tool records, outputs and QC rows modulo the
  sample name.
- E1 QC and provenance (Ludmil, 25 Sep): `extracted_low_score_mapped_reads`, `extracted_reads`,
  `extract_low_score_no_as_records` and `extract_low_score_mode` in `pks.qc.summary.tsv` (with an
  add-up check); `aligner_bwa_mem`, `bwa_min_score`, `extract_low_score_mapped`, `extract_min_as`
  in `provenance/input_header.tsv`; the two parameters and a per-state count in `RUN_REPORT.txt`.
- Two separate optional checks on the retained pks-island reads (Ludmil, 25 Sep report
  section 4, V8a and V8b), both off by default and deliberately never reported as one
  number. `--pks_locus_specificity` (V8a, `Modules/pks_specificity.nf`) asks whether a
  supporting read really belongs to the pks/clb locus, re-challenging it against a
  pks-NEGATIVE panel -- pks-negative Enterobacterales, NRPS/PKS homolog genomes and
  non-colibactin MIBiG clusters, gut and oral backgrounds, UniVec, optionally GRCh38 --
  and reporting SURVIVES / NOT_DISTINGUISHABLE_FROM_CROSS_MAPPING / UNDERPOWERED with
  `arms_used`, since Arm A (KrakenUniq) is optional and an Arm-B-only verdict says so.
  `--pks_carrier_attribution` (V8b, `Modules/pks_carrier.nf`) asks which organism carries
  it, against a pks-POSITIVE carrier panel (pks+ *E. coli* beyond IHE3034, pks+
  *Klebsiella*, *Citrobacter*, *Enterobacter*), reporting CARRIER_ATTRIBUTED /
  CARRIER_GENUS_ONLY / CARRIER_AMBIGUOUS / CARRIER_UNDETERMINED with strain and genus
  together. Both panels are built by `scripts/build_specificity_panels.py` from public
  sequence, with membership screened against the canonical island in opposite directions
  and every accepted and rejected candidate recorded. New per-sample outputs
  (`locus_specificity/`, `carrier_attribution/`), two column groups in
  `pks.master_summary.tsv`, two cohort-report sections, both panels tracked as F15
  dependencies. See `docs/running/specificity.md`.
- Long-read input, `--read_type short|ont|hifi` (owner decision, 29 Sep; wgs and rna only):
  fastplong read QC (`--long_read_min_length`, QC metric names kept), map-ont/map-hifi depletion
  (splice/splice:hq for the RNA human pass) against `--hg38_db_long` / `--rna_splice_db_long` /
  `--t2t_phix_db_long`, minimap2 profiling against the IHE3034 FASTA with featureCounts `-L`.
  Long-read WGS is untiered; targeted assembly (on by default) is switched off with a log
  line and refused only when `--tumor_targeted_assembly true` is on the command line;
  `--pks_taxa` / `--pks_community_taxa` are refused (Bracken's read-length model is short-read). `pks_align_env` gains minimap2 2.31;
  new `conda_envs/fastplong_env.yml` (fastplong 0.7.1). See `docs/running/long_reads.md`.
- RNA-seq lane, `--sample_type rna` (Ludmil, 26 Sep mode memo section 4): the human depletion
  pass is splice-aware (`minimap2 -ax splice:sr` against `--rna_splice_db`, `-j --rna_junctions_bed`
  when given; `--hg38_db` not required); T2T+phiX and pangenome passes stay `-ax sr`.
  `--rna_library_type polyA|total|unknown`. Profiling unchanged; `rna_evidence.tsv`, no tier, no
  assembly, no MAGs. Presets per mode live in `lib/HostDepletion.groovy`, executed by
  `tests/check_depletion_plan.sh`.
- WES lane, `--sample_type wes` (Ludmil, 26 Sep mode memo section 3): standard depletion and
  profiling, `wes_evidence.tsv`, no tier, no targeted assembly, bowtie2/both only. WES QC columns:
  input duplicate fraction (submitter's 0x400), on/off-target fraction with `--wes_capture_bed`
  (digest in provenance), pks-island duplicate fraction by position collapse. RNA adds
  `clb_transcripts_detected` (`--rna_transcript_min_reads`, default 3) and a residual-human
  fraction on a capped sample of the final non-host reads (`docs/running/wes_rna.md`).
- Sample-type routing: `--sample_type metagenome | wgs | wes | rna` (required, no default).
- Evidence tiers from read count, *clb* gene count and island breadth together, replacing a
  two-tier scheme whose lower boundary had been rejected during threshold analysis.
- Targeted island reassembly for tumour data using MEGAHIT and metaSPAdes independently, with
  an assembler-agreement call.
- Genome-resolved analysis for metagenomes: binning, completeness, taxonomy, annotation and
  *clb* protein search.
- Community context: prophage inventory, producer-versus-neighbour interactions, island
  mobility. Published to `community_context/`, separate from the draft genomes.
- Strain typing on assembled output: sequence type, clonal complex and phylogroup against a
  45,761-genome panel. Refuses to report a type from an incomplete allele profile.
- Optional taxonomic prefilter with translated-homology rescue, off by default and resolved
  automatically for metagenomes when a Kraken database is supplied.
- `scripts/build_master_summary.py`, joining every stage into one row per sample.
- Population-scale *clb* nucleotide and protein profile models; see `ref/hmm/PROVENANCE.md`.
- Continuous integration, `CITATION.cff`, `CITATIONS.md`.
- `tests/test_version_consistency.py` (GAMEPLAN rc1 release mechanics): fails if the
  `nextflow.config` manifest version, `CITATION.cff`, the top CHANGELOG heading, the README
  badge, the workflow figure, `--help` or `RUN_REPORT.txt` disagree. RUN_REPORT already carries
  the git commit (and a dirty-tree flag). The workflow figure's legend now reads the manifest
  (it said v0.0.2 while the manifest said 0.0.2dev).

- QC terminal-state schema (GAMEPLAN rc1): `terminal_state` in `cohort/qc/pks.qc.summary.tsv`
  and `pks.master_summary.tsv`, exactly one of `completed_positive`, `completed_negative`,
  `completed_untiered`, `failed_at_<stage>` (first missing stage in pipeline order) or
  `excluded_<reason>`. `classifyPksReadEvidence` emits an `evidence_outcome` QC fragment.

- E3 input-alignment provenance (Ludmil, 25 Sep): for BAM/CRAM input, `extractReads` records
  each header's `@SQ` fingerprint (md5 over ordered SN:LN), M5-set fingerprint, contig count,
  build guess and `@PG` chain (ID/PN/VN along PP, never CL) in
  `by_sample/<s>/provenance/input_header.tsv`; `cohort/provenance/` joins them and a
  counts-only summary that `RUN_REPORT.txt` prints as INPUT ALIGNMENTS
  (`scripts/header_provenance.py`, adapted from the E3 validation script).

- Metagenome paired assembly and abundance (Ludmil, 26 Sep mode memo section 5): the MAG
  lane rebuilds R1/R2 + singletons from the host-depleted stream by mate suffix
  (`recover_recruited_mates.py --all-reads`, external sort, every read kept) and feeds MEGAHIT
  `-1/-2/-r` and read-to-contig alignment as pairs; `genomes/assembly/mate_pairing.tsv` carries
  the A7 mate-tag diagnostics (also added to targeted recruitment stats; duplicate mate slots
  are kept as singletons instead of overwriting). `pks_reads_per_million` (per million
  host-depleted reads) and `host_depleted_reads` join metagenome `read_evidence.tsv`, the
  cohort report and the master summary; raw `pks_reads` kept.

- `tests/check_mode_execution.sh` (26 Sep mode memo section 6): one real execution per mode --
  wgs, wes, rna (splice:sr + junction BED on a GRCh38 slice), metagenome, wgs+ont, rna+ont --
  asserting the branch, outputs, disabled features, depletion preset and terminal state.
- `process.beforeScript` puts the conda base first on PATH so each task's activated env wins
  over tools earlier on the launching shell's PATH (a user minimap2 2.24 had been shadowing
  the pinned 2.31 in every host-depletion pass on TSCC).
- Divergence-from-IHE3034 QC (Ludmil, 25 Sep A5): `pksProfiler_align` runs
  `scripts/ihe3034_divergence_qc.py` on the counts BAM and publishes
  `alignment/ihe3034_divergence.qc.tsv` (a QC fragment) and `alignment/ihe3034_divergence_by_gene.tsv`.
  Reads measured are the primary MAPQ >= 40 alignments featureCounts assigns to a clb gene;
  divergence is NM / alignment columns. `pks.qc.summary.tsv` gains `divergence_ihe3034_reads`,
  `nm_ihe3034_median`, `divergence_ihe3034_median/_q1/_q3`, `divergence_ihe3034_frac_gt_2pct/_5pct`
  and `clb_read_length_median/_q1/_q3/_min/_max`, and RUN_REPORT describes them. QC only: no
  count, tier or call changes, and it is not a specificity or identity measure.

**Changed**
- `targetedPksRecruit` streams Bowtie2's SAM straight into `recover_recruited_mates.py --sam -`
  instead of writing it uncompressed to scratch and reading it back (A8, Ludmil 25 Sep).
  Recruitment IDs and stats are byte-identical (`tests/test_recruitment_streaming.py`, and the
  wgs_assembly mode run before and after).
- The bamCoverage track is `alignment/coverage_50bp_bins.bedGraph` (was `coverage.bedgraph`),
  with `--binSize 50` stated (A9, Ludmil 25 Sep): mean depth per 50 bp bin over the whole
  IHE3034 genome, for the circos figure only; `coverage.txt` is per-base depth over the island,
  which breadth and tiers use. `docs/output.md` says how they differ.
- `tests/run_checks.sh`: each validation check gets its own Nextflow log and also accepts the
  expected message from that log's ERROR line, because Nextflow intermittently exits after
  logging the error but before printing `ERROR ~` (seen on the base commit too). The
  non-zero exit is still required.
- **Extraction from bwa mem alignments also takes low-score mapped reads** (Ludmil, 25 Sep E1):
  `extractReads` adds primary mapped records with `AS:i` < `--extract_min_as` (30) to the FLAG 4
  reads, in the same single decode. bwa leaves a read unmapped only below `-T` (default 30);
  GDC aligned TCGA with `-T 0`, so about a third of genuine *pks* reads were "mapped" and never
  extracted. `--extract_low_score_mapped auto` (default) turns this on per sample when the
  header's `@PG` chain runs `bwa mem`, off otherwise (STAR, minimap2, no `@PG`); `on`/`off`
  force it. `extracted_unmapped_reads` keeps its meaning (FLAG 4 only); the extraction total in
  the attrition chain is now `extracted_reads`. Validated in `validation/E1` (synthetic: 99.9-100%
  of hidden reads recovered; 20-sample TCGA pilot: tumours gained evidence, matched normals no
  *clb* read). See `docs/running/bam.md`.
- **Read-level taxonomy is not carrier evidence** (Ludmil, 26 Sep mode memo section 5):
  `prefilter/diamond_rescue_taxonomy.tsv` column `organism` is renamed `read_taxon` (the
  cohort report still reads the old name). README, `docs/running/taxonomy.md`,
  `docs/prefilter.md`, `docs/output.md`, the cohort report's species and DIAMOND-rescue panels
  and the Bracken barplot titles call Kraken/Bracken/DIAMOND-rescue output "read-level
  taxonomic support"; carrier language is kept for contig and genome-bin evidence.
  `tests/test_taxonomy_wording.py` scans the docs and source for carrier language beside
  read-level taxonomy, and for NM/divergence text claiming specificity or uniqueness.
- **F21, real execution in CI**: the execution checks now run on every push to `main` and
  nightly, not only on manual dispatch, with the pipeline's conda environments cached
  between runs and the outputs retained as build artefacts. Two new execution checks:
  `tests/check_resume_reference_change.sh` (change a reference's content, `-resume`, and
  assert that exactly the tasks carrying its F15 digest re-ran while the rest stayed
  cached -- the first execution proof of the F15 mechanism) and
  `tests/check_v8_lanes_execution.sh` (one real run of each V8 lane on synthetic reads
  whose correct answer the fixture knows).
- **Targeted recruitment specificity and classes** (Ludmil, 25 Sep A1/A2): recruitment Bowtie2
  drops `-k 1` (best-alignment reporting, so MAPQ is meaningful) and every read needs MAPQ >=
  `--recruit_min_mapq` (20). Reads are classed `island_internal` (query coverage >=
  `--recruit_min_query_coverage`, 0.80), `boundary_spanning` (crosses an island boundary, or
  reaches it within `--recruit_boundary_tolerance_bp` (10) with the clip pointing out; no coverage
  floor) or `flank_only`; a 150 bp read with a 60 bp internal local match is now rejected.
  `pks_recruitment.tsv` reports reads/fragments per class and `island_supporting_fragments`
  (internal + boundary); `final_pks_evidence.tsv`, the master summary and the report carry them
  plus `assembly_input_diagnostic`, which counts island-supporting fragments only. Flank-only
  fragments still reach the assemblers (`--recruit_include_flank_only`). Structural confirmation
  input only; primary evidence and positivity are unchanged.
- **Host depletion defaults to `--frag=no`** (`--minimap2_frag_mode no`; E0, owner decision 29 Sep):
  every short read is judged on its own, so depletion no longer depends on read order / input route.
  On all 100 TCGA CRC positives the old default kept a route-dependent read set in 46 samples and
  `--frag=no` in none; no pks read count or tier changed in any of 200 runs. `--minimap2_frag_mode default`
  reproduces v0.0.1.
- Mode-specific evidence schemas (Ludmil, 26 Sep mode memo sections 3, 4, 6):
  `classify_tumor_pks_evidence.py --mode {wgs,metagenome,wes,rna}` writes `read_evidence.tsv`
  (wgs tiered, metagenome untiered), `wes_evidence.tsv` or `rna_evidence.tsv`. Every file has a
  `mode` column (second); `pks.master_summary.tsv` and `pks_cohort_report.tsv` gain `mode` and
  `evidence_label`, and WES/RNA rows are shown as their own evidence, never as a genomic tier.
- `--sample_type` values renamed (owner decision on Ludmil's 26 Sep mode memo):
  `tumor_wgs` → `wgs`, `tumor_wes` → `wes`, `tumor_rna` → `rna`; `metagenome` unchanged. The old
  values are refused before any task runs with `--sample_type tumor_wgs was renamed to wgs`.
  The tier threshold parameters keep their `tumor_*_min_*` names and apply to `wgs` only.
- minimap2 pinned to **2.31** for every mode (`conda_envs/minimap2_env.yml`,
  `conda_envs/targeted_alignment_env.yml`; was 2.28), so WGS, WES, RNA and metagenome runs
  deplete with one aligner version and the RNA lane's `-ax splice:sr` preset is available.
- Island breadth tier thresholds converted for the change of measurement axis: `broad_island`
  0.10 -> **0.075**, `extensive_island` 0.20 -> **0.15**. v0.0.1 measured breadth from a 50 bp
  binned coverage track; v0.0.2 measures per base, which reads a median 0.794x lower on the same
  data. This is a unit conversion preserving v0.0.1's effective stringency, not a recalibration.
  On 420 TCGA tumours it returns 8 boundary samples to the tier v0.0.1 gave them. Specificity is
  unaffected either way: 0 of 420 matched normals are positive at breadth >= 0.05 and above.
- All MAG output consolidated under `mags/`; `pks_summary/mag/` is gone.
- Profile HMM defaults now point at the population-scale models. The v0.0.1 benchmark model
  is retained and selectable with `--hmm_model`.
- `nextflow.config` carries a real manifest, so the Nextflow version requirement is enforced
  rather than documented.

**Removed**
- Container profiles (`docker`, `singularity`, `podman`, `shifter`, `charliecloud`, `arm`).
  No process declared a container image, so they could not work; they now fail loudly as
  unknown profiles instead of silently running on the host.
- nf-core template placeholders and the unused `input` parameter.

**Fixed**
- **V8a counted reads it never measured as lost to a decoy** (V8 real-data pilot, 5 Oct).
  minimap2 `-x sr` pairs two consecutive records named `X/1` and `X/2` and strips the
  suffix; island reads come from a coordinate-sorted BAM whose names end in `/1` `/2`, so
  adjacent mates had no Arm B row under their read_clb_gene name and were reported as
  non-survivors absorbed by `unknown_missing_from_arm_b` (and, on synthetic reads, turned
  a 100%-retention sample into NOT_DISTINGUISHABLE at R=0.66). Every V8a and V8b minimap2
  call now runs `--frag=no`. `specificity_decision.py` defines retention (and its CI)
  over measured reads, reports unmeasured reads by kind, keeps them out of the absorber
  tally, adds `retention_if_unmeasured_lost`, and returns `INCOMPLETE_MEASUREMENT` above
  `--pks_specificity_max_unmeasured_fraction` (0.05). The V8 lane fixture now writes
  adjacent `/1` `/2` mates so `tests/check_v8_lanes_execution.sh` catches a regression.
- **V8b CARRIER_UNDETERMINED conflated "no evidence" with "evidence that does not
  separate carriers".** It now means fewer than the minimum reads *aligned* to the panel;
  many aligned reads that all tie across genera are `CARRIER_AMBIGUOUS`. The summary adds
  `n_mapped`, `n_off_map` and `tied_genus_sets`, and an alignment to a contig absent from
  the carrier map (index and map from different builds) fails the task instead of being
  counted as unmapped.
- `summarize_tumor_pks_contigs.py` crashed (`dict contains fields not in fieldnames`) on any
  read_evidence column it did not repeat -- the same failure the `mode` column caused before.
  It now ignores extra evidence columns (found by the V2 columns).
- `tests/check_v8_lanes_execution.sh` runs `Modules/pks_specificity.nf` and `pks_carrier.nf`
  standalone, which loads no `nextflow.config`: conda was off, so the V8 tasks ran whatever
  minimap2/python3 the launching shell had (from a clean PATH: `minimap2: command not
  found`). Its config now enables conda with the same cache and PATH fix; `Modules/lib` links
  to `../lib` so the standalone modules can render the tool-version record.
- `scoreCompetitiveRealignment` (V8a) ran `python3` in `minimap2_env.yml`, which has no python,
  so it used whichever python the launching shell had. It now runs in `pks_align_env`
  (python 3.10, samtools 1.21), like the verdict step after it. Found by the tool-version
  record's construction test.
- An unaligned BAM/CRAM (no `@SQ` header lines) failed extraction with `samtools quickcheck`'s
  bare "had no targets in header". `extractReads` now says it is unaligned and to run it as
  FASTQ; documented in the README input section and `docs/running/bam.md`.
- `plotPKSTaxa` (`--pks_taxa`) and `plotCommunityTaxa` (`--pks_community_taxa`) both published
  `by_sample/<s>/figures/taxa_barplots.pdf`, so a run with both kept whichever finished last.
  They now publish `pks_island_taxa_barplots.pdf` and `community_taxa_barplots.pdf`, with
  report manifest kinds `figure_taxa_pks_island` / `figure_taxa_community`.
- Boolean flags were truthy when set to `false`, because Nextflow passes `--flag false` as a
  string. `--enable_mags false` would have enabled MAGs.
- `plotPKSTaxa` and `pksMAG` were each invoked twice, which Nextflow rejects at runtime. Both
  now use distinct aliases.
- `params.prefilter_dir` was undefined, so prefilter output would have been published to a
  directory named `null`.
- geNomad prophage calls on tumour contigs were computed and then discarded.
- `params.pks_reference_fasta` and `params.pks_recruit_index` interpolated before
  `params.pks_genome` was defined, yielding `null.fna`.

## v0.0.1 — 2026-09-01

First immutable release of the benchmarked workflow. See `RELEASE_NOTES_v0.0.1.md`, which
includes the disposition of a 16-point external review.
