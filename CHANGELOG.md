# Changelog

## v0.0.2dev — unreleased

Extends the immutable v0.0.1 release with island reassembly, genome-resolved analysis and
strain typing. Full detail in `RELEASE_NOTES_v0.0.2-dev.md`.

**Added**
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
  profiling, `wes_evidence.tsv`, no tier, no targeted assembly, bowtie2/both only. Duplicate-rate
  and capture/off-target QC are a documented TODO (`docs/running/wes_rna.md`).
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

**Changed**
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
