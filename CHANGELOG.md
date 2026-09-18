# Changelog

## v0.0.2dev — unreleased

Extends the immutable v0.0.1 release with island reassembly, genome-resolved analysis and
strain typing. Full detail in `RELEASE_NOTES_v0.0.2-dev.md`.

**Added**
- Sample-type routing: `--sample_type auto | metagenome | tumor_wgs | tumor_wes | tumor_rna`.
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

**Changed**
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
