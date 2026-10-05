# Two optional checks on the island reads: V8a and V8b

Both run on the reads the alignment lane already retained, never on the raw input, and
neither re-runs primary detection. Both are off by default. They answer **different
questions** and are deliberately kept apart — separate panels, separate scripts, separate
outputs, separate vocabularies.

|  | `--pks_locus_specificity` (V8a) | `--pks_carrier_attribution` (V8b) |
|---|---|---|
| Question | Does this read really belong to the pks/clb locus? | Which organism is the likely carrier? |
| Panel | **pks-negative**: Enterobacterales, NRPS/PKS homolog clusters, gut/oral backgrounds, GRCh38, UniVec | **pks-positive**: pks+ *E. coli* beyond IHE3034, pks+ *Klebsiella*, *Citrobacter*, *Enterobacter* |
| Outcome | `SURVIVES` / `NOT_DISTINGUISHABLE_FROM_CROSS_MAPPING` / `UNDERPOWERED` / `INCOMPLETE_MEASUREMENT` | `CARRIER_ATTRIBUTED` / `CARRIER_GENUS_ONLY` / `CARRIER_AMBIGUOUS` / `CARRIER_UNDETERMINED` |
| Published to | `by_sample/<sample>/locus_specificity/` | `by_sample/<sample>/carrier_attribution/` |
| Module | `Modules/pks_specificity.nf` | `Modules/pks_carrier.nf` |

**They are never one number.** A read that fits another carrier better than IHE3034 is
evidence about the organism, not against the locus: every V8b panel member carries the
island. A read that loses to a decoy is evidence about the locus, and says nothing about
any organism: no V8a panel member carries the island.

## Why V8a exists

`MAPQ >= 40` against IHE3034 shows confident mapping to *that* reference, not uniqueness
against the wider microbial universe. `clbB/clbC/clbH/clbI/clbJ/clbK/clbN/clbO` carry NRPS
and PKS megasynthase domains shared with unrelated bacterial secondary-metabolite clusters
— the same domain-homology failure already seen once at the assembled-bin level. V8a
re-challenges each supporting read against a panel built to contain exactly those
look-alikes, and reports what fraction survives.

Two independent arms may decide:

* **Arm B** (always): dual-index competitive realignment. One concatenated true+decoy
  index and two lone-family indices must independently agree before a read is
  `confirmed_true`; exact ties are `ambiguous`, never forced. See
  `scripts/competitive_realign_score.py`.
* **Arm A** (optional, `--pks_specificity_arm_a_dir`): a KrakenUniq reclassification
  against the same taxon set. Without it the lane runs Arm B alone and records
  `arms_used=arm_b_only`, so a one-arm verdict can never be read as a two-arm one.

**Measured versus lost.** Retention is survivors over the reads the challenge actually
*measured*. A supporting read with no Arm B row, or with no alignment to either side of
the panel in any architecture (or, with Arm A, no Arm A row), is **unmeasured**: it is
counted by kind (`n_unmeasured_missing_from_arm_b`, `n_unmeasured_unaligned`,
`n_unmeasured_missing_from_arm_a`), never named as an absorbing decoy, and
`retention_if_unmeasured_lost` gives the fail-closed bound beside the estimate. More than
`--pks_specificity_max_unmeasured_fraction` (default 0.05) unmeasured makes the verdict
`INCOMPLETE_MEASUREMENT`, because a pass or fail would describe a different read set.
Every V8 alignment runs minimap2 with `--frag=no`: under `-x sr` it otherwise pairs two
consecutive records named `X/1` and `X/2` and strips the suffix, which is how a third of
a coordinate-sorted island FASTQ once lost its names (V8 pilot, 5 Oct 2026).

## Why V8b exists, and how far it can go

The island is mobile and is carried at high identity by several genera. V8b asks which of
them best explains a sample's reads. It reports strain-level and genus-level attribution
**together**, because short reads over a highly conserved island usually cannot separate
strains; `CARRIER_GENUS_ONLY` is the honest answer in that case and the script says so
rather than picking one. An alignment-score margin — like NM — is a divergence metric
relative to the panel, never proof of carrier identity.

`CARRIER_UNDETERMINED` means fewer than `--pks_carrier_min_attributed_reads` reads
aligned to the panel at all. A sample with plenty of aligned reads that do not separate
the carriers — including every read tying across genera, so that nothing is attributed —
is `CARRIER_AMBIGUOUS`, and `tied_genus_sets` says which genera its reads could not tell
apart. A read aligned only to contigs the carrier map does not name means the index and
the map are from different builds, and the task fails rather than calling it unmapped.

## Building the panels

Both panels are built from public sequence by one script, and membership is screened from
sequence rather than trusted from the accession list: each candidate assembly is aligned
against the canonical clb island and must fall **below** the coverage ceiling (V8a) or
**at or above** the floor (V8b). Every candidate, accepted and rejected, lands in
`<panel>_pks_screen.tsv`.

```bash
# V8a: pks-negative decoys, plus UniVec, optionally GRCh38, plus non-colibactin
# NRPS/PKS clusters from a MIBiG release
python scripts/build_specificity_panels.py --panel v8a \
    --members ref/specificity_panels/v8a_locus_decoy_panel.tsv \
    --outdir <panel>/v8a --island-fasta island.fna --minimap2 $(command -v minimap2) \
    --univec --mibig-gbk https://dl.secondarymetabolites.org/mibig/mibig_gbk_4.0.tar.gz \
    --true-fasta indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic.fna

# V8b: pks-positive carriers
python scripts/build_specificity_panels.py --panel v8b \
    --members ref/specificity_panels/v8b_carrier_panel.tsv \
    --outdir <panel>/v8b --island-fasta island.fna --minimap2 $(command -v minimap2)
```

The island FASTA is the canonical interval cut from the pipeline's own reference:

```bash
samtools faidx indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic.fna \
    NC_017628.1:2193827-2244594 > island.fna
```

GRCh38 is opt-in (`--grch38 <fasta>`): it adds roughly 3 GB of sequence and a ~7 GB index.

## Running the lanes

```bash
nextflow run main.nf ... \
    --pks_locus_specificity true \
    --pks_specificity_concatenated_index <panel>/v8a/v8a_locus_concatenated.mmi \
    --pks_specificity_true_index         <panel>/v8a/v8a_true_index.mmi \
    --pks_specificity_decoy_index        <panel>/v8a/v8a_panel.mmi \
    --pks_specificity_true_contigs       <panel>/v8a/v8a_true_contigs.txt \
    --pks_specificity_decoy_contigs      <panel>/v8a/v8a_decoy_contigs.txt \
    --pks_specificity_island_length      50768

nextflow run main.nf ... \
    --pks_carrier_attribution true \
    --pks_carrier_panel_index <panel>/v8b/v8b_panel.mmi \
    --pks_carrier_map         <panel>/v8b/v8b_carrier_map.tsv
```

Both lanes need alignment profiling (`--profiling_method bowtie2` or `both`) and short
reads. V8a additionally needs `--sample_type wgs` or `metagenome`: it recomputes the
sample's tier from the surviving reads, and WES and RNA write a different evidence schema.

Each lane can also be run on its own against a cohort's retained outputs, without
re-running the pipeline:

```bash
nextflow run Modules/pks_specificity.nf --sample sheet.csv --sample_dir out/ ...
nextflow run Modules/pks_carrier.nf     --sample sheet.csv --sample_dir out/ ...
```

## Where the results appear

* `by_sample/<sample>/locus_specificity/locus_specificity_summary.tsv` and
  `…_per_read.tsv`, plus the Arm B score table and the alignments
* `by_sample/<sample>/carrier_attribution/carrier_attribution.tsv` and
  `carrier_per_read.tsv`
* `pks.master_summary.tsv` gains two column groups: `locus_specificity*` and `carrier_*`
* The cohort report gains one section per lane

## Thresholds

Every threshold is a parameter, none is hardcoded, and all of them are placeholders until
calibrated against real positive and negative controls. The defaults are in
`Modules/pks_specificity.nf` and `Modules/pks_carrier.nf`; the reasoning behind each is in
`scripts/specificity_decision.py` and `scripts/carrier_attribution.py`.
