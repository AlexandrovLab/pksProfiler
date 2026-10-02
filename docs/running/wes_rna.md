# Exome (WES) and RNA-seq input

Ludmil's 26 Sep mode memo: WGS is the calibrated tumour mode. WES and RNA-seq get lanes of
their own, with their own output schema, and **neither is ever given a WGS evidence tier or
targeted assembly**. Their results are supplementary evidence of a different kind, not a
second estimate of the WGS presence call, and should not be pooled with WGS prevalence.

| | `--sample_type wes` | `--sample_type rna` |
|---|---|---|
| Host depletion | standard (`minimap2 -ax sr`) | splice-aware human pass (see below) |
| Profiling | Bowtie2 + featureCounts, as WGS | Bowtie2 + featureCounts, as WGS |
| Output | `by_sample/<sample>/wes_evidence.tsv` | `by_sample/<sample>/rna_evidence.tsv` |
| Tier | none | none |
| Targeted assembly, MAGs | no | no |
| Read as | exploratory off-target DNA evidence | expressed *clb* sequence |

Both require `--profiling_method bowtie2` or `both` (the evidence file comes from the
alignment lane). `--tumor_targeted_assembly` is on by default for WGS; for these modes it is
inert and the run logs that it will not assemble. `--prefilter_mode auto` resolves to `off`.

## WES

```bash
nextflow run main.nf -profile local \
  --sample samples.csv --input_data_type bam \
  --sample_type wes \
  --wes_capture_bed /refs/capture_kit_targets.bed \
  --hg38_db      /refs/human-GRC-db.mmi \
  --t2t_phix_db  /refs/human-GCA-phix-db.mmi \
  --outdir       results
```

```text
sample   mode  clb_wes_reads  wes_clb_genes_detected  wes_island_breadth_1x  wes_island_breadth_2x  wes_island_breadth_3x
```

Exome capture enriches human exons; any *clb* read is off-target by construction, so read
counts depend on capture kit and off-target rate as much as on the tissue. A WES positive can
be meaningful; **a WES negative is not biological absence.** If a WES positive/negative call
is ever needed it must be calibrated on WES independently -- the WGS thresholds do not apply.

**WES QC (memo section 3).** `cohort/qc/pks.qc.summary.tsv` reports total,
extracted-unmapped, post-fastp and host-depleted read counts for WES as for WGS, and for a
`wes` run (only) six more columns, just before `status`:

| Column | Definition |
|---|---|
| `wes_input_primary_mapped` | primary mapped records in the input alignment (`samtools view -F 2308`) |
| `wes_input_duplicate_fraction` | of those, the fraction flagged `0x400` -- the **submitter's** duplicate marking; pksProfiler does not re-mark the alignment. An input nobody ran a marker on reads `0.000000`; check the `@PG` chain in `provenance/input_header.tsv` |
| `wes_capture_on_target_fraction` | of those, the fraction overlapping `--wes_capture_bed` (`samtools view -c -L`) |
| `wes_capture_off_target_fraction` | 1 − on-target; `capture_bed_not_provided` when no BED was given |
| `wes_pks_island_reads` | primary island reads in the profiling BAM (MAPQ ≥ 40, `pks_contig:pks_start_1based-pks_end_1based`) |
| `wes_pks_island_duplicate_fraction` | of those, the fraction sharing (contig, strand, unclipped 5′ position) with an earlier read |

The two duplicate rates are different quantities and are never combined. The input's
`0x400` flags describe the library; they do not survive extraction to FASTQ, so the pks
reads are collapsed by position on the single-end profiling alignment instead (mates are
separate records, so a duplicated pair counts twice). Position collapse cannot tell a PCR
duplicate from two fragments that start at the same base, so at high island depth it
overestimates duplication; off-target WES island depth is low, where the two coincide.

The input counts come from one extra decode of the input's primary mapped records
(`wesInputQC`, `Modules/mode_qc.nf`), teed into `samtools flagstat` and the BED count; it is
a separate task from extraction, which reads the unmapped records only. **Alignment input
only**: a FASTQ WES sheet has no alignment to count, and its three input columns read `NA`
(the pks-island columns still come from the profiling BAM). The BED's content digest is
recorded in the provenance record (`runs/provenance/<session>.json`, `dependencies.wes_capture_bed`),
and the path in `RUN_REPORT.txt` under DATABASES. Use the kit's target (not bait) BED, with
the input alignment's contig names.

The master summary carries the same six columns for a `wes` run.

## RNA-seq

```bash
nextflow run main.nf -profile local \
  --sample samples.csv --input_data_type bam \
  --sample_type rna \
  --rna_splice_db     /refs/grch38_splice_sr/GRCh38.d1.vd1.splice_sr.mmi \
  --rna_junctions_bed /refs/grch38_splice_sr/gencode.v36.junctions.bed \
  --rna_library_type  polyA \
  --t2t_phix_db       /refs/human-GCA-phix-db.mmi \
  --outdir            results
```

**Splice-aware host depletion.** `minimap2 -ax sr` is not splice-aware: a human read that
spans an exon-exon junction aligns in two pieces or not at all, and survives depletion as if
it were non-human. For `rna` the first (human) pass is therefore

```text
minimap2 -ax splice:sr [--frag=no] [-j <rna_junctions_bed>] --secondary=no <rna_splice_db>
```

and `--hg38_db` is not required (it is not used). The T2T+phiX and optional pangenome passes
stay `-ax sr`: there is no junction annotation for T2T-CHM13, phiX or the pangenome
haplotypes, and those passes exist to catch residual human and spike-in reads GRCh38 missed.
`--minimap2_frag_mode` applies to every call. The QC metric names are unchanged
(`reads_after_hg38` is the reads left after the spliced human pass). The table of presets per
sample type is `lib/HostDepletion.groovy`; `tests/check_depletion_plan.sh` executes it.

| Parameter | Required | |
|---|---|---|
| `--rna_splice_db` | yes | a minimap2 index built from GRCh38 **with `-x splice:sr`** (an index carries its preset's k-mer settings), or a FASTA |
| `--rna_junctions_bed` | no | BED12 of annotated junctions, passed as `-j`; improves alignment across short exons |
| `--rna_library_type` | no | `polyA`, `total` or `unknown` (default); recorded in `rna_evidence.tsv` and `RUN_REPORT.txt` |

Getting the references (public GRCh38.d1.vd1 + GENCODE v36; no sample data involved, ~12 GB).
Download the prebuilt files and verify their MD5s:

```bash
bash scripts/download_rna_references.sh /path/to/dbs/rna_grch38_splice_sr_gencode_v36
```

They are shared at <https://drive.google.com/drive/folders/10np5NSeAPRpHz1a22drGybs4enTP-OdR>
(index built with minimap2 2.31). Or build them yourself:

```bash
minimap2 -x splice:sr -d GRCh38.d1.vd1.splice_sr.mmi GRCh38.d1.vd1.fa
paftools.js gff2bed gencode.v36.annotation.gtf > gencode.v36.junctions.bed
```

GRCh38.d1.vd1 is the GDC reference that TCGA's harmonised data use, and GENCODE v36 matches GDC's RNA
annotation. On TSCC the built files are in
`/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/refs_public/rna_grch38_splice_sr_gencode_v36/`.

**Output.**

```text
sample   mode  clb_rna_reads  rna_clb_genes_detected  rna_expression_breadth_1x  rna_expression_breadth_2x  rna_expression_breadth_3x  rna_library_type  clb_transcripts_detected  rna_residual_human_fraction  rna_residual_checked_reads
```

**RNA QC (memo section 4).**

`clb_transcripts_detected` is the number of *clb* genes (clbA-clbS) with at least
`--rna_transcript_min_reads` (default 3) featureCounts-assigned RNA reads. It is
gene-level: the IHE3034 annotation (`ref/annotations/IHE3034.clbA-clbS.gff`) defines genes,
not operons or transcription units, so each gene is its own unit. `rna_clb_genes_detected`
is the same count at ≥ 1 read; the higher floor keeps a single stray read from counting as
a transcript.

The **residual-human** check measures human reads that survived every depletion pass
(`rnaResidualHuman`, `Modules/mode_qc.nf`; `scripts/rna_residual_human.py`):

1. a systematic sample of the final host-depleted reads -- every *k*-th read, *k* chosen so
   at most `--rna_residual_human_sample_reads` (default 200,000) are taken;
2. realigned to `--rna_splice_db` with the human pass's own preset and junctions
   (`minimap2 -ax splice:sr [-j <rna_junctions_bed>]`) **plus `--rna_residual_human_args`**
   (default `-n 2 -m 20 -s 30`: fewer minimizers, a lower chaining score and a lower
   alignment score than the preset's floors);
3. `rna_residual_human_fraction` = reads with a primary alignment / reads sampled.

Re-running the depletion command on its own output would find nothing by construction; a
human read survives depletion by scoring under the preset's floors (divergent, short-exon,
junction-heavy), and the relaxed check is the search for exactly those. On simulated 100 bp
reads against a 200 kb random reference (minimap2 2.31, `-x splice:sr` index), the default
preset placed 403 / 209 / 68 of 500 reads at 8 / 12 / 16 % divergence and the relaxed floors
441 / 303 / 150; neither placed any of 5,000 random 100 bp reads on a 20 Mb random reference
(`-n 1 -m 15 -s 25` placed one, which is why the default stops short of it). The check only
measures: the reads that go on to profiling are unchanged, and nothing is removed.
Long reads (`--read_type ont|hifi`) use the long-read splice preset without the relaxation
(the floors are on a different scale there), so for them it is a consistency check.

`pks.qc.summary.tsv` carries `rna_residual_checked_reads`, `rna_residual_human_reads` and
`rna_residual_human_fraction` for an `rna` run, and checks that the sample was drawn from
exactly the final host-depleted count. The master summary and the cohort report (row
tooltip, and `--table` columns when a run has RNA rows) carry `clb_transcripts_detected`
and `rna_residual_human_fraction`.

Interpret RNA only as **expressed *clb* sequence**. Counts reflect expression and library
capture (poly(A) selection under-represents bacterial transcripts, which are not
polyadenylated), gene breadth is expression breadth rather than island completeness, and an
RNA negative does not establish genomic absence. Profiling is unchanged (Bowtie2 +
featureCounts against the IHE3034 genome), with no tier, no targeted assembly and no MAGs.
