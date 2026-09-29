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

**TODO (memo section 3, WES QC).** The QC table reports total, extracted-unmapped,
post-fastp and host-depleted read counts for WES exactly as for WGS. Not yet reported:
the **duplicate rate** (WES duplicate behaviour must be checked separately from WGS) and the
**capture kit / off-target fraction**. Both need inputs the pipeline does not currently
take (a capture BED, and duplicate marking on the full alignment rather than on the unmapped
reads it extracts), so they are left open rather than approximated.

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

Building the references (public GRCh38 + GENCODE; no sample data involved):

```bash
minimap2 -x splice:sr -d GRCh38.d1.vd1.splice_sr.mmi GRCh38.d1.vd1.fa
paftools.js gff2bed gencode.v36.annotation.gtf > gencode.v36.junctions.bed
```

On TSCC these are being built under
`/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/ref_indices/grch38_splice_sr/`.

**Output.**

```text
sample   mode  clb_rna_reads  rna_clb_genes_detected  rna_expression_breadth_1x  rna_expression_breadth_2x  rna_expression_breadth_3x  rna_library_type
```

Interpret RNA only as **expressed *clb* sequence**. Counts reflect expression and library
capture (poly(A) selection under-represents bacterial transcripts, which are not
polyadenylated), gene breadth is expression breadth rather than island completeness, and an
RNA negative does not establish genomic absence. Profiling is unchanged (Bowtie2 +
featureCounts against the IHE3034 genome), with no tier, no targeted assembly and no MAGs.

Not yet reported for RNA (memo section 4 QC): residual-human reads after depletion as a
separate metric, and `clb_transcripts_detected`, which needs a transcript model of the island.
