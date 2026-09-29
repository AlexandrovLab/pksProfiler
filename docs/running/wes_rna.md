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
