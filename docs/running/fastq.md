# FASTQ mode

Use this mode when each sample has either one gzip-compressed FASTQ file or paired, gzip-compressed R1 and R2 FASTQ files. A sample sheet may contain single- and paired-FASTQ samples together.

## Sample sheet

For paired FASTQs, create a CSV with `patient`, `fastq1`, and `fastq2`:

```csv
patient,fastq1,fastq2
sample1,/absolute/path/to/sample1_R1.fastq.gz,/absolute/path/to/sample1_R2.fastq.gz
sample2,/absolute/path/to/sample2_R1.fastq.gz,/absolute/path/to/sample2_R2.fastq.gz
```

For one FASTQ per sample, `fastq2` may be omitted:

```csv
patient,fastq1
sample1,/absolute/path/to/sample1.unmapped.fastq.gz
sample2,/absolute/path/to/sample2.unmapped.fastq.gz
```

Alternatively, retain the `fastq2` column and leave its value empty for single-FASTQ samples.

Sample identifiers must be unique and may contain letters, numbers, periods, underscores, and hyphens; the first character must be alphanumeric. Paths should be absolute on a cluster. Editable templates are available for [paired FASTQs](../../examples/sample_sheets/fastq.csv) and a [single FASTQ](../../examples/sample_sheets/fastq_single.csv).

## Run alignment profiling

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.fastq.csv \
    --input_data_type fastq \
    --hg38_db /references/pksProfiler/human-GRC-db.mmi \
    --t2t_phix_db /references/pksProfiler/human-GCA-phix-db.mmi \
    --outdir results_fastq
```

`bowtie2` is the default profiling method.

## Outputs to check

- `results_fastq/pks_summary/gene_counts/pks.gene.counts.align.txt`
- `results_fastq/pks_summary/coverage_plots/<sample>.pks.circos.pdf` for samples with aligned *pks* reads
- `results_fastq/pks_summary/qc/pks.qc.summary.tsv` for per-sample read attrition and *pks* QC
- per-sample alignment, count, and coverage files under `results_fastq/pks_per_sample/`

Filtered and host-depleted FASTQs are not copied to the results directory by default. Add `--save_intermediates true` to publish them under `results_fastq/unmapped_reads/` and `results_fastq/host_depleted_reads/`.

## Resume

```bash
nextflow run main.nf -resume [the same options]
```

[Return to the main README](../../README.md)
