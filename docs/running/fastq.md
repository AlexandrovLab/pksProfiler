# Paired FASTQ mode

Use this mode when each sample has paired, gzip-compressed R1 and R2 FASTQ files.

## Sample sheet

Create a CSV with the columns `patient`, `fastq1`, and `fastq2`:

```csv
patient,fastq1,fastq2
sample1,/absolute/path/to/sample1_R1.fastq.gz,/absolute/path/to/sample1_R2.fastq.gz
sample2,/absolute/path/to/sample2_R1.fastq.gz,/absolute/path/to/sample2_R2.fastq.gz
```

Sample identifiers must be unique, and paths should be absolute on a cluster. An editable template is available at [`examples/sample_sheets/fastq.csv`](../../examples/sample_sheets/fastq.csv).

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
- filtered and host-depleted reads under `results_fastq/host_depleted_reads/`
- per-sample alignment, count, and coverage files under `results_fastq/pks_per_sample/`

## Resume

```bash
nextflow run main.nf -resume [the same options]
```

[Return to the main README](../../README.md)
