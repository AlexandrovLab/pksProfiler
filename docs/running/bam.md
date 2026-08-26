# BAM mode

Use this mode when each sample is stored in a BAM file. pksProfiler extracts every primary record flagged as unmapped, regardless of whether its mate is mapped.

## Sample sheet

Create a CSV with the columns `patient` and `bam`:

```csv
patient,bam
sample1,/absolute/path/to/sample1.bam
sample2,/absolute/path/to/sample2.bam
```

Sample identifiers must be unique, and paths should be absolute on a cluster. An editable template is available at [`examples/sample_sheets/bam.csv`](../../examples/sample_sheets/bam.csv).

## Run alignment profiling

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.bam.csv \
    --input_data_type bam \
    --hg38_db /references/pksProfiler/human-GRC-db.mmi \
    --t2t_phix_db /references/pksProfiler/human-GCA-phix-db.mmi \
    --outdir results_bam
```

`bowtie2` is the default profiling method. You may add `--profiling_method bowtie2` explicitly, but it is not required.

## Outputs to check

- `results_bam/pks_summary/gene_counts/pks.gene.counts.align.txt`
- `results_bam/pks_summary/coverage_plots/<sample>.pks.circos.pdf` for samples with aligned *pks* reads
- per-sample alignment, count, and coverage files under `results_bam/pks_per_sample/`

Valid samples without qualifying reads remain in the combined count matrix with zeros.

## Resume

Repeat the same command with `-resume`:

```bash
nextflow run main.nf -resume [the same options]
```

[Return to the main README](../../README.md)
