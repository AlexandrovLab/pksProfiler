# BAM and CRAM mode

Use this mode when each sample is stored in a BAM or CRAM file. pksProfiler extracts every primary record flagged as unmapped into one FASTQ stream, regardless of whether its mate is mapped.

## Sample sheet

Create a CSV with the columns `patient` and `alignment`:

```csv
patient,alignment
sample1,/absolute/path/to/sample1.bam
sample2,/absolute/path/to/sample2.cram
```

Sample identifiers must be unique and may contain letters, numbers, periods, underscores, and hyphens; the first character must be alphanumeric. Paths should be absolute on a cluster. An editable template is available at [`examples/sample_sheets/bam.csv`](../../examples/sample_sheets/bam.csv).

## Run alignment profiling

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.bam.csv \
    --input_data_type auto \
    --hg38_db /references/pksProfiler/human-GRC-db.mmi \
    --t2t_phix_db /references/pksProfiler/human-GCA-phix-db.mmi \
    --outdir results_bam
```

`bowtie2` is the default profiling method. You may add `--profiling_method bowtie2` explicitly, but it is not required.

CRAM references may be embedded or resolved by HTSlib. If that lookup is unavailable, add `--cram_reference /absolute/path/to/reference.fa`. When supplied, the pipeline validates the CRAM sequence names, lengths, and MD5 values against the FASTA before extraction. Use `--input_data_type cram` to reject any accidental BAM input. Existing `patient,bam` sample sheets remain supported.

## Outputs to check

- `results_bam/pks_summary/gene_counts/pks.gene.counts.align.txt`
- `results_bam/pks_summary/coverage_plots/<sample>.pks.circos.pdf` for samples with aligned *pks* reads
- `results_bam/pks_summary/qc/pks.qc.summary.tsv` for per-sample read attrition and *pks* QC
- per-sample alignment, count, and coverage files under `results_bam/pks_per_sample/`

Valid samples without qualifying reads remain in the combined count matrix with zeros.
Extracted and host-depleted FASTQs are not copied to the results directory by default. Add `--save_intermediates true` if they need to be retained outside the Nextflow work directory.

For performance, each BAM/CRAM is decoded once. The QC summary counts extracted primary-unmapped reads but does not perform an additional full-file scan to count every primary alignment record.

## Resume

Repeat the same command with `-resume`:

```bash
nextflow run main.nf -resume [the same options]
```

[Return to the main README](../../README.md)
