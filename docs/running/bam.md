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

**The input must be an alignment against a reference.** Extraction keeps the records an
aligner flagged unmapped, and checks each file with `samtools quickcheck`, which rejects a
BAM/CRAM whose header has no `@SQ` lines. An *unaligned* BAM (from `samtools import`, or a
basecaller's uBAM, common for ONT/PacBio) therefore fails in `extractReads` with
`has no @SQ lines in its header: it is an unaligned bam`. Convert it to FASTQ and use FASTQ
mode instead:

```bash
samtools fastq -@ 4 sample.unaligned.bam | gzip > sample.fastq.gz   # then --input_data_type fastq
```

(For paired uBAMs, `samtools fastq -1 R1.fastq.gz -2 R2.fastq.gz -0 /dev/null -s /dev/null`.)

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

CRAM input requires `--cram_reference /absolute/path/to/reference.fa`; the pipeline validates the CRAM sequence names, lengths, and MD5 values against the FASTA before extraction. Use `--input_data_type cram` to reject any accidental BAM input. Existing `patient,bam` sample sheets remain supported.

## Outputs to check

- `results_bam/cohort/gene_counts/pks.gene.counts.align.txt`
- `results_bam/by_sample/<sample>/figures/circos.pdf` for samples with aligned *pks* reads
- `results_bam/cohort/qc/pks.qc.summary.tsv` for per-sample read attrition and *pks* QC
- per-sample alignment, count, and coverage files under `results_bam/by_sample/<sample>/alignment/`

Valid samples without qualifying reads remain in the combined count matrix with zeros.
Extracted and host-depleted FASTQs are not copied to the results directory by default. Add `--save_intermediates true` if they need to be retained outside the Nextflow work directory.

For performance, each BAM/CRAM is decoded once. The QC summary counts extracted primary-unmapped reads but does not perform an additional full-file scan to count every primary alignment record.

## Resume

Repeat the same command with `-resume`:

```bash
nextflow run main.nf -resume [the same options]
```

[Return to the main README](../../README.md)

## Reading the gene count table

```text
results/cohort/gene_counts/pks.gene.counts.align.txt   one row per gene, one column per sample
results/by_sample/<sample>/figures/circos.pdf  circular map of coverage over the island
results/cohort/qc/pks.qc.summary.tsv                   how many reads survived each stage
```

**How to read it.** A sample with zero counts stays in the table as a row of zeros — it does not
disappear, so you can tell "tested and negative" from "never ran". In the QC summary,
`reads_clb_genes_align` is the read count that matters and `num_clb_genes_align` is how many of the
19 genes were seen at all.

> One gene with a couple of reads is not a positive result. Short stretches of *clb* genes resemble
> sequence found elsewhere, so a small number of reads can appear in a sample with no island.
> That is the entire reason for step 2.
