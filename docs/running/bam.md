# BAM and CRAM mode

Use this mode when each sample is stored in a BAM or CRAM file. pksProfiler extracts every primary record flagged as unmapped into one FASTQ stream, regardless of whether its mate is mapped. For a `bwa mem` alignment it also extracts primary mapped records with a low alignment score; see [Low-score mapped reads](#low-score-mapped-reads-bwa-mem--t-0) below.

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

## Low-score mapped reads (bwa mem -T 0)

`bwa mem` reports a read as unmapped (FLAG 4) only when its best local alignment scores below
`-T`, which is 30 by default. With `-T 0`, bwa instead reports the best short seed hit, however
poor, and marks the read mapped. GDC aligned the TCGA BAMs with `bwa mem -T 0` (it is in their
`@PG CL`); Mutographs used `-T 30`. So in a TCGA BAM about a third of genuine *pks* reads are
"mapped", mostly clipped, at MAPQ 0-9, and extracting FLAG 4 records alone never sees them.

Those are exactly the reads that `-T 30` would have left unmapped: primary mapped records with
`AS:i` < 30. Extraction adds them:

| Parameter | Default | Meaning |
|---|---|---|
| `--extract_low_score_mapped` | `auto` | `auto`: on when the input's `@PG` chain runs `bwa mem` (or `bwa-mem2 mem`), off otherwise. `on` / `off` force it |
| `--extract_min_as` | `30` | integer >= 1; primary mapped records with `AS:i` below it are extracted (bwa's default `-T`) |

- It is decided per sample from the header, before any record is decoded, and logged
  (`extractReads <sample>: low-score mapped extraction on_auto_bwa_mem ...` in the run log).
  `auto` stays off for STAR, minimap2 and other aligners, whose `AS` is on a different scale,
  and for a bwa `@PG` without a `bwa mem` command line (bwa aln writes no `AS` anyway).
- The input is still decoded once. Secondary and supplementary records are never extracted.
  A mapped primary record without an `AS` tag is never added; it is counted
  (`extract_low_score_no_as_records`).
- Read names, `/1` `/2` suffixes and orientation are the same as for FLAG 4 reads
  (`samtools fastq -N`; a reverse-strand record is reverse-complemented back).
- To detect `bwa mem` the header's `@PG CL` is read in memory: only the program name, the
  subcommand and the `-T` value are looked at, and no part of the command line is written.
- On a `-T 30` BAM it is on too, and adds nothing: bwa already left every such read unmapped.

Validation (projects/pksprofiler/validation/E1): on reads simulated through GDC's own aligner
settings the rule recovered 99.9-100% of the hidden *pks* reads; on a 20-sample TCGA pilot,
tumours gained *clb* evidence (7 of 10 moved up one tier) and the matched normals gained no
*clb* read.

QC columns (`cohort/qc/pks.qc.summary.tsv`):

| Column | Meaning |
|---|---|
| `extracted_unmapped_reads` | FLAG 4 reads extracted, as before |
| `extracted_low_score_mapped_reads` | primary mapped reads with `AS` < `--extract_min_as` (0 when off) |
| `extracted_reads` | everything extraction wrote: the two above, summed. This is what enters fastp (`filter_input_reads`) and the extraction step of the attrition chain |
| `extract_low_score_no_as_records` | mapped primary records with no `AS` tag, not added (only when on) |
| `extract_low_score_mode` | `on_auto_bwa_mem`, `off_auto_not_bwa_mem`, `on_requested` or `off_requested` |

The extracted FASTQ keeps its name, `<sample>.UNMAPPED.fastq.gz`, although with this on it also
holds the low-score mapped reads.

## Outputs to check

- `results_bam/cohort/gene_counts/pks.gene.counts.align.txt`
- `results_bam/by_sample/<sample>/figures/circos.pdf` for samples with aligned *pks* reads
- `results_bam/cohort/qc/pks.qc.summary.tsv` for per-sample read attrition and *pks* QC
- per-sample alignment, count, and coverage files under `results_bam/by_sample/<sample>/alignment/`

Valid samples without qualifying reads remain in the combined count matrix with zeros.
Extracted and host-depleted FASTQs are not copied to the results directory by default. Add `--save_intermediates true` if they need to be retained outside the Nextflow work directory.

For performance, each BAM/CRAM is decoded once. The QC summary counts the extracted reads (FLAG 4 and, for bwa mem input, low-score mapped) but does not perform an additional full-file scan to count every primary alignment record.

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
