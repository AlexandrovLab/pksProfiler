# HMM and combined profiling modes

The DNA HMM method can detect divergent or fragmented *clb* sequences. It is complementary to alignment, so its counts do not need to equal the Bowtie2 counts.

Both BAM and paired FASTQ sample sheets are supported. Follow the [BAM](bam.md) or [FASTQ](fastq.md) guide to prepare the appropriate CSV.

## HMM-only mode

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.csv \
    --input_data_type fastq \
    --profiling_method hmm \
    --hg38_db /references/pksProfiler/human-GRC-db.mmi \
    --t2t_phix_db /references/pksProfiler/human-GCA-phix-db.mmi \
    --outdir results_hmm
```

The combined output is:

```text
results_hmm/cohort/gene_counts/pks.gene.counts.hmm.txt
```

## Alignment plus HMM mode

Change the profiling method to `both`:

```bash
--profiling_method both
```

This produces both alignment and HMM count matrices and produces coverage plots from the alignment results.

## HMM chunking

For a large input, allow the allocated HMM CPUs to scan chunks in parallel:

```bash
--hmm_chunking true
```

Complete example:

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.csv \
    --input_data_type fastq \
    --profiling_method hmm \
    --hmm_chunking true \
    --hg38_db /references/pksProfiler/human-GRC-db.mmi \
    --t2t_phix_db /references/pksProfiler/human-GCA-phix-db.mmi \
    --outdir results_hmm_chunked
```

The default HMM E-value threshold is `1e-10`. See the parameter table in the [main README](../../README.md#adding-capabilities) before changing it.

[Return to the main README](../../README.md)

## Alignment versus profile models

`--profiling_method bowtie2 | hmm | both`

**Alignment** (the default) lines each read up against a single reference *E. coli* genome and
counts how many land on each gene. A read spanning two genes is credited to the one it covers
most, rather than being thrown away. Strict, and easy to explain.

**Profile models** (`hmm`) compare reads against a statistical description of what each *clb* gene
looks like across 2,868 different colibactin-producing genomes, instead of one reference. Better
at spotting versions that have drifted from the reference, and fragments too short to align
confidently.

They are not meant to give identical numbers. `both` runs them side by side, which is the honest
option when you care about divergent strains.

Models ship with the repository; see [`ref/hmm/PROVENANCE.md`](../../ref/hmm/PROVENANCE.md) for how they
were built and validated. They carry no pre-set score cutoffs, so detection depends on
`--hmm_evalue` and `--hmm_protein_evalue`.
