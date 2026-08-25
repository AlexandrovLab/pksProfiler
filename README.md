# pksProfiler

**pksProfiler** is a Nextflow pipeline for detecting, quantifying, and visualizing the *pks* (polyketide synthase) pathogenicity island in sequencing data. It reports gene-level signals across the 19-gene *clbA–clbS* cluster.

## Key features

- Host-read depletion and read filtering
- Alignment-based *clb* profiling
- HMM-based detection of divergent or fragmented *clb* sequences
- Per-sample and combined gene-count tables
- Coverage visualization across the *pks* island
- Optional taxonomic profiling of reads mapped to the *pks* island
- Optional parallel HMM execution through FASTA chunking
- Support for BAM and paired-end FASTQ input

## Workflow

The main workflow consists of:

1. Input validation and read preparation
2. Host/reference read depletion
3. Alignment-based and/or HMM-based *pks* profiling
4. Coverage plotting and count-table generation
5. Optional KrakenUniq/Bracken taxonomic profiling

## Workflow schematic

<img src="workflow_logo/v1.png" width="800" alt="pksProfiler workflow">

## Requirements

- Linux
- [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- Conda or Mamba
- Git
- Access to the required human-reference indexes
- A KrakenUniq/Bracken database only when taxonomic profiling is enabled

The individual bioinformatics tools are installed through the Conda environment specifications included in `conda_envs/`.

## Installation

Clone the repository:

```bash
git clone https://github.com/ammalabbasi/pksProfiler.git
cd pksProfiler
```

Confirm that Nextflow is available:

```bash
nextflow -version
```

Validate the pipeline source:

```bash
nextflow lint main.nf
```

## Reference inputs

The pipeline includes the *pks*-positive *Escherichia coli* reference indexes, the 19-gene annotation, and the DNA HMM database used for *clb* profiling.

Users must provide the host-depletion indexes:

- `--hg38_db`: path or prefix for the GRCh38 index
- `--t2t_phix_db`: path or prefix for the T2T/phiX index

The paths must be accessible from the compute nodes used to run the pipeline.

For optional taxonomic profiling, users must also provide:

- `--kraken_db`: KrakenUniq/Bracken database directory
- `--bracken_read_length`: positive integer matching a read length supported by the selected Bracken database

## Input sample sheets

The sample sheet is supplied with:

```text
--sample samples.csv
```

### BAM input

Use the columns `patient` and `bam`:

```csv
patient,bam
sample1,/path/to/sample1.bam
sample2,/path/to/sample2.bam
```

Run with:

```text
--input_data_type bam
```

### Paired-end FASTQ input

Use the columns `patient`, `fastq1`, and `fastq2`:

```csv
patient,fastq1,fastq2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

Run with:

```text
--input_data_type fastq
```

Sample identifiers should be unique.

## Running pksProfiler

### Basic run

```bash
nextflow run main.nf \
    -profile tscc \
    --sample samples.csv \
    --input_data_type fastq \
    --profiling_method both \
    --hg38_db /path/to/grch38/index \
    --t2t_phix_db /path/to/t2t_phix/index \
    --outdir results
```

For execution outside TSCC, omit `-profile tscc` or select another configured execution profile.

### Profiling methods

The `--profiling_method` parameter accepts:

- `bowtie2`: alignment-based profiling only
- `hmm`: HMM-based profiling only
- `both`: run both methods

Example:

```bash
--profiling_method hmm
```

### Optional HMM chunking

HMM chunking is disabled by default.

Enable it with:

```bash
--hmm_chunking true
```

When enabled and more than one CPU is allocated, the merged FASTA is split into chunks and processed by parallel single-threaded `nhmmscan` jobs.

Example:

```bash
nextflow run main.nf \
    -profile tscc \
    --sample samples.csv \
    --input_data_type fastq \
    --profiling_method hmm \
    --hmm_chunking true \
    --hg38_db /path/to/grch38/index \
    --t2t_phix_db /path/to/t2t_phix/index \
    --outdir results_hmm_chunked
```

Chunked and non-chunked execution should be validated on the same sample before large production runs.

### Optional taxonomic profiling

Taxonomic profiling is disabled by default.

Enable it with:

```bash
--pks_taxa true
```

A complete example is:

```bash
nextflow run main.nf \
    -profile tscc \
    --sample samples.csv \
    --input_data_type fastq \
    --profiling_method bowtie2 \
    --pks_taxa true \
    --kraken_db /path/to/kraken_bracken_database \
    --bracken_read_length 150 \
    --hg38_db /path/to/grch38/index \
    --t2t_phix_db /path/to/t2t_phix/index \
    --outdir results_with_taxonomy
```

The Bracken read length must be supported by the selected database.

### Resuming a run

Use Nextflow’s cache to resume an interrupted run:

```bash
nextflow run main.nf \
    -resume \
    -profile tscc \
    --sample samples.csv \
    --input_data_type fastq \
    --profiling_method both \
    --hg38_db /path/to/grch38/index \
    --t2t_phix_db /path/to/t2t_phix/index \
    --outdir results
```

Use the same parameters and work directory when resuming.

## Output structure

By default, outputs are written under `results/`. Use `--outdir` to select another location.

The main output directories are:

```text
results/
├── unmapped_reads/
├── host_depleted_reads/
├── pks_per_sample/
└── pks_summary/
```

Outputs include:

- Per-sample alignment count tables
- Per-sample HMM count tables
- Combined *clb* count matrices
- HMM hit read identifiers and filtered FASTA files
- Coverage tracks and *pks* island plots
- Optional KrakenUniq and Bracken reports
- Optional taxonomic summary plots

HMM count tables report all 19 genes, `clbA` through `clbS`. Genes without qualifying hits are reported with a count of zero.

## Validation

Check the pipeline syntax:

```bash
nextflow lint main.nf
```

Inspect resolved TSCC configuration:

```bash
nextflow config -profile tscc -flat
```

Before a large run, test the pipeline on a small sample and verify:

- expected output files are produced;
- count tables contain all 19 *clb* genes;
- valid zero-read samples produce zero-count outputs;
- corrupt inputs terminate with an error;
- `-resume` reuses completed tasks;
- chunked and non-chunked HMM results are equivalent.

## License

This project is distributed under the BSD 2-Clause License. See [LICENSE](LICENSE).

## Citation

A manuscript citation will be added when available.