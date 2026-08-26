# pksProfiler

[![License: BSD 2-Clause](https://img.shields.io/badge/License-BSD_2--Clause-blue.svg)](LICENSE)

**pksProfiler** is a Nextflow pipeline for detecting, quantifying, and visualizing the *pks* (polyketide synthase) pathogenicity island in BAM or paired-end FASTQ data. It reports counts for each of the 19 *clbA–clbS* genes and can optionally identify species supported by *pks*-aligned reads and summarize their *clb* gene support.

## Quick start

```bash
git clone https://github.com/ammalabbasi/pksProfiler.git
cd pksProfiler

nextflow run main.nf \
    -profile conda \
    --sample examples/sample_sheets/fastq.csv \
    --input_data_type fastq \
    --profiling_method bowtie2 \
    --hg38_db /path/to/human-GRC-db.mmi \
    --t2t_phix_db /path/to/human-GCA-phix-db.mmi \
    --outdir results
```

The example sheet shows the required format, but its placeholder FASTQ paths must be replaced. See [Running pksProfiler](#running-pksprofiler) for complete BAM, FASTQ, HMM-chunked, and taxonomy-enabled examples.

## Installation

### Requirements

- Linux
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) 24.10 or newer
- Java 17 or newer, as required by Nextflow
- Conda or Mamba
- Git
- Storage for the input data, Nextflow work directory, and references

Pipeline tools are installed from the specifications in [`conda_envs/`](conda_envs). After installing Nextflow, clone and check the pipeline:

```bash
git clone https://github.com/ammalabbasi/pksProfiler.git
cd pksProfiler
nextflow -version
nextflow lint main.nf
```

### Reference inputs

The repository contains the *pks*-positive *E. coli* Bowtie2 index, the *clbA–clbS* annotation, and the DNA HMM database.

Two Minimap2 indexes are required for host/reference depletion:

| Parameter | Required index | Download |
|---|---|---|
| `--hg38_db` | GRCh38 (`human-GRC-db.mmi`) | [pksProfiler reference indexes (Google Drive)](https://drive.google.com/drive/folders/10np5NSeAPRpHz1a22drGybs4enTP-OdR?usp=share_link) |
| `--t2t_phix_db` | T2T-CHM13 plus phiX (`human-GCA-phix-db.mmi`) | [pksProfiler reference indexes (Google Drive)](https://drive.google.com/drive/folders/10np5NSeAPRpHz1a22drGybs4enTP-OdR?usp=share_link) |

Download both `.mmi` files from the shared folder, place them in a stable reference directory, and pass their full paths to the pipeline. For example:

```bash
--hg38_db /path/to/pksProfiler_reference_indexes/human-GRC-db.mmi \
--t2t_phix_db /path/to/pksProfiler_reference_indexes/human-GCA-phix-db.mmi
```

Files must be readable from every compute node.

Taxonomic profiling additionally requires a KrakenUniq-compatible database with Bracken files. The pipeline was validated with the 8 August 2023 Microbial database. It is very large: `database.kdb` is approximately 535 GB, in addition to the companion archive and extracted files.

Download both required files into one directory:

```bash
mkdir -p /path/to/krakenuniq_2023
cd /path/to/krakenuniq_2023

wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/database.kdb
wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/kuniq_microbialdb_minus_kdb.20230808.tgz
tar -xzf kuniq_microbialdb_minus_kdb.20230808.tgz
```

The directory should contain files such as `database.kdb`, `database.idx`, `taxDB`, `seqid2taxid.map`, and `database150mers.kmer_distrib`. Pass the **directory**, not one of its files:

```bash
--kraken_db /path/to/krakenuniq_2023 \
--bracken_read_length 150
```

Other KrakenUniq collections and their direct AWS links are listed in the [KrakenUniq section of the AWS index collection](https://benlangmead.github.io/aws-indexes/k2/#krakenuniq). Download both the `.kdb` and `.tar.gz` links for the same collection. Do **not** use a Kraken 2-only database.

## What the workflow does

<img src="workflow_logo/v1.png" width="800" alt="pksProfiler workflow">

1. **Prepare reads:** extract primary unmapped reads from BAM input or accept paired FASTQs, then filter reads with fastp.
2. **Deplete host/reference reads:** remove reads matching the supplied GRCh38 and T2T/phiX Minimap2 indexes.
3. **Profile the island:** quantify *clb* genes with Bowtie2, DNA HMMs, or both.
4. **Summarize results:** generate 19-gene count matrices and an island coverage plot for samples with aligned *pks* reads.
5. **Optionally assign taxonomy:** run KrakenUniq and Bracken only on reads aligned to the *pks* island, producing species-by-*clb*-gene support and taxonomic summaries.

Alignment and HMM profiling are complementary and are not expected to produce identical counts. Alignment provides stringent reference-based evidence; HMM profiling can recover more divergent or fragmented matches.

## Running pksProfiler

### Main parameters

| Parameter | Values/default | Required | Description |
|---|---|---:|---|
| `--sample` | CSV path | Yes | Sample sheet described below |
| `--input_data_type` | `bam` (default), `fastq` | Yes | Selects the sample-sheet format |
| `--profiling_method` | `bowtie2` (default), `hmm`, `both` | No | Profiling method(s) to run |
| `--hg38_db` | `.mmi` path | Yes | GRCh38 Minimap2 index |
| `--t2t_phix_db` | `.mmi` path | Yes | T2T/phiX Minimap2 index |
| `--outdir` | `results` | No | Output directory |
| `--hmm_evalue` | `1e-10` | No | Positive HMM E-value threshold |
| `--hmm_chunking` | `false` | No | Parallelize HMM scanning across chunks |
| `--pks_taxa` | off | No | Enable taxonomy; specify the flag without a value |
| `--kraken_db` | directory | With taxonomy | KrakenUniq/Bracken database directory |
| `--bracken_read_length` | positive integer | With taxonomy | Read length supported by the Bracken database |

Advanced reference and output parameters are defined near the top of [`main.nf`](main.nf).

### BAM input

The CSV must contain `patient` and `bam`. Each BAM may contain aligned and unmapped records; pksProfiler extracts primary records flagged as unmapped.

```csv
patient,bam
sample1,/data/sample1.bam
sample2,/data/sample2.bam
```

An editable example is available at [`examples/sample_sheets/bam.csv`](examples/sample_sheets/bam.csv).

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.bam.csv \
    --input_data_type bam \
    --profiling_method bowtie2 \
    --hg38_db /references/human-GRC-db.mmi \
    --t2t_phix_db /references/human-GCA-phix-db.mmi \
    --outdir results_bam
```

### Paired-end FASTQ input

The CSV must contain `patient`, `fastq1`, and `fastq2`.

```csv
patient,fastq1,fastq2
sample1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz
sample2,/data/sample2_R1.fastq.gz,/data/sample2_R2.fastq.gz
```

An editable example is available at [`examples/sample_sheets/fastq.csv`](examples/sample_sheets/fastq.csv).

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.fastq.csv \
    --input_data_type fastq \
    --profiling_method bowtie2 \
    --hg38_db /references/human-GRC-db.mmi \
    --t2t_phix_db /references/human-GCA-phix-db.mmi \
    --outdir results_fastq
```

Sample identifiers must be unique. Input paths should be absolute on a cluster.

### HMM chunking

For large inputs, enable chunked HMM execution:

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.fastq.csv \
    --input_data_type fastq \
    --profiling_method hmm \
    --hmm_chunking true \
    --hg38_db /references/human-GRC-db.mmi \
    --t2t_phix_db /references/human-GCA-phix-db.mmi \
    --outdir results_hmm
```

### Taxonomic profiling

Taxonomy is optional and depends on alignment profiling, so use `bowtie2` or `both`. KrakenUniq is run only on reads that overlap the *pks* island.

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.fastq.csv \
    --input_data_type fastq \
    --profiling_method both \
    --pks_taxa \
    --kraken_db /references/krakenuniq_database \
    --bracken_read_length 150 \
    --hg38_db /references/human-GRC-db.mmi \
    --t2t_phix_db /references/human-GCA-phix-db.mmi \
    --outdir results_taxonomy
```

Omit `--pks_taxa` to disable taxonomy; do not write `--pks_taxa false`.

### Resume an interrupted run

Repeat the same command with `-resume`, keeping the same work directory and parameters:

```bash
nextflow run main.nf -resume [the same pipeline options]
```

## Example results

Small synthetic examples illustrate the expected count-table format:

- [positive alignment result](examples/results/synthetic_positive/pks.gene.counts.align.txt)
- [positive HMM result](examples/results/synthetic_positive/pks.gene.counts.hmm.txt)
- [negative alignment result](examples/results/synthetic_negative/pks.gene.counts.align.txt)
- [negative HMM result](examples/results/synthetic_negative/pks.gene.counts.hmm.txt)

A positive sample contains one row per *clb* gene and one count column per sample. For example:

```text
Gene    synthetic_pks_positive
clbA    1
clbB    38
clbC    10
clbD    3
clbE    1
clbF    4
clbG    5
clbH    19
clbI    12
clbJ    20
clbK    20
clbL    5
clbM    5
clbN    16
clbO    9
clbP    5
clbQ    2
clbR    1
clbS    2
```

A successfully processed negative sample is retained with explicit zeros:

```text
Gene    synthetic_pks_negative
clbA    0
clbB    0
...     ...
clbS    0
```

These are output-format demonstrations, not a bundled end-to-end test dataset.

### Example coverage plot

The alignment workflow produces a circular view of coverage across the 19-gene island for a positive sample:

![Synthetic pks-positive coverage plot](examples/plots/synthetic_positive/synthetic_pks_positive.pks.circos.png)

[Download the example PDF](examples/plots/synthetic_positive/synthetic_pks_positive.pks.circos.pdf).

A sample with no aligned *pks* reads is retained in the count matrix with zero values, but no empty coverage PDF is generated.

## Output files

```text
results/
├── unmapped_reads/
├── host_depleted_reads/
├── pks_per_sample/
└── pks_summary/
    ├── gene_counts/
    │   ├── pks.gene.counts.align.txt
    │   └── pks.gene.counts.hmm.txt
    ├── coverage_plots/
    │   └── <sample>.pks.circos.pdf
    └── taxonomy/                    # only with --pks_taxa
        ├── pks.clb_species_support.tsv
        ├── bracken.genus.mpa.report.txt
        ├── bracken.species.mpa.report.txt
        └── plots/
            └── <sample>.pks_island_taxa_barplots.pdf
```

Combined count matrices always contain all 19 genes. Valid samples without qualifying signal receive zero counts. Coverage plots are produced only for samples with aligned *pks* reads.

The optional `pks.clb_species_support.tsv` contains species as rows and *clb* genes as columns. It is based on direct KrakenUniq classifications of *pks*-aligned reads. Bracken reports contain re-estimated genus/species abundances and may be empty when support is below the conservative threshold.

## Running on HPC systems

| Profile | Scheduler/use |
|---|---|
| `conda` | Workstation or one allocated compute node |
| `mamba` | Local execution with Mamba |
| `tscc` | UC San Diego TSCC (Slurm) |
| `slurm` | Generic Slurm cluster |
| `biowulf` | NIH Biowulf |
| `pbspro` | PBS Pro cluster |
| `lsf` | IBM Spectrum LSF cluster |
| `sge` | Sun/Oracle Grid Engine cluster |

For a generic cluster, keep site accounts, partitions, and QOS settings in a local file:

```groovy
// site.config
process {
    queue = 'my_partition'
    clusterOptions = '--account=my_account'
}
```

```bash
nextflow run main.nf \
    -profile slurm \
    -c site.config \
    [pipeline options]
```

See the [HPC execution guide](docs/hpc.md) for scheduler examples, Biowulf notes, driver-job guidance, and resource monitoring.

### Runtime and resource reports

```bash
RUN_TAG=$(date +%Y%m%d_%H%M%S)

nextflow run main.nf \
    [pipeline options] \
    -with-trace "run.${RUN_TAG}.trace.txt" \
    -with-report "run.${RUN_TAG}.report.html" \
    -with-timeline "run.${RUN_TAG}.timeline.html"
```

## License

This project is distributed under the [BSD 2-Clause License](LICENSE).

## Citation

A manuscript citation will be added when available.
