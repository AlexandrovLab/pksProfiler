# Parameter reference

## Main options

| Parameter | Values/default | Required | Description |
|---|---|---:|---|
| `--sample` | CSV path | Yes | BAM or FASTQ sample sheet |
| `--input_data_type` | `bam` (default), `fastq` | No | Selects the sample-sheet format |
| `--profiling_method` | `bowtie2` (default), `hmm`, `both` | No | Profiling method(s) to run |
| `--hg38_db` | `.mmi` path | Yes | GRCh38 Minimap2 index |
| `--t2t_phix_db` | `.mmi` path | Yes | T2T/phiX Minimap2 index |
| `--outdir` | `results` | No | Output directory |
| `--hmm_evalue` | `1e-10` | No | Positive HMM E-value threshold |
| `--hmm_chunking` | `false` | No | Parallelize HMM scanning across chunks |
| `--pks_taxa` | off | No | Enables taxonomy when present |
| `--kraken_db` | directory | With taxonomy | KrakenUniq/Bracken database directory |
| `--bracken_read_length` | positive integer | With taxonomy | Read length supported by the Bracken database |

## Input-sheet columns

| Input type | Required columns |
|---|---|
| BAM | `patient,bam` |
| Paired FASTQ | `patient,fastq1,fastq2` |

Sample identifiers must be unique. Paths should be absolute when running on a cluster.

## Internal reference parameters

The default *pks* reference, annotation, cytoband, adapter, HMM, environment, and output paths are defined near the top of [`main.nf`](../../main.nf). Most users should not change them.

[Return to the main README](../../README.md)
