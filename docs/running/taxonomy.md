# Taxonomy mode

Taxonomy is an optional addition to alignment profiling. KrakenUniq is run only on reads aligned to the *pks* island.

## Before running

1. Install the KrakenUniq/Bracken database described in the main [installation guide](../../README.md#4-optional-install-the-taxonomy-database).
2. Confirm the database directory contains `database.kdb`, `database.idx`, `taxDB`, and the Bracken distribution matching the chosen read length.
3. Prepare a [BAM](bam.md) or [paired FASTQ](fastq.md) sample sheet.

## Run taxonomy

```bash
nextflow run main.nf \
    -profile conda \
    --sample samples.csv \
    --input_data_type fastq \
    --profiling_method bowtie2 \
    --pks_taxa \
    --kraken_db /references/krakenuniq_2023 \
    --bracken_read_length 150 \
    --hg38_db /references/pksProfiler/human-GRC-db.mmi \
    --t2t_phix_db /references/pksProfiler/human-GCA-phix-db.mmi \
    --outdir results_taxonomy
```

Pass the KrakenUniq database **directory**, not `database.kdb`. The value of `--bracken_read_length` must match an available file such as `database150mers.kmer_distrib`.

To run HMM profiling in the same execution, use:

```bash
--profiling_method both
```

Taxonomy cannot be used with `--profiling_method hmm` alone because it depends on *pks*-aligned reads.

## Outputs to check

```text
results_taxonomy/pks_summary/taxonomy/
├── pks.clb_species_support.tsv
├── bracken.genus.mpa.report.txt
├── bracken.species.mpa.report.txt
└── plots/
    └── <sample>.pks_island_taxa_barplots.pdf
```

`pks.clb_species_support.tsv` contains species as rows and *clb* genes as columns. Bracken reports contain re-estimated genus and species abundances.

To disable taxonomy, omit `--pks_taxa`; do not write `--pks_taxa false`.

[Return to the main README](../../README.md)
