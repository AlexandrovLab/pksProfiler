# Taxonomy mode

Taxonomy is an optional addition to alignment profiling. KrakenUniq is run only on reads aligned to the *pks* island.

## What it measures

Everything this mode writes is **read-level taxonomic support** (Ludmil, 26 Sep mode memo
section 5). KrakenUniq assigns each read to the taxon its k-mers most resemble, and Bracken
re-estimates abundance from those assignments. A species in these tables is a taxon some
*clb*-aligned reads were assigned to — not the organism shown to carry the island. Short island
reads are shared across hosts of a mobile element, so the same reads can support several species.
Carrier-genome language is reserved for contig or genome-bin evidence that links *pks* to a
reconstructed genome ([MAGs](mags.md), [assembly](assembly.md)).

## Before running

1. Install the KrakenUniq/Bracken database described in the main [installation guide](../../README.md#databases).
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
results_taxonomy/cohort/taxonomy/
├── pks.clb_species_support.tsv
├── bracken.genus.mpa.report.txt
└── bracken.species.mpa.report.txt

results_taxonomy/by_sample/<sample>/taxonomy/
├── bracken.{G,S}.report.txt                --pks_taxa (island-aligned reads)
├── bracken.{G,S}.krakenreport.txt, .mpa.krakenreport.txt
├── community_bracken.{G,S}.report.txt      --pks_community_taxa (whole non-host report)
└── community_bracken.{G,S}.krakenreport.txt

results_taxonomy/by_sample/<sample>/figures/
├── pks_island_taxa_barplots.pdf    --pks_taxa (Bracken on island-aligned reads)
└── community_taxa_barplots.pdf     --pks_community_taxa (Bracken on the whole non-host report)
```

`pks.clb_species_support.tsv` contains species as rows and *clb* genes as columns: how many reads aligned to each gene were assigned to that species (read-level taxonomic support). Bracken reports contain re-estimated genus and species abundances, also read-level.

To disable taxonomy, omit `--pks_taxa`; do not write `--pks_taxa false`.

[Return to the main README](../../README.md)
