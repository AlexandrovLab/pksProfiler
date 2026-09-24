# Island reassembly (tumour data)

Counting reads tells you how much *clb* sequence is present. Reassembling those reads tells you
whether they form a real, connected island or scattered fragments that merely resemble one.

This page covers `--sample_type tumor_wgs`: evidence tiering and targeted reassembly. Prophage
and neighbouring-gene mobility context is a metagenome-only capability -- see
[MAGs](mags.md) and [community context](community_context.md).

## Evidence tiering

This runs automatically for tumour data. Set `--sample_type` and you also get a tier per sample.

```bash
nextflow run main.nf -profile local \
  --sample samples.csv \
  --input_data_type bam \
  --sample_type tumor_wgs \
  --hg38_db      /refs/human-GRC-db.mmi \
  --t2t_phix_db  /refs/human-GCA-phix-db.mmi \
  --outdir       results
```

**What you get**

```text
results/by_sample/<sample>/read_evidence.tsv
```

One row, with the tier plus the three numbers behind it:

```text
sample      read_evidence      pks_reads  clb_genes_detected  island_breadth_1x
TUMOR_01    extensive_island   1497       15                  0.731800
```

Because you set `--sample_type tumor_wgs`, samples that clear the gate also go straight into
step 3 — targeted assembly is on by default.

## Targeted reassembly

Counting reads tells you *how much*. Assembly tells you whether those reads form a real, connected
island or are scattered fragments that merely resemble one.

pksProfiler collects the island reads, reassembles them **twice using two independent methods**,
lines both results up against the reference island, and tells you whether the two agree. Two
attempts rather than one because agreement is itself evidence — a real island reconstructs the
same way twice.

This is already running if you did step 2. To be explicit:

```bash
nextflow run main.nf -profile local \
  --sample samples.csv \
  --input_data_type bam \
  --sample_type tumor_wgs \
  --tumor_targeted_assembly true \
  --hg38_db      /refs/human-GRC-db.mmi \
  --t2t_phix_db  /refs/human-GCA-phix-db.mmi \
  --outdir       results
```

**What you get**, per qualifying sample:

```text
results/by_sample/<sample>/contigs/
├── evidence/          the tier and the numbers behind it
├── recruitment/       the reads pulled out for assembly
├── megahit/           contigs from the first assembly method
├── metaspades/        contigs from the second
└── final_evidence/
    ├── <sample>.final_pks_evidence.tsv
    └── <sample>.raw_depth_megahit_metaspades_IHE3034.svg
```

The SVG is the figure to look at: reference genes along the top, your assembled contigs beneath,
read depth underneath that. `final_pks_evidence.tsv` ends in one of:

| Call | Meaning |
|---|---|
| `complete_island` | essentially the whole island reconstructed |
| `near_complete_island` | most of it, small gaps |
| `partial_island` | a real but incomplete stretch |
| `fragmented_island` | pieces that do not join up |
| `limited_contig_support` | assembly barely produced anything |
| `read_evidence_only` | reads were there, contigs were not |

**Only `broad_island` and `extensive_island` samples are assembled.** Below that there are too few
reads for assembly to succeed even when an island is genuinely present, so a failure would be
telling you about read depth, not biology. To widen it anyway:

```bash
  --tumor_contig_tiers multi_gene,broad_island,extensive_island
```

---

Back to the [README](../../README.md).
