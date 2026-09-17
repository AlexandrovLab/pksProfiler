# Output reference

A run writes three trees and nothing else:

```text
results/
├── by_sample/<sample>/     everything this one sample produced
├── cohort/                 everything that spans samples
└── runs/                   Nextflow's own report, timeline and trace
```

Before v0.0.2 a sample's artefacts were split across six top-level trees with
inconsistent nesting, and finding one sample's results meant knowing which stage had
produced each file. Now the sample directory is the unit: if it is not in
`by_sample/<sample>/`, it is not about that sample.

Files inside `by_sample/<sample>/` carry no `<sample>.` prefix — the directory already
says which sample they are. Every published table still has a `sample` column, so a file
stays self-describing if you copy it somewhere else.

## by_sample/&lt;sample&gt;/

```text
by_sample/<sample>/
├── read_evidence.tsv       evidence tier, reads, genes, breadth      (step 2-3)
├── counts.txt              reads per clb gene
├── coverage.txt            island coverage
├── alignment/              sorted.bam, .bai, .sam, coverage.bedgraph, alignment.qc.tsv
├── hmm/                    nhmmscan tblout, hmm_counts.tsv, qc       --profiling_method hmm
├── contigs/                                                          (step 4, tumour)
│   ├── recruitment/        recruited read counts
│   ├── megahit/            contigs + canonical-reference PAF
│   ├── metaspades/         contigs + canonical-reference PAF
│   └── final_evidence/     final_pks_evidence.tsv
├── genomes/                                                          (step 5, --enable_mags)
│   ├── assembly/  bins/  checkm2/  gtdbtk/  annotation/
│   ├── pks_mag_summary.tsv      one row per draft genome
│   └── mag_status.tsv
├── community/                                                        (step 6)
│   ├── prophages/          geNomad virus summaries
│   ├── genomic_context/    per-genome neighbouring-gene context
│   ├── annotation/         Prokka + clb HMM on assembled contigs
│   ├── community_prophage_inventory.tsv
│   ├── pks_community_interactions.tsv
│   ├── pks_island_mobility.tsv
│   └── contig_pks_context.tsv / _summary.tsv / _mobility.tsv
├── strain/                 ST, clonal complex, phylogroup            (step 7)
├── taxonomy/               bracken reports, island reads             --pks_taxa
├── prefilter/              kept and rescued reads                    --prefilter_mode
├── figures/                circos.pdf, taxa_barplots.pdf, contig_validation.svg
└── intermediates/          extracted / filtered / host-depleted FASTQ  --save_intermediates
```

Both sample types write to the same paths. There is no `tumor_wgs/` or `metagenome/`
level: `--sample_type` decides which directories exist, not where they sit.

## cohort/

```text
cohort/
├── gene_counts/            pks.gene.counts.align.txt, pks.gene.counts.hmm.txt
├── qc/                     pks.qc.summary.tsv — read attrition through every stage
├── taxonomy/               combined bracken tables, clb taxonomy support
├── tumor_eligibility/      per-sample tumour screen
└── pks.master_summary.tsv  written by build_master_summary.py (see below)
```

Start with `cohort/qc/pks.qc.summary.tsv` when a result looks wrong — it shows how many
reads survived extraction, fastp and each host-depletion pass, so you can see where a
sample lost its reads.

## runs/

`report.html`, `timeline.html` and `trace.txt`, written by Nextflow itself and
configured in `nextflow.config`. They overwrite on re-run: a resumed run that reused a
report name previously died with `FileAlreadyExistsException` after every task had
already succeeded.

`trace.txt` carries per-process CPU time and peak RSS, which is where runtime and memory
figures should come from.

## One table with everything

`scripts/build_master_summary.py` joins every stage into one row per sample:

```bash
python3 scripts/build_master_summary.py --results results \
  --output results/cohort/pks.master_summary.tsv
```

Column groups appear only when the stage that produces them ran:

| Group | Columns |
|---|---|
| read profiling | `input_reads`, `reads_after_*`, `num_clb_genes_*`, `reads_clb_genes_*` |
| evidence tier | `read_evidence`, `pks_reads`, `clb_genes_detected`, `island_breadth_1x` |
| contig reassembly | `final_structural_evidence`, `assembler_agreement` |
| draft genomes | `mag_bins_total`, `mag_bins_pks_positive`, `pks_mag_taxonomy`, `pks_mag_completeness`, `pks_mag_clb_genes` |
| community context | `community_prophages_total`, `community_pks_producers`, `community_neighbours_assessed`, `community_neighbours_with_prophage`, `island_nearby_integrase`, `island_nearby_trna`, `island_in_prophage` |
| strain typing | `typing_status`, `ST`, `clonal_complex`, `phylogroup` |

Where a sample has several draft genomes, the reported one is the *pks*-positive genome
carrying the most *clb* genes. A sample that was profiled and found negative keeps its
row, with `NA` in the later groups.

It only reads published output, so it is safe to re-run at any time.

## Primary result files

Column-level specifications are still being written. Until then, the pages under
`docs/running/` describe the files each mode produces, and every table carries a header
row.

| File | Produced by | One row per |
|---|---|---|
| `cohort/gene_counts/pks.gene.counts.*.txt` | every run | *clb* gene |
| `cohort/qc/pks.qc.summary.tsv` | every run | sample |
| `cohort/pks.master_summary.tsv` | `build_master_summary.py` | sample — every stage joined |
| `by_sample/*/read_evidence.tsv` | `--sample_type tumor_wgs` or `metagenome` | sample |
| `by_sample/*/contigs/final_evidence/final_pks_evidence.tsv` | targeted reassembly | sample |
| `by_sample/*/genomes/pks_mag_summary.tsv` | `--enable_mags` | draft genome |
| `by_sample/*/community/community_prophage_inventory.tsv` | `--enable_mags` | prophage |
| `by_sample/*/community/pks_community_interactions.tsv` | `--enable_mags` | producer × neighbour pair |
| `by_sample/*/community/pks_island_mobility.tsv` | `--enable_mags` | draft genome |
| `by_sample/*/community/contig_pks_context.tsv` | `--tumor_full_contig_context` | *clb* gene on a contig |
| `by_sample/*/strain/strain_types.tsv` | assembled output | typed unit |

## Reading results from before v0.0.2

Runs published under the old layout (`pks_per_sample/`, `pks_summary/`, `mags/`,
`community_context/`, `strain_typing/`, `prefilter/`) are not migrated. Keep them as
they are; `build_master_summary.py` reads the new layout only.

---

Back to the [README](../README.md).
