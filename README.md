# pksProfiler

[![CI](https://github.com/AlexandrovLab/pksProfiler/actions/workflows/ci.yml/badge.svg)](https://github.com/AlexandrovLab/pksProfiler/actions/workflows/ci.yml)
[![Version](https://img.shields.io/badge/version-0.0.2dev-orange.svg)](CHANGELOG.md)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.10-23aa62.svg)](https://www.nextflow.io/)
[![License: BSD 2-Clause](https://img.shields.io/badge/License-BSD_2--Clause-blue.svg)](LICENSE)

**pksProfiler** searches sequencing data for the *pks* island — a cluster of 19 bacterial genes
(*clbA* through *clbS*) that together produce **colibactin**, a molecule that damages human DNA and
leaves a recognisable mutation pattern in colorectal tumours.

You give it BAM, CRAM or FASTQ files. It tells you whether the island is there, **how convincing
the evidence is**, and — depending on your data — which bacterium is carrying it, what strain that
bacterium is, and whether the island looks able to move between bacteria.

<img src="workflow_logo/v2-animated.svg" width="1000" alt="pksProfiler workflow: a shared trunk of read extraction, quality filtering and human-read removal, splitting by sample type into a tumour lane that scores evidence and reassembles the island, and a metagenome lane that assembles and bins the community, with strain typing on assembled output from both">

*Left: metagenomes, where the question is which organism carries the island. Right: tumour
genomes, where the question is how strong the evidence is and whether the island reassembles.
Both share the trunk at the top.*

## Pipeline summary

1. Extract reads, discard low-quality ones, remove everything matching the human genome.
2. **Quantify the island** — counts for each of the 19 *clb* genes, island coverage, a coverage plot. *Always runs.*
3. **Score the evidence** — an [evidence tier](#evidence-tiers) per sample from read count, gene count and coverage breadth together. *Always runs for tumour and metagenome data.*
4. **Reassemble the island** from its own reads, using two independent assemblers, and report whether they agree. *Tumour data.*
5. **Recover draft genomes** from a mixed sample and identify which one carries the island. *Metagenomes.*
6. **Describe the neighbourhood** — prophages, mobility markers, and each neighbouring organism's DNA-damage-response genes. Written to `by_sample/<sample>/community/`, separate from the genomes themselves.
7. **Type the strain** — sequence type, clonal complex and phylogroup against a 45,761-genome panel. *Whenever assembly produced contigs.*

> **Typing needs a complete allele profile, and tumour-derived assemblies rarely give one.**
> MLST reports a type only from exact allele calls; a partial profile is reported as
> `insufficient_loci`, never guessed at. On fragmented, low-coverage tumour contigs the usual
> outcome is `insufficient_loci` rather than a sequence type — in the v0.0.2 test set, no
> tumour-derived assembly was typeable.

Steps 1–3 and 7 need no flags. The rest you switch on — see [Adding capabilities](#adding-capabilities).

## Evidence tiers

Read-level evidence is reported as a tier, not a yes/no call. All three criteria must be met, and a
sample gets the highest tier it satisfies.

| Tier | *clb* reads | *clb* genes | Island breadth at ≥1× | Reassembled by default? |
|---|---:|---:|---:|:--:|
| `extensive_island` | ≥100 | ≥10 | ≥15% | yes |
| `broad_island` | ≥30 | ≥8 | ≥7.5% | yes |
| `multi_gene` | ≥5 | ≥3 | ≥1% | no |
| `localized_indeterminate` | ≥1, but fails one of the above | | | no |
| `negative` | 0 | — | — | no |

*Breadth* is the fraction of the island's 50,767 bp covered by at least one read. Every threshold
is a parameter, e.g. `--tumor_broad_island_min_pks_reads 50`.

> **Why 15% and 7.5% and not 20% and 10%.** v0.0.1 measured breadth from a binned coverage
> track — 50 bp bins, where a bin counted as covered if any single base in it was. v0.0.2
> measures per base with `samtools depth`, which is the more honest quantity and reads a median
> **0.794×** lower on the same data. The v0.0.1 thresholds of 20% and 10% were set on the binned
> axis, so carrying them across unchanged would have made v0.0.2 stricter by accident rather than
> by choice. 0.20 × 0.794 = 0.159 and 0.10 × 0.794 = 0.079; a grid search against the v0.0.1 tier
> assignments lands on the same pair. These are a **unit conversion, not a recalibration** — they
> preserve v0.0.1's effective stringency on the new axis and are no better validated than the
> numbers they replace.

> **These tiers rank how strong the evidence is. They are not a validated positive/negative test.**
>
> They were derived from cohorts in which every sample already carried at least one *pks* read, so
> what they measure directly is how well designated positives are retained. They cannot tell you
> where true positive stops and background mapping begins: that boundary was never fitted.
>
> There is one piece of specificity evidence. Across **420 matched TCGA tumour/normal pairs**, no
> normal was called positive at island breadth ≥ 0.05 or above — 0 of 420, at every threshold at or
> above that point. That bounds the false-positive rate in matched normal tissue; it does not
> validate the tier boundaries themselves, which sit well above it.
>
> A small nonzero breadth is **indeterminate**, not positive. No tier here is a clinical result.

## Quick start

Requires Linux, [Nextflow](https://www.nextflow.io/docs/latest/install.html) 24.10 or newer
(enforced by the pipeline manifest), **Java 17+ on every node that runs a task**, and Conda or
Mamba. **There is no container support** — no process declares a container image, so
`-profile singularity` and friends are deliberately absent rather than present-and-broken.
Per-tool environments are built from `conda_envs/` on the first run.

Run a fixed release directly, which is the reproducible form and needs no clone:

```bash
nextflow run AlexandrovLab/pksProfiler -r v0.0.1 --help
```

Or clone to work on it. `main` is a moving target, so record the commit if you clone:

```bash
git clone https://github.com/AlexandrovLab/pksProfiler.git
cd pksProfiler
git rev-parse --short HEAD    # note this alongside your results
```

A sample sheet needs a unique `patient` column plus your files:

```csv
patient,bam
TUMOR_01,/data/TUMOR_01.bam
```

`fastq1`/`fastq2` and `cram`/`cram_reference` columns work the same way; with
`--input_data_type auto` each row is detected independently, so one sheet may mix them. Templates
are in [`examples/sample_sheets/`](examples/sample_sheets) — **replace the placeholder paths with
your own files.**

The minimum run, which quantifies the island and scores the evidence:

```bash
nextflow run main.nf -profile local \
  --sample       samples.csv \
  --input_data_type bam \
  --sample_type  tumor_wgs \
  --hg38_db      /refs/human-GRC-db.mmi \
  --t2t_phix_db  /refs/human-GCA-phix-db.mmi \
  --outdir       results
```

Start with `results/cohort/qc/pks.qc.summary.tsv`, which shows how many reads survived each
stage — it is the fastest way to see where a sample lost its reads.

## Adding capabilities

Each row adds flags to the command above. Nothing is replaced.

| To also get | Add | Details |
|---|---|---|
| island reassembly, two assemblers | `--sample_type tumor_wgs` (already above) | [assembly](docs/running/assembly.md) |
| prophage and mobility context | `--tumor_full_contig_context true --genomad_db <dir>` | [assembly](docs/running/assembly.md) |
| which organism carries it | `--sample_type metagenome --enable_mags true` plus `--gtdbtk_db --checkm2_db --genomad_db` | [MAGs](docs/running/mags.md) |
| community and neighbour context | automatic with `--enable_mags` | [community context](docs/running/community_context.md) |
| strain type and phylogroup | automatic; disable with `--enable_strain_typing false` | [strain typing](docs/running/strain_typing.md) |
| species abundances | `--pks_taxa true` or `--pks_community_taxa true`, plus `--kraken_db --bracken_read_length` | [taxonomy](docs/running/taxonomy.md) |
| profile-model search instead of alignment | `--profiling_method hmm` or `both` | [HMM modes](docs/running/hmm.md) |
| faster profiling on metagenomes | `--prefilter_mode balanced --kraken_db <dir>` | [prefilter](docs/prefilter.md) |

## Databases

Only the first two are ever required.

| Flag | Needed for |
|---|---|
| `--hg38_db`, `--t2t_phix_db` | removing human reads — every run |
| `--pangenome_db` | optional extra human-read removal pass |
| `--kraken_db` | species abundances, and the optional prefilter |
| `--genomad_db` | prophage detection |
| `--gtdbtk_db`, `--checkm2_db` | naming and quality-scoring draft genomes |

## Output

```text
results/
├── by_sample/<sample>/     everything this one sample produced
│   ├── read_evidence.tsv   evidence tier, reads, genes, breadth
│   ├── counts.txt          reads per clb gene
│   ├── alignment/          bam, coverage, per-sample QC
│   ├── contigs/            recruitment, both assemblies, final evidence
│   ├── genomes/            draft genomes, completeness, species, annotation
│   ├── community/          prophages, flanking genes, producer-neighbour tables
│   ├── strain/             sequence type, clonal complex, phylogroup
│   └── figures/            coverage plot, taxa barplots, contig validation
├── cohort/                 gene counts, QC summary, taxonomy, master summary
└── runs/                   Nextflow report, timeline and trace
```

The sample directory is the unit: if it is not under `by_sample/<sample>/`, it is not about
that sample. Both sample types write the same paths — `--sample_type` decides which
directories exist, not where they sit.

A sample with no *pks* reads stays in every table as a row of zeros rather than disappearing, so
"tested and negative" is distinguishable from "never ran". Full layout and per-file notes are in
the [output reference](docs/output.md).

### One table with everything

Each stage writes its own results, which is awkward to read across. This joins them into one row
per sample — whatever ran, with `NA` where a stage did not:

```bash
python3 scripts/build_master_summary.py --results results \
  --output results/cohort/pks.master_summary.tsv
```

Columns appear only for stages that actually ran, so the width of the table tells you what the run
did. With everything enabled you get read counts, the evidence tier and breadth, the structural
call from reassembly, the *pks*-positive genome's species and completeness, prophage and neighbour
counts, island mobility flags, and the sequence type and phylogroup — about 30 columns. It reads
only published output, so it is safe to re-run at any time without re-running the pipeline.

### Example output

| Sample | Alignment counts | HMM counts |
|---|---|---|
| Synthetic positive | [table](examples/results/synthetic_positive/pks.gene.counts.align.txt) | [table](examples/results/synthetic_positive/pks.gene.counts.hmm.txt) |
| Synthetic negative | [table](examples/results/synthetic_negative/pks.gene.counts.align.txt) | [table](examples/results/synthetic_negative/pks.gene.counts.hmm.txt) |

![Coverage plot for a pks-positive sample](examples/plots/synthetic_positive/synthetic_pks_positive.pks.circos.png)

## Running on HPC

Swap `-profile local` for your scheduler and add `-resume` so an interrupted run continues rather
than restarting. Profiles for TSCC, Slurm, PBS Pro, LSF and SGE are in [`conf/`](conf).

Two things that catch people out: **Java 17+ must be present on the compute nodes**, not just the
login node, and concurrent runs each need their own launch directory or they will collide on
Nextflow's session lock. See the [HPC guide](docs/hpc.md).

## Documentation

| Guide | |
|---|---|
| [BAM and CRAM mode](docs/running/bam.md) | alignment profiling from aligned input |
| [FASTQ mode](docs/running/fastq.md) | single or paired FASTQ |
| [HMM and combined modes](docs/running/hmm.md) | profile-model search, and how it differs from alignment |
| [Island reassembly](docs/running/assembly.md) | evidence tiers, targeted assembly, mobility context |
| [Genome-resolved analysis](docs/running/mags.md) | recovering draft genomes from a mixed sample |
| [Community context](docs/running/community_context.md) | prophages, island mobility, neighbouring organisms |
| [Strain typing](docs/running/strain_typing.md) | sequence type, clonal complex, phylogroup |
| [Taxonomy](docs/running/taxonomy.md) | species abundances |
| [Read prefilter](docs/prefilter.md) | optional narrowing before profiling |
| [Output reference](docs/output.md) | what every directory and file contains |
| [Glossary](docs/glossary.md) | reads, contigs, MAGs, breadth, tiers and the rest |
| [HPC guide](docs/hpc.md) | schedulers, resources, known pitfalls |
| [Reference models](ref/hmm/PROVENANCE.md) | how the *clb* profile models were built and validated |
| [Tool citations](CITATIONS.md) | every tool the pipeline calls, with DOIs |

The test suite runs on every push. To run it yourself: `bash tests/run_checks.sh` for the
full set, or `python3 -m unittest discover -s tests` for the fast checks that need no Nextflow.

## Citation

A manuscript is in preparation. Until it is available, cite the software itself — GitHub and
Zenodo both read [`CITATION.cff`](CITATION.cff), so the "Cite this repository" button on the
repository page produces a correctly formatted reference. Record the release tag or commit you
ran.

Please also cite the underlying tools your run used; they are listed with DOIs in
[`CITATIONS.md`](CITATIONS.md).

The biological premise — that colibactin-producing *E. coli* leave a distinctive mutational
signature in colorectal cancer — rests on Pleguezuelos-Manzano *et al.* and
Dziubańska-Kusibab *et al.*, both *Nature* 2020.

## License

[BSD 2-Clause](LICENSE).
