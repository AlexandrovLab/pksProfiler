# pksProfiler v1.0.0

pksProfiler finds and characterises the *pks* island — the 19 *clb* genes that produce
colibactin, a bacterial toxin that damages human DNA — in sequencing data. One pipeline
covers gut metagenomes, whole-genome, whole-exome and RNA sequencing, short reads or long.

## What every run does, in every assay

- **Takes FASTQ, BAM or CRAM.** Aligned input has its reads extracted first; FASTQ enters
  directly. BAM and CRAM are detected per file.
- **Removes human reads** against GRCh38 and T2T/phiX, splice-aware for RNA, and reports
  how many reads survive each step.
- **Counts reads on all 19 *clb* genes**, `clbA` through `clbS`, by alignment to the
  canonical IHE3034 island. Nucleotide HMMs can run instead of, or alongside, alignment —
  the two are always reported as separate numbers, never merged.
- **Measures how much of the island is covered**, as the fraction of its 50,768 bp seen at
  1×, 2× and 3×.
- **Writes a per-sample evidence file** whose fields suit the assay, plus one cohort table
  joining every sample and every stage that ran.
- **Records what it did**: every parameter, reference, database fingerprint and tool
  version, so a run can be reproduced or audited later.

## What each assay adds

- **Gut metagenomes** — *clb* reads per million host-depleted reads, so cohorts of
  different depth can be compared. Optionally assembles the community (MEGAHIT), groups
  contigs into genomes (MetaBAT2), assesses their quality (CheckM2) and taxonomy (GTDB-Tk),
  and reports which recovered genome carries the island, with its prophage and
  flanking-gene context and seven-locus *E. coli* strain type.
- **Whole-genome sequencing** — short reads get an evidence tier: `extensive_island`,
  `broad_island`, `multi_gene`, `localized_indeterminate` or `negative`, each requiring a
  read count, a gene count and an island breadth together. Qualifying samples then have the
  island reassembled from their own reads with two assemblers independently (MEGAHIT and
  metaSPAdes), and agreement between them is reported.
- **Whole-exome sequencing** — off-target *clb* reads are counted and reported as
  exploratory evidence, with capture QC when a target BED is supplied.
- **RNA-seq** — expressed *clb* sequence and expression breadth, with the number of genes
  passing a configurable read threshold.

## Optional analyses, switched on per run

- Island-read taxonomy, and whole-community taxonomy (KrakenUniq/Bracken)
- Locus-detection specificity against a *pks*-negative decoy panel
- Carrier attribution against a *pks*-positive panel, reporting ambiguity at strain and
  genus level
- A prefilter that keeps the bacterial fraction and rescues *clb*-like reads from what it
  discards

## What you get out

Everything lands under one results directory.

| | |
|---|---|
| `cohort/pks_cohort_report.html` | **The report to open first.** Self-contained HTML: every sample with its evidence, the measurements behind each call, whichever optional analyses ran, and a closing section stating where each number came from |
| `cohort/pks.master_summary.tsv` | One row per sample, one column per measurement, joining every stage that ran |
| `cohort/qc/pks.qc.summary.tsv` | Read attrition step by step, which stages ran, and how each sample ended |
| `cohort/gene_counts/pks.gene.counts.align.txt` and `.hmm.txt` | One *clb* gene per row, one sample per column — alignment and HMM kept apart |
| `by_sample/<sample>/` | That sample's evidence file, alignments, coverage plots, assembled contigs, recovered genomes and optional-analysis outputs |
| `RUN_REPORT.txt` | Parameters, references, database fingerprints, software versions, and what failed |

Where WGS reassembly ran, the report links a per-sample figure showing the reference gene
track, the contigs each assembler produced, and read depth across the island, so a call can
be inspected rather than taken on trust.

## Running it

Linux, Nextflow 24.10 or newer, Java 17, and Conda or Mamba. Every process uses a pinned
environment built on first run. You supply a sample sheet and the human-depletion indexes;
the bacterial reference and *clb* models ship with the pipeline. The README has a quick
start that runs in one command once the references are in place.

## What CI checks, and what it does not

Continuous integration on GitHub runs the unit suite and compiles the workflow against
the Nextflow version this release declares as its floor. It does **not** run the
end-to-end execution checks: that job has never passed on a GitHub-hosted image, for
reasons specific to that image rather than to the pipeline, and it is kept available to
run on demand rather than left failing in front of the release.

Those end-to-end checks are run and recorded on the cluster instead, under the pipeline's
default strict tool-version enforcement, and the release gate verifies that record
against a content digest of the tagged tree, so it cannot be a result from an earlier
state of the code. One further check, the metagenome genome-binning run, needs about
100 GB of reference data and is recorded as skipped under a signed waiver that the gate
prints on every tag build.

The short version: the pipeline's behaviour is verified where it actually runs, and the
release says so rather than implying a green tick covers more than it does.

## Reading the results honestly

- **A missing result is not a negative.** Check the QC table first: a sample that failed a
  stage, or had too few reads to be informative, is not the same as one with no island.
- **Tiers are evidence categories, not probabilities**, and not a clinical result.
- **Depth decides what a negative means.** Shallow sequencing can miss a real island,
  and metagenomes need substantially more depth than tumour genomes before absence means
  anything. The tools to quantify this per cohort are documented separately.
- **Taxonomy at read level is support, not attribution.** Finding *clb* reads classified to
  a species is not the same as showing that species carries the island; that needs the
  genome-resolved or carrier-attribution analyses.
- **Do not pool raw prevalence across cohorts** without checking that aligner, duplicate
  handling and upstream filtering match, because those differ between public datasets and
  can move prevalence on their own.
