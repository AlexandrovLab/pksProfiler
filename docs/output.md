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
├── read_evidence.tsv       evidence tier, reads, genes, breadth      (step 2-3; wgs, metagenome)
├── wes_evidence.tsv        WES off-target DNA evidence, not tiered   (--sample_type wes)
├── rna_evidence.tsv        RNA expression evidence, not tiered       (--sample_type rna)
├── counts.txt              reads per clb gene
├── coverage.txt            island coverage
├── alignment/              sorted.bam, .bai, coverage.bedgraph, alignment.qc.tsv,
│                           ihe3034_divergence.qc.tsv, ihe3034_divergence_by_gene.tsv
│                           (A5: divergence from IHE3034 of the clb-assigned reads; QC only)
│                           coverage.bedgraph is **raw per-bin depth** over reads that
│                           aligned to the island reference at MAPQ >= 40. It is not a
│                           library-normalised abundance and is not comparable between
│                           samples; it drives the circos figure only. Gene counts are
│                           the primary read-level evidence.
├── hmm/                    nhmmscan tblout, hmm_counts.tsv, qc       --profiling_method hmm
├── contigs/                                                          (step 4, tumour)
│   ├── recruitment/        pks_recruitment.tsv: reads/fragments per recruitment class,
│   │                       island_supporting_fragments (island_internal + boundary_spanning)
│   ├── megahit/            contigs + canonical-reference PAF
│   ├── metaspades/         contigs + canonical-reference PAF
│   └── final_evidence/     final_pks_evidence.tsv
├── genomes/                                                          (step 5, --enable_mags)
│   ├── assembly/  bins/  checkm2/  gtdbtk/  annotation/
│   │   assembly/mate_pairing.tsv: pairs rebuilt after host depletion (see below)
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
├── provenance/             input_header.tsv: the input BAM/CRAM's reference and aligner
│                           (E3; alignment input only, see below)
├── figures/                circos.pdf, contig_validation.svg,
│                           pks_island_taxa_barplots.pdf  (--pks_taxa: Bracken on island reads)
│                           community_taxa_barplots.pdf   (--pks_community_taxa: whole community)
└── intermediates/          extracted / filtered / host-depleted FASTQ  --save_intermediates
```

Both sample types write to the same paths. There is no `wgs/` or `metagenome/`
level: `--sample_type` decides which directories exist, not where they sit.

## cohort/

```text
cohort/
├── gene_counts/            pks.gene.counts.align.txt, pks.gene.counts.hmm.txt
├── qc/                     pks.qc.summary.tsv — read attrition through every stage
├── provenance/             input_headers.tsv (every input's header provenance, one row per
│                           sample) and input_header_summary.tsv (counts only)   BAM/CRAM input
├── taxonomy/               combined bracken tables, clb taxonomy support
└── pks.master_summary.tsv  written by build_master_summary.py (see below)
```

Start with `cohort/qc/pks.qc.summary.tsv` when a result looks wrong — it shows how many
reads survived extraction, fastp and each host-depletion pass, so you can see where a
sample lost its reads.

### Input alignment provenance (E3)

For BAM/CRAM input, `extractReads` reads the header (`@SQ` and `@PG` lines only, no read
record) and writes `by_sample/<sample>/provenance/input_header.tsv`:

| Column | Meaning |
|---|---|
| `sq_fingerprint` | md5 over the ordered `SN:LN` lines. Same reference in the same order, same value; any other build differs |
| `m5_set_fingerprint` | md5 over the sorted `M5` values; `NA` unless every `@SQ` has one. Insensitive to contig renaming (`chr1` vs `1`) and order |
| `n_contigs`, `n_contigs_with_m5` | `@SQ` count, and how many carry `M5` |
| `chr1_name`, `chr1_length`, `build_guess` | `GRCh38` (248,956,422), `GRCh37` (249,250,621), `T2T-CHM13` (248,387,328), `unknown`, or `no_chr1` |
| `n_pg`, `pg_chain` | `@PG` records along their `PP` links as `ID\|PN\|VN`; `CL` is never recorded |
| `pg_programs` | the same chain as `PN:VN` |
| `aligner_bwa_mem` | E1: `yes` when a `@PG` command line runs `bwa mem` / `bwa-mem2 mem`, `unconfirmed` for a bwa `@PG` without one, else `no`. The `CL` is read in memory for this and for `-T`; nothing else of it is kept |
| `bwa_min_score` | the `-T` of that `bwa mem` (`default` = 30), `NA` without one |
| `extract_low_score_mapped` | what `--extract_low_score_mapped` resolved to for this sample: `on_auto_bwa_mem`, `off_auto_not_bwa_mem`, `on_requested`, `off_requested` |
| `extract_min_as` | the `AS` threshold in force, `NA` when off |

`cohort/provenance/input_header_summary.tsv` counts distinct fingerprints (with build and
contig count), distinct `PN:VN` chains and inputs per `extract_low_score_mapped` state, and
`RUN_REPORT.txt` prints it as **INPUT ALIGNMENTS**. It carries no sample identifier, `@PG` ID or path.

## runs/

`report.html`, `timeline.html` and `trace.txt`, written by Nextflow itself and
configured in `nextflow.config`. They overwrite on re-run: a resumed run that reused a
report name previously died with `FileAlreadyExistsException` after every task had
already succeeded.

`trace.txt` carries per-process CPU time and peak RSS, which is where runtime and memory
figures should come from.

## One table with everything

`cohort/pks.master_summary.tsv`, written automatically at the end of every run by the
`masterSummary` process, which wraps `scripts/build_master_summary.py` and joins every
stage into one row per sample. Every file it reads (MAG, community/prophage context,
strain typing, and the align and targeted-assembly outputs `pks_cohort_report.html`
also reads) is handed to the task as a declared input, with a manifest naming each by
sample and kind; neither it nor the cohort report reads the published tree. So a slower
stage that has not finished publishing, or a file left over from an older run at the
same `--outdir`, cannot change either table.

It is still safe to run by hand against a finished (or partly finished) results
directory, exactly as before -- it only reads published output and does not require
re-running the pipeline:

```bash
python3 scripts/build_master_summary.py --results results \
  --output results/cohort/pks.master_summary.tsv
```

Column groups appear only when the stage that produces them ran:

| Group | Columns |
|---|---|
| sample mode | `mode` (the sample type) and `evidence_label` (what the row's evidence may be read as); always second and third |
| terminal state | `terminal_state`, copied from the QC summary (fourth) |
| read profiling | `input_alignment_records`, `total_primary_reads`, `extracted_unmapped_reads`, `filter_input_reads`, `reads_after_*`, `num_clb_genes_*`, `reads_clb_genes_*` |
| evidence tier | `read_evidence`, `pks_reads_per_million` and `host_depleted_reads` (metagenome), `pks_reads`, `clb_genes_detected`, `island_breadth_1x` |
| WES evidence | `clb_wes_reads`, `wes_clb_genes_detected`, `wes_island_breadth_1x` |
| RNA evidence | `clb_rna_reads`, `rna_clb_genes_detected`, `rna_expression_breadth_1x`, `rna_library_type` |
| contig reassembly | `final_structural_evidence`, `assembler_agreement` |
| draft genomes | `mag_bins_total`, `mag_bins_pks_positive`, `pks_mag_taxonomy`, `pks_mag_completeness`, `pks_mag_clb_genes` |
| community context | `community_prophages_total`, `community_pks_producers`, `community_neighbours_assessed`, `community_neighbours_with_prophage`, `island_nearby_integrase`, `island_nearby_trna`, `island_in_prophage` |
| strain typing | `typing_status`, `ST`, `clonal_complex`, `phylogroup` |

A draft genome counts as *pks*-positive when its own assembly aligns across the canonical
island at a positive locus tier (`multi_gene`, `broad_island` or `extensive_island`, from
`locus_tier` in `pks_mag_summary.tsv`). HMM domain hits alone do not make a genome positive.
`mag_bins_total` counts every recovered bin. Where a sample has several positive genomes,
the reported one has the strongest locus tier, then the most canonical-locus *clb* genes
(`pks_mag_clb_genes` is that count), then the greatest island breadth. Contigs that were never
binned are reported separately (`mag_unbinned_*`) and are never counted as a genome. A sample that was profiled and found negative keeps its
row, with `NA` in the later groups.

### Evidence by mode

Ludmil, 26 Sep mode memo: every mode measures the same three things -- *clb*-assigned reads,
*clb* genes with at least one read, and breadth over the canonical island at >=1x/2x/3x --
and every mode names them for what they can mean. Each file's second column is `mode`.

| Mode | File | Columns after `sample`, `mode` | Tier? | Read as |
|---|---|---|---|---|
| `wgs` | `read_evidence.tsv` | `read_evidence`, `pks_reads`, `clb_genes_detected`, `island_breadth_{1,2,3}x` | yes (short reads) | genomic *pks* evidence |
| `metagenome` | `read_evidence.tsv` | same, plus `host_depleted_reads` and **`pks_reads_per_million`**; `read_evidence` is `untiered` | no | microbial read evidence; read the rate first |
| `wes` | `wes_evidence.tsv` | `clb_wes_reads`, `wes_clb_genes_detected`, `wes_island_breadth_{1,2,3}x` | no | exploratory off-target DNA evidence; a negative is not absence |
| `rna` | `rna_evidence.tsv` | `clb_rna_reads`, `rna_clb_genes_detected`, `rna_expression_breadth_{1,2,3}x`, `rna_library_type` | no | expressed *clb* sequence; breadth is expression breadth, not island completeness; a negative is not genomic absence |

**Metagenome abundance** (26 Sep mode memo section 5). `pks_reads_per_million` is
`pks_reads * 1e6 / host_depleted_reads`, where `host_depleted_reads` is the final
host-depleted count (`reads_after_pangenome` when that pass ran, else
`reads_after_t2t_phix`) -- host depletion, not the prefiltered profiling input, so the rate
does not move with `--prefilter_mode`. It is the primary metagenome quantity because stool
cohorts differ greatly in depth; the raw `pks_reads` stays beside it. `NA` for other modes
(and when nothing survived depletion). It is in `read_evidence.tsv`, the cohort report TSV
and `pks.master_summary.tsv`.

**Metagenome read pairs** (section 5, A7). Host depletion evaluates reads independently;
the MAG lane then rebuilds R1/R2 plus singletons from the mate suffixes of the complete
host-depleted stream (`scripts/recover_recruited_mates.py --all-reads`) and gives MEGAHIT
`-1/-2` plus `-r`, and read-to-contig alignment the pairs as pairs. Every read is kept
(pairs + singletons == the host-depleted stream). `genomes/assembly/mate_pairing.tsv`
reports `input_reads`, `paired_fragments`, `singleton_reads` and the mate-tag diagnostics
`mate_tag_1`, `mate_tag_2`, `fragments_without_mate_tag`, `orphan_reads`,
`duplicate_mate_slot_overwrites` (the last warns when > 0 and fails the task above 0.1% of
reads). Profiling is unchanged. The MAG lane assembles with MEGAHIT only; metaSPAdes runs
in the WGS targeted lane, not here.

The cohort report TSV (`cohort/pks_cohort_report.tsv`) keeps one set of measurement columns
for every mode, so it carries `mode` and `evidence_label` beside them; its `read_evidence`
column holds a tier only for `wgs` rows (`rna_expression` / `wes_offtarget` otherwise).

## Primary result files

Every table carries a header row. The three tables produced by every run are specified
below; the pages under `docs/running/` describe the files each optional mode adds.

| File | Produced by | One row per |
|---|---|---|
| `cohort/gene_counts/pks.gene.counts.*.txt` | every run | *clb* gene |
| `cohort/qc/pks.qc.summary.tsv` | every run | sample |
| `cohort/pks.master_summary.tsv` | every run | sample — every stage joined |
| `by_sample/*/read_evidence.tsv` | `--sample_type wgs` or `metagenome` | sample |
| `by_sample/*/wes_evidence.tsv` | `--sample_type wes` | sample |
| `by_sample/*/rna_evidence.tsv` | `--sample_type rna` | sample |
| `by_sample/*/contigs/final_evidence/final_pks_evidence.tsv` | targeted reassembly | sample |

### `cohort/qc/pks.qc.summary.tsv`

Read attrition through every stage, one row per sample. Counts are read counts and
decrease monotonically down the depletion chain; a violation is a fatal error, not a
warning.

| Column | Meaning |
|---|---|
| `Sample` | sample identifier, as given in the sheet's `patient` column |
| `input_alignment_records` | records in the input alignment, from its index — secondary and supplementary included. Absent for FASTQ input and when the alignment has no index |
| `total_primary_reads` | exact count of primary records in the input alignment. Only written with `--exact_input_counts`, which costs a second decode of every input |
| `extracted_unmapped_reads` | FLAG 4 (unmapped) primary reads extraction wrote out |
| `extracted_low_score_mapped_reads` | E1: primary *mapped* reads with `AS:i` below `--extract_min_as` (30) that extraction also wrote; 0 when `--extract_low_score_mapped` resolved off ([BAM and CRAM mode](running/bam.md#low-score-mapped-reads-bwa-mem--t-0)) |
| `extracted_reads` | everything extraction wrote: `extracted_unmapped_reads` + `extracted_low_score_mapped_reads` (checked). This, not the FLAG 4 count, is what enters fastp and the extraction step of the attrition chain |
| `extract_low_score_no_as_records` | mapped primary records without an `AS` tag, never added; only when low-score extraction was on |
| `extract_low_score_mode` | `on_auto_bwa_mem`, `off_auto_not_bwa_mem`, `on_requested` or `off_requested` |
| `filter_input_reads` | reads entering fastp. For alignment input this equals the extracted count; for FASTQ input it is the library |
| `taxonomy_status` | why a sample is or is not in the species table: `no_pks_reads`, `below_rank_threshold`, `no_species_identified`, `species_identified`. `NA` when taxonomy did not run |
| `clb_species_reported` | how many species that sample contributed to `pks.clb_species_support.tsv` |
| `reads_mapped_ihe3034_mapq40` | reads that mapped to the reference at all, at MAPQ ≥ 40 (named `reads_mapped_ihe3034` before rc1) — the denominator the *clb* counts sit inside |
| `unmapped_reads` | unmapped primary records extracted (`-f 4 -F 2304`) |
| `reads_after_fastp` | surviving adapter and quality filtering |
| `reads_after_hg38` | surviving GRCh38 depletion |
| `reads_after_t2t_phix` | surviving T2T + PhiX depletion |
| `reads_after_pangenome` | surviving optional pangenome depletion; `NA` unless `--pangenome_db` |
| `num_clb_genes_align` | *clb* genes with at least one assigned read, alignment profiling (0–19) |
| `reads_clb_genes_align` | reads assigned to *clb* genes, alignment profiling |
| `divergence_ihe3034_reads` | reads the divergence QC measured: primary alignments at MAPQ ≥ 40 that featureCounts' largest-overlap rule assigns to one *clb* gene (A5). 0 when none; the columns below are then `NA` |
| `nm_ihe3034_median` | median edit distance (`NM`) of those reads to IHE3034 |
| `divergence_ihe3034_median`, `_q1`, `_q3` | median and quartiles of divergence from IHE3034, `NM` / alignment columns (M/=/X + I + D); 0 is identical to IHE3034 |
| `divergence_ihe3034_frac_gt_2pct`, `_frac_gt_5pct` | fraction of those reads more than 2% / 5% divergent from IHE3034 |
| `clb_read_length_median`, `_q1`, `_q3`, `_min`, `_max` | length of those reads, clipped bases included |
| `num_clb_genes_hmm` | as above, HMM profiling; `NA` unless `--profiling_method hmm` or `both` |
| `reads_clb_genes_hmm` | as above, HMM profiling |
| `status` | `complete`, `incomplete` (sample ran partially), or `no_qc_produced` (sample entered the run but emitted no QC) |
| `terminal_state` | how the sample ended; exactly one value per sample, below |

**Reading the divergence columns (Ludmil, 25 Sep A5).** They describe how far the
clb-supporting reads sit from one reference, IHE3034, and nothing more. MAPQ ≥ 40 means a read
placed confidently on IHE3034, not that no other microbe would take it, so a low divergence is
not evidence of uniqueness, organism identity or carrier identity, and a high one is not
evidence against the call. Specificity is the competitive-mapping analysis
(`Modules/pks_specificity.nf`, run standalone), stated relative to the panel it tested. Soft-clipped bases are outside
the alignment, so a local-mode alignment can understate divergence; long reads carry their
sequencing error in `NM`. `by_sample/<sample>/alignment/ihe3034_divergence_by_gene.tsv` gives
reads, median `NM` and median divergence per *clb* gene. Computing them changes no count, tier
or call.

A sample that failed mid-run still gets a row. The table is built from whatever
fragments exist, so one sample hitting a time cap no longer costs the cohort its QC —
check `status` before treating an `NA` as a biological zero.

**`terminal_state`** (GAMEPLAN v1.0.0-rc1). Every expected sample ends in exactly one of:

| Value | Meaning |
|---|---|
| `completed_positive` | every stage ran; tiered short-read `wgs` at a formal tier (`multi_gene`, `broad_island`, `extensive_island`) |
| `completed_negative` | every stage ran; tiered short-read `wgs` with 0 *clb* reads or `localized_indeterminate` (not a formal positive tier) |
| `completed_untiered` | every stage ran; the evidence carries no call: `wes`, `rna`, `metagenome`, long-read `wgs`, or an HMM-only run with no evidence step |
| `failed_at_<stage>` | the first stage, in pipeline order, the sample has no QC for: `extraction`, `read_qc`, `host_depletion`, `prefilter`, `diamond_rescue`, `candidate_merge`, `alignment_profiling`, `hmm_profiling`, `evidence`, `taxonomy` |
| `excluded_<reason>` | named in `build_qc_summary.py --excluded-samples` (`sample<TAB>reason`) |

`main.nf` passes the stages the run should have produced (`extraction` only for BAM/CRAM,
the prefilter trio only with `--prefilter_mode balanced`, `taxonomy` only with `--pks_taxa`).
The positive/negative split is the evidence tier, written by `classifyPksReadEvidence` as
an `evidence_outcome` QC fragment; the tier in `read_evidence.tsv` stays authoritative.
Not part of the terminal state: targeted assembly, MAGs, community context and strain
typing, which have their own status columns. The pipeline itself excludes no sample yet;
`--excluded-samples` is for a curated exclusion list passed to a hand run.

### `cohort/gene_counts/pks.gene.counts.*.txt`

One row per *clb* gene, one column per sample, plus a leading `Gene` column. Values are
read counts assigned to that gene in that sample. `.align.txt` comes from Bowtie2 plus
featureCounts, `.hmm.txt` from `nhmmscan`. Genes are always all 19 of clbA–clbS, in
order, whether or not any sample carried reads for them.
| `by_sample/*/genomes/pks_mag_summary.tsv` | `--enable_mags` | draft genome |
| `by_sample/*/community/community_prophage_inventory.tsv` | `--enable_mags` | prophage |
| `by_sample/*/community/pks_community_interactions.tsv` | `--enable_mags` | producer × neighbour pair |
| `by_sample/*/community/pks_island_mobility.tsv` | `--enable_mags` | draft genome |
| `by_sample/*/strain/strain_types.tsv` | assembled output | typed unit |

## Reading results from before v0.0.2

Runs published under the old layout (`pks_per_sample/`, `pks_summary/`, `mags/`,
`community_context/`, `strain_typing/`, `prefilter/`) are not migrated. Keep them as
they are; `build_master_summary.py` reads the new layout only.

---

Back to the [README](../README.md).
