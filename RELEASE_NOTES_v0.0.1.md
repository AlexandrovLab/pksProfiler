# pksProfiler v0.0.1

This is the first immutable development release of the benchmarked pksProfiler workflow. It is intended for reproducible method evaluation. It is **not** yet a calibrated pan-cancer presence/absence caller: nonzero counts alone must not be interpreted as a high-confidence biological positive.

## Benchmark evidence

The release candidate and pre-CRAM snapshot were run on the same metagenomic sample (`ERR525841`) with Bowtie2, HMM, taxonomy, and cohort summaries enabled under the TSCC profile.

| Metric | Pre-CRAM snapshot | v0.0.1 candidate | Result |
|---|---:|---:|---|
| Driver wall time | 11:28:53 | 8:51:29 | 22.85% faster |
| HMM process | 7:18:46 | 5:49:14 | 20.40% faster |
| fastp process | 1:59:58 | 1:19:46 | 33.51% faster |
| Bowtie2/featureCounts process | 24:46 | 19:44 | 20.30% faster |
| Alignment-assigned reads | 318,639 | 328,911 | +10,272 (+3.22%) |
| HMM, taxonomy, QC outputs | baseline | identical | preserved |

The alignment-count increase is expected: v0.0.1 uses `featureCounts --largestOverlap`, assigning an overlapping read to the gene with greatest aligned overlap rather than silently discarding it. All 19 `clbA`–`clbS` genes remain detected. Both benchmark workflows completed with exit code 0, and 16/16 release regression/validator tests pass.

## Tongwu Zhang review disposition

Statuses deliberately distinguish completed work from remaining scientific development.

| # | Concern | Status | Disposition |
|---:|---|---|---|
| 1 | Overlapping-gene assignment | Substantially addressed | Alignment uses `--largestOverlap`; taxonomy uses greatest overlap. Exact-tie/all-boundary fixtures still need expansion. |
| 2 | Paired templates and duplicates | Open | Evidence is counted per read. Template counts, duplicate-collapsed loci, and paired Bowtie2/featureCounts semantics are not implemented. |
| 3 | Formal biological call model | Partial | Stage and method QC are reported. Calibrated classes, breadth/depth thresholds, `call_reason`, and contamination-aware calls are absent. |
| 4 | Batch contamination/index hopping | Open | No batch-aware locus/sequence comparison or contamination score exists. Sparse signals remain indeterminate. |
| 5 | Specificity/sensitivity beyond one strain | Partial | HMM is complementary evidence and human-pangenome depletion is optional. Bowtie2 still uses one IHE3034 reference without microbial decoys. |
| 6 | CRAM/reference handling | Substantially addressed | `alignment`/legacy `bam`, auto/BAM/CRAM/FASTQ, SN/LN/M5 validation, and match/mismatch fixtures are present. Reference digest injection remains. |
| 7 | T2T/phiX contract | Partial | Documentation warns index composition cannot be inferred from its name. The combined parameter remains and does not independently verify phiX. |
| 8 | Complete stage QC | Partial | Extraction, fastp, depletion, alignment/HMM evidence and conservation checks are published. Full alignment categories, checksums, manifests, and explicit QC-fail calls remain. |
| 9 | Immutable reproducibility | Partial | This annotated tag freezes tested source and primary tool versions are pinned. Complete reference/input digests and containers remain. |
| 10 | Release test suite | Partial | Sixteen tests cover core references, summaries, QC, CRAM, bedGraph, taxonomy joins, and duplicate IDs. Boundary ties, equivalence, contamination and concurrency fixtures remain. |
| 11 | Taxonomy correctness | Addressed for benchmarked path | Real KrakenUniq/Bracken genus/species, plotting, and species-by-clb attribution completed. Variable trimmed lengths versus one Bracken `-r` still need study. |
| 12 | CPU/HPC accounting | Partial | Portable stage profiles and HMM chunking work; HMM time fell 20.40%. Some simultaneous samtools commands can still oversubscribe. |
| 13 | Storage amplification | Partial | Preparation intermediates are opt-in, but SAM/per-base coverage and work retention remain large; each benchmark directory was about 183 GiB. |
| 14 | Strict deterministic summaries | Partial | Exactly 19 unique ordered genes are required; malformed/duplicate/missing data fail. Atomic promotion, checksums, locks, and `RUN_COMPLETE.tsv` remain. |
| 15 | Hard-coded references | Open | `plotPKS.R` still hard-codes IHE3034 coordinates/contig, and several thresholds are not fully parameterized. |
| 16 | Permissions/temp paths | Open | Portable HPC profiles exist, but ACL/umask enforcement and complete temp/cache control remain untested. |

## Additional confirmed fixes

- Scientific-notation bedGraph is accepted; empty, negative, NaN, and Inf values are rejected.
- CRAM-reference mismatches and duplicate sample identifiers fail before analysis.
- Strict validators prevent malformed cohort matrices from being silently published.
- Final host-depleted FASTQ can be published with `--save_intermediates true`; automatic checksums remain future work.

## Interpretation limits

- Do not equate any nonzero count with a positive sample.
- Treat sparse/two-locus evidence as partial or indeterminate until contamination and template evidence are assessed.
- Bowtie2 and HMM are complementary; their counts are not interchangeable.
- Short-read taxonomy is context, not proof that a species carries an intact pks island.

## Deferred production gates

Before pan-cancer biological calling: implement template/duplicate metrics, calibrated calls, contamination controls, multi-reference/decoy benchmarking, phiX validation, complete provenance/checksums, boundary/tie fixtures, and removal of hard-coded reference assumptions.
