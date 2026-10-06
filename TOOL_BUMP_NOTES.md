# Tool bump notes (feature/tool-bump, 2026-10-06)

Every pinned conda tool is moved to its latest bioconda release (checked against
`api.anaconda.org/package/bioconda/<pkg>` on 2026-10-06). Base: `release/v1.0.0-rc1` @ f505b4e.
This file lists what changed, what broke, what was fixed, and what can change counts. Read it before
using a bumped run beside a pinned one.

## Pins

| env | tool | old | new | note |
|---|---|---|---|---|
| minimap2, pks_align, samtools | samtools (+htslib) | 1.21 / 1.22 | 1.24 | |
| pks_align, targeted_recruit | bowtie2 | 2.5.4 | 2.5.5 | |
| fastp | fastp | 0.24.0 | 1.3.7 | major |
| pks_align | deeptools | 3.5.6 | 4.0.0 | major; **python 3.10 -> 3.12** (4.0.0 ships py312/py313 builds only) |
| pks_hmm | hmmer / seqtk / seqkit | 3.3.2 / 1.4 / 2.6.1 | 3.4 / 1.5 / 2.14.0 | |
| krakenUniq_bracken | bracken | 2.9 | 3.1p1 | prints `v3.0.1`, see below |
| targeted_spades | spades | 4.2.0 | 4.3.0 | |
| metabat2 | metabat2 | 2.17 | **==2.18** | `=2.18` fuzzy-matches the untagged `2.18_23_gc869c52` git snapshot, so the pin is exact |
| mlst | mlst | 2.23.0 | 2.35.0 | |
| prokka | prokka | 1.14.6 | 1.15.6 | |
| prefilter, checkm2 | diamond | 2.1.10 / 2.1.11 | 2.2.8 | checkm2 then resolves to bioconda build `pyhdfd78af_2` (2026-10-04) |

Already latest, unchanged: minimap2 2.31, subread 2.1.1, bedtools 2.31.1, fastplong 0.7.1, krakenuniq
1.0.4, krakentools 1.2.1, megahit 1.2.9, checkm2 1.1.0, gtdbtk 2.7.2, genomad 1.12.0. Every env was
dry-run solved (`conda create --dry-run`, libmamba, private package cache), then built and run.

## Breaking changes found, and the fix

1. **bracken 3.1p1 reports itself as `Bracken v3.0.1`.** The upstream v3.1 tag never updated the VERSION
   string in the `bracken` script, and bioconda's 3.1p1 builds that tag. Strict `--tool_version_check`
   would have stopped every Bracken task (exit 78, `VERSION_MISMATCH`), and RUN_REPORT would have flagged
   `INSTALLED_MISMATCH`. Fix: `scripts/tool_version_aliases.tsv` records this one known self-report
   (package, conda version, printed version, source). `record_tool_versions.sh` and `finalize_provenance.py`
   both read it. An alias covers exactly that conda version: a stale 2.9 binary under the 3.1p1 pin, or 3.0.1
   under another pin, is still a mismatch. A unit test fails if an alias names a version no env pins any more.
2. **deeptools 4.0.0 needs python >= 3.12.** `pks_align_env` moves from python 3.10 to 3.12 (pandas 2.3.3,
   numpy 2.5). Its scripts ran under 3.12 in every execution check.
3. **metabat2 `=2.18` resolves to a git snapshot**, so it is pinned `==2.18`. MetaBAT2 (2.17 and 2.18) aborts
   on `--version`, so it was always recorded `version_unread`. It is now read from `metabat2 -h` and the
   no-argument usage of `jgi_summarize_bam_contig_depths`, and both record `2.18 ok`.
4. `tests/test_tool_versions.py` fixtures that read the real `minimap2_env.yml` expect samtools 1.24.

No call site needed a CLI change. Checked against each tool's release notes and source: every flag
the pipeline passes still exists with the same meaning (fastp `-l --adapter_fasta --cut_tail -i --stdin -w
--compression --json --html -o`; bamCoverage `-b -o --binSize --normalizeUsing None --outFileFormat
bedgraph --numberOfProcessors`, now `deeptools.bamCoverage2:main` on the Rust core with the same argparse;
diamond `makedb`, `blastx --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore`; bracken
`-d -i -o -w -r -l -t`; bowtie2, samtools, nhmmscan/hmmsearch, seqkit `split2`/`grep`, seqtk `seq -a`,
metaspades.py, metabat2, prokka, mlst `--scheme ecoli_achtman_4 --quiet`). The fastp JSON writer is
unchanged, so the awk read of `summary.{before,after}_filtering.total_reads` still works. HMMER 3.4 reads
the repo's 3.3.2-pressed models.

## Behaviour changes that can move counts

**fastp 0.24.0 -> 1.3.7** (source diff of `main.cpp`, `evaluator.cpp`, `adaptertrimmer.cpp`, `filter.cpp`)
- *Poly-G auto-enable:* fastp switches poly-G trimming on when the **first read name** looks like a
  two-colour instrument. Old prefixes: `@NS @NB @NDX @A0`. New: `@FS @MN @SH @NS @NB @NDX @VL @VH @A @NA @LH`,
  so any name starting with `A` now counts. It applies only to file input (`-i`): the single-file path,
  which is every BAM/CRAM sample. `--stdin` (paired FASTQ) never evaluates.
- *Adapter trimming tolerates one indel* when the exact/Hamming match fails. This applies to `--adapters`,
  and to the adapter auto-detected for single-file input.
- *SE adapter auto-detection rewritten* (0.26.0: known-adapter scan over the first reads). It now also
  requires the input to be a regular file. Nextflow's staged symlinks pass, because the check is `stat`.
- *Adapter dimers* (reads <= 2 bp after trimming) are a new failure category. They were already removed by
  `-l 45`, so retention is unchanged. JSON `adapter_trimmed_reads` now also counts `--adapter_fasta` trims.
- Quality-window and complexity filters became SIMD/integer arithmetic. Same results.

Synthetic A/B (60,000 101-bp reads from the public IHE3034 genome: 30% adapter read-through, 40% of those with
a 1-base adapter indel, 10% poly-G tails, 15% low-quality tails, 3% N-rich; pipeline command; JSON only):

| first read name | input | adapters | reads kept new/old | bases kept new/old |
|---|---|---|---|---|
| `HWI-ST...` | file | auto | 1.0000 | 1.0000 |
| `HWI-ST...` | file | `--adapter_fasta` | 0.9963 | 0.9954 |
| any | `--stdin` | auto | 1.0000 | 1.0000 |
| any | `--stdin` | `--adapter_fasta` | 0.9964 | 0.9954 |
| `AB12CDXX...` (A, not A0) | file | auto | **0.9894** | **0.9714** |
| `A00123...` (A0) | file | auto | 1.0000 | 1.0000 |

So with the default (no `--adapters`), retention is unchanged unless the first read name newly matches the
two-colour list. When it does, poly-G tails are trimmed: here 1.1% fewer reads and 2.9% fewer bases, but this
synthetic set has 10% poly-G-tailed reads. Real 4-colour HiSeq reads (TCGA) rarely end in long G runs, so the
real effect should be far smaller. It is what `fastp_retained_fraction` in the 10-sample comparison measures.

**samtools 1.21/1.22 -> 1.24** (release notes 1.21.1-1.24). Default CRAM *output* is 3.1 (1.22); the pipeline
writes no CRAM, only the test fixtures do. HTSlib no longer fetches CRAM references from EBI (1.22); the
pipeline requires `--cram_reference`. `samtools fastq` warns on coordinate-sorted input (1.22; stderr only).
The `samtools depth` overlap-removal fix (#2170) applies to pairs, and the profiling BAM holds unpaired reads.
`view --subsample` seed: unused. No count change expected.

**bracken 2.9 -> 3.1p1.** The default `-t` went from 0 to 10; the pipeline always passes `-t 2`, so no effect.
`est_abundance.py` treats blank rank codes as `-`. Bioconda's `fix_bug.patch` stops the `-w` (kraken-style)
report from descending below the estimation level. On three synthetic KrakenUniq reports over the real
`krakenUniq_8_8_2023` taxonomy and its `database100mers.kmer_distrib`, the abundance tables (`-o`) and the
`kreport2mpa` tables were identical. Only the G-level `-w` report's "reads at this node" column changed for
genus rows (now the estimate, as for leaves).

**diamond 2.1.10/2.1.11 -> 2.2.8.** Removed options are not ones the pipeline uses. 2.1.12 fixed a parse
error on FASTQ quality lines starting with `@` (the blastx query here is FASTQ). Synthetic blastx A/B
(20,000 reads vs the clb reference proteins): byte-identical output.

**seqkit 2.9** fixed ID parsing for headers like `id\tx y`. This touches only the HMM lane's
qualifying-read FASTA, not counts. **bowtie2 2.5.5** adds `-d`; defaults are unchanged. **SPAdes 4.3.0**:
RAM and error-code changes; `--phred-offset 33` is still passed. **prokka 1.15.6**: new sprot/HAMAP
databases (product names can change), no bundled binaries. **mlst**: databases frozen at PubMLST
2024-12-31 (2.33.0); `ecoli_achtman_4` still exists; default output unchanged.

## Synthetic A/B through the whole pipeline (release-rc1 vs tool-bump)

`tests/check_e2e_execution.sh` on both branches: 7 runs (paired FASTQ, BAM, CRAM, HMM, 5-sample cohort,
single-end, bwa -T 0 low-score), 18 sample-runs. Base used the existing rc1 env cache; tool-bump used a
fresh private cache. `pks.qc.summary.tsv`, `pks.master_summary.tsv`, every `read_evidence.tsv` and every
`counts.txt` were compared cell by cell: **0 differing cells**. That covers filter_input_reads,
reads_after_fastp, host depletion, pks_reads, clb genes, breadth, tier and the HMM counts. Every published
file of the cohort run was also compared byte for byte. Identical: recruitment tables, MEGAHIT and
metaSPAdes contigs, final evidence. Different, as expected: BAM/BAI and recruitment FASTQ.gz (headers and
compression; their counts are identical), circos PDFs, RUN_REPORT (versions). The deepTools bedGraph is
identical per 50-bp bin (max difference 0). deepTools 4 just no longer splits a trailing zero run at the old
5-Mb chunk boundary (70 -> 69 lines).

Public-genome A/B on the IHE3034 assembly in `indices/` (the MAG/typing tools that no execution check runs
on a real bin):
- **mlst** 2.23.0 vs 2.35.0: identical (`ecoli_achtman_4`, ST95, same seven alleles).
- **prokka** 1.14.6 vs 1.15.6 (`--metagenome`, as `prokkaAnnotate`): 4,807 vs 4,805 CDS, 56 CDS coordinates
  differ; rRNA/tRNA/tmRNA identical. **hmmsearch** (3.3.2 vs 3.4, clb protein HMM, E <= 1e-5) on each one's
  proteins finds the same 19/19 clb genes with the same number of hit proteins per gene.
- **CheckM2** 1.1.0 with diamond 2.1.11 vs 2.2.8: identical (completeness 100.0, contamination 0.87, same
  model). Either pin now resolves to bioconda build `pyhdfd78af_2` (uploaded 2026-10-04: pandas 3,
  scikit-learn 1.9.1, tensorflow 2.19, keras 3). Both warn `InconsistentVersionWarning` (models pickled with
  scikit-learn 1.6.1). The prediction is unaffected here, but it is worth watching on real bins.
- **DIAMOND blastx** 2.1.10 vs 2.2.8 (prefilter call, 20,000 synthetic reads): byte-identical output.

Harness findings, which apply to anyone running these envs outside Nextflow's driver (the driver already
sets both): the pinned prokka env's GNU parallel 20260722 cannot exec without `SHELL` set (prokka then
dies at the IS blastp step). The bumped env's parallel 20260922 works without it. checkm2 imports from
`~/.local` unless `PYTHONNOUSERSITE=1`, and needs a short `TMPDIR` (`AF_UNIX path too long` on a long
Lustre path).

## Checks run on tool-bump (synthetic / public data, private conda cache)
`tests/run_checks.sh` (1287 unit tests OK; also OK under the bumped py3.12 `pks_align_env`;
22/22 workflow-construction checks) and `check_e2e_execution`, `check_mode_execution`,
`check_v8_lanes_execution`, `check_resume_reference_change`, `check_malformed_input_execution`,
`check_metagenome_mag_execution`: all pass. Every run is under strict `--tool_version_check`, with no
mismatch.

## Not yet measured
Real-sample effect, mainly fastp's poly-G rule on TCGA read names, and bowtie2/featureCounts on real reads.
The 10-sample pinned-vs-updated comparison is written for Ammal:
`claude-code/projects/pksprofiler/validation/TOOL_BUMP/`.
