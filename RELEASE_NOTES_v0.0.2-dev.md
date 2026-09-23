# pksProfiler v0.0.2 (in development)

Extends the immutable v0.0.1 release with sample-type routing, tumour contig analysis gated
on evidence tier, and genome-resolved MAG reconstruction with prophage and neighbouring-gene
context. Ported from the tested v3 integration tree onto v0.0.1 under version
control.

## What is new

**Sample-type routing.** `--sample_type auto | metagenome | tumor_wgs | tumor_wes | tumor_rna`.
Routing is explicit; the metagenome MAG path is branched off the **complete** host-depleted
read stream before profiling, so genome-resolved analysis sees every host-depleted read.

**Tumour contig analysis** (`--sample_type tumor_wgs`). `classifyTumorPksEvidence` assigns each
sample an evidence tier from read count, clb genes detected at >=1 assigned read, and canonical
island breadth at >=1x. Samples in a qualifying tier are re-joined to the complete host-depleted
stream and passed to targeted assembly: recruit -> MEGAHIT **and** metaSPAdes -> align contigs to
the canonical reference -> `summarizeTargetedPksEvidence`. Read pairing is preserved throughout.

**MAG reconstruction** (`--enable_mags`, metagenome only). MEGAHIT assembly, MetaBAT2 binning,
CheckM2 completeness, GTDB-Tk classification, Prokka annotation, clb protein HMM search,
geNomad prophage detection, community prophage summary, and genomic-context extraction.

**Prophage and neighbouring-gene context on tumour contigs** (`--tumor_full_contig_context`).
Prokka annotation, clb protein HMM search, geNomad prophage detection, and contig context for
assembled tumour contigs.

## Evidence tiers

The frozen three-tier system from `tier_threshold_development_20260909`. All three criteria must
hold for a tier; classification proceeds from the highest tier downward.

| Tier | pks reads | clb genes (>=1 read) | island breadth >=1x | qualifies for contig analysis |
|---|---:|---:|---:|:--:|
| `negative` | 0 | — | — | no |
| `localized_indeterminate` | >=1 | fails a multi_gene criterion | | no |
| `multi_gene` | >=5 | >=3 | >=1% | no (classifies only) |
| `broad_island` | >=30 | >=8 | >=7.5% | yes |
| `extensive_island` | >=100 | >=10 | >=15% | yes |

`localized_indeterminate` is **not** a formal positive tier.

Contig analysis runs on **`broad_island` and `extensive_island` only** (`--tumor_contig_tiers`,
default `broad_island,extensive_island`). `multi_gene` samples are classified and reported but
not assembled: the depth-matched control sweep showed targeted assembly is underpowered at that
depth, so an assembly failure there would reflect power, not biology. Widen with
`--tumor_contig_tiers multi_gene,broad_island,extensive_island` when that is the question.

### Changed from v3-working

v3-working used a two-tier `moderate`/`strong` scheme whose `moderate` boundary was
**11 reads / 3 genes / 1%** — that is candidate `multi_gene_B`, which the threshold analysis
evaluated and **rejected** in favour of 5 reads. Its `broad_island` tier (30/8/10%) had no
representation at all. v0.0.2 replaces both with the frozen three-tier system:

- `moderate` (11/3/1%) -> `multi_gene` (**5**/3/1%)
- *(absent)* -> `broad_island` (30/8/10%, as v0.0.2 first introduced it)
- `strong` (100/10/20%) -> `extensive_island` (100/10/20%, as v0.0.2 first introduced it --
  matching `strong`'s percentage, not carried over unchanged from it)

**Superseded.** Both breadth percentages above were later recalibrated -- `broad_island`
10% -> **7.5%**, `extensive_island` 20% -> **15%** -- for a change of measurement axis: v0.0.1
measured breadth from a 50 bp-binned coverage track, v0.0.2 measures per base, which reads a
median 0.794x lower on the same data (see `CHANGELOG.md`). The **Evidence tiers** table above
carries the current, frozen values (7.5%/15%); this section is a historical note on where the
read-count and gene-count thresholds came from, not a second source for the breadth values.

Thresholds are validated as non-decreasing across tiers, since the highest-tier-down cascade is
otherwise unreachable.

### Bug fixed in the port

`tumorEligibilityStatus` and `tumorContigContext` in `pks_mag.nf` reference
`params.tumor_min_clb_genes`, `params.tumor_min_pks_reads` and `params.tumor_min_nonhost_reads`,
which were **defined nowhere** in v3-working — a null-pointer failure latent only because
`--tumor_full_contig_context` defaults to false. They are now defined and aligned with the
lowest positive tier (3 genes, 5 reads); the non-host read floor defaults to 0 (disabled),
because no analysis supports a specific value. The tier classifier is the authoritative gate.

Evidence-file parsing was also hardened: the gate now resolves `read_evidence` by column name
and fails explicitly on an empty or malformed file, instead of indexing line 2 field 2 blindly.

## Taxonomic prefilter and clb homology rescue

Ported, **off by default** (`--prefilter_mode off | balanced`). With `off`, host-depleted reads
go straight to profiling exactly as in v0.0.1.

With `balanced`:

1. `krakenPrefilter` — `krakenuniq` against `--kraken_db`, then `route_kraken_reads.py` splits
   reads by `--primary_taxid` (default **91347, Enterobacterales** — the order, not the
   Enterobacteriaceae family) into `primary` and `non_target`.
2. `buildClbDiamondDb` → `diamondRescue` — `diamond makedb` over
   `ref/clb_reference_proteins.faa` (all 19 clb proteins, clbS included), then
   **`diamond blastx`** across the *non-target* reads at `--diamond_evalue 1e-5`,
   `--diamond_min_aa 25`, `--diamond_min_query_coverage 0.5`. Translated search, so divergent
   homologs still hit; there is no `blastn` step anywhere in the pipeline.
3. `mergePksCandidates` — rescued reads rejoin `primary` to become `PROFILING_READS`.

`--diamond_rescue false` with `balanced` is **refused**: without the rescue the mode degenerates
into a blunt taxonomic filter that discards clb-carrying reads placed outside the target taxon.

The prefilter narrows **only what gets profiled**. MAG assembly and the tumour re-join both
draw on `MAPPED_READS`, so neither is affected.

### Two taxonomy questions, two flags

| Flag | Input | Question answered |
|---|---|---|
| `--pks_taxa` | reads extracted from the island (`extractPksIslandReads` → `Bracken`) | which taxa do the island reads come from? |
| `--pks_community_taxa` | the whole non-host krakenuniq report (`sampleBracken`) | what is in the sample overall? |

`--pks_taxa` is v0.0.1's behaviour, unchanged, and works with the prefilter on or off.
`--pks_community_taxa` is new and requires `--prefilter_mode balanced`, since that is what
produces the krakenuniq report. v3-working had replaced the former with the latter; v0.0.2
keeps both.

## Divergent-read rescue: 844-protein DIAMOND database

The prefilter's DIAMOND database was built from `ref/clb_reference_proteins.faa` — **19 single
reference proteins**, one per gene. It now builds from
`ref/clb_reference_proteins_exact_v1.faa`: the **844 unique exact_v1 training proteins**, all 19
genes represented (clbB 166, clbJ 124, clbK 90 …). Same machinery, materially better recall on
divergent reads. The 19-sequence file is retained and selectable via `--clb_protein_fasta`.

`select_diamond_reads.py` discards `sseqid` and filters purely on e-value, alignment length and
query coverage, so the different header shape needs no parsing change.

**Read-level profiling stays on `nhmmscan` against the nucleotide profile.** Protein profiles
are more sensitive to divergent sequence, but HMMER has no protein-HMM-vs-nucleotide search, and
translating 150 bp reads gives 50 aa fragments scored against profiles up to 3,206 aa — with no
frameshift tolerance and no calibrated per-domain cutoffs. `nhmmscan` is fragment-native. The
protein profiles already do the divergent-sequence work where it is statistically sound:
`hmmsearchClb` on Prokka proteins from MAG bins, and `hmmsearchTumorContigs` on tumour contigs.

## Prophage, nearby-gene and community context

### Bug fixed: tumour prophage calls were computed and discarded

`tumorContigContext` declared `path(virusSummary)` as an input and staged the file, but its
script never referenced it — `genomadTumorContigs` ran geNomad on tumour contigs and the result
reached the work directory and then nothing. A clb island sitting inside a predicted provirus on
a tumour contig was not flagged as such. The metagenome path integrated prophage correctly all
along, via `communityProphageSummary`; the tumour path was the asymmetric one.

`extract_genomic_context.py` now takes `--genomad` and emits four columns per clb hit:
`in_prophage`, `prophage_id`, `prophage_virus_score`, `prophage_taxonomy`. Provirus intervals
are parsed with the same conventions as `build_community_prophage.py` — `topology == provirus`,
interval from `coordinates`, host contig by stripping the `|provirus_...` suffix from `seq_name`.

It also emits **`nearby_trna`**, which the tumour mobility table previously lacked while the
metagenome one had it. This matters specifically: the canonical pks island integrates at a tRNA
locus, so a tRNA in the flank is a targeted signal, not generic annotation. tRNA features need a
separate GFF pass because `parse_prokka_gff` keeps CDS only.

The tumour mobility table now carries all four new columns alongside `has_integrase` and
`has_transposase`.

### Context window is now a parameter

`--context_window_bp` (default 50000) replaces the hardcoded script default, and is wired into
both lanes. Previously the flanking window could not be changed from the command line.

### What the metagenome community analysis produces

Unchanged, and worth stating because it is the reason assembly must never be prefiltered.
`communityProphageSummary` runs over **all** bins, not only pks+ ones, and emits:

- `community_prophage_inventory.tsv` — per prophage: host taxonomy, coordinates, virus score,
  viral taxonomy, `has_recA`, `has_lexA`, `clbS_like`, evidence level.
- `pks_community_interactions.tsv` — producer bin x recipient bin, with the recipient's prophage
  count, recA, lexA and clbS_like, a susceptibility hypothesis and an interpretation.
- `pks_island_mobility.tsv` — distinct clb genes, provisional class, biosynthetic clb detected,
  clbP/clbS, contig span, `nearby_integrase`, `nearby_trna`, mobility interpretation.

recA and lexA are the SOS switch: colibactin damages DNA, SOS induces, prophage excises — so
recording them in the *recipient* is what makes the interaction table a mechanistic hypothesis
rather than a co-occurrence list. `clbS_like` in a neighbour marks self-resistance.

**Narrowing to Enterobacterales before assembly would delete the non-Enterobacterales bins that
`pks_community_interactions` is entirely about.** That is the strongest reason the assembly
stream is never prefiltered.

## Prefilter: sample-type-aware default

`--prefilter_mode` gains **`auto`**, now the default:

| sample_type | `--kraken_db` given | resolves to |
|---|---|---|
| `metagenome` | yes | `balanced` |
| `metagenome` | no | `off` |
| anything else | either | `off` |

Metagenomes carry far more non-target background, so the narrowing pays there; tumour WGS keeps
v0.0.1 behaviour unless asked. `auto` never turns a previously-working command into a failure,
because it only selects `balanced` when the required database is present. The resolution is
logged. An explicit `off` or `balanced` always wins.

**The prefilter has never been gated on sample type** and still is not — it narrows
`PROFILING_READS` for any sample type. What is exempt is assembly, in both lanes.

### The read-stream invariant

| Consumer | Stream |
|---|---|
| metagenome MAG assembly | `MAPPED_READS` — complete |
| tumour targeted assembly | `MAPPED_READS` — complete |
| alignment profiling | `PROFILING_READS` |
| HMM profiling | `PROFILING_READS` |

Profiling may narrow; assembly never may. The tumour lane decides eligibility from profiled
reads and then assembles from the complete stream. `tests/test_read_stream_invariant.py` guards
this, including that the MAG branch stays *positioned* above the prefilter block — moving it
below would silently inherit the narrowed stream while the code still read `MAPPED_READS`.

With `balanced` active alongside any assembly path, a `log.warn` states that read-level and
assembly-level clb evidence are not computed on the same reads, since the two numbers land in
the same results directory and can legitimately disagree.

### Read evidence is now computed for metagenomes too

`classifyTumorPksEvidence` is renamed **`classifyPksReadEvidence`** — it was never
tumour-specific — and publishes to `pks_summary/read_evidence/<sample>/` rather than under
`targeted_assembly/`. It runs for `tumor_wgs` and `metagenome`.

For metagenomes the tier is an **annotation only**. `pksMAG` is not tier-conditioned, and a test
asserts it. The thresholds were derived from TCGA tumour breadth distributions and are not
calibrated for metagenome read depth; gating genome-resolved analysis on a reference-based
nucleotide tier would filter out exactly the divergent carriers that lane exists to find. The
contig-count gate remains the metagenome gate.

## Strain typing: ST, clonal complex, phylogroup

New `Modules/pks_typing.nf`, on by default (`--enable_strain_typing`). Runs on **assembled units
only** — MAG bins individually, and the tumour contig set as one unit. Reads are never typed.

```
assembly -> mlst (ecoli_achtman_4) -> ST
         -> assign_strain_type.py  -> clonal complex, phylogroup, panel context
         -> strainTypeSummary      -> one table per sample
```

Output under `results/strain_typing/<sample>/`.

**The reference panel** is `ref/typing/st_phylogroup.tsv`, derived from
`pksProfiler_reference_build/cohorts/genome_roles.tsv` — 45,761 *E. coli* genomes reduced to
2,614 STs, each with clonal complex, majority phylogroup, a phylogroup-consistency figure, genome
count, and the fraction of that ST which is pks-positive.

**It refuses to guess.** Below `--typing_min_assembly_bp` (500 kb) MLST is not attempted; a
profile with fewer than seven loci is reported `insufficient_loci` with no ST, not a partial
call. Seven loci with no matching ST is `novel_allele_combination`.

**Panel-informed flags.** pks+ genomes in the panel are 92% phylogroup **B2** (4,386/4,760),
concentrated in ST73 (1,618) and ST95 (823). A typed non-B2 assignment, an ST with no pks+
genome in the panel, or an ST whose phylogroup is inconsistent (<0.9) each raise a note rather
than passing silently.

`tests/test_strain_typing.py` covers all of it, including that a partial profile produces no ST
and that typing is wired to assemblies rather than reads.

## Not ported from v3-working

Nothing outstanding — `pks_prefilter`, `pks_mag` and `pks_targeted` are all in.

## Coordinates

v0.0.1's island coordinates are retained (`pks_shift = 2193827`, `pks_island_len = 50767`).
v3-working had shifted to a 0-based `2193826 / 50768`. Both are internally consistent — the
classifier derives its breadth denominator from these same params — and the 1 bp difference is
0.002% of island length, far below any tier breadth threshold. Retaining v0.0.1's values avoids
changing released behaviour.

## Databases required

`--enable_mags` and `--tumor_enable_mags` each require `--gtdbtk_db`, `--checkm2_db` **and**
`--genomad_db`; `pksMAG` runs all three unconditionally. `--tumor_full_contig_context` requires
`--genomad_db`. `scripts/download_mag_databases.sbatch` and
`scripts/download_genomad_database.sbatch` fetch them.

## Profile HMMs: the exact_v1 population build

Both shipped profile sets are now the population-scale **`exact_v1`** build from
`pksProfiler_reference_build`, replacing the v3-working protein model whose provenance
matched no known build.

| File | Used by | Tool | Indexes shipped |
|---|---|---|---|
| `ref/hmm/clb_population_protein_exact_v1.hmm` | `hmmsearchClb`, tumour contig HMM search | `hmmsearch` | no — not needed |
| `ref/hmm/clb_population_dna_exact_v1.hmm` | `pksProfilerHMM` | `nhmmscan` | yes — required |

SHA256 verified against `models/exact_v1/SHA256SUMS`. 29 MB added; the protein indexes are
omitted because `hmmsearch` reads the plain `.hmm`, saving ~8 MB over shipping everything.
A validation guard fails early if `--hmm_model` is set to an unpressed profile.

### Why exact_v1 rather than the original build

Both derive from the same 2,868 lineage-split training genomes and differ only in how
near-identical sequences were collapsed. `exact_v1` removes only **100%-identical** duplicates,
after excluding sequences outside 80–120% of each gene's modal length. The original clustered at
**99% identity** — far too aggressive on an island this conserved, where 99% for clbA (244 aa)
tolerates roughly two amino acids.

| | raw sequences | retained |
|---|---:|---:|
| exact_v1 | 52,865 | **844** |
| original (99% clustered) | 52,865 | **65** |

Protein training depth: exact_v1 median `NSEQ` 23, minimum 7, no single-sequence models. The
99% build had median 2 and built clbG, clbL, clbM, clbO and clbP from **one sequence** each;
its DNA counterpart had eight such genes. A profile trained on a single sequence models no
positional variation, which removes the reason to use an HMM. Exact deduplication still
suppresses abundant ST73/ST95 clones, since identical copies collapse either way.

`tests/test_hmm_models.py` pins the checksums, the 19 genes in both alphabets, and asserts no
gene has `NSEQ` below 2 — a regression guard against reintroducing the clustered build.

### v0.0.1 benchmark retained

`ref/hmm/clb_all_dna.hmm` and its indexes are untouched and still in the tree. The
reference-build README treats it as an immutable benchmark until sensitivity, specificity,
locus-reconstruction and MAG-fragment tests are complete. Reproduce v0.0.1 numbers with
`--hmm_model ref/hmm/clb_all_dna.hmm`.

### Still uncalibrated

Neither profile set carries GA/TC/NC cutoffs, so detection rests entirely on `--hmm_evalue`
(1e-10) and `--hmm_protein_evalue` (1e-5). The 964 validation and 928 test genomes remain
reserved and unused. Two audit-table anomalies to resolve before calibration: **clbJ** lost 687
sequences to the length filter, far more than any other gene, and **clbK** yielded only 1,731
raw sequences against ~2,840 for every other gene — roughly 1,100 training genomes produced no
clbK call at all. Full detail in `ref/hmm/PROVENANCE.md`.

## Interpretation limits

Unchanged from v0.0.1, and reinforced: the tier system was derived from **positive-only** cohorts
in which every sample carries at least one pks read. It ranks evidence strength and measures
retention of designated positives. It does **not** establish the boundary between a true positive
and background mapping, and provides no specificity or false-positive rate. NIH-designated
negatives or matched controls are required before any tier is used as a clinical or biological
call. Small nonzero breadth remains indeterminate, not a positive.

## Fixed after the initial port

`params.pks_reference_fasta` and `params.pks_recruit_index` were initially declared before
`params.pks_genome`, so the GString interpolated null and every tumour run failed with
`Canonical PKS reference not found: null.fna`. They are now declared immediately after
`params.pks_genome` and derive from it. Existence guards for the reference FASTA and the
Bowtie2 recruitment index run when the targeted path is enabled.

**Boolean flags were truthy when set to false.** Nextflow passes `--flag false` as the
*String* `"false"`, which is truthy in Groovy, so `--enable_mags false` would have enabled MAGs
and `--diamond_rescue false` was silently ignored. All six v0.0.2 booleans are now coerced with
`toString().toBoolean()` at the top of the workflow — before the MAG branch, which is the
earliest consumer — and only the coerced locals are tested. `tests/test_tumor_tiers.py` asserts
no flag is read directly after that point.

**Nextflow 26 strict-DSL notes.** The 26.04.6 parser rejects `instanceof` in a ternary condition
and will not call a local closure as a function, so the coercion is inlined per flag rather than
sharing a helper. Unit tests are text assertions and cannot catch a compilation failure — always
run an actual `nextflow run` after editing `main.nf`.

## Verification

113 unit tests pass (`python3 -m unittest discover -s tests`). Nextflow 26.04.6 parses `main.nf`
and all routing, tier, prefilter, taxonomy and database gates were confirmed to fire. `nextflow lint` does not exist
in the TSCC Nextflow build and is no longer invoked by `tests/run_checks.sh`.
No pipeline run has yet been executed against real data on this tree.
