# Island reassembly (tumour data)

Counting reads tells you how much *clb* sequence is present. Reassembling those reads tells you
whether they form a real, connected island or scattered fragments that merely resemble one.

This page covers `--sample_type wgs`: evidence tiering and targeted reassembly. Prophage
and neighbouring-gene mobility context is a metagenome-only capability -- see
[MAGs](mags.md) and [community context](community_context.md).

## Evidence tiering

This runs automatically for tumour data. Set `--sample_type` and you also get a tier per sample.

```bash
nextflow run main.nf -profile local \
  --sample samples.csv \
  --input_data_type bam \
  --sample_type wgs \
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
sample      mode  read_evidence      pks_reads  clb_genes_detected  island_breadth_1x
TUMOR_01    wgs   extensive_island   1497       15                  0.731800
```

Because you set `--sample_type wgs`, samples that clear the gate also go straight into
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
  --sample_type wgs \
  --tumor_targeted_assembly true \
  --hg38_db      /refs/human-GRC-db.mmi \
  --t2t_phix_db  /refs/human-GCA-phix-db.mmi \
  --outdir       results
```

**What you get**, per qualifying sample:

```text
results/by_sample/<sample>/contigs/
├── evidence/          the tier and the numbers behind it
├── recruitment/       the reads pulled out for assembly, and pks_recruitment.tsv (classes below)
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

### Which reads are recruited

Recruitment aligns the sample's non-host reads (Bowtie2 `--very-sensitive-local`, ordinary
best-alignment reporting) against the canonical IHE3034 island plus 5 kb of flank either side
(`indices/pks_island/`, whose FASTA header records `island=2193827-2244594 flank=5000`). The
island therefore sits at bases 5,001-55,768 of that reference. Every read must first reach
MAPQ >= `--recruit_min_mapq`, `--pks_recruit_min_aligned` aligned bases and
`--pks_recruit_min_identity`; it is then put in one of three classes (Ludmil, 25 Sep A1/A2):

| Class | Rule |
|---|---|
| `island_internal` | inside the island, and query coverage (aligned read bases / read length) >= `--recruit_min_query_coverage` |
| `boundary_spanning` | crosses an island boundary by more than `--recruit_boundary_tolerance_bp` on each side, **or** lies in the island with its outer end within the tolerance of a boundary and a soft clip longer than the tolerance pointing out of the island. No query-coverage floor: a real junction read is clipped where the host locus differs from IHE3034 |
| `flank_only` | in a flank (touching the island by at most the tolerance), query coverage >= `--recruit_min_query_coverage`. A flank read clipped *into* the island is flank-only: that is what a pks-negative relative sharing the flank looks like |
| rejected | everything else, e.g. a 150 bp read with only a 60 bp local match (coverage 0.40) |

A fragment takes the best class of its mates (island_internal > boundary_spanning > flank_only).
`island_supporting_fragments` = island_internal + boundary_spanning. The ±5 kb flanks can be
shared by pks-negative bacteria, so **only island-supporting fragments count as island support**
and drive `assembly_input_diagnostic` in `final_pks_evidence.tsv`:

| `assembly_input_diagnostic` | Meaning |
|---|---|
| `no_island_supporting_fragments` | nothing island-derived reached the assemblers (flank-only fragments may have) |
| `island_fragments_no_contig` | island-supporting fragments went in, no island contig came out |
| `island_contigs_recovered` | at least one contig covers island sequence |
| `recruitment_classes_unavailable` | recruitment stats written before the classes existed |

Flank-only fragments still go to the assemblers for context unless
`--recruit_include_flank_only false`. `recruitment/pks_recruitment.tsv` reports reads and
fragments per class (`reads_island_internal`, `fragments_boundary_spanning`, ...), the rejection
reasons (`reads_rejected_low_mapq`, `_short_alignment`, `_low_identity`, `_low_query_coverage`),
`island_supporting_fragments`, and every threshold used.

| Parameter | Default | |
|---|---|---|
| `--pks_recruit_min_aligned` | 60 | aligned bases (M/=/X) a recruitment alignment needs |
| `--pks_recruit_min_identity` | 0.90 | (aligned - NM) / aligned |
| `--recruit_min_query_coverage` | 0.80 | island_internal and flank_only reads |
| `--recruit_boundary_tolerance_bp` | 10 | how close to a boundary counts as reaching it; must be < half of `--pks_recruit_min_aligned` |
| `--recruit_min_mapq` | 20 | every recruited read |
| `--recruit_include_flank_only` | true | pass flank-only fragments to the assemblers |
| `--pks_recruit_flank_bp` | 5000 | flank the recruitment index was built with; checked against its FASTA header |

These change only what goes into structural confirmation; the primary evidence tier and
positivity are untouched.

**Only `broad_island` and `extensive_island` samples are assembled.** Below that there are too few
reads for assembly to succeed even when an island is genuinely present, so a failure would be
telling you about read depth, not biology. To widen it anyway:

```bash
  --tumor_contig_tiers multi_gene,broad_island,extensive_island
```

---

Back to the [README](../../README.md).
