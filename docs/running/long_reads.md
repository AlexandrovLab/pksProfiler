# Long-read input (ONT, PacBio HiFi)

`--read_type short | ont | hifi` (default `short`). Long reads are accepted for
`--sample_type wgs` and `rna` only; `metagenome` and `wes` refuse anything but `short` before
a task is submitted ("long-read input is supported for wgs and rna only").

What changes for `ont` / `hifi`:

| Step | Short reads | Long reads |
|---|---|---|
| Read QC | fastp (`-l 45 --cut_tail`) | **fastplong** (`-l --long_read_min_length`, default 200) |
| Human pass, wgs | `-ax sr` vs `--hg38_db` | `-ax map-ont` / `-ax map-hifi` vs `--hg38_db_long` |
| Human pass, rna | `-ax splice:sr` vs `--rna_splice_db` | `-ax splice` (ont) / `-ax splice:hq` (hifi) vs `--rna_splice_db_long`, same optional `-j --rna_junctions_bed` |
| T2T + phiX pass | `-ax sr` vs `--t2t_phix_db` | `-ax map-ont` / `-ax map-hifi` vs `--t2t_phix_db_long` |
| Pangenome passes | `-ax sr` | `-ax map-ont` / `-ax map-hifi` (build the pangenome `.mmi` files with that preset, or they keep their short-read k/w) |
| `--minimap2_frag_mode` | applied to every call | not applied: `--frag` only changes how minimap2 treats short paired reads |
| *clb* profiling | Bowtie2 vs the IHE3034 index | **minimap2 `-ax map-ont` / `-ax map-hifi`** vs the IHE3034 genome FASTA (`--pks_reference_fasta`), same MAPQ filter (`--pks_align_min_mapq`), featureCounts `-L` |
| WGS evidence | tiered | **untiered** (`read_evidence=untiered`, "WGS evidence (not tiered: long-read)") |
| RNA evidence | `rna_evidence.tsv` | the same schema |
| Targeted assembly | wgs | **refused** for wgs (pass `--tumor_targeted_assembly false`; it is on by default) |

Why separate `*_long` indexes: a minimap2 `.mmi` stores the k-mer and window size of the
preset it was built with, and those are not overridden by `-x` at mapping time. A short-read
`sr` index would silently map long reads with short-read seeds. Pass a FASTA (minimap2 then
indexes it with the right preset in each task) or a prebuilt index:

```bash
minimap2 -x map-ont  -d GRCh38.map-ont.mmi   GRCh38.fa      # --hg38_db_long, ont
minimap2 -x map-hifi -d GRCh38.map-hifi.mmi  GRCh38.fa      # --hg38_db_long, hifi
minimap2 -x splice   -d GRCh38.d1.vd1.splice.mmi GRCh38.d1.vd1.fa   # --rna_splice_db_long, ont
minimap2 -x splice:hq -d GRCh38.splice_hq.mmi GRCh38.fa     # --rna_splice_db_long, hifi
```

On TSCC, the ONT splice index is being built at
`/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/ref_indices/grch38_splice_long/GRCh38.d1.vd1.splice.mmi`.

QC names are unchanged so the QC table and its contracts hold: for long reads,
`filter_input_reads` and `reads_after_fastp` are fastplong's before/after read counts.

Why WGS long reads are untiered: the tier thresholds were calibrated on short-read Bowtie2
alignments (read counts of ~150 bp reads, breadth from their depth). A long read covers far
more of the island per read, so the same numbers would not mean the same thing. Calibrate
before tiering.

Long-read sample sheets must not set `fastq2`: fastplong takes one stream per sample.
Not covered yet: direct-RNA ONT libraries may want `-uf -k14`, and `--pks_taxa`'s Bracken
read-length model is a short-read concept.
