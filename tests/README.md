# Lightweight regression checks

Run the repository checks with:

```bash
bash tests/run_checks.sh
```

The checks do not run bioinformatics workloads or require the large host and
taxonomy databases. They verify:

- Nextflow syntax for the main workflow and included modules;
- the 19-gene GFF and HMM reference contracts;
- presence of the complete bundled Bowtie2 index;
- mate-safe taxonomy joins;
- combined species-support table structure;
- positive and negative example-table structure; and
- rejection of duplicate sample identifiers before task execution.

These checks complement, rather than replace, end-to-end runs on representative
BAM, FASTQ, HMM, and taxonomy inputs.

## Real execution per mode

`tests/check_mode_execution.sh` (Ludmil, 26 Sep mode memo section 6) runs the pipeline for
real, once per mode, on tiny inputs built from public references
(`tests/fixtures/generate_mode_fixtures.py`: IHE3034 reads, a 30 kb GRCh38 slice with a
three-exon transcript at real GT/AG sites and its BED12, ONT-like long reads simulated from
both): `wgs`, `wes`, `rna` (`-ax splice:sr -j` against a splice:sr index built on the slice),
`metagenome`, `wgs --read_type ont` and `rna --read_type ont`, plus `wgs_assembly` (Ludmil,
25 Sep A1/A2): an extensive_island sample with flank-only, boundary-spanning and chimeric
(60 bp island + 90 bp random) pairs, so targeted recruitment and MEGAHIT/metaSPAdes run for
real and every recruitment class in `pks_recruitment.tsv` is counted.
`tests/fixtures/assert_mode_outputs.py` checks which evidence file each mode wrote, the
counts, the host-depletion preset the task really ran, that tiers/assembly/MAGs are off
where they should be, `terminal_state`, and the master summary. About 8-9 minutes with a
warm conda cache:

```bash
CONDA_CACHE_DIR=/path/to/cache bash tests/check_mode_execution.sh   # MODES="rna wgs_ont" to subset
```

The metagenome paired-assembly path (`--enable_mags`) needs the GTDB-Tk, CheckM2 and geNomad
databases and is exercised by `tests/check_metagenome_mag_execution.sh`.
