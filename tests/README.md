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

## The V8a and V8b lanes, run for real

`tests/check_v8_lanes_execution.sh` runs each of Ludmil's two 25 Sep section-4 checks once,
end to end, on fixtures cut from the repository's own IHE3034 reference
(`tests/fixtures/generate_v8_fixtures.py`) — so the answer is known by construction:

- **V8a, locus-detection specificity.** A sample mixing reads from a 20 kb island slice with
  reads from a 12 %-divergent NRPS/PKS "homolog" slice taken from a different part of the
  island. Island reads must come out `confirmed_true`, homolog reads `confirmed_decoy`, and
  the sample's retention must land between the two extremes — a lane that confirmed or
  rejected everything would otherwise pass. The lane is run without Arm A, so it must also
  report `arms_used=arm_b_only`, and it must refuse outright when no panel is given.
- **V8b, carrier attribution.** The same island as two carriers 3 % apart, one *E. coli* and
  one *Klebsiella*. Reads taken from the *Klebsiella* copy must be attributed to that
  carrier, and the sample call must name *Klebsiella*.

`tests/fixtures/assert_v8_outputs.py` also checks that neither lane's output carries the
other's vocabulary — Ludmil's "should not be conflated", enforced on real output rather than
on source. About two minutes with a warm conda cache:

```bash
CONDA_CACHE_DIR=/path/to/cache bash tests/check_v8_lanes_execution.sh
```

## Resume after a reference changes (F15, executed)

`tests/check_resume_reference_change.sh` is the one test that executes the F15 dependency-
digest mechanism rather than reading it as source. Three real runs on the e2e fixtures: a
clean run, a `-resume` with nothing changed (every task must be `CACHED`, or nothing below
means anything), and a `-resume` after the clb annotation's **content** changes. Every
process carrying the `pks_annotation` digest must then re-run, and the upstream processes
that do not carry it must still be cached — re-running everything would satisfy "it re-ran"
while proving the digest does nothing.

```bash
CONDA_CACHE_DIR=/path/to/cache bash tests/check_resume_reference_change.sh
```

Both scripts, plus `check_e2e_execution.sh`, `check_malformed_input_execution.sh` and
`check_mode_execution.sh`, now run in CI on every push to `main` and nightly (F21), with the
conda environments cached between runs and the outputs kept as build artefacts.
