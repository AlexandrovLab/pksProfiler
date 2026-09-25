#!/usr/bin/env bash
# Run one real, non-preview Nextflow execution end to end and check its real output.
#
# Why this exists. Ludmil's audit, "not yet confirmed" item u-3: every existing test is
# either (a) a Python unit test that pattern-matches .nf/.py source as text, or (b)
# `nextflow run -preview` (tests/check_workflow_construction.sh), which wires the DSL2
# channel graph and runs zero tasks. Neither has ever run bowtie2, minimap2, fastp or
# samtools on real data and looked at the output -- which is exactly how F01 got through
# both tiers at once: samtools fastq silently dropped a category of reads because of
# which routing flags were and weren't passed (see tests/test_extract_read_routing.py's
# docstring), and nothing here would have noticed until a real cohort run did.
#
# This script:
#   1. builds tiny, real fixtures (tests/fixtures/generate_e2e_fixtures.py) -- paired
#      FASTQ reads sliced out of the repo's already-shipped GCF_000025745.1 reference at
#      four clbA-clbS genes ("positive"), paired FASTQ reads of unrelated pseudorandom
#      sequence ("negative"), and two small synthetic host references;
#   2. builds REAL minimap2 .mmi indices from those synthetic host references
#      (minimap2 -d), so --hg38_db/--t2t_phix_db point at real indices, not stand-ins;
#   3. runs main.nf for real (no -preview) on both samples, with --sample_type tumor_wgs
#      --profiling_method bowtie2, through a real conda environment;
#   4. asserts on the real output values (tests/fixtures/assert_e2e_outputs.py) --
#      the read-evidence tier, the featureCounts gene matrix, the cohort QC counts, and
#      -- Ludmil's u-1 ask -- that BOTH mates of every fragment survive extraction and
#      host depletion, checked against the real host-depleted FASTQ a real run produced;
#   5. packages the same two read sets as BAM (aligned against the synthetic host
#      reference, so every read comes out unmapped with a real @SQ header -- see the
#      comment above the BAM-building loop below for why) and runs the pipeline again
#      with --input_data_type bam, so extractReads' own `samtools fastq` invocation --
#      the exact process F01 is about -- is exercised for real too, not just the
#      paired-FASTQ path that skips it entirely.
#   6. converts the same two BAMs to CRAM (`samtools view -C -T`, same synthetic hg38
#      reference) and runs the pipeline a third time with --input_data_type cram
#      --cram_reference, so extractReads' CRAM branch (--cram_reference required,
#      validate_cram_reference.py's @SQ/MD5 check, the `samtools fastq --reference`
#      decode) runs for real too -- not just BAM's referenceless decode path.
#   7. runs the pipeline a fourth time on the original paired-FASTQ sheet with
#      --sample_type metagenome --profiling_method hmm (tumor_wgs rejects hmm-alone;
#      see main.nf), so pksProfilerHMM/hmm_best_hit.py -- shipped since v0.0.1 but never
#      executed by any test in this repo -- run for real on both the positive fixture
#      (real nhmmscan hits, not just the bowtie2 lane's) and the negative one (a real
#      zero-hit nhmmscan run, a different code path than "no reads reached the process
#      at all").
#   8. runs the pipeline a fifth time on a 5-sample cohort (positive, negative, plus
#      three new fixtures spanning broad_island, extensive_island and
#      localized_indeterminate -- see BROAD_WINDOWS/EXTENSIVE_WINDOWS/
#      BORDERLINE_WINDOWS in generate_e2e_fixtures.py), so the cohort-level reducers
#      (masterTableAlign, masterQCSummary, cohortReport, masterSummary) aggregate more
#      than two rows at least once, with membership checked exactly (every sample
#      appears, none duplicated) rather than just "both of the original two are in
#      there somewhere".
#   9. runs the pipeline a sixth time on single-end FASTQ (fastq1 only, no fastq2) --
#      every run above is paired-end. No new fixtures: one FASTQ half each of the
#      existing broad/negative fixtures is already a real, independently-validated
#      single-end read set (see the comment above this run below for the exact
#      numbers, which differ from their paired-end namesakes since single-end halves
#      the bases covered).
#
# This is slower than the other two test tiers (real conda envs, real alignment) but
# still small: at most ~110 reads per sample, tiny synthetic references, no
# assembly/MAG/taxonomy lane (the default flags don't trigger any of those for any of
# these fixtures -- see the comment above LOCI in generate_e2e_fixtures.py).
#
#     CONDA_CACHE_DIR=/path/to/cache bash tests/check_e2e_execution.sh
#
# CONDA_CACHE_DIR defaults to the lab's shared cache, which already has every conda env
# this run needs pre-built (samtools/fastp/minimap2/bowtie2/subread/deeptools/R) --
# Nextflow's conda envs are hashed from file CONTENT, not path, and conda_envs/*.yml is
# byte-identical to the tree that cache was built from, so pointing at it reuses those
# environments instead of re-solving them. Rebuilding from scratch on first use is a
# few minutes per env; consider that budget before ignoring this variable.

set -uo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONDA_CACHE_DIR=${CONDA_CACHE_DIR:-/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache}

if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: nextflow is required for e2e execution checks." >&2
    exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is required for e2e execution checks." >&2
    exit 1
fi

# minimap2 is needed once, up front, to build the two tiny synthetic host indices --
# not by Nextflow, which builds its own copy inside the pipeline's minimap2 conda env.
# Preference order: an explicit override, then whatever the shared conda cache already
# built for conda_envs/minimap2_env.yml (exact version match with the run below), then
# anything on PATH.
resolve_minimap2() {
    if [[ -n "${MINIMAP2_BIN:-}" ]]; then
        echo "$MINIMAP2_BIN"; return 0
    fi
    if [[ -d "$CONDA_CACHE_DIR" ]]; then
        local hit
        hit=$(grep -Fls "conda_envs/minimap2_env.yml" "$CONDA_CACHE_DIR"/env-*/conda-meta/history 2>/dev/null | head -1)
        if [[ -n "$hit" ]]; then
            local env_bin="${hit%/conda-meta/history}/bin/minimap2"
            if [[ -x "$env_bin" ]]; then
                echo "$env_bin"; return 0
            fi
        fi
    fi
    command -v minimap2 || true
}

MINIMAP2_BIN=$(resolve_minimap2)
if [[ -z "$MINIMAP2_BIN" ]]; then
    echo "ERROR: no minimap2 found. Set MINIMAP2_BIN, build the minimap2_env conda" >&2
    echo "       env under CONDA_CACHE_DIR, or put minimap2 on PATH." >&2
    exit 1
fi

# conda_envs/minimap2_env.yml also installs samtools (needed to sort/index the BAM
# fixture below); prefer the copy right beside the resolved minimap2 over a second,
# possibly different, lookup.
MINIMAP2_ENV_SAMTOOLS="$(dirname "$MINIMAP2_BIN")/samtools"
if [[ ! -x "$MINIMAP2_ENV_SAMTOOLS" ]]; then
    MINIMAP2_ENV_SAMTOOLS=$(command -v samtools || true)
fi
if [[ -z "$MINIMAP2_ENV_SAMTOOLS" ]]; then
    echo "ERROR: no samtools found (needed to build the BAM fixtures). Put samtools on PATH." >&2
    exit 1
fi

work=$(mktemp -d "${TMPDIR:-/tmp}/pks-e2e.XXXXXX")
echo "work dir: $work (kept on failure for inspection)"
cleanup() { [[ "${KEEP_WORK:-0}" == "1" ]] || rm -rf "$work"; }
trap cleanup EXIT

fixtures="$work/fixtures"
mkdir -p "$fixtures"

echo "Building fixtures (synthetic reads + synthetic host references)"
if ! python3 "$repo_dir/tests/fixtures/generate_e2e_fixtures.py" --out-dir "$fixtures" > "$work/fixtures.log" 2>&1; then
    cat "$work/fixtures.log" >&2
    echo "ERROR: fixture generation failed." >&2
    exit 1
fi
cat "$work/fixtures.log"

echo "Building real minimap2 indices for the synthetic host references ($MINIMAP2_BIN)"
if ! "$MINIMAP2_BIN" -d "$fixtures/hg38.mmi" "$fixtures/synthetic_hg38.fa" > "$work/mm2_hg38.log" 2>&1; then
    cat "$work/mm2_hg38.log" >&2
    echo "ERROR: could not build the synthetic hg38 index." >&2
    exit 1
fi
if ! "$MINIMAP2_BIN" -d "$fixtures/t2t.mmi" "$fixtures/synthetic_t2t_phix.fa" > "$work/mm2_t2t.log" 2>&1; then
    cat "$work/mm2_t2t.log" >&2
    echo "ERROR: could not build the synthetic t2t+phix index." >&2
    exit 1
fi

POS_SAMPLE=e2e_positive
NEG_SAMPLE=e2e_negative

printf 'patient,fastq1,fastq2\n' > "$work/sheet.csv"
printf '%s,%s,%s\n' "$POS_SAMPLE" "$fixtures/positive_R1.fastq.gz" "$fixtures/positive_R2.fastq.gz" >> "$work/sheet.csv"
printf '%s,%s,%s\n' "$NEG_SAMPLE" "$fixtures/negative_R1.fastq.gz" "$fixtures/negative_R2.fastq.gz" >> "$work/sheet.csv"

# Per-process cpu/memory floors in conf/base.config were sized for real cohorts (up to
# 16 cpus/task, 64 GB/task) and are irrelevant to an 8-read fixture. withLabel selectors
# take precedence over a plain `process {}` block regardless of file order (Nextflow
# config precedence: process < withLabel < withName), so a generic override here would
# silently lose to conf/base.config's per-label blocks; matching labels one for one is
# what actually takes effect. `executor.cpus`/`executor.memory` are set high enough that
# Nextflow's local executor -- sized from the JVM's available-processors/available-memory
# counts, which some sandboxed shells under-report regardless of the node's real
# capacity -- never refuses to admit a task that asks for more than this.
#
# targeted_recruit/targeted_assembly are u-3 follow-up additions: the original
# positive/negative fixtures were both outside --tumor_contig_tiers (default
# "broad_island,extensive_island"), so targetedPksAssembly (bowtie2 recruit + MEGAHIT,
# --tumor_targeted_assembly defaults true) never ran and its real memory floor --
# conf/base.config's 64 GB for targeted_assembly -- was never exercised by this suite.
# The broad_island/extensive_island cohort fixtures below land inside that tier on
# purpose (that's the whole point of scaling the cohort), so they do trigger it now,
# and 64 GB > this config's 32 GB executor cap failed the run outright before these two
# labels were added.
cat > "$work/e2e_local.config" <<'EOF'
executor {
    cpus   = 32
    memory = '32 GB'
}
process {
    withLabel:extract_reads     { cpus = 1 }
    withLabel:filter_reads      { cpus = 1 }
    withLabel:map_reads         { cpus = 1 }
    withLabel:pks_align         { cpus = 1 }
    withLabel:pks_hmm           { cpus = 1 }
    withLabel:process_low       { cpus = 1 }
    withLabel:sample_stage      { cpus = 1 }
    withLabel:targeted_recruit  { cpus = 1; memory = 1.GB }
    withLabel:targeted_assembly { cpus = 1; memory = 1.GB }
}
EOF

echo "Running the real pipeline (a real execution, not a channel-graph-only preview): tumor_wgs, paired FASTQ, bowtie2"
nf_log="$work/nextflow.log"
(
    cd "$work" && nextflow run "$repo_dir/main.nf" \
        -work-dir "$work/work" \
        -c "$work/e2e_local.config" \
        --sample "$work/sheet.csv" \
        --input_data_type fastq \
        --sample_type tumor_wgs \
        --profiling_method bowtie2 \
        --hg38_db "$fixtures/hg38.mmi" \
        --t2t_phix_db "$fixtures/t2t.mmi" \
        --save_intermediates true \
        --conda_cache_dir "$CONDA_CACHE_DIR" \
        --outdir "$work/results"
) > "$nf_log" 2>&1
nf_status=$?

if [[ $nf_status -ne 0 ]]; then
    echo "ERROR: the real pipeline run failed (exit $nf_status). Last 60 lines:" >&2
    tail -60 "$nf_log" >&2
    echo "Full log: $nf_log" >&2
    KEEP_WORK=1
    exit 1
fi

echo "Pipeline finished. Checking real output values."
if ! python3 "$repo_dir/tests/fixtures/assert_e2e_outputs.py" \
        --results "$work/results" \
        --positive-sample "$POS_SAMPLE" \
        --negative-sample "$NEG_SAMPLE"; then
    echo "ERROR: e2e output assertions failed. Results kept at: $work/results" >&2
    KEEP_WORK=1
    exit 1
fi

# ---------------------------------------------------------------------------------
# Second input form: the same two read sets, packaged as BAM instead of FASTQ, through
# extractReads (Modules/extract_reads.nf). This is the one process whose exact
# `samtools fastq` invocation F01 is about (see tests/test_extract_read_routing.py's
# docstring) -- the FASTQ run above never calls it at all, so it is the one part of
# Ludmil's u-3 wish list (BAM/CRAM/single FASTQ/paired FASTQ/HMM/cohort aggregation)
# that a paired-FASTQ-only run cannot cover by construction. It also exercises a
# different mate-suffix code path than the FASTQ run above: extractReads' `samtools
# fastq -N` derives /1 and /2 from each record's FLAG bits, where filterReads' own
# awk-based tagging (for direct FASTQ input) reads it from the QNAME instead.
#
# The BAM is built by aligning the reads against the synthetic hg38 reference with
# minimap2, exactly as mapReads does for host depletion -- not `samtools import`,
# which writes an unaligned BAM with no @SQ header at all ("had no targets in
# header", which is what extractReads' own `samtools quickcheck` correctly refuses).
# A real BAM/CRAM input to this pipeline is always a whole-genome alignment with a
# real reference dictionary; a header-less synthetic stand-in would not be
# representative of that even if quickcheck were relaxed to accept it. The reads
# are E. coli sequence aligned against unrelated synthetic "human" sequence, so
# every one comes out unmapped (which is exactly the read category extractReads
# keeps) with a real reference header attached.
# ---------------------------------------------------------------------------------
echo
echo "Building BAM fixtures (aligned-and-unmapped, like a real WGS BAM's non-host reads)"
for name in positive negative; do
    if ! "$MINIMAP2_BIN" -ax sr "$fixtures/hg38.mmi" \
            "$fixtures/${name}_R1.fastq.gz" "$fixtures/${name}_R2.fastq.gz" 2> "$work/mm2_bam_${name}.log" \
            | "$MINIMAP2_ENV_SAMTOOLS" sort -o "$fixtures/${name}.bam" -; then
        cat "$work/mm2_bam_${name}.log" >&2
        echo "ERROR: could not build the ${name} BAM fixture." >&2
        exit 1
    fi
    "$MINIMAP2_ENV_SAMTOOLS" index "$fixtures/${name}.bam"
done

POS_BAM_SAMPLE=e2e_bam_positive
NEG_BAM_SAMPLE=e2e_bam_negative

printf 'patient,alignment\n' > "$work/bam_sheet.csv"
printf '%s,%s\n' "$POS_BAM_SAMPLE" "$fixtures/positive.bam" >> "$work/bam_sheet.csv"
printf '%s,%s\n' "$NEG_BAM_SAMPLE" "$fixtures/negative.bam" >> "$work/bam_sheet.csv"

echo "Running the real pipeline on BAM input: tumor_wgs, bowtie2"
nf_bam_log="$work/nextflow_bam.log"
(
    cd "$work" && nextflow run "$repo_dir/main.nf" \
        -work-dir "$work/work_bam" \
        -c "$work/e2e_local.config" \
        --sample "$work/bam_sheet.csv" \
        --input_data_type bam \
        --sample_type tumor_wgs \
        --profiling_method bowtie2 \
        --hg38_db "$fixtures/hg38.mmi" \
        --t2t_phix_db "$fixtures/t2t.mmi" \
        --save_intermediates true \
        --conda_cache_dir "$CONDA_CACHE_DIR" \
        --outdir "$work/results_bam"
) > "$nf_bam_log" 2>&1
nf_bam_status=$?

if [[ $nf_bam_status -ne 0 ]]; then
    echo "ERROR: the real BAM-input pipeline run failed (exit $nf_bam_status). Last 60 lines:" >&2
    tail -60 "$nf_bam_log" >&2
    echo "Full log: $nf_bam_log" >&2
    KEEP_WORK=1
    exit 1
fi

echo "BAM-input pipeline finished. Checking real output values."
if ! python3 "$repo_dir/tests/fixtures/assert_e2e_outputs.py" \
        --results "$work/results_bam" \
        --positive-sample "$POS_BAM_SAMPLE" \
        --negative-sample "$NEG_BAM_SAMPLE" \
        --expect-extracted-unmapped-reads 8; then
    echo "ERROR: BAM-input e2e output assertions failed. Results kept at: $work/results_bam" >&2
    KEEP_WORK=1
    exit 1
fi

# ---------------------------------------------------------------------------------
# Third input form: CRAM. Same two BAMs, converted to CRAM against the same synthetic
# hg38 reference used to build them (`samtools view -C -T`) -- so the CRAM's @SQ M5s
# are computed from, and therefore match, the exact FASTA passed as --cram_reference,
# the same way a real CRAM and its reference dictionary agree. extractReads.nf refuses
# CRAM input without --cram_reference (Ludmil, revised report finding 10) and then runs
# scripts/validate_cram_reference.py, which would fail loudly on any accidental
# mismatch here -- so a passing run is a genuine confirmation the reference contract
# holds, not just that decode succeeded.
# ---------------------------------------------------------------------------------
echo
echo "Building CRAM fixtures (same alignments as the BAM fixtures, same synthetic hg38 reference)"
for name in positive negative; do
    if ! "$MINIMAP2_ENV_SAMTOOLS" view -C -T "$fixtures/synthetic_hg38.fa" \
            -o "$fixtures/${name}.cram" "$fixtures/${name}.bam" 2> "$work/cram_${name}.log"; then
        cat "$work/cram_${name}.log" >&2
        echo "ERROR: could not build the ${name} CRAM fixture." >&2
        exit 1
    fi
    "$MINIMAP2_ENV_SAMTOOLS" index "$fixtures/${name}.cram"
done

POS_CRAM_SAMPLE=e2e_cram_positive
NEG_CRAM_SAMPLE=e2e_cram_negative

printf 'patient,alignment\n' > "$work/cram_sheet.csv"
printf '%s,%s\n' "$POS_CRAM_SAMPLE" "$fixtures/positive.cram" >> "$work/cram_sheet.csv"
printf '%s,%s\n' "$NEG_CRAM_SAMPLE" "$fixtures/negative.cram" >> "$work/cram_sheet.csv"

echo "Running the real pipeline on CRAM input: tumor_wgs, bowtie2, --cram_reference"
nf_cram_log="$work/nextflow_cram.log"
(
    cd "$work" && nextflow run "$repo_dir/main.nf" \
        -work-dir "$work/work_cram" \
        -c "$work/e2e_local.config" \
        --sample "$work/cram_sheet.csv" \
        --input_data_type cram \
        --cram_reference "$fixtures/synthetic_hg38.fa" \
        --sample_type tumor_wgs \
        --profiling_method bowtie2 \
        --hg38_db "$fixtures/hg38.mmi" \
        --t2t_phix_db "$fixtures/t2t.mmi" \
        --save_intermediates true \
        --conda_cache_dir "$CONDA_CACHE_DIR" \
        --outdir "$work/results_cram"
) > "$nf_cram_log" 2>&1
nf_cram_status=$?

if [[ $nf_cram_status -ne 0 ]]; then
    echo "ERROR: the real CRAM-input pipeline run failed (exit $nf_cram_status). Last 60 lines:" >&2
    tail -60 "$nf_cram_log" >&2
    echo "Full log: $nf_cram_log" >&2
    KEEP_WORK=1
    exit 1
fi

echo "CRAM-input pipeline finished. Checking real output values."
if ! python3 "$repo_dir/tests/fixtures/assert_e2e_outputs.py" \
        --results "$work/results_cram" \
        --positive-sample "$POS_CRAM_SAMPLE" \
        --negative-sample "$NEG_CRAM_SAMPLE" \
        --expect-extracted-unmapped-reads 8; then
    echo "ERROR: CRAM-input e2e output assertions failed. Results kept at: $work/results_cram" >&2
    KEEP_WORK=1
    exit 1
fi

# ---------------------------------------------------------------------------------
# Fourth run: --profiling_method hmm alone, on the original paired-FASTQ sheet.
# --sample_type tumor_wgs rejects --profiling_method hmm outright (main.nf: "requires
# --profiling_method bowtie2 or both for read-level clb screening"), so this is
# --sample_type metagenome instead -- the same sample sheet format, just a different
# profiling lane. Reuses $work/sheet.csv's positive/negative FASTQ pairs so the exact
# same real reads that produce a real bowtie2 hit also produce (or don't) a real
# nhmmscan hit, under the population DNA HMM this repo already ships pressed
# (ref/hmm/clb_population_dna_exact_v1.hmm) -- 19 small profiles, real hits in well
# under a second even without --hmm_chunking.
# ---------------------------------------------------------------------------------
echo
echo "Running the real pipeline with --profiling_method hmm (sample_type metagenome)"
nf_hmm_log="$work/nextflow_hmm.log"
(
    cd "$work" && nextflow run "$repo_dir/main.nf" \
        -work-dir "$work/work_hmm" \
        -c "$work/e2e_local.config" \
        --sample "$work/sheet.csv" \
        --input_data_type fastq \
        --sample_type metagenome \
        --profiling_method hmm \
        --hg38_db "$fixtures/hg38.mmi" \
        --t2t_phix_db "$fixtures/t2t.mmi" \
        --save_intermediates true \
        --conda_cache_dir "$CONDA_CACHE_DIR" \
        --outdir "$work/results_hmm"
) > "$nf_hmm_log" 2>&1
nf_hmm_status=$?

if [[ $nf_hmm_status -ne 0 ]]; then
    echo "ERROR: the real HMM-profiling pipeline run failed (exit $nf_hmm_status). Last 60 lines:" >&2
    tail -60 "$nf_hmm_log" >&2
    echo "Full log: $nf_hmm_log" >&2
    KEEP_WORK=1
    exit 1
fi

echo "HMM-profiling pipeline finished. Checking real output values."
if ! python3 "$repo_dir/tests/fixtures/assert_e2e_outputs.py" \
        --results "$work/results_hmm" \
        --positive-sample "$POS_SAMPLE" \
        --negative-sample "$NEG_SAMPLE" \
        --profiling-method hmm; then
    echo "ERROR: HMM-profiling e2e output assertions failed. Results kept at: $work/results_hmm" >&2
    KEEP_WORK=1
    exit 1
fi

# ---------------------------------------------------------------------------------
# Fifth run: a 5-sample cohort (positive, negative, and three new fixtures spanning
# broad_island, extensive_island and localized_indeterminate -- see
# BROAD_WINDOWS/EXTENSIVE_WINDOWS/BORDERLINE_WINDOWS in generate_e2e_fixtures.py) so
# the cohort-level reducers aggregate more than two rows at least once. A 2-sample
# cohort cannot show a bug in sorting, de-duplication, or a boundary condition that
# only appears at N>2; this can.
# ---------------------------------------------------------------------------------
echo
echo "Running the real pipeline on a 5-sample cohort: tumor_wgs, bowtie2"
BROAD_SAMPLE=e2e_broad
EXTENSIVE_SAMPLE=e2e_extensive
BORDERLINE_SAMPLE=e2e_borderline

printf 'patient,fastq1,fastq2\n' > "$work/cohort_sheet.csv"
printf '%s,%s,%s\n' "$POS_SAMPLE" "$fixtures/positive_R1.fastq.gz" "$fixtures/positive_R2.fastq.gz" >> "$work/cohort_sheet.csv"
printf '%s,%s,%s\n' "$NEG_SAMPLE" "$fixtures/negative_R1.fastq.gz" "$fixtures/negative_R2.fastq.gz" >> "$work/cohort_sheet.csv"
printf '%s,%s,%s\n' "$BROAD_SAMPLE" "$fixtures/broad_R1.fastq.gz" "$fixtures/broad_R2.fastq.gz" >> "$work/cohort_sheet.csv"
printf '%s,%s,%s\n' "$EXTENSIVE_SAMPLE" "$fixtures/extensive_R1.fastq.gz" "$fixtures/extensive_R2.fastq.gz" >> "$work/cohort_sheet.csv"
printf '%s,%s,%s\n' "$BORDERLINE_SAMPLE" "$fixtures/borderline_R1.fastq.gz" "$fixtures/borderline_R2.fastq.gz" >> "$work/cohort_sheet.csv"

nf_cohort_log="$work/nextflow_cohort.log"
(
    cd "$work" && nextflow run "$repo_dir/main.nf" \
        -work-dir "$work/work_cohort" \
        -c "$work/e2e_local.config" \
        --sample "$work/cohort_sheet.csv" \
        --input_data_type fastq \
        --sample_type tumor_wgs \
        --profiling_method bowtie2 \
        --hg38_db "$fixtures/hg38.mmi" \
        --t2t_phix_db "$fixtures/t2t.mmi" \
        --save_intermediates true \
        --conda_cache_dir "$CONDA_CACHE_DIR" \
        --outdir "$work/results_cohort"
) > "$nf_cohort_log" 2>&1
nf_cohort_status=$?

if [[ $nf_cohort_status -ne 0 ]]; then
    echo "ERROR: the real 5-sample cohort pipeline run failed (exit $nf_cohort_status). Last 60 lines:" >&2
    tail -60 "$nf_cohort_log" >&2
    echo "Full log: $nf_cohort_log" >&2
    KEEP_WORK=1
    exit 1
fi

echo "5-sample cohort pipeline finished. Checking real output values."
if ! python3 "$repo_dir/tests/fixtures/assert_e2e_outputs.py" \
        --results "$work/results_cohort" \
        --positive-sample "$POS_SAMPLE" \
        --negative-sample "$NEG_SAMPLE" \
        --cohort-samples "${BROAD_SAMPLE}:broad_island,${EXTENSIVE_SAMPLE}:extensive_island,${BORDERLINE_SAMPLE}:localized_indeterminate"; then
    echo "ERROR: 5-sample cohort e2e output assertions failed. Results kept at: $work/results_cohort" >&2
    KEEP_WORK=1
    exit 1
fi

# ---------------------------------------------------------------------------------
# Sixth run: single-end FASTQ (fastq1 only, no fastq2). Ludmil's u-3 wish list names
# this explicitly and every real run above is paired-end. main.nf's fastq branch only
# adds fastq2 to the reads list `if (fastq2)` is non-empty (main.nf's sample_sheet_fastq
# mapping), and preflight.py only requires a fastq1 column -- fastq2 is checked for
# existence/non-emptiness like any other optional column, never required. Single-end is
# a real, already-supported input mode, not one invented for this test.
#
# No new fixtures needed: filterReads' fastp lane runs single-input files through `fastp
# -i` directly (Modules/filter_reads.nf's `input_list.size() == 1` branch), so any one
# of the existing paired-end fixture's own FASTQ halves is already a valid, independently
# real 150 bp single-end read set --
#   - positive: broad_R1.fastq.gz alone -- the 17 BROAD_WINDOWS mate1 reads (see
#     generate_e2e_fixtures.py), each individually confirmed by the same standalone
#     bowtie2 --very-sensitive --no-unal | samtools view -q40 run that validated the
#     paired broad_island fixture (every mate maps uniquely, MAPQ 42, CIGAR 150M) --
#     but taken alone (no mate2), covering half the bases the paired run does: 17 * 150
#     = 2550 bp of the 50,768 bp island = 5.02% breadth. That clears multi_gene's floor
#     (>=5 reads, >=3 genes, >=1% breadth) but not broad_island's (>=7.5% breadth), so
#     single-end lands on multi_gene here -- a real, different tier than its paired-end
#     namesake, not the same result recomputed.
#   - negative: negative_R1.fastq.gz alone -- 4 reads of unrelated pseudorandom sequence.
#
# Neither fastp nor host depletion mate-suffixes a single-input read (there is no mate
# to tag), so assert_e2e_outputs.py's --single-end flag swaps the paired-end u-1
# mate-survival check for a single-end read-identity one, and the expected reads/genes/
# breadth are overridden to match this fixture instead of the paired-end constants.
# ---------------------------------------------------------------------------------
echo
echo "Running the real pipeline on single-end FASTQ input (fastq1 only, no fastq2)"
POS_SE_SAMPLE=e2e_single_end_positive
NEG_SE_SAMPLE=e2e_single_end_negative

printf 'patient,fastq1\n' > "$work/single_end_sheet.csv"
printf '%s,%s\n' "$POS_SE_SAMPLE" "$fixtures/broad_R1.fastq.gz" >> "$work/single_end_sheet.csv"
printf '%s,%s\n' "$NEG_SE_SAMPLE" "$fixtures/negative_R1.fastq.gz" >> "$work/single_end_sheet.csv"

nf_se_log="$work/nextflow_single_end.log"
(
    cd "$work" && nextflow run "$repo_dir/main.nf" \
        -work-dir "$work/work_single_end" \
        -c "$work/e2e_local.config" \
        --sample "$work/single_end_sheet.csv" \
        --input_data_type fastq \
        --sample_type tumor_wgs \
        --profiling_method bowtie2 \
        --hg38_db "$fixtures/hg38.mmi" \
        --t2t_phix_db "$fixtures/t2t.mmi" \
        --save_intermediates true \
        --conda_cache_dir "$CONDA_CACHE_DIR" \
        --outdir "$work/results_single_end"
) > "$nf_se_log" 2>&1
nf_se_status=$?

if [[ $nf_se_status -ne 0 ]]; then
    echo "ERROR: the real single-end FASTQ pipeline run failed (exit $nf_se_status). Last 60 lines:" >&2
    tail -60 "$nf_se_log" >&2
    echo "Full log: $nf_se_log" >&2
    KEEP_WORK=1
    exit 1
fi

echo "Single-end FASTQ pipeline finished. Checking real output values."
if ! python3 "$repo_dir/tests/fixtures/assert_e2e_outputs.py" \
        --results "$work/results_single_end" \
        --positive-sample "$POS_SE_SAMPLE" \
        --negative-sample "$NEG_SE_SAMPLE" \
        --single-end \
        --positive-source-fastq "$fixtures/broad_R1.fastq.gz" \
        --negative-source-fastq "$fixtures/negative_R1.fastq.gz" \
        --positive-expected-reads 17 \
        --positive-expected-genes 8 \
        --positive-gene-names clbS,clbQ,clbD,clbA,clbG,clbF,clbM,clbN \
        --positive-min-breadth 0.03 \
        --positive-max-breadth 0.075 \
        --negative-expected-reads 4; then
    echo "ERROR: single-end FASTQ e2e output assertions failed. Results kept at: $work/results_single_end" >&2
    KEEP_WORK=1
    exit 1
fi

echo
echo "e2e execution check passed (paired FASTQ, BAM, CRAM, single-end FASTQ, HMM profiling, and a 5-sample cohort; positive/negative/broad/extensive/borderline fixtures)"
