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
#
# This is slower than the other two test tiers (real conda envs, real alignment) but
# still small: 8 reads per sample, tiny synthetic references, no assembly/MAG/taxonomy
# lane (the default flags don't trigger any of those for either fixture -- see the
# comment above LOCI in generate_e2e_fixtures.py).
#
# Still open against Ludmil's full wish list (BAM, CRAM, single FASTQ, paired FASTQ,
# HMM, cohort aggregation): CRAM input, --profiling_method hmm, and a real multi-sample
# cohort large enough to say something about aggregation beyond "two rows landed in one
# table". See the commit message for why these were left for a follow-up.
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
# 16 cpus/task) and are irrelevant to an 8-read fixture. withLabel selectors take
# precedence over a plain `process {}` block regardless of file order (Nextflow config
# precedence: process < withLabel < withName), so a generic override here would
# silently lose to conf/base.config's per-label blocks; matching labels one for one is
# what actually takes effect. `executor.cpus` is set high enough that Nextflow's local
# executor -- sized from the JVM's available-processors count, which some sandboxed
# shells report as 1 regardless of the node's real core count -- never refuses to admit
# a task that asks for more cpus than that.
cat > "$work/e2e_local.config" <<'EOF'
executor {
    cpus   = 32
    memory = '32 GB'
}
process {
    withLabel:extract_reads  { cpus = 1 }
    withLabel:filter_reads   { cpus = 1 }
    withLabel:map_reads      { cpus = 1 }
    withLabel:pks_align      { cpus = 1 }
    withLabel:pks_hmm        { cpus = 1 }
    withLabel:process_low    { cpus = 1 }
    withLabel:sample_stage   { cpus = 1 }
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

echo
echo "e2e execution check passed (paired FASTQ and BAM input forms, positive and negative fixtures)"
