#!/usr/bin/env bash
# Run a real, non-preview Nextflow execution on deliberately malformed BAM input, and
# check that the pipeline fails visibly and informatively rather than hanging,
# crashing obscurely, or silently producing wrong output.
#
# F21 (Ludmil's audit): "...zero, positive and malformed inputs." Every other real
# execution check in this suite (tests/check_e2e_execution.sh,
# tests/check_metagenome_mag_execution.sh) proves the pipeline computes the right
# answer on good input; this one proves the opposite direction.
#
# Three samples, one cohort, one real run:
#   - a healthy BAM (real reads, real synthetic host references -- reuses
#     tests/fixtures/generate_e2e_fixtures.py's "positive" fixture, built the same way
#     check_e2e_execution.sh's own BAM-input run builds it);
#   - truncated.bam -- a real BAM cut off partway through. Its leading BGZF magic is
#     intact, so extractReads.nf's htsfile-based format detection passes it through to
#     `samtools quickcheck -v`, which is the check that must (and empirically does,
#     samtools 1.21) fail it with "was missing EOF block when one should be present.";
#   - bad_magic.bam -- a real BAM with its first 4 bytes overwritten, so `htsfile`
#     reports "unknown data" and extractReads.nf's own explicit branch fires: "input is
#     neither BAM nor CRAM according to htsfile: ...".
#     (see tests/fixtures/generate_malformed_input_fixtures.py for both.)
#
# `params.sample_failure_strategy` defaults to 'ignore' (main.nf) for exactly this
# lane (extractReads/filterReads/mapReads/pksProfiler_align/pksProfiler_hmm, conf/
# base.config's withName selector, Ludmil's revised report finding 8) so one sample's
# bad input must not take the whole cohort down -- this run proves the healthy sample
# still completes correctly, the two malformed samples are visibly recorded as failed
# (not silently dropped, and not reported as though they succeeded), and the pipeline's
# own RUN_REPORT.txt WHAT FAILED section points at each failed task's real, specific
# .command.err diagnostic. This is also the live half of Ludmil's finding-8 acceptance
# test that tests/test_sample_stage_failure_policy.py's own docstring says needs a live
# Nextflow run: "force one sample to fail ... other samples must complete and cohort
# reducers must still run."
#
# A hang is a real failure mode this guards against too: the whole run is wrapped in
# `timeout`, and a timeout (exit 124) is reported as its own distinct error rather than
# being confused with an ordinary nonzero pipeline exit.
#
#     CONDA_CACHE_DIR=/path/to/cache bash tests/check_malformed_input_execution.sh
#
# See tests/check_e2e_execution.sh's own header for what CONDA_CACHE_DIR does and why.

set -uo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONDA_CACHE_DIR=${CONDA_CACHE_DIR:-/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache}
NEXTFLOW_TIMEOUT=${NEXTFLOW_TIMEOUT:-600}  # seconds; a real hang must not block CI/CD forever

if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: nextflow is required for e2e execution checks." >&2
    exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is required for e2e execution checks." >&2
    exit 1
fi

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

MINIMAP2_ENV_SAMTOOLS="$(dirname "$MINIMAP2_BIN")/samtools"
if [[ ! -x "$MINIMAP2_ENV_SAMTOOLS" ]]; then
    MINIMAP2_ENV_SAMTOOLS=$(command -v samtools || true)
fi
if [[ -z "$MINIMAP2_ENV_SAMTOOLS" ]]; then
    echo "ERROR: no samtools found (needed to build the BAM fixtures). Put samtools on PATH." >&2
    exit 1
fi

work=$(mktemp -d "${TMPDIR:-/tmp}/pks-malformed.XXXXXX")
echo "work dir: $work (kept on failure for inspection)"
cleanup() { [[ "${KEEP_WORK:-0}" == "1" ]] || rm -rf "$work"; }
trap cleanup EXIT

fixtures="$work/fixtures"
mkdir -p "$fixtures"

echo "Building the healthy fixture (real reads, real synthetic host references)"
if ! python3 "$repo_dir/tests/fixtures/generate_e2e_fixtures.py" --out-dir "$fixtures" > "$work/fixtures.log" 2>&1; then
    cat "$work/fixtures.log" >&2
    echo "ERROR: fixture generation failed." >&2
    exit 1
fi
cat "$work/fixtures.log"

echo "Building the malformed fixtures (truncated.bam, bad_magic.bam)"
if ! python3 "$repo_dir/tests/fixtures/generate_malformed_input_fixtures.py" \
        --out-dir "$fixtures" --samtools-bin "$MINIMAP2_ENV_SAMTOOLS" > "$work/malformed.log" 2>&1; then
    cat "$work/malformed.log" >&2
    echo "ERROR: malformed fixture generation failed." >&2
    exit 1
fi
cat "$work/malformed.log"

echo "Building a real minimap2 index for the synthetic hg38 reference ($MINIMAP2_BIN)"
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

echo "Building the healthy BAM fixture (aligned-and-unmapped, like check_e2e_execution.sh's)"
if ! "$MINIMAP2_BIN" -ax sr "$fixtures/hg38.mmi" \
        "$fixtures/positive_R1.fastq.gz" "$fixtures/positive_R2.fastq.gz" 2> "$work/mm2_bam.log" \
        | "$MINIMAP2_ENV_SAMTOOLS" sort -o "$fixtures/healthy.bam" -; then
    cat "$work/mm2_bam.log" >&2
    echo "ERROR: could not build the healthy BAM fixture." >&2
    exit 1
fi
"$MINIMAP2_ENV_SAMTOOLS" index "$fixtures/healthy.bam"

HEALTHY_SAMPLE=malformed_check_healthy
TRUNCATED_SAMPLE=malformed_check_truncated
BAD_MAGIC_SAMPLE=malformed_check_bad_magic

printf 'patient,alignment\n' > "$work/sheet.csv"
printf '%s,%s\n' "$HEALTHY_SAMPLE" "$fixtures/healthy.bam" >> "$work/sheet.csv"
printf '%s,%s\n' "$TRUNCATED_SAMPLE" "$fixtures/truncated.bam" >> "$work/sheet.csv"
printf '%s,%s\n' "$BAD_MAGIC_SAMPLE" "$fixtures/bad_magic.bam" >> "$work/sheet.csv"

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

echo "Running the real pipeline (one healthy BAM, two malformed BAMs; default"
echo "params.sample_failure_strategy=ignore), bounded by a ${NEXTFLOW_TIMEOUT}s timeout"
nf_log="$work/nextflow.log"
nf_work_dir="$work/work"
(
    cd "$work" && timeout "$NEXTFLOW_TIMEOUT" nextflow run "$repo_dir/main.nf" \
        -work-dir "$nf_work_dir" \
        -c "$work/e2e_local.config" \
        --sample "$work/sheet.csv" \
        --input_data_type bam \
        --sample_type tumor_wgs \
        --profiling_method bowtie2 \
        --hg38_db "$fixtures/hg38.mmi" \
        --t2t_phix_db "$fixtures/t2t.mmi" \
        --save_intermediates true \
        --conda_cache_dir "$CONDA_CACHE_DIR" \
        --outdir "$work/results"
) > "$nf_log" 2>&1
nf_status=$?

if [[ $nf_status -eq 124 ]]; then
    echo "ERROR: the pipeline HUNG (killed after ${NEXTFLOW_TIMEOUT}s) on malformed input" \
         "-- this is exactly the failure mode this check exists to rule out. Last 60 lines:" >&2
    tail -60 "$nf_log" >&2
    echo "Full log: $nf_log" >&2
    KEEP_WORK=1
    exit 1
fi

if [[ $nf_status -ne 0 ]]; then
    # Not necessarily wrong -- a whole-run failure is still not a hang and not silent --
    # but params.sample_failure_strategy=ignore means the expected shape of a healthy
    # sample alongside malformed ones is a completed run with the failure recorded, not
    # a failed run. Surface it clearly rather than assume either interpretation.
    echo "NOTE: the pipeline exited nonzero ($nf_status) rather than completing with the"
    echo "      malformed samples recorded as failed. This may still be correct (an"
    echo "      obscure hang or crash is what this check rules out, not a nonzero exit"
    echo "      per se) but check the log before assuming the assertions below are"
    echo "      meaningful against \$work/results. Last 60 lines:"
    tail -60 "$nf_log"
fi

echo "Pipeline finished (exit $nf_status). Checking real output values."
if ! python3 "$repo_dir/tests/fixtures/assert_malformed_input_outputs.py" \
        --results "$work/results" \
        --work-dir "$nf_work_dir" \
        --healthy-sample "$HEALTHY_SAMPLE" \
        --truncated-sample "$TRUNCATED_SAMPLE" \
        --bad-magic-sample "$BAD_MAGIC_SAMPLE"; then
    echo "ERROR: malformed-input assertions failed. Results kept at: $work/results" >&2
    echo "       Nextflow work dir kept at: $nf_work_dir" >&2
    KEEP_WORK=1
    exit 1
fi

echo
echo "malformed-input execution check passed (healthy BAM unaffected; truncated BAM and" \
     "bad-magic BAM both fail with their real, specific diagnostics recorded)"
