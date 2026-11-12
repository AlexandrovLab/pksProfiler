#!/usr/bin/env bash
# v3 review item 18: does --prefilter_mode move pks_reads_per_million?
#
# The comment in main.nf used to say the rate "does not move with --prefilter_mode",
# reasoning that the denominator is taken at host depletion, upstream of the prefilter.
# That is true of the DENOMINATOR and this check confirms it. The review's point is about
# the other half: under `balanced`, PROFILING_READS is the narrowed CANDIDATE_OUT.reads,
# so pks_reads -- the NUMERATOR -- is counted on a prefiltered stream.
#
# main.nf already refuses `balanced` without --diamond_rescue true, on the grounds that
# "without the rescue, clb reads outside taxid <primary> are discarded". So the design
# intends the numerator to survive. This runs the same fixture, same sample, twice,
# differing in exactly one flag, and checks whether it does.
#
# WHAT IS ASSERTED
#   host_depleted_reads   IDENTICAL between modes   -- the denominator is pre-prefilter
#   pks_reads             reported either way; a difference is the finding, not a failure
#   pks_reads_per_million follows from the two above
#
# A difference in pks_reads does NOT fail this check: it is the measurement item 18 asks
# for, and the honest outcome either way. The check fails only if a run fails, if the
# denominator moves (which would contradict the pipeline's own design), or if neither
# mode finds any pks read at all -- a 0-vs-0 comparison proves nothing and must not be
# read as parity.
#
# Public references and synthetic reads only (tests/fixtures/generate_mode_fixtures.py).
#
#     CONDA_CACHE_DIR=/path/to/cache KRAKEN_DB=/path/to/krakenUniq_8_8_2023 \
#         bash tests/check_prefilter_rate_parity.sh
#
# KRAKEN_DB is required: `balanced` is the whole point of the comparison and needs a real
# krakenuniq database. On TSCC:
#   KRAKEN_DB=/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023

PKS_TOOL_VERSION_CHECK="${PKS_TOOL_VERSION_CHECK:-strict}"
set -uo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONDA_CACHE_DIR=${CONDA_CACHE_DIR:-/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache}
GRCH38=${GRCH38:-/tscc/projects/ps-lalexandrov/shared/Reference_Genomes/GRCh38.d1.vd1/GRCh38.d1.vd1.fa}
KRAKEN_DB=${KRAKEN_DB:-}
MINIMAP2_BIN=${MINIMAP2_BIN:-/tscc/projects/ps-lalexandrov/amabbasi/anaconda3/envs/minimap2_env/bin/minimap2}

for tool in nextflow python3; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool is required." >&2; exit 1; }
done
if [[ -z "$KRAKEN_DB" ]]; then
    echo "ERROR: KRAKEN_DB must point at a krakenuniq database; --prefilter_mode balanced" >&2
    echo "       cannot run without one, and that mode is the subject of this check." >&2
    exit 1
fi
[[ -d "$KRAKEN_DB" ]] || { echo "ERROR: KRAKEN_DB is not a directory: $KRAKEN_DB" >&2; exit 1; }
[[ -x "$MINIMAP2_BIN" ]] || { echo "ERROR: minimap2 not executable: $MINIMAP2_BIN" >&2; exit 1; }

work=$(mktemp -d "${TMPDIR:-/tmp}/pks-prefilter-parity-XXXXXX")
KEEP_WORK=0
cleanup() { (( KEEP_WORK )) && echo "work kept at $work" >&2 || rm -rf "$work"; }
trap cleanup EXIT

fixtures="$work/fixtures"
echo "=== building fixtures (public references, synthetic reads)"
python3 "$repo_dir/tests/fixtures/generate_mode_fixtures.py" \
    --out-dir "$fixtures" --grch38 "$GRCH38" > "$work/fixtures.log" 2>&1 || {
        echo "ERROR: fixture generation failed. Last 30 lines:" >&2
        tail -30 "$work/fixtures.log" >&2; KEEP_WORK=1; exit 1; }

# The generator writes FASTAs; the .mmi indexes are built here, exactly as
# check_mode_execution.sh does it. Without this the run reaches preflight and stops with
# "--hg38_db does not exist", which is what happened on the first two attempts.
echo "=== building host indexes ($MINIMAP2_BIN)"
"$MINIMAP2_BIN" -d "$fixtures/hg38.mmi" "$fixtures/synthetic_hg38.fa" 2> "$work/mm2.log" &&
"$MINIMAP2_BIN" -d "$fixtures/t2t.mmi" "$fixtures/synthetic_t2t_phix.fa" 2>> "$work/mm2.log" || {
    cat "$work/mm2.log" >&2; echo "ERROR: index build failed." >&2; KEEP_WORK=1; exit 1; }

cat > "$work/local.config" <<'EOF'
executor {
    name   = 'local'
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

sample=prefilter_parity
printf 'patient,fastq1,fastq2\n%s,%s,%s\n' "$sample" \
    "$fixtures/positive_R1.fastq.gz" "$fixtures/positive_R2.fastq.gz" > "$work/sheet.csv"

run_one() {
    local mode=$1; shift
    local started=$SECONDS
    echo
    echo "=== --prefilter_mode $mode"
    (
        cd "$work" && nextflow run "$repo_dir/main.nf" \
            --tool_version_check "$PKS_TOOL_VERSION_CHECK" \
            -work-dir "$work/work_$mode" \
            -c "$work/local.config" \
            --sample "$work/sheet.csv" \
            --input_data_type fastq \
            --sample_type metagenome \
            --profiling_method bowtie2 \
            --conda_cache_dir "$CONDA_CACHE_DIR" \
            --outdir "$work/results_$mode" \
            --hg38_db "$fixtures/hg38.mmi" \
            --t2t_phix_db "$fixtures/t2t.mmi" \
            --prefilter_mode "$mode" \
            "$@"
    ) > "$work/nextflow_$mode.log" 2>&1
    local status=$?
    if [[ $status -ne 0 ]]; then
        echo "ERROR: --prefilter_mode $mode run failed (exit $status). Last 40 lines:" >&2
        tail -40 "$work/nextflow_$mode.log" >&2
        KEEP_WORK=1; return 1
    fi
    echo "    completed in $((SECONDS - started))s"
}

run_one off || exit 1
run_one balanced --kraken_db "$KRAKEN_DB" --diamond_rescue true || exit 1

echo
python3 "$repo_dir/tests/fixtures/assert_prefilter_rate_parity.py" \
    --sample "$sample" \
    --off "$work/results_off" \
    --balanced "$work/results_balanced" \
    --out "${PREFILTER_PARITY_OUT:-$repo_dir/validation/item18_prefilter_rate_parity.tsv}"
status=$?
(( status )) && KEEP_WORK=1
exit $status
