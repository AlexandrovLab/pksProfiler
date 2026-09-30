#!/usr/bin/env bash
# F15, proved by execution: change a reference file's CONTENT, resume, and the tasks that
# depend on it must re-run.
#
# Why this exists. F15 is the finding that Nextflow hashes a task from its declared inputs
# and its rendered script, so a reference, annotation or helper that appears only as a
# path string inside the script block is neither -- rebuild the index in place or edit the
# helper, and every task hash is unchanged, `-resume` reuses the old work, and the results
# no longer correspond to the current references. lib/Provenance.groovy fixes it by
# interpolating content/metadata digests into the affected scripts.
#
# tests/test_dependency_digests.py checks that mechanism as source text. It cannot check
# that Nextflow's resume actually behaves the way the mechanism intends, because that is a
# property of a real run's task hashing, not of the file. Ludmil's F21 asks for exactly
# this kind of execution test ("resume/reference-change tests"), and it is the one part of
# the F15 fix nothing has ever executed.
#
# Three real runs on the tiny e2e fixtures:
#
#   run 1  a clean run; every task executes
#   run 2  -resume, nothing changed; every task must be CACHED. Without this control,
#          run 3's "the task re-ran" would prove nothing -- a resume that never caches
#          anything would pass the same assertion.
#   run 3  the clb annotation GFF's CONTENT is changed (a comment line appended, so the
#          parsed annotation is identical and only the bytes differ, which is precisely
#          the case a path string cannot see), then -resume. Every process carrying the
#          pks_annotation digest must re-run, and the upstream processes that do not carry
#          it -- extraction, read QC, host depletion -- must still be CACHED. Re-running
#          everything would also satisfy "it re-ran", and would mean the digest is doing
#          nothing; this pins both directions.
#
#     CONDA_CACHE_DIR=/path/to/cache bash tests/check_resume_reference_change.sh
#
# Same fixtures, same conda environments and the same CONDA_CACHE_DIR contract as
# tests/check_e2e_execution.sh -- see that script's header for what the cache buys.

set -uo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONDA_CACHE_DIR=${CONDA_CACHE_DIR:-/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache}

for tool in nextflow python3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool is required for the resume/reference-change check." >&2
        exit 1
    fi
done

resolve_minimap2() {
    if [[ -n "${MINIMAP2_BIN:-}" ]]; then echo "$MINIMAP2_BIN"; return 0; fi
    if [[ -d "$CONDA_CACHE_DIR" ]]; then
        local hit
        hit=$(grep -Fls "conda_envs/minimap2_env.yml" "$CONDA_CACHE_DIR"/env-*/conda-meta/history 2>/dev/null | head -1)
        if [[ -n "$hit" ]]; then
            local env_bin="${hit%/conda-meta/history}/bin/minimap2"
            [[ -x "$env_bin" ]] && { echo "$env_bin"; return 0; }
        fi
    fi
    command -v minimap2 || true
}

MINIMAP2_BIN=$(resolve_minimap2)
if [[ -z "$MINIMAP2_BIN" ]]; then
    echo "ERROR: no minimap2 found. Set MINIMAP2_BIN or put minimap2 on PATH." >&2
    exit 1
fi

work=$(mktemp -d "${TMPDIR:-/tmp}/pks-resume.XXXXXX")
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

for name in hg38 t2t; do
    src="synthetic_hg38.fa"; [[ $name == t2t ]] && src="synthetic_t2t_phix.fa"
    if ! "$MINIMAP2_BIN" -d "$fixtures/$name.mmi" "$fixtures/$src" > "$work/mm2_$name.log" 2>&1; then
        cat "$work/mm2_$name.log" >&2
        echo "ERROR: could not build the synthetic $name index." >&2
        exit 1
    fi
done

# The reference under test: a private copy of the shipped clb annotation, so the run reads
# a file this script may edit without ever touching the repo's own ref/ tree.
annotation="$work/IHE3034.clbA-clbS.gff"
cp "$repo_dir/ref/annotations/IHE3034.clbA-clbS.gff" "$annotation"

printf 'patient,fastq1,fastq2\n' > "$work/sheet.csv"
printf 'resume_positive,%s,%s\n' "$fixtures/positive_R1.fastq.gz" "$fixtures/positive_R2.fastq.gz" >> "$work/sheet.csv"

cat > "$work/local.config" <<'EOF'
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

run_pipeline() {   # run_pipeline <label> <extra nextflow args...>
    local label=$1; shift
    local log="$work/nextflow_${label}.log"
    (
        cd "$work" && nextflow run "$repo_dir/main.nf" \
            -work-dir "$work/work" \
            -c "$work/local.config" \
            -with-trace "$work/trace_${label}.txt" \
            --sample "$work/sheet.csv" \
            --input_data_type fastq \
            --sample_type wgs \
            --profiling_method bowtie2 \
            --hg38_db "$fixtures/hg38.mmi" \
            --t2t_phix_db "$fixtures/t2t.mmi" \
            --pks_genome_annotation "$annotation" \
            --conda_cache_dir "$CONDA_CACHE_DIR" \
            --outdir "$work/results_${label}" \
            "$@"
    ) > "$log" 2>&1
    local status=$?
    if [[ $status -ne 0 ]]; then
        echo "ERROR: pipeline run '$label' failed (exit $status). Last 60 lines:" >&2
        tail -60 "$log" >&2
        KEEP_WORK=1
        exit 1
    fi
}

echo
echo "Run 1: clean run"
run_pipeline first

echo "Run 2: -resume with nothing changed (the control)"
run_pipeline unchanged -resume

echo "Run 3: -resume after changing the clb annotation's content"
# Appending a comment line: the parsed annotation is byte-for-byte equivalent in meaning,
# and only Provenance.code()'s content digest can see the difference. A path string, an
# mtime check or a "same size" heuristic all miss it.
printf '# content change for tests/check_resume_reference_change.sh\n' >> "$annotation"
run_pipeline changed -resume

echo
echo "Checking which tasks re-ran"
if ! python3 "$repo_dir/tests/fixtures/assert_resume_reference_change.py" \
        --unchanged-trace "$work/trace_unchanged.txt" \
        --changed-trace "$work/trace_changed.txt" \
        --repo "$repo_dir"; then
    echo "ERROR: resume/reference-change assertions failed. Traces kept at: $work" >&2
    KEEP_WORK=1
    exit 1
fi

echo
echo "Resume/reference-change check passed."
