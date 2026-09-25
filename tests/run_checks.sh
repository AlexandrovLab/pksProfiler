#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

echo "Running lightweight regression tests"
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s "$repo_dir/tests" -p 'test_*.py' -v

if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: Nextflow is required for workflow validation checks." >&2
    exit 1
fi

echo "Checking duplicate-sample rejection"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/pksprofiler-checks.XXXXXX")
trap 'rm -rf "$test_dir"' EXIT

printf '%s\n' \
    'patient,bam' \
    'duplicate,/does/not/exist/first.bam' \
    'duplicate,/does/not/exist/second.bam' \
    > "$test_dir/duplicate.csv"

set +e
validation_output=$(
    cd "$test_dir"
    nextflow run "$repo_dir/main.nf" \
        -work-dir "$test_dir/work" \
        --sample "$test_dir/duplicate.csv" \
        --input_data_type bam \
        --hg38_db unused-for-validation \
        --t2t_phix_db unused-for-validation \
        --outdir "$test_dir/results" \
        2>&1
)
validation_status=$?
set -e

if [[ "$validation_status" -eq 0 ]] ||
   ! grep -Fq "duplicate patient value 'duplicate'" <<< "$validation_output"; then
    printf '%s\n' "$validation_output" >&2
    echo "ERROR: Duplicate-sample validation regression failed." >&2
    exit 1
fi

echo "All checks passed"

# preflight.py now rejects a run outright on ANY problem, including a merely-missing
# reference file -- and main.nf's own parameter validation (sample_type, tier,
# mag-routing) only runs after preflight succeeds. The checks below exist to reach
# that later, main.nf-level validation, so their fastq1/reference paths must actually
# exist; only the specific flag each check passes should be what's being tested. A
# literal nonexistent placeholder here silently produces a preflight-level failure
# instead, and the expected message never has a chance to appear.
printf 'placeholder' > "$test_dir/placeholder.fastq.gz"
printf 'placeholder' > "$test_dir/placeholder-hg38.mmi"
printf 'placeholder' > "$test_dir/placeholder-t2t.mmi"
mkdir -p "$test_dir/placeholder-gtdbtk-db" "$test_dir/placeholder-genomad-db"
printf 'placeholder' > "$test_dir/placeholder-checkm2-db"

printf '%s\n' \
    'patient,fastq1' \
    "sample1,$test_dir/placeholder.fastq.gz" \
    > "$test_dir/single.csv"

check_validation() {
    local name=$1
    local expected=$2
    shift 2

    local output
    local status
    set +e
    output=$(
        cd "$test_dir"
        nextflow run "$repo_dir/main.nf" \
            -work-dir "$test_dir/work-$name" \
            --sample "$test_dir/single.csv" \
            --input_data_type fastq \
            --hg38_db "$test_dir/placeholder-hg38.mmi" \
            --t2t_phix_db "$test_dir/placeholder-t2t.mmi" \
            --outdir "$test_dir/results-$name" \
            "$@" \
            2>&1
    )
    status=$?
    set -e

    if [[ "$status" -eq 0 ]] || ! grep -Fq -- "$expected" <<< "$output"; then
        printf '%s\n' "$output" >&2
        echo "ERROR: $name validation regression failed." >&2
        exit 1
    fi
}

echo "Checking MAG routing validation"
check_validation \
    invalid-sample-type \
    'Unknown --sample_type: invalid' \
    --sample_type invalid
check_validation \
    tumor-mag \
    '--enable_mags requires --sample_type metagenome' \
    --sample_type tumor_wgs --enable_mags true \
    --gtdbtk_db "$test_dir/placeholder-gtdbtk-db" \
    --checkm2_db "$test_dir/placeholder-checkm2-db" \
    --genomad_db "$test_dir/placeholder-genomad-db"
check_validation \
    missing-gtdbtk \
    '--gtdbtk_db is not set' \
    --sample_type metagenome --enable_mags true
check_validation \
    missing-genomad \
    '--genomad_db is not set' \
    --sample_type metagenome --enable_mags true \
    --gtdbtk_db "$test_dir/placeholder-gtdbtk-db" \
    --checkm2_db "$test_dir/placeholder-checkm2-db"
echo "All checks passed"

echo "Checking tumour contig tier validation"
check_validation \
    unknown-tier \
    'Unknown tier in --tumor_contig_tiers: bogus' \
    --sample_type tumor_wgs --tumor_contig_tiers bogus
check_validation \
    non-monotonic-tiers \
    'Tier thresholds must be non-decreasing' \
    --sample_type tumor_wgs --tumor_broad_island_min_pks_reads 1
check_validation \
    tumor-hmm-only \
    '--sample_type tumor_wgs requires --profiling_method bowtie2 or both' \
    --sample_type tumor_wgs --profiling_method hmm

echo
echo "Building the workflow graph for every flag combination"
# Not piped: a pipeline's exit status is the last command's, so `| tail` would swallow
# the failure this check exists to report.
bash "$repo_dir/tests/check_workflow_construction.sh"

echo "All checks passed"
