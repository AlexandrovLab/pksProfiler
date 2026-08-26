#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

echo "Running lightweight regression tests"
PYTHONDONTWRITEBYTECODE=1 python3 \
    "$repo_dir/tests/test_regressions.py"

if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: Nextflow is required for workflow linting." >&2
    exit 1
fi

echo "Linting the Nextflow workflow"
nextflow lint "$repo_dir/main.nf"

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
   ! grep -Fq 'Duplicate patient value(s): duplicate' <<< "$validation_output"; then
    printf '%s\n' "$validation_output" >&2
    echo "ERROR: Duplicate-sample validation regression failed." >&2
    exit 1
fi

echo "All checks passed"
