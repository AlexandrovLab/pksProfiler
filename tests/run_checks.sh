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
# Owner decision on Ludmil's 26 Sep mode memo: the old names fail, naming the new one.
check_validation \
    renamed-tumor-wgs \
    '--sample_type tumor_wgs was renamed to wgs' \
    --sample_type tumor_wgs
check_validation \
    renamed-tumor-rna \
    '--sample_type tumor_rna was renamed to rna' \
    --sample_type tumor_rna
check_validation \
    tumor-mag \
    '--enable_mags requires --sample_type metagenome' \
    --sample_type wgs --enable_mags true \
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

echo "Checking targeted-recruitment rule validation (Ludmil, 25 Sep A1/A2)"
check_validation \
    recruit-query-coverage \
    '--recruit_min_query_coverage must be a fraction in (0,1]' \
    --sample_type wgs --recruit_min_query_coverage 1.5
check_validation \
    recruit-mapq \
    '--recruit_min_mapq must be a non-negative integer' \
    --sample_type wgs --recruit_min_mapq high
check_validation \
    recruit-tolerance \
    '--recruit_boundary_tolerance_bp must be less than half of --pks_recruit_min_aligned' \
    --sample_type wgs --recruit_boundary_tolerance_bp 30
check_validation \
    recruit-flank-mismatch \
    'recruitment classes would be misplaced' \
    --sample_type wgs --pks_recruit_flank_bp 4000

echo "Checking tumour contig tier validation"
check_validation \
    unknown-tier \
    'Unknown tier in --tumor_contig_tiers: bogus' \
    --sample_type wgs --tumor_contig_tiers bogus
check_validation \
    non-monotonic-tiers \
    'Tier thresholds must be non-decreasing' \
    --sample_type wgs --tumor_broad_island_min_pks_reads 1
check_validation \
    tumor-hmm-only \
    '--sample_type wgs requires --profiling_method bowtie2 or both' \
    --sample_type wgs --profiling_method hmm
check_validation \
    wes-hmm-only \
    '--sample_type wes requires --profiling_method bowtie2 or both' \
    --sample_type wes --profiling_method hmm
# Ludmil, 26 Sep mode memo section 4: rna needs its splice index, not --hg38_db.
check_validation \
    rna-missing-splice-db \
    '--rna_splice_db is not set' \
    --sample_type rna
check_validation \
    rna-bad-library-type \
    '--rna_library_type must be polyA, total or unknown' \
    --sample_type rna --rna_splice_db "$test_dir/placeholder-hg38.mmi" --rna_library_type ribo
# Owner decision, 29 Sep: long reads for wgs and rna only; long-read wgs cannot assemble.
check_validation \
    long-read-metagenome \
    'long-read input is supported for wgs and rna only' \
    --sample_type metagenome --read_type ont \
    --hg38_db_long "$test_dir/placeholder-hg38.mmi" --t2t_phix_db_long "$test_dir/placeholder-t2t.mmi"
check_validation \
    unknown-read-type \
    '--read_type must be one of short, ont, hifi' \
    --sample_type wgs --read_type nanopore
# rc1 pass, item 1: targeted assembly defaults on and is switched off for long reads;
# only an explicit request on the command line is refused.
check_validation \
    long-read-wgs-assembly \
    '--tumor_targeted_assembly true is not supported with --read_type hifi' \
    --sample_type wgs --read_type hifi --tumor_targeted_assembly true \
    --hg38_db_long "$test_dir/placeholder-hg38.mmi" --t2t_phix_db_long "$test_dir/placeholder-t2t.mmi"
check_validation \
    long-read-missing-index \
    '--hg38_db_long is not set' \
    --sample_type wgs --read_type ont
check_validation \
    long-read-pks-taxa \
    "--pks_taxa / --pks_community_taxa are not supported with --read_type ont: Bracken's read-length model" \
    --sample_type wgs --read_type ont --pks_taxa true \
    --hg38_db_long "$test_dir/placeholder-hg38.mmi" --t2t_phix_db_long "$test_dir/placeholder-t2t.mmi"

echo
echo "Running the host-depletion plan (lib/HostDepletion.groovy)"
bash "$repo_dir/tests/check_depletion_plan.sh"

echo
echo "Building the workflow graph for every flag combination"
# Not piped: a pipeline's exit status is the last command's, so `| tail` would swallow
# the failure this check exists to report.
bash "$repo_dir/tests/check_workflow_construction.sh"

echo "All checks passed"
