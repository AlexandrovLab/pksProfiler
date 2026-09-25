#!/usr/bin/env bash
# Build the workflow graph for each flag combination. Run nothing.
#
# Why this exists. On 21 September two runs died on defects that every unit test passed
# and that `nextflow run --help` did not see, because both scripts compiled fine:
#
#   * F10 gave Bracken a second output, which made `Bracken(...)` a multi-channel object;
#     the call site still did `.set{}` and the next `.map` failed with "Multi-channel
#     output cannot be applied to operator map". That is EVERY --pks_taxa run, and it
#     died 11 seconds in.
#   * cohortReport wrote its report into a `cohort/` subdirectory while declaring the
#     output at the task root. That one died 13 minutes in, after every upstream task
#     had already succeeded.
#
# `nextflow run -preview` wires every channel and executes no task, so the first class of
# failure surfaces in about three seconds without a scheduler, without data, and without
# touching a single patient file. Verified: with the Bracken fix reverted, the taxonomy
# case below fails with exactly the error the real run produced.
#
# The sample sheet points at a one-byte file. Preflight checks that inputs exist and are
# non-empty, which that satisfies, and -preview never reads it. The databases are checked
# for existence only -- no record is opened.
#
#     bash tests/check_workflow_construction.sh
#
# DB_ROOT may be overridden. A combination whose database is absent is SKIPPED and named
# in the summary; it is not silently treated as a pass.

set -uo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DB_ROOT=${DB_ROOT:-/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs}

if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: nextflow is required for workflow-construction checks." >&2
    exit 1
fi

HG38="$DB_ROOT/human-GRC-db.mmi"
T2T="$DB_ROOT/human-GCA-phix-db.mmi"
KRAKEN="$DB_ROOT/krakenUniq_8_8_2023"
GENOMAD="$DB_ROOT/genomad_db_v1.9"
GTDBTK="$DB_ROOT/gtdbtk_r232"
# CheckM2 takes --database_path, which is the .dmnd file, not the directory.
CHECKM2="$DB_ROOT/checkm2_v1.1.0/CheckM2_database/uniref100.KO.1.dmnd"

for required in "$HG38" "$T2T"; do
    if [[ ! -e "$required" ]]; then
        echo "ERROR: $required not found. Set DB_ROOT to a directory holding the" >&2
        echo "       host-depletion references; every combination needs them." >&2
        exit 1
    fi
done

work=$(mktemp -d "${TMPDIR:-/tmp}/pks-construction.XXXXXX")
trap 'rm -rf "$work"' EXIT
printf 'x' > "$work/input.bam"
printf 'patient,bam\nCONSTRUCTION_CHECK,%s/input.bam\n' "$work" > "$work/sheet.csv"

pass=0; fail=0; skip=0
failed_names=(); skipped_names=()

construct() {
    local name=$1; shift
    local n=$((pass + fail + skip + 1))
    local out status
    out=$(cd "$work" && nextflow run "$repo_dir/main.nf" -preview \
            -work-dir "$work/w$n" \
            --sample "$work/sheet.csv" \
            --input_data_type bam \
            --hg38_db "$HG38" \
            --t2t_phix_db "$T2T" \
            --outdir "$work/r$n" \
            "$@" 2>&1)
    status=$?
    if [[ $status -eq 0 ]] && ! grep -q '\[FAILED\]' <<< "$out"; then
        printf '  ok    %s\n' "$name"
        pass=$((pass + 1))
    else
        printf '  FAIL  %s\n' "$name"
        grep -E '^\[ERROR\]|Multi-channel|already been used|No such variable|Missing output|Unknown method|Cannot invoke' <<< "$out" \
            | sed 's/^/          /' | head -4
        fail=$((fail + 1)); failed_names+=("$name")
    fi
}

need() {   # need <name> <path...> -- skip, loudly, if any is missing
    local name=$1; shift
    for p in "$@"; do
        if [[ ! -e "$p" ]]; then
            printf '  skip  %s (missing %s)\n' "$name" "$(basename "$p")"
            skip=$((skip + 1)); skipped_names+=("$name")
            return 1
        fi
    done
    return 0
}

echo "Building the workflow graph for each flag combination (no tasks run)"

construct "tumour, alignment only"        --sample_type tumor_wgs --profiling_method bowtie2
# tumor_wgs + hmm is rejected by design ("requires bowtie2 or both"); that rejection
# is asserted in run_checks.sh, so the HMM wiring is exercised on the metagenome side.
construct "tumour, both methods"          --sample_type tumor_wgs --profiling_method both
construct "tumour, targeted assembly"     --sample_type tumor_wgs --profiling_method bowtie2 \
                                          --tumor_targeted_assembly true
construct "tumour, widened tier gate"     --sample_type tumor_wgs --profiling_method bowtie2 \
                                          --tumor_targeted_assembly true \
                                          --tumor_contig_tiers multi_gene,broad_island,extensive_island

# The one that broke on 21 September.
if need "tumour, taxonomy" "$KRAKEN"; then
    construct "tumour, taxonomy" --sample_type tumor_wgs --profiling_method bowtie2 \
              --pks_taxa true --kraken_db "$KRAKEN" --bracken_read_length 150
fi

if need "tumour, contig context" "$GENOMAD"; then
    construct "tumour, contig context" --sample_type tumor_wgs --profiling_method bowtie2 \
              --genomad_db "$GENOMAD"
fi

if need "tumour, everything on" "$KRAKEN" "$GENOMAD"; then
    construct "tumour, everything on" --sample_type tumor_wgs --profiling_method both \
              --tumor_targeted_assembly true \
              --tumor_contig_tiers multi_gene,broad_island,extensive_island \
              --pks_taxa true --kraken_db "$KRAKEN" --bracken_read_length 150 \
              --genomad_db "$GENOMAD"
fi

construct "metagenome, profiling only"    --sample_type metagenome --profiling_method bowtie2
construct "metagenome, HMM only"          --sample_type metagenome --profiling_method hmm
construct "metagenome, both methods"      --sample_type metagenome --profiling_method both

if need "metagenome, MAGs" "$GTDBTK" "$CHECKM2" "$GENOMAD"; then
    construct "metagenome, MAGs" --sample_type metagenome --profiling_method bowtie2 \
              --enable_mags true --gtdbtk_db "$GTDBTK" --checkm2_db "$CHECKM2" \
              --genomad_db "$GENOMAD"
fi

if need "metagenome, prefilter" "$KRAKEN"; then
    construct "metagenome, prefilter" --sample_type metagenome --profiling_method bowtie2 \
              --prefilter_mode balanced --kraken_db "$KRAKEN"
fi

echo
echo "construction checks: $pass passed, $fail failed, $skip skipped"
if (( skip )); then
    printf 'skipped (database not present): %s\n' "${skipped_names[*]}"
fi
if (( fail )); then
    printf 'FAILED: %s\n' "${failed_names[*]}" >&2
    exit 1
fi
