#!/usr/bin/env bash
# One real execution of each V8 lane, on synthetic reads whose correct answer is known.
#
# Ludmil, 25 Sep report section 4, asks two different questions and says they must not be
# conflated. The unit tests check each decision rule on hand-written SAM records; this
# runs both lanes for real -- real minimap2, real Nextflow channel wiring through each
# module's standalone entry point, real decision scripts -- on fixtures cut from the
# repository's own IHE3034 reference (tests/fixtures/generate_v8_fixtures.py):
#
#   V8a  locus-detection specificity, pks-NEGATIVE panel.
#        A sample mixing reads from the true island slice with reads from a 12%-divergent
#        NRPS/PKS "homolog" slice. Island reads must come out confirmed_true and homolog
#        reads confirmed_decoy. This is the domain-homology failure mode V8a exists for,
#        reproduced at read level.
#
#   V8b  carrier attribution, pks-POSITIVE panel.
#        The same island as two carriers, one E. coli and one Klebsiella copy 3% apart.
#        Reads taken from the Klebsiella copy must be attributed to the Klebsiella
#        carrier, and the sample-level call must name Klebsiella -- which V8a, whose panel
#        contains no carriers at all, could never say.
#
#     CONDA_CACHE_DIR=/path/to/cache bash tests/check_v8_lanes_execution.sh
#
# Same CONDA_CACHE_DIR contract as tests/check_e2e_execution.sh. Nothing here reads sample
# data: every base comes from the shipped public reference or from a seeded RNG.

set -uo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONDA_CACHE_DIR=${CONDA_CACHE_DIR:-/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache}

for tool in nextflow python3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool is required for the V8 lane execution check." >&2
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

work=$(mktemp -d "${TMPDIR:-/tmp}/pks-v8.XXXXXX")
echo "work dir: $work (kept on failure for inspection)"
cleanup() { [[ "${KEEP_WORK:-0}" == "1" ]] || rm -rf "$work"; }
trap cleanup EXIT

fixtures="$work/fixtures"
echo "Building V8 fixtures (slices of the shipped IHE3034 reference)"
if ! python3 "$repo_dir/tests/fixtures/generate_v8_fixtures.py" --out-dir "$fixtures" \
        > "$work/fixtures.log" 2>&1; then
    cat "$work/fixtures.log" >&2
    echo "ERROR: V8 fixture generation failed." >&2
    exit 1
fi
cat "$work/fixtures.log"

echo "Building real minimap2 indices for both panels ($MINIMAP2_BIN)"
for name in v8a_true v8a_decoy v8a_concatenated v8b_carriers; do
    if ! "$MINIMAP2_BIN" -x sr -d "$fixtures/$name.mmi" "$fixtures/$name.fna" \
            > "$work/mm2_$name.log" 2>&1; then
        cat "$work/mm2_$name.log" >&2
        echo "ERROR: could not build the $name index." >&2
        exit 1
    fi
done

cat > "$work/local.config" <<'EOF'
executor {
    cpus   = 8
    memory = '8 GB'
}
process {
    withLabel:sample_stage   { cpus = 1 }
    withLabel:process_medium { cpus = 1; memory = 2.GB }
}
EOF

# ---------------------------------------------------------------------------------
# V8a: locus-detection specificity, through the module's own standalone entry point.
# Arm A is not given, so the lane must run Arm B alone and say so (arms_used).
# ---------------------------------------------------------------------------------
echo
echo "Running the V8a locus-specificity lane for real"
v8a_log="$work/nextflow_v8a.log"
(
    cd "$work" && nextflow run "$repo_dir/Modules/pks_specificity.nf" \
        -work-dir "$work/work_v8a" \
        -c "$work/local.config" \
        --sample "$fixtures/v8a_sheet.csv" \
        --sample_dir "$work/v8a_out" \
        --pks_specificity_true_index "$fixtures/v8a_true.mmi" \
        --pks_specificity_decoy_index "$fixtures/v8a_decoy.mmi" \
        --pks_specificity_concatenated_index "$fixtures/v8a_concatenated.mmi" \
        --pks_specificity_true_contigs "$fixtures/v8a_true_contigs.txt" \
        --pks_specificity_decoy_contigs "$fixtures/v8a_decoy_contigs.txt" \
        --pks_specificity_island_length 50768 \
        --pks_specificity_min_supporting_reads 10 \
        --conda_cache_dir "$CONDA_CACHE_DIR"
) > "$v8a_log" 2>&1
if [[ $? -ne 0 ]]; then
    echo "ERROR: the V8a lane failed. Last 60 lines:" >&2
    tail -60 "$v8a_log" >&2
    KEEP_WORK=1
    exit 1
fi

# The lane must refuse to run without its panel, naming what is missing -- a panel path
# silently defaulting to null would produce an empty check that looks like a pass.
echo "Checking the V8a lane refuses to run without its panel"
refusal=$(cd "$work" && nextflow run "$repo_dir/Modules/pks_specificity.nf" \
    -work-dir "$work/work_v8a_refuse" \
    --sample "$fixtures/v8a_sheet.csv" --sample_dir "$work/v8a_refuse" 2>&1)
if [[ $? -eq 0 ]] || ! grep -Fq -- "--pks_specificity_concatenated_index" <<< "$refusal"; then
    printf '%s\n' "$refusal" >&2
    echo "ERROR: the V8a lane did not refuse a run with no panel." >&2
    KEEP_WORK=1
    exit 1
fi

# ---------------------------------------------------------------------------------
# V8b: carrier attribution, the other question, the other panel.
# ---------------------------------------------------------------------------------
echo
echo "Running the V8b carrier-attribution lane for real"
v8b_log="$work/nextflow_v8b.log"
(
    cd "$work" && nextflow run "$repo_dir/Modules/pks_carrier.nf" \
        -work-dir "$work/work_v8b" \
        -c "$work/local.config" \
        --sample "$fixtures/v8b_sheet.csv" \
        --sample_dir "$work/v8b_out" \
        --pks_carrier_panel_index "$fixtures/v8b_carriers.mmi" \
        --pks_carrier_map "$fixtures/v8b_carrier_map.tsv" \
        --pks_carrier_min_attributed_reads 10 \
        --conda_cache_dir "$CONDA_CACHE_DIR"
) > "$v8b_log" 2>&1
if [[ $? -ne 0 ]]; then
    echo "ERROR: the V8b lane failed. Last 60 lines:" >&2
    tail -60 "$v8b_log" >&2
    KEEP_WORK=1
    exit 1
fi

echo
echo "Checking the real outputs"
if ! python3 "$repo_dir/tests/fixtures/assert_v8_outputs.py" \
        --fixtures "$fixtures" \
        --v8a-arm-b "$work/v8a_out/v8a_sample/locus_specificity/arm_b.tsv" \
        --v8a-summary "$work/v8a_out/v8a_sample/locus_specificity/locus_specificity_summary.tsv" \
        --v8a-per-read "$work/v8a_out/v8a_sample/locus_specificity/locus_specificity_per_read.tsv" \
        --v8b-summary "$work/v8b_out/v8b_sample/carrier_attribution/carrier_attribution.tsv" \
        --v8b-per-read "$work/v8b_out/v8b_sample/carrier_attribution/carrier_per_read.tsv"; then
    echo "ERROR: V8 output assertions failed. Outputs kept at: $work" >&2
    KEEP_WORK=1
    exit 1
fi

echo
echo "V8 lane execution check passed."
