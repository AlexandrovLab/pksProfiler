#!/usr/bin/env bash
# One real execution per mode (Ludmil, 26 Sep mode memo section 6).
#
# "Add one real execution test per modality. Verify which branches run, which outputs
# appear, and which tiers/assembly features are intentionally disabled." Six real,
# non-preview Nextflow runs on tiny inputs built from public references
# (tests/fixtures/generate_mode_fixtures.py):
#
#   wgs         paired pks reads, bowtie2; read_evidence tiered (multi_gene), completed_positive
#   wgs_assembly  extensive_island reads plus flank-only, boundary-spanning and chimeric
#               (60 bp island + 90 bp random) pairs; targeted recruitment (Bowtie2
#               best-alignment, MAPQ + query-coverage + boundary rule) and MEGAHIT/metaSPAdes
#               run for real, and every recruitment class is counted (Ludmil, 25 Sep A1/A2)
#   wes         same reads; wes_evidence.tsv only, no tier, no assembly
#   rna         pks reads + reads across an exon-exon junction on a GRCh38 slice; human pass
#               minimap2 -ax splice:sr -j <junction BED> against a splice:sr index built here;
#               every junction read depleted, rna_evidence.tsv
#   metagenome  same paired reads; read_evidence untiered with pks_reads_per_million
#               (the paired MAG assembly path needs the GTDB-Tk/CheckM2/geNomad databases:
#               tests/check_metagenome_mag_execution.sh runs it)
#   wgs + ont   3 kb pks and human long reads; fastplong, map-ont depletion against a FASTA,
#               minimap2 profiling; untiered, targeted assembly switched off
#   rna + ont   3 kb pks reads + spliced human long reads; -ax splice -j depletion
#
# tests/fixtures/assert_mode_outputs.py checks each run's real output: the evidence file
# of that mode and no other, the counts, the preset the host pass really ran with (from
# the rendered task script), what was disabled, terminal_state and the master summary.
#
#     CONDA_CACHE_DIR=/path/to/cache bash tests/check_mode_execution.sh
#
# CONDA_CACHE_DIR: as tests/check_e2e_execution.sh. The first run in a fresh cache builds
# fastplong_env and pks_align_env (minimap2 2.31), a few minutes each.

set -uo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONDA_CACHE_DIR=${CONDA_CACHE_DIR:-/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache}
GRCH38=${GRCH38:-/tscc/projects/ps-lalexandrov/shared/Reference_Genomes/GRCh38.d1.vd1/GRCh38.d1.vd1.fa}
MODES=${MODES:-"wgs wgs_assembly wes rna metagenome wgs_ont rna_ont"}

for tool in nextflow python3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool is required for the mode execution checks." >&2
        exit 1
    fi
done
if [[ ! -f "$GRCH38.fai" ]]; then
    echo "ERROR: $GRCH38.fai not found (set GRCH38 to a public GRCh38 FASTA with a .fai)." >&2
    exit 1
fi

# minimap2 builds the fixture indexes; same lookup as check_e2e_execution.sh.
resolve_minimap2() {
    if [[ -n "${MINIMAP2_BIN:-}" ]]; then echo "$MINIMAP2_BIN"; return 0; fi
    local hit env_bin
    for hit in $(grep -Fls "conda_envs/minimap2_env.yml" "$CONDA_CACHE_DIR"/env-*/conda-meta/history 2>/dev/null); do
        env_bin="${hit%/conda-meta/history}/bin/minimap2"
        if [[ -x "$env_bin" ]] && "$env_bin" --version 2>/dev/null | grep -q '^2\.31'; then
            echo "$env_bin"; return 0
        fi
    done
    command -v minimap2 || true
}
MINIMAP2_BIN=$(resolve_minimap2)
if [[ -z "$MINIMAP2_BIN" ]]; then
    echo "ERROR: no minimap2 found (set MINIMAP2_BIN)." >&2
    exit 1
fi

work=$(mktemp -d "${TMPDIR:-/tmp}/pks-modes.XXXXXX")
echo "work dir: $work (kept on failure for inspection)"
cleanup() { [[ "${KEEP_WORK:-0}" == "1" ]] || rm -rf "$work"; }
trap cleanup EXIT
fixtures="$work/fixtures"

echo "Building fixtures (public references and reads simulated from them)"
if ! python3 "$repo_dir/tests/fixtures/generate_mode_fixtures.py" --out-dir "$fixtures" \
        --grch38 "$GRCH38" > "$work/fixtures.log" 2>&1; then
    cat "$work/fixtures.log" >&2; echo "ERROR: fixture generation failed." >&2; exit 1
fi

echo "Building indexes ($MINIMAP2_BIN)"
"$MINIMAP2_BIN" -d "$fixtures/hg38.mmi" "$fixtures/synthetic_hg38.fa" 2> "$work/mm2.log" &&
"$MINIMAP2_BIN" -d "$fixtures/t2t.mmi" "$fixtures/synthetic_t2t_phix.fa" 2>> "$work/mm2.log" &&
"$MINIMAP2_BIN" -x splice:sr -d "$fixtures/human_splice_sr.mmi" "$fixtures/human_slice.fa" 2>> "$work/mm2.log" || {
    cat "$work/mm2.log" >&2; echo "ERROR: index build failed." >&2; exit 1; }

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

paired_sheet() { printf 'patient,fastq1,fastq2\n%s,%s,%s\n' "$1" "$fixtures/positive_R1.fastq.gz" "$fixtures/positive_R2.fastq.gz"; }
single_sheet() { printf 'patient,fastq1\n%s,%s\n' "$1" "$2"; }

failed=()
passed=()
run_mode() {
    local mode=$1; shift
    local sample="mode_${mode}"
    local started=$SECONDS
    case "$mode" in
        wgs|wes|metagenome) paired_sheet "$sample" > "$work/$mode.csv" ;;
        wgs_assembly) printf 'patient,fastq1,fastq2\n%s,%s,%s\n' "$sample" \
                          "$fixtures/recruitment_R1.fastq.gz" "$fixtures/recruitment_R2.fastq.gz" > "$work/$mode.csv" ;;
        rna)     single_sheet "$sample" "$fixtures/rna_short.fastq.gz" > "$work/$mode.csv" ;;
        wgs_ont) single_sheet "$sample" "$fixtures/wgs_ont.fastq.gz" > "$work/$mode.csv" ;;
        rna_ont) single_sheet "$sample" "$fixtures/rna_ont.fastq.gz" > "$work/$mode.csv" ;;
    esac
    echo
    echo "=== $mode: nextflow run $*"
    (
        cd "$work" && nextflow run "$repo_dir/main.nf" \
            -work-dir "$work/work_$mode" \
            -c "$work/local.config" \
            --sample "$work/$mode.csv" \
            --input_data_type fastq \
            --conda_cache_dir "$CONDA_CACHE_DIR" \
            --outdir "$work/results_$mode" \
            "$@"
    ) > "$work/nextflow_$mode.log" 2>&1
    local status=$?
    if [[ $status -ne 0 ]]; then
        echo "ERROR: $mode run failed (exit $status). Last 40 lines:" >&2
        tail -40 "$work/nextflow_$mode.log" >&2
        failed+=("$mode"); KEEP_WORK=1; return
    fi
    if python3 "$repo_dir/tests/fixtures/assert_mode_outputs.py" --mode "$mode" \
            --results "$work/results_$mode" --work-dir "$work/work_$mode" --sample "$sample" \
            --fixtures "$fixtures/mode_fixtures.json" --log "$work/nextflow_$mode.log"; then
        passed+=("$mode ($((SECONDS - started))s)")
    else
        failed+=("$mode"); KEEP_WORK=1
    fi
}

for mode in $MODES; do
    case "$mode" in
        wgs)        run_mode wgs --sample_type wgs --profiling_method bowtie2 \
                        --hg38_db "$fixtures/hg38.mmi" --t2t_phix_db "$fixtures/t2t.mmi" ;;
        wgs_assembly) run_mode wgs_assembly --sample_type wgs --profiling_method bowtie2 \
                        --hg38_db "$fixtures/hg38.mmi" --t2t_phix_db "$fixtures/t2t.mmi" ;;
        wes)        run_mode wes --sample_type wes --profiling_method bowtie2 \
                        --hg38_db "$fixtures/hg38.mmi" --t2t_phix_db "$fixtures/t2t.mmi" ;;
        rna)        run_mode rna --sample_type rna --profiling_method bowtie2 \
                        --rna_splice_db "$fixtures/human_splice_sr.mmi" \
                        --rna_junctions_bed "$fixtures/junctions.bed" --rna_library_type polyA \
                        --t2t_phix_db "$fixtures/t2t.mmi" ;;
        metagenome) run_mode metagenome --sample_type metagenome --profiling_method bowtie2 \
                        --prefilter_mode off \
                        --hg38_db "$fixtures/hg38.mmi" --t2t_phix_db "$fixtures/t2t.mmi" ;;
        wgs_ont)    run_mode wgs_ont --sample_type wgs --read_type ont --profiling_method bowtie2 \
                        --hg38_db_long "$fixtures/human_slice.fa" \
                        --t2t_phix_db_long "$fixtures/synthetic_t2t_phix.fa" ;;
        rna_ont)    run_mode rna_ont --sample_type rna --read_type ont --profiling_method bowtie2 \
                        --rna_splice_db_long "$fixtures/human_slice.fa" \
                        --rna_junctions_bed "$fixtures/junctions.bed" \
                        --t2t_phix_db_long "$fixtures/synthetic_t2t_phix.fa" ;;
        *) echo "ERROR: unknown mode $mode" >&2; failed+=("$mode") ;;
    esac
done

echo
echo "mode execution checks: ${#passed[@]} passed, ${#failed[@]} failed (total ${SECONDS}s)"
for p in "${passed[@]}"; do echo "  passed: $p"; done
if (( ${#failed[@]} )); then
    printf '  FAILED: %s\n' "${failed[@]}" >&2
    echo "work kept at $work" >&2
    exit 1
fi
