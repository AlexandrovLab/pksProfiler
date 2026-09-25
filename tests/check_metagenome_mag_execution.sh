#!/usr/bin/env bash
# Run one real, non-preview Nextflow execution of the metagenome --enable_mags lane, on
# a fixture engineered to recover ZERO real bins, and check its real per-sample output.
#
# Why this exists. Found via a real (non-mocked) execution this session: build_mag_summary.py's
# per-sample summary table (genomes/pks_mag_summary.tsv) requires a real CheckM2 report and a
# real GTDB-Tk summary before it will report anything for a sample -- and checkm2Predict/
# gtdbtkClassify only ever run on real bins (metabat2Bin.out.bins, optional: true). A sample
# with zero real bins therefore never had either report, and Modules/pks_mag.nf's
# summary_input_ch used a strict join on both -- so that sample was silently dropped from
# pks_mag_summary.tsv entirely, even when it had a real, positive call from U1's own "unbinned
# pool" (metabat2Bin's --unbinned output, added earlier this session for exactly this "MetaBAT2
# recovered nothing" case). U1's fix computed the right evidence; this session's fix
# (stubBinlessMagQuality) is what actually gets it reported. No unit test catches this: it
# requires an assembly small enough that metabat2 --unbinned pools everything (real tool
# behaviour, not something you can fake with a stand-in database), which only a real execution
# can produce.
#
# This script:
#   1. builds a real fixture (tests/fixtures/generate_metagenome_mag_fixtures.py): paired reads
#      tiled at real, overlapping depth across the complete pks/clb island (all 19 genes) of the
#      repo's own shipped GCF_000025745.1 (IHE3034) reference, sized so the resulting assembly is
#      well under MetaBAT2's 200 kb --minClsSize floor -- confirmed in real testing to make
#      metabat2 --unbinned pool the whole assembly rather than recover a bin;
#   2. builds REAL minimap2 indices for two tiny synthetic host references (same technique as
#      tests/check_e2e_execution.sh);
#   3. runs main.nf for real (no -preview): --sample_type metagenome --profiling_method bowtie2
#      --enable_mags true, against the real production GTDB-Tk/CheckM2/geNomad databases (no
#      small stand-in exists for these -- same reasoning as tests/check_workflow_construction.sh's
#      need()/skip pattern being wrong for an execution test, not just tests/
#      check_tumor_mag_execution.sh's now-removed rationale for the same point);
#   4. asserts on real output (tests/fixtures/assert_metagenome_mag_outputs.py): mag_status.tsv
#      shows zero real bins recovered (confirming the fixture actually landed in the scenario
#      this test exists to cover, not just assuming it), and -- the regression this test
#      guards -- genomes/pks_mag_summary.tsv still gets written, with a real unbinned-pool row
#      carrying real clb gene evidence, and no real checkm2/gtdbtk directories (since no real
#      bin ever existed for either tool to run on).
#
# This is slower than the fast unit-test tier (real conda environments: megahit, metabat2,
# checkm2, gtdbtk, prokka, hmmer, geNomad; a real, if tiny, CheckM2/GTDB-Tk database read for
# the samples this test does NOT expect to hit, since neither actually runs here) but still
# small: one sample, ~1010 read pairs, no multi-sample cohort. Deliberately excluded from the
# fast tier and from CI's default jobs -- see .github/workflows/ci.yml.
#
#     CONDA_CACHE_DIR=/path/to/cache DB_ROOT=/path/to/dbs bash tests/check_metagenome_mag_execution.sh

set -uo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONDA_CACHE_DIR=${CONDA_CACHE_DIR:-/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache}
DB_ROOT=${DB_ROOT:-/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs}

if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: nextflow is required for metagenome MAG execution checks." >&2
    exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is required for metagenome MAG execution checks." >&2
    exit 1
fi

GTDBTK_DB="$DB_ROOT/gtdbtk_r232"
CHECKM2_DB="$DB_ROOT/checkm2_v1.1.0/CheckM2_database/uniref100.KO.1.dmnd"
GENOMAD_DB="$DB_ROOT/genomad_db_v1.9"
for required in "$GTDBTK_DB" "$CHECKM2_DB" "$GENOMAD_DB"; do
    if [[ ! -e "$required" ]]; then
        echo "ERROR: $required not found. This test needs the real GTDB-Tk, CheckM2 and" >&2
        echo "       geNomad databases. Set DB_ROOT to a directory holding all three." >&2
        exit 1
    fi
done

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

work=$(mktemp -d "${TMPDIR:-/tmp}/pks-metagenome-mag.XXXXXX")
echo "work dir: $work (kept on failure for inspection)"
cleanup() { [[ "${KEEP_WORK:-0}" == "1" ]] || rm -rf "$work"; }
trap cleanup EXIT

fixtures="$work/fixtures"
mkdir -p "$fixtures"

echo "Building fixtures: synthetic host refs (base E2E generator, reused for its host-ref builder only)"
if ! python3 "$repo_dir/tests/fixtures/generate_e2e_fixtures.py" --out-dir "$fixtures" --negative-pairs 1 > "$work/fixtures_base.log" 2>&1; then
    cat "$work/fixtures_base.log" >&2
    echo "ERROR: base fixture generation failed." >&2
    exit 1
fi
cat "$work/fixtures_base.log"

echo "Building fixtures: real-depth metagenome MAG reads (pks island, sized to stay unbinned)"
if ! python3 "$repo_dir/tests/fixtures/generate_metagenome_mag_fixtures.py" --out-dir "$fixtures" > "$work/fixtures_mag.log" 2>&1; then
    cat "$work/fixtures_mag.log" >&2
    echo "ERROR: metagenome MAG fixture generation failed." >&2
    exit 1
fi
cat "$work/fixtures_mag.log"

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

MAG_SAMPLE=e2e_metagenome_mag

printf 'patient,fastq1,fastq2\n' > "$work/sheet.csv"
printf '%s,%s,%s\n' "$MAG_SAMPLE" "$fixtures/metagenome_mag_R1.fastq.gz" "$fixtures/metagenome_mag_R2.fastq.gz" >> "$work/sheet.csv"

cat > "$work/metagenome_mag_local.config" <<'EOF'
executor {
    cpus   = 16
    memory = '64 GB'
}
process {
    withLabel:extract_reads   { cpus = 1 }
    withLabel:filter_reads    { cpus = 1 }
    withLabel:map_reads       { cpus = 1 }
    withLabel:pks_align       { cpus = 1 }
    withLabel:pks_hmm         { cpus = 1 }
    withLabel:process_low     { cpus = 1 }
    withLabel:sample_stage    { cpus = 1 }
    withLabel:mag_assembly    { cpus = 2 }
    withLabel:mag_binning     { cpus = 1 }
    withLabel:mag_gtdbtk      { cpus = 4 }
    withLabel:mag_prophage    { cpus = 2 }
    withLabel:mag_hmm         { cpus = 1 }
}
EOF

echo "Running the real pipeline: metagenome, bowtie2, --enable_mags true"
nf_log="$work/nextflow.log"
(
    cd "$work" && nextflow run "$repo_dir/main.nf" \
        -work-dir "$work/work" \
        -c "$work/metagenome_mag_local.config" \
        --sample "$work/sheet.csv" \
        --input_data_type fastq \
        --sample_type metagenome \
        --profiling_method bowtie2 \
        --hg38_db "$fixtures/hg38.mmi" \
        --t2t_phix_db "$fixtures/t2t.mmi" \
        --enable_mags true \
        --gtdbtk_db "$GTDBTK_DB" \
        --checkm2_db "$CHECKM2_DB" \
        --genomad_db "$GENOMAD_DB" \
        --save_intermediates true \
        --conda_cache_dir "$CONDA_CACHE_DIR" \
        --outdir "$work/results"
) > "$nf_log" 2>&1
nf_status=$?

if [[ $nf_status -ne 0 ]]; then
    echo "ERROR: the real pipeline run failed (exit $nf_status). Last 100 lines:" >&2
    tail -100 "$nf_log" >&2
    echo "Full log: $nf_log" >&2
    KEEP_WORK=1
    exit 1
fi

echo "Pipeline finished. Checking real output values."
if ! python3 "$repo_dir/tests/fixtures/assert_metagenome_mag_outputs.py" \
        --results "$work/results" \
        --mag-sample "$MAG_SAMPLE"; then
    echo "ERROR: metagenome MAG assertions failed. Results kept at: $work/results" >&2
    KEEP_WORK=1
    exit 1
fi

echo
echo "metagenome MAG (zero-bin / unbinned-pool) execution check passed"
