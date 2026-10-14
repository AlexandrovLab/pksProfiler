#!/usr/bin/env bash
# rc2 regression: the v1.0.0-rc1 Mayo metagenome run (--enable_mags) crashed ~37 h in with
#   Invalid method invocation 'call' with arguments: [<sample>, null, [bin.18.locus_evidence.tsv, ...]]
# after gtdbtkClassify failed (exit 1, errorStrategy 'ignore') for some samples. This drives
# pksMAG's real channel logic (Modules/pks_mag.nf: magQualityFallback, magSummaryInput,
# magCommunityInput) with fabricated channels in that exact shape -- seconds, local executor,
# no conda, no data -- and checks that every sample is kept (or, on the raw rc1 wiring,
# skipped with a warning) and that the failed tool is named per bin.
#
#     bash tests/check_mag_channel_regression.sh
set -uo pipefail
repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
command -v nextflow >/dev/null 2>&1 || { echo "ERROR: nextflow is required" >&2; exit 1; }

work=$(mktemp -d "${TMPDIR:-/tmp}/pks-mag-channel.XXXXXX")
cleanup() { [[ "${KEEP_WORK:-0}" == "1" ]] || rm -rf "$work"; }
trap cleanup EXIT

# Nextflow loads lib/ (ToolVersions etc., which Modules/pks_mag.nf needs to compile) from the
# launched script's own directory, so launch from a scratch project that links the repo's lib/.
mkdir -p "$work/proj"
ln -s "$repo_dir/lib" "$work/proj/lib"
sed "s#'../../Modules/#'$repo_dir/Modules/#" "$repo_dir/tests/fixtures/mag_channel_regression.nf" > "$work/proj/main.nf"
( cd "$work" && nextflow -log "$work/nextflow.log" run "$work/proj/main.nf" \
      -work-dir "$work/work" --outdir "$work/results" ) > "$work/run.log" 2>&1
status=$?
if [[ $status -ne 0 ]]; then
    tail -60 "$work/run.log" >&2
    echo "ERROR: channel regression run failed (exit $status) -- the rc1 crash shape is back?" >&2
    KEEP_WORK=1; exit 1
fi

rows="$work/results/channel_rows.tsv"
fail=0
expect() {   # expect <description> <exact line>
    if ! grep -qxF -- "$2" "$rows"; then echo "FAIL: $1"$'\n'"  missing line: $2" >&2; fail=1; else echo "ok: $1"; fi
}
expect_absent() {
    if grep -q -- "$2" "$rows"; then echo "FAIL: $1"$'\n'"  unexpected: $(grep -- "$2" "$rows")" >&2; fail=1; else echo "ok: $1"; fi
}
T=$'\t'
expect "sample with both tools reaches magSummaryTable" "summary${T}s_ok${T}quality_report.tsv${T}gtdbtk.bac120.summary.tsv${T}2${T}1${T}2"
expect "GTDB-Tk failure: sample still reaches magSummaryTable, with the stand-in" "summary${T}s_gtdb_fail${T}quality_report.tsv${T}failed.gtdbtk.bac120.summary.tsv${T}2${T}0${T}2"
expect "CheckM2 failure: sample still reaches magSummaryTable, with the stand-in" "summary${T}s_checkm_fail${T}failed.checkm2_quality_report.tsv${T}gtdbtk.bac120.summary.tsv${T}1${T}0${T}1"
expect "binless sample unchanged (binless stub path)" "summary${T}s_binless${T}stub.checkm2.tsv${T}stub.gtdbtk.tsv${T}0${T}0${T}1"
expect "GTDB-Tk stand-in names every bin gtdbtk_failed" "gtdbtk_stub${T}s_gtdb_fail${T}bin.18${T}gtdbtk_failed;bin.19${T}gtdbtk_failed"
expect "CheckM2 stand-in is header-only (every field NA), marked by its Status column" "checkm2_stub${T}s_checkm_fail${T}"
expect "community input keeps the GTDB-Tk-failed sample via the stand-in" "community${T}s_gtdb_fail${T}failed.gtdbtk.bac120.summary.tsv${T}2${T}2${T}2${T}2"
expect "community input for a healthy sample" "community${T}s_ok${T}gtdbtk.bac120.summary.tsv${T}2${T}2${T}2${T}2"
expect "raw rc1 wiring: healthy samples still flow" "community_raw${T}s_ok${T}6"
expect_absent "raw rc1 wiring: the left-orphaned GTDB-Tk-failed row is skipped, not passed on" "community_raw${T}s_gtdb_fail"
expect_absent "no sentinel row escapes" "__stubBinlessMagQuality"
if ! grep -q "communityProphageSummary skipped for s_gtdb_fail" "$work/nextflow.log"; then
    echo "FAIL: the skipped orphan row was not reported with a warning" >&2; fail=1
else echo "ok: the skipped orphan row is reported with a warning"; fi
if ! grep -q "gtdbtk produced no output for s_gtdb_fail" "$work/nextflow.log"; then
    echo "FAIL: the GTDB-Tk failure was not reported with a warning" >&2; fail=1
else echo "ok: the GTDB-Tk failure is reported with a warning"; fi

if [[ $fail -ne 0 ]]; then
    echo "--- channel_rows.tsv ---" >&2; cat "$rows" >&2; KEEP_WORK=1; exit 1
fi
echo "MAG channel regression (failed GTDB-Tk / CheckM2) check passed"
