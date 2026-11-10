#!/usr/bin/env bash
# Produce tests/LOCAL_EXECUTION_CHECKS.tsv: the execution-check results a GitHub runner
# cannot produce, recorded against the exact commit being released.
#
# rc2 review 2.2. The tag build runs the end-to-end, malformed-input, mode, V8-lane and
# resume checks on the runner, on the exact tagged commit. Two things cannot go there:
#
#   tests/check_metagenome_mag_execution.sh   needs GTDB-Tk r232 (~100 GB), CheckM2 and
#                                             geNomad. No stand-in database can simulate
#                                             the tool behaviour the check exists for.
#   tests/assert_required_coverage.py         a ZERO-skip unit run needs every real tool,
#                                             R included; the runner has no R.
#
# So the release gate verifies a recorded result instead, and
# tests/verify_local_execution_record.py refuses it unless the SOURCE DIGEST in the record
# matches the tree being released -- a content digest of main.nf, nextflow.config, lib/,
# Modules/, scripts/, conf/, conda_envs/, ref/ and tests/, excluding this record itself.
# Committing the record therefore does not invalidate it, so: run this on the tree you are
# about to tag, commit the file it writes (no amend needed), and tag that. Change anything
# under those paths afterwards and the record stops being accepted, which is the point.
#
#   RUN IT ON A COMPUTE NODE, NOT A LOGIN NODE. On TSCC: platinum + hcp-ddp302 (the
#   condo nodes may not mount the restricted Lustre tree).
#
#   DB_ROOT=/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs \
#   CONDA_CACHE_DIR=/path/to/private/conda_cache \
#       bash tests/record_local_execution_checks.sh
#
# Environment:
#   DB_ROOT          where the GTDB-Tk / CheckM2 / geNomad databases live   (required)
#   CONDA_CACHE_DIR  a conda cache OUTSIDE the work dir                     (required)
#   SIGNED_OFF       the sign-off line; defaults to the git user, which is enough only
#                    if that is really you
#   OUT              the record to write (default tests/LOCAL_EXECUTION_CHECKS.tsv)
#
# It does NOT fail the shell on a failing check: a failure is recorded as `fail`, and the
# verifier is what refuses to release it. Recording a failure honestly is the point.
set -uo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OUT=${OUT:-"$repo_dir/tests/LOCAL_EXECUTION_CHECKS.tsv"}
LOG_DIR=${LOG_DIR:-"$repo_dir/.local_execution_logs"}

if [[ -z "${DB_ROOT:-}" ]]; then
    echo "ERROR: set DB_ROOT to the directory holding the GTDB-Tk, CheckM2 and geNomad databases" >&2
    exit 2
fi
# SKIP_CHECKS: space-separated check paths to record as `skipped` instead of running.
# SKIP_REASON (required with it) goes into the record as a `skip_waiver:` header line,
# under the same sign-off as everything else -- the verifier refuses a skip without it,
# and refuses any skip of the unit gate. A skip is a signed decision, never a default.
if [[ -n "${SKIP_CHECKS:-}" && -z "${SKIP_REASON:-}" ]]; then
    echo "ERROR: SKIP_CHECKS without SKIP_REASON. Skipping a release check is a signed" >&2
    echo "       decision: say why, and the reason goes into the record." >&2
    exit 2
fi
if [[ " ${SKIP_CHECKS:-} " == *" tests/assert_required_coverage.py "* ]]; then
    echo "ERROR: the unit gate cannot be skipped: it runs anywhere Python runs." >&2
    exit 2
fi

if [[ -z "${CONDA_CACHE_DIR:-}" ]]; then
    echo "ERROR: set CONDA_CACHE_DIR to a conda cache outside the work directory" >&2
    exit 2
fi

commit=$(git -C "$repo_dir" rev-parse HEAD 2>/dev/null || echo "unknown")
dirty=$(git -C "$repo_dir" status --porcelain 2>/dev/null)
if [[ -n "$dirty" ]]; then
    echo "WARNING: the working tree is dirty, so the commit recorded below is not exactly" >&2
    echo "         what ran. The source digest is, and it is what the gate checks." >&2
fi

mkdir -p "$LOG_DIR"
rows=()

waiver_lines=()

run_check() {
    local name=$1 label=$2
    shift 2
    if [[ " ${SKIP_CHECKS:-} " == *" $name "* ]]; then
        echo "==> $label"
        echo "    skipped  ($SKIP_REASON)"
        rows+=("$name"$'\t'"skipped"$'\t'"0"$'\t'"0"$'\t'"-")
        # since=<commit> binds the reason to the tree it describes. verify_local_
        # execution_record.py recomputes what has changed inside the source digest since
        # that commit and expires the waiver if the code the check exercises moved. Set
        # SKIP_SINCE to the commit the reason is about -- normally the one at which the
        # check last passed, NOT HEAD, or the binding asserts nothing.
        if [[ -z "${SKIP_SINCE:-}" ]]; then
            echo "ERROR: SKIP_CHECKS needs SKIP_SINCE too: the commit whose state" >&2
            echo "       SKIP_REASON describes, usually where the check last passed." >&2
            exit 2
        fi
        waiver_lines+=("# skip_waiver: $name: since=$SKIP_SINCE $SKIP_REASON")
        return 0
    fi
    local log="$LOG_DIR/$(basename "$label").log"
    echo "==> $label"
    local start=$SECONDS
    "$@" > "$log" 2>&1
    local status=$?
    local seconds=$((SECONDS - start))
    local verdict=pass
    (( status == 0 )) || verdict=fail
    # How many assertions/tests the check reported. A check that ran and asserted nothing
    # is not evidence, and the verifier refuses a zero.
    # Two bugs lived here. `grep -c` prints 0 AND exits 1 when nothing matches, so a
    # trailing `|| echo 0` fired as well and wrote "0\n0" -- an embedded newline that made
    # the whole record unparseable, which the verifier then rejected as malformed rather
    # than reporting the real status. And the pattern was anchored at column 0, while
    # check_metagenome_mag_execution.sh prints its assertions indented ("  ok  ..."), so a
    # passing run recorded 0 assertions, which the verifier refuses by design.
    local count
    count=$(grep -cE '^[[:space:]]*(ok|PASS|OK:|\[ok\])' "$log" 2>/dev/null) || count=0
    count=${count//[^0-9]/}
    count=${count:-0}
    if [[ "$name" == "tests/assert_required_coverage.py" ]]; then
        count=$(sed -n 's/^tests run *\([0-9][0-9]*\).*/\1/p' "$log" | tail -1)
        count=${count:-0}
    fi
    rows+=("$name"$'\t'"$verdict"$'\t'"$count"$'\t'"$seconds"$'\t'"${log#$repo_dir/}")
    echo "    $verdict  ($count assertions, ${seconds}s)  -> ${log#$repo_dir/}"
}

run_check "tests/assert_required_coverage.py" "assert_required_coverage" \
    python3 "$repo_dir/tests/assert_required_coverage.py" -q

run_check "tests/check_metagenome_mag_execution.sh" "check_metagenome_mag_execution" \
    env DB_ROOT="$DB_ROOT" CONDA_CACHE_DIR="$CONDA_CACHE_DIR" \
        bash "$repo_dir/tests/check_metagenome_mag_execution.sh"

signed_off=${SIGNED_OFF:-"$(git -C "$repo_dir" config user.name) <$(git -C "$repo_dir" config user.email)> ran these on $(hostname) and confirms the results above"}

# Computed AFTER the checks ran, over the tree they ran against. The record itself is
# excluded from it, so committing the record does not invalidate it.
source_digest=$(python3 "$repo_dir/tests/verify_local_execution_record.py" --print-source-digest)

{
    echo "# The execution checks a GitHub runner cannot run, recorded against the commit"
    echo "# being released. Verified by tests/verify_local_execution_record.py, which the"
    echo "# release-gate CI job runs on every v* tag and which REFUSES a record whose"
    echo "# commit is not the commit being released."
    echo "#"
    echo "# Produced by tests/record_local_execution_checks.sh. Do not hand-edit."
    echo "#"
    echo "# commit: $commit"
    echo "# source_digest: $source_digest"
    echo "# recorded_on: $(date -u +%Y-%m-%d)"
    echo "# recorded_by: $(git -C "$repo_dir" config user.name) <$(git -C "$repo_dir" config user.email)>"
    echo "# host: $(hostname)"
    echo "# signed_off: $signed_off"
    for w in ${waiver_lines[@]+"${waiver_lines[@]}"}; do echo "$w"; done
    echo "#"
    echo "# DB_ROOT=$DB_ROOT"
    echo "# CONDA_CACHE_DIR=$CONDA_CACHE_DIR"
    echo "#"
    printf 'check\tstatus\ttests\tduration_s\tlog\n'
    printf '%s\n' "${rows[@]}"
} > "$OUT"

echo
echo "Wrote $OUT"
echo
python3 "$repo_dir/tests/verify_local_execution_record.py" --record "$OUT"
verify=$?
echo
if (( verify == 0 )); then
    echo "Now commit it, as the last commit before the tag:"
    echo "    git add ${OUT#$repo_dir/} && git commit -m 'Record the local execution checks for the release'"
    echo "No amend is needed: the record excludes itself from the source digest."
fi
exit "$verify"
