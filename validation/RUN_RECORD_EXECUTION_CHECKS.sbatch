#!/bin/bash
#SBATCH --job-name=record_checks
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Write tests/LOCAL_EXECUTION_CHECKS.tsv for the tree as it stands: the execution
# evidence a GitHub runner cannot produce, recorded against this exact source digest.
#
# v3 review item 16. The record currently carries a waiver whose text is false for the
# tree it signs -- the MAG module changed inside the digest scope after the waived run.
# tests/check_metagenome_mag_execution.sh has since PASSED on this tree (job 14254810,
# 1h32m, all assertions green), so the honest fix is to record that pass and delete the
# waiver rather than reword it.
#
# This re-runs the checks; it does not read the earlier result. Budget a couple of hours:
# the MAG check alone took 1h32m, and the zero-skip unit run needs every real tool.
#
# It must NOT run on a login node -- the script says so itself.

set -uo pipefail

# ~/.bashrc puts tools/minimap2 (2.24) ahead of conda and exports LD_LIBRARY_PATH
# (ascat's libcrypto), which breaks /usr/bin/python3's hashlib. Re-exec clean.
if [[ "${RECORD_CLEAN_ENV:-}" != "1" ]]; then
    exec env -i \
        RECORD_CLEAN_ENV=1 REPO="${REPO:-}" \
        HOME="$HOME" USER="$USER" LOGNAME="${LOGNAME:-$USER}" \
        PATH=/usr/local/bin:/usr/bin:/bin \
        TERM="${TERM:-dumb}" \
        $(env | grep -E '^SLURM_' | sed 's/^/ /') \
        bash "$0" "$@"
fi

# Override with REPO=<path> sbatch ... to record a different worktree. The release
# record is digest-bound, so after merging to main it must be regenerated in the
# worktree that has main checked out, not here.
REPO=${REPO:-/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/Tool/latest/.worktrees/rc4}
CONDA_BASE=/tscc/projects/ps-lalexandrov/amabbasi/anaconda3

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
export JAVA_HOME="$CONDA_BASE/envs/nf-env/lib/jvm"
export PATH="$JAVA_HOME/bin:$PATH"
export NXF_OPTS="-Xms1g -Xmx8g"

export DB_ROOT=/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs
export CONDA_CACHE_DIR=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache
export GRCH38=/tscc/projects/ps-lalexandrov/shared/Reference_Genomes/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
export TMPDIR="${SLURM_TMPDIR:-/scratch/$USER/job_${SLURM_JOB_ID:-$$}}"
mkdir -p "$TMPDIR" 2>/dev/null || export TMPDIR=/tmp

# nextflow by name, no python3 on PATH that could shadow a task's pinned one.
mkdir -p "$TMPDIR/bin"
ln -sf "$CONDA_BASE/envs/nf-env/bin/nextflow" "$TMPDIR/bin/nextflow"
export PATH="$TMPDIR/bin:$PATH"

# assert_required_coverage.py runs the unit suite and REFUSES to call it green if a
# consequential test skipped. Eight test_low_score_extraction tests and three
# test_mate_loss_audit tests execute a real samtools; under a bare env -i PATH they all
# skip, and the first attempt was rejected with "16 skipped, budget is 1".
#
# pks_env supplies samtools 1.22, Rscript, bedtools and featureCounts. It goes at the END
# of PATH on purpose: python3 must keep resolving to /usr/bin/python3, whose user
# site-packages carry pandas and matplotlib, and a pks_env python ahead of it would not.
# The tests look for samtools "beside python OR on PATH", so being on PATH at all is
# enough. Nextflow prepends each task's own environment when it activates it, so this
# sits behind the pinned envs rather than shadowing them -- which is why the MAG check
# still passed under --tool_version_check strict.
export PATH="$PATH:$CONDA_BASE/envs/pks_env/bin"

cd "$REPO" || exit 1
echo "host    : $(hostname)"
echo "commit  : $(git rev-parse --short HEAD)"
echo "started : $(date -Is)"
echo

bash tests/record_local_execution_checks.sh
status=$?

echo
echo "finished: $(date -Is)"
if [[ $status -eq 0 ]]; then
    echo "Review tests/LOCAL_EXECUTION_CHECKS.tsv, then commit it. A failing check is"
    echo "recorded as 'fail' rather than hidden -- read it before committing."
fi
exit $status
