#!/bin/bash
#SBATCH --job-name=item16_mag
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

# v3 item 16: re-run the MAG execution check on THIS tree, so the release record can
# carry a real PASS instead of a waiver.
#
# What the review found: the current waiver states, in substance, that nothing the MAG
# check exercises changed since the earlier run. That is false for the exact tree the
# record signs -- the recorded source digest matches the released tree, that tree
# contains the 6 Nov locus-evidence join fix in Modules/pks_mag.nf, and Modules/ is
# inside the digest scope. A verified false statement in a signed release artifact.
#
# Item 17's machinery (already committed) stops this recurring: a waiver must now name
# the commit it was reasoned against and expires when digest-scope files change. But it
# does not repair the existing one, and the verifier will now reject it. Re-running is
# the clean fix -- evidence instead of a better-worded excuse.
#
# The check is already a SINGLE-SAMPLE test: one sample, ~1010 read pairs, an assembly
# deliberately sized under MetaBAT2's 200 kb floor so --unbinned pools the whole thing.
# That zero-bin case is the regression it guards (genomes/pks_mag_summary.tsv used to be
# silently skipped for such a sample). Reads are tiled from the repo's own shipped
# GCF_000025745.1 (IHE3034) reference across all 19 clb genes, against two tiny synthetic
# host references -- public reference sequence and synthetic reads only, no sample data.
#
# Cost is dominated by the real databases it must read: GTDB-Tk r232 (94 GB), CheckM2
# and geNomad. 128 GB and 8 CPUs covers the mag_gtdbtk label (4 cpus) with headroom.
#
# On success, record the result with tests/record_local_execution_checks.sh so the
# release record carries the PASS, and drop the old waiver.

set -euo pipefail

# ~/.bashrc puts tools/minimap2 (2.24) ahead of conda and exports LD_LIBRARY_PATH
# (ascat's libcrypto), LIBRARY_PATH, CPATH and PYTHONPATH; sbatch copies the submitting
# environment in, and Nextflow copies it again into every task. Re-exec under env -i
# keeping only SLURM_* and the basics.
if [[ "${ITEM16_CLEAN_ENV:-}" != "1" ]]; then
    exec env -i \
        ITEM16_CLEAN_ENV=1 \
        HOME="$HOME" USER="$USER" LOGNAME="${LOGNAME:-$USER}" \
        PATH=/usr/local/bin:/usr/bin:/bin \
        TERM="${TERM:-dumb}" \
        $(env | grep -E '^SLURM_' | sed 's/^/ /') \
        bash "$0" "$@"
fi

REPO=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/Tool/latest/.worktrees/rc4
CONDA_BASE=/tscc/projects/ps-lalexandrov/amabbasi/anaconda3

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

export JAVA_HOME="$CONDA_BASE/envs/nf-env/lib/jvm"
export PATH="$JAVA_HOME/bin:$PATH"
export NXF_OPTS="-Xms1g -Xmx8g"

export CONDA_CACHE_DIR=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache
export DB_ROOT=/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs
export GRCH38=/tscc/projects/ps-lalexandrov/shared/Reference_Genomes/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
export TMPDIR="${SLURM_TMPDIR:-/scratch/$USER/job_${SLURM_JOB_ID:-$$}}"
mkdir -p "$TMPDIR" 2>/dev/null || export TMPDIR=/tmp

# The check calls `nextflow` and `python3` by name.
mkdir -p "$TMPDIR/bin"
ln -sf "$CONDA_BASE/envs/nf-env/bin/nextflow" "$TMPDIR/bin/nextflow"
ln -sf /usr/bin/python3 "$TMPDIR/bin/python3"
export PATH="$TMPDIR/bin:$PATH"

echo "host     : $(hostname)"
echo "commit   : $(cd "$REPO" && git rev-parse --short HEAD)"
echo "nextflow : $(nextflow -v 2>&1 | head -1)"
echo "DB_ROOT  : $DB_ROOT"
echo "started  : $(date -Is)"
echo

bash "$REPO/tests/check_metagenome_mag_execution.sh"
status=$?

echo
echo "finished : $(date -Is)"
if [[ $status -eq 0 ]]; then
    echo "PASS. Record it so the release record carries the result rather than a waiver:"
    echo "  cd $REPO && bash tests/record_local_execution_checks.sh"
else
    echo "FAILED (exit $status) -- do NOT record a pass."
fi
exit $status
