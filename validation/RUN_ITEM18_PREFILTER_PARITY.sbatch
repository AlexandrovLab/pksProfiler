#!/bin/bash
#SBATCH --job-name=item18_parity
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# v3 item 18: run the same synthetic metagenome fixture twice, once with
# --prefilter_mode off and once with balanced, and compare the rate.
#
# The denominator (host_depleted_reads) is taken at host depletion, upstream of the
# prefilter, so it must be identical -- the check asserts that. The question is the
# numerator: `balanced` narrows PROFILING_READS to CANDIDATE_OUT.reads, so pks_reads is
# counted on a prefiltered stream. main.nf already refuses balanced without
# --diamond_rescue true because "clb reads outside taxid <primary> are discarded"
# without it; this measures whether the rescue in fact keeps the numerator whole.
#
# 16 CPUs / 72 GB matches the prefilter_kraken label (16 / 64 GB) with headroom for the
# driver. The krakenuniq database is 544 GB on Lustre and is memory-mapped rather than
# preloaded, so the kraken step is I/O-bound; the fixture's reads are negligible.
#
# Synthetic reads over public references only. No sample data is read.

set -euo pipefail

# ~/.bashrc prepends tools/minimap2 (2.24) ahead of conda and exports LD_LIBRARY_PATH
# (ascat's libcrypto), LIBRARY_PATH, CPATH and PYTHONPATH; sbatch copies the submitting
# environment into the job, and Nextflow copies it again into every task. A stale
# CONDA_SHLVL also makes conda.sh skip adding condabin. Re-exec under env -i keeping only
# SLURM_* and the basics, exactly as the CMPipeline drivers do.
if [[ "${ITEM18_CLEAN_ENV:-}" != "1" ]]; then
    exec env -i \
        ITEM18_CLEAN_ENV=1 \
        HOME="$HOME" USER="$USER" LOGNAME="${LOGNAME:-$USER}" \
        PATH=/usr/local/bin:/usr/bin:/bin \
        TERM="${TERM:-dumb}" \
        $(env | grep -E '^SLURM_' | sed 's/^/ /') \
        bash "$0" "$@"
fi

REPO=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/Tool/latest/.worktrees/rc4
CONDA_BASE=/tscc/projects/ps-lalexandrov/amabbasi/anaconda3
DBS=/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Nextflow by absolute path, and its own JRE: compute nodes default to Java 1.8, which
# Nextflow 24.10 refuses.
export JAVA_HOME="$CONDA_BASE/envs/nf-env/lib/jvm"
export PATH="$JAVA_HOME/bin:$PATH"
NEXTFLOW="$CONDA_BASE/envs/nf-env/bin/nextflow"
export NXF_OPTS="-Xms1g -Xmx8g"

export CONDA_CACHE_DIR=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache
export KRAKEN_DB="$DBS/krakenUniq_8_8_2023"
export GRCH38=/tscc/projects/ps-lalexandrov/shared/Reference_Genomes/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
export TMPDIR="${SLURM_TMPDIR:-/scratch/$USER/job_${SLURM_JOB_ID:-$$}}"
mkdir -p "$TMPDIR" 2>/dev/null || export TMPDIR=/tmp
export PREFILTER_PARITY_OUT="$REPO/validation/item18_prefilter_rate_parity.tsv"

# The check script calls `nextflow` and `python3` by name; give it the right ones.
mkdir -p "$TMPDIR/bin"
ln -sf "$NEXTFLOW" "$TMPDIR/bin/nextflow"
ln -sf /usr/bin/python3 "$TMPDIR/bin/python3"
export PATH="$TMPDIR/bin:$PATH"

echo "host       : $(hostname)"
echo "nextflow   : $("$NEXTFLOW" -v 2>&1 | head -1)"
echo "java       : $(java -version 2>&1 | head -1)"
echo "kraken db  : $KRAKEN_DB"
echo "started    : $(date -Is)"
echo

bash "$REPO/tests/check_prefilter_rate_parity.sh"
status=$?

echo
echo "finished   : $(date -Is)"
echo "result     : $PREFILTER_PARITY_OUT"
exit $status
