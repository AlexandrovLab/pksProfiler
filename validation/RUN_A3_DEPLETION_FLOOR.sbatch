#!/bin/bash
#SBATCH --job-name=a3_depletion_floor
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Tier A 3 -- benchmark the host-depletion removal floor.
#
# Ludmil's v3 review: the pipeline keeps only minimap2-unmapped reads
# (`samtools fastq -f 4`), so any read minimap2 aligns at all is discarded, with no
# pksProfiler-specific benchmark of MAPQ, identity, aligned fraction or alignment
# score behind that choice.
#
# This job simulates reads of known origin, pushes them through the pipeline's own
# minimap2 depletion command against the real human index, and tabulates what the
# floor would remove, per class, per metric.
#
# ONLY public reference sequence and reads simulated from it. No sample file is
# opened, and the output is counts and distributions.
#
# Classes:
#   pks_identical / _2pct / _5pct  reads from the pks island; any removal here is a
#                                  FALSE removal and a direct sensitivity loss
#   related_bacteria               E. coli IHE3034 backbone, the carrier organism
#   low_complexity                 homopolymer / short-repeat, the spurious class
#   true_human                     GRCh38 primary assembly; SHOULD be removed
#
# platinum + hcp-ddp302 because the human reference FASTA lives on the restricted
# Lustre tree, which condo nodes may not mount.

set -euo pipefail

# A clean PATH: ~/.bashrc prepends tools/minimap2 (2.24) ahead of conda, and exports
# LD_LIBRARY_PATH / PYTHONPATH that break unrelated binaries. Everything here is
# called by absolute path, so the minimal PATH is all that is needed.
export PATH=/usr/local/bin:/usr/bin:/bin
unset LD_LIBRARY_PATH LIBRARY_PATH CPATH PYTHONPATH || true

REPO=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/Tool/latest/.worktrees/rc4
DBS=/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs

# minimap2 2.28 -- the version the pipeline pins, NOT the 2.24 on the login PATH.
MINIMAP2=/tscc/projects/ps-lalexandrov/amabbasi/anaconda3/envs/minimap2_env/bin/minimap2

# The two depletion stages EVERY run applies, in pipeline order (map_reads.nf: the
# human pass at line 85, then the T2T+phiX pass at line 103). A run with
# --pangenome_db set applies the same rule once more per pangenome index; there are
# 96 of them in shared/pangenome_individual_mmi, hence 98 stages in that configuration.
HG38_INDEX=$DBS/human-GRC-db.mmi
T2T_INDEX=$DBS/human-GCA-phix-db.mmi
TOTAL_PIPELINE_STAGES=98

HUMAN_FASTA=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/tools/vep_cache/Homo_sapiens.GRCh38.dna.primary_assembly.fa
PKS_FASTA=$REPO/indices/pks_island/pks_island.fna
BACTERIA_FASTA=$REPO/indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic.fna

OUTDIR=$REPO/validation
SCRATCH=${SLURM_TMPDIR:-/scratch/$USER/job_$SLURM_JOB_ID}
mkdir -p "$OUTDIR"
mkdir -p "$SCRATCH" 2>/dev/null || SCRATCH=$OUTDIR

echo "host        : $(hostname)"
echo "minimap2    : $($MINIMAP2 --version)"
echo "stage 1     : $HG38_INDEX"
echo "stage 2     : $T2T_INDEX"
echo "started     : $(date -Is)"
echo

for f in "$HG38_INDEX" "$T2T_INDEX" "$HUMAN_FASTA" "$PKS_FASTA" "$BACTERIA_FASTA" "$MINIMAP2"; do
    [ -r "$f" ] || { echo "MISSING OR UNREADABLE: $f" >&2; exit 1; }
done

/usr/bin/python3 -E "$REPO/scripts/a3_depletion_floor_benchmark.py" \
    --index          "$HG38_INDEX" \
    --index          "$T2T_INDEX" \
    --human-fasta    "$HUMAN_FASTA" \
    --pks-fasta      "$PKS_FASTA" \
    --bacteria-fasta "$BACTERIA_FASTA" \
    --preset sr \
    --reads 5000 \
    --total-pipeline-stages "$TOTAL_PIPELINE_STAGES" \
    --threads "${SLURM_CPUS_PER_TASK:-2}" \
    --minimap2 "$MINIMAP2" \
    --tmpdir "$SCRATCH" \
    --out "$OUTDIR/a3_depletion_floor.tsv"

echo
echo "finished    : $(date -Is)"
echo "result      : $OUTDIR/a3_depletion_floor.tsv"
