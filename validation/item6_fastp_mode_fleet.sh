#!/bin/bash
# Item 6 fleet: the single-end vs paired fastp comparison across 50 public samples.
#
#   bash item6_fastp_mode_fleet.sh select    pick 50 pairs across cohorts -> manifest
#   bash item6_fastp_mode_fleet.sh submit    one SLURM job per selected pair
#   bash item6_fastp_mode_fleet.sh status    queue + parts progress
#   bash item6_fastp_mode_fleet.sh collect   one row per sample + the verdict numbers
#
# Samples come from the public infant-metagenome corpus (SRA accessions fetched from
# ENA/SRA; public data, fine to process directly). Five pairs per cohort, spread over
# ten cohorts, because adapter read-through is a property of the library prep: one
# cohort cannot clear the others.
#
# Per sample it runs validation/item6_fastp_mode_comparison.sh -- arm A is the
# pipeline's exact single-end invocation, arm B native paired with orphan retention,
# same pinned fastp 1.3.7 -- and keeps only the counts.

set -euo pipefail

REPO=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/Tool/latest/.worktrees/rc4
CORPUS=/tscc/lustre/restricted/alexandrov-ddn/users/shared/infant_metagenomes
COMPARE=$REPO/validation/item6_fastp_mode_comparison.sh
MANIFEST=$REPO/validation/item6_fleet_manifest.tsv
PARTS=$REPO/validation/item6_parts
MERGED=$REPO/validation/item6_fastp_mode.50samples.tsv

COHORTS=(Asnicar_2017 Backhed_2015 Baumann_Dudenhoeffer_2018 Bonham_2023 Brooks_2017
         Busi_2021 Casaburi_2021 Ferretti_2018 Gasparrini_2019 Gibson_2016)
PER_COHORT=5
MIN_BYTES=100000000     # <100 MB: too shallow to say anything
MAX_BYTES=3000000000    # >3 GB: wastes queue time without adding information

case "${1:-}" in
  select)
    : > "$MANIFEST"
    total=0
    for cohort in "${COHORTS[@]}"; do
        [ -d "$CORPUS/$cohort" ] || { echo "  $cohort: missing, skipped"; continue; }
        picked=0
        # Deterministic spread: name-sorted, stepped through the size-filtered list.
        mapfile -t ones < <(find "$CORPUS/$cohort" -maxdepth 1 -name "*_1.fastq.gz" \
                                -size +"$((MIN_BYTES / 1024))"k -size -"$((MAX_BYTES / 1024))"k \
                                2>/dev/null | sort)
        [ "${#ones[@]}" -eq 0 ] && { echo "  $cohort: no pairs in size window"; continue; }
        step=$(( ${#ones[@]} / PER_COHORT )); [ "$step" -lt 1 ] && step=1
        for (( i = 0; i < ${#ones[@]} && picked < PER_COHORT; i += step )); do
            r1=${ones[$i]}
            r2=${r1%_1.fastq.gz}_2.fastq.gz
            [ -r "$r2" ] || continue
            printf '%s\t%s\t%s\n' "$cohort" "$r1" "$r2" >> "$MANIFEST"
            picked=$((picked + 1)); total=$((total + 1))
        done
        echo "  $cohort: $picked"
    done
    echo "selected $total pairs -> $MANIFEST"
    ;;

  submit)
    [ -s "$MANIFEST" ] || { echo "run 'select' first"; exit 1; }
    mkdir -p "$PARTS"
    submitted=0 skipped=0
    while IFS=$'\t' read -r cohort r1 r2; do
        acc=$(basename "$r1" _1.fastq.gz)
        part=$PARTS/${cohort}.${acc}.tsv
        [ -s "$part" ] && { skipped=$((skipped + 1)); continue; }
        sbatch --account=ddp302 --partition=platinum --qos=hcp-ddp302 \
            --cpus-per-task=2 --mem=8G --time=03:00:00 \
            --job-name=item6_fastp --output="$PARTS/%x_%j.log" \
            --wrap "W=/scratch/\$USER/job_\$SLURM_JOB_ID/item6; THREADS=2 bash '$COMPARE' '$r1' '$r2' \"\$W\" '$part' && rm -rf \"\$W\"" >/dev/null
        submitted=$((submitted + 1))
    done < "$MANIFEST"
    echo "submitted $submitted, skipped $skipped already done"
    ;;

  status)
    echo "queued/running: $(squeue -u "$USER" -n item6_fastp -h 2>/dev/null | wc -l)"
    echo "parts done:     $(ls "$PARTS"/*.tsv 2>/dev/null | wc -l) of $(wc -l < "$MANIFEST" 2>/dev/null || echo '?')"
    ;;

  collect)
    ls "$PARTS"/*.tsv >/dev/null 2>&1 || { echo "no parts yet"; exit 1; }
    {
      printf 'cohort\tsample\treads_in\tse_surviving\tpe_surviving\tse_adapter_reads\tpe_adapter_reads\tonly_se\tonly_pe\tdisagree_fraction\tpe_adapter_fraction\n'
      for part in "$PARTS"/*.tsv; do
        base=$(basename "$part" .tsv); cohort=${base%%.*}; acc=${base#*.}
        awk -F'\t' -v cohort="$cohort" -v acc="$acc" '
            $1 == "reads_in"                 { reads=$2 }
            $1 == "reads_surviving"          { se=$2; pe=$3 }
            $1 == "adapter_trimmed_reads"    { sea=$2; pea=$3 }
            $1 == "survivors_only_this_mode" { ose=$2; ope=$3 }
            END { if (reads > 0) printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.6f\t%.6f\n",
                      cohort, acc, reads, se, pe, sea, pea, ose, ope,
                      (ose + ope) / reads, pea / reads }' "$part"
      done
    } > "$MERGED"
    echo "collected $(( $(wc -l < "$MERGED") - 1 )) samples -> $MERGED"
    echo
    awk -F'\t' 'NR > 1 {
            n++; d=$10; a=$11
            if (d > dmax) { dmax=d; dwho=$1 }
            if (a > amax) { amax=a; awho=$1 }
            dsum+=d; asum+=a
        }
        END {
            printf "  samples: %d\n", n
            printf "  survivor disagreement: mean %.6f, worst %.6f (%s)\n", dsum/n, dmax, dwho
            printf "  PE adapter-trimmed fraction: mean %.6f, worst %.6f (%s)\n", asum/n, amax, awho
        }' "$MERGED"
    ;;

  *)
    echo "usage: bash item6_fastp_mode_fleet.sh [select|submit|status|collect]" ;;
esac
