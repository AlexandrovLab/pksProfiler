#!/bin/bash
# Item 6 (Ludmil v3): does single-end fastp actually differ from paired fastp
# on a real paired metagenome sample?
#
# Arm A replicates Modules/filter_reads.nf exactly: mate-suffix the names with
# the module's own awk, concatenate R1+R2, run fastp single-end over --stdin.
# Arm B is the prescribed fix: native paired fastp with orphan retention.
# Same binary (the pipeline-pinned fastp 1.3.7 from the conda cache), same
# -l 45 --cut_tail, same threads.
#
# Sample: a PUBLIC SRA accession (per the 2026-10-27 rule, public ENA/SRA data
# is fine to process directly). Everything written here is counts.
#
# Usage: bash item6_fastp_mode_comparison.sh R1.fastq.gz R2.fastq.gz WORKDIR OUT_TSV

set -euo pipefail

R1=${1:?R1.fastq.gz}
R2=${2:?R2.fastq.gz}
WORK=${3:?scratch workdir}
OUT=${4:?output tsv}

FASTP=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/pksProfiler/conda_cache/env-d15612242667a9471affc562777eeb05/bin/fastp
THREADS=${THREADS:-2}
mkdir -p "$WORK"

suffix() {  # the exact awk from Modules/filter_reads.nf
    gzip -dc "$1" | awk -v mate="$2" '
        NR % 4 == 1 {
            split($0, header_parts, " ")
            read_id=header_parts[1]
            if (read_id !~ /\/[12]$/) { read_id=read_id "/" mate }
            $0=read_id substr($0, length(header_parts[1]) + 1)
        }
        { print }'
}

json_num() { # first "key": N after an optional section marker
    awk -v section="$2" -v key="$3" '
        section != "" && $0 ~ "\"" section "\"" { inside=1; next }
        (section == "" || inside) && $0 ~ "\"" key "\"" {
            value=$0; gsub(/[^0-9]/, "", value); print value; exit
        }' "$1"
}

echo "fastp: $($FASTP --version 2>&1)"
echo "sample: $(basename "$R1" | sed 's/_1.*//') (public SRA)"
echo

echo "== arm A: single-end stream, exactly as filter_reads.nf =="
{ suffix "$R1" 1; suffix "$R2" 2; } | \
    "$FASTP" -l 45 --cut_tail --stdin -w "$THREADS" --compression 4 \
        --json "$WORK/se.json" --html "$WORK/se.html" -o "$WORK/se.fastq.gz" \
        2> "$WORK/se.log"

echo "== arm B: native paired, with orphan retention =="
"$FASTP" -l 45 --cut_tail -w "$THREADS" --compression 4 \
    -i "$R1" -I "$R2" \
    -o "$WORK/pe_R1.fastq.gz" -O "$WORK/pe_R2.fastq.gz" \
    --unpaired1 "$WORK/pe_u1.fastq.gz" --unpaired2 "$WORK/pe_u2.fastq.gz" \
    --json "$WORK/pe.json" --html "$WORK/pe.html" \
    2> "$WORK/pe.log"

echo "== survivor name sets (disk-based) =="
names() { gzip -dc "$1" | awk -v add="$2" 'NR % 4 == 1 { n=substr($1, 2); if (add != "" && n !~ /\/[12]$/) n = n add; print n }'; }
export LC_ALL=C
names "$WORK/se.fastq.gz" ""      | sort -T "$WORK" > "$WORK/se.names"
{ names "$WORK/pe_R1.fastq.gz" "/1"; names "$WORK/pe_R2.fastq.gz" "/2"
  names "$WORK/pe_u1.fastq.gz" "/1"; names "$WORK/pe_u2.fastq.gz" "/2"; } | sort -T "$WORK" > "$WORK/pe.names"

ONLY_SE=$(comm -23 "$WORK/se.names" "$WORK/pe.names" | wc -l)
ONLY_PE=$(comm -13 "$WORK/se.names" "$WORK/pe.names" | wc -l)
BOTH=$(comm -12 "$WORK/se.names" "$WORK/pe.names" | wc -l)

SE_IN=$(json_num "$WORK/se.json" before_filtering total_reads)
SE_OUT=$(json_num "$WORK/se.json" after_filtering total_reads)
SE_ADAPT_READS=$(json_num "$WORK/se.json" adapter_cutting adapter_trimmed_reads)
SE_ADAPT_BASES=$(json_num "$WORK/se.json" adapter_cutting adapter_trimmed_bases)
PE_IN=$(json_num "$WORK/pe.json" before_filtering total_reads)
PE_OUT=$(json_num "$WORK/pe.json" after_filtering total_reads)
PE_ADAPT_READS=$(json_num "$WORK/pe.json" adapter_cutting adapter_trimmed_reads)
PE_ADAPT_BASES=$(json_num "$WORK/pe.json" adapter_cutting adapter_trimmed_bases)
PE_ORPHANS=$(( $(gzip -dc "$WORK/pe_u1.fastq.gz" | wc -l) / 4 + $(gzip -dc "$WORK/pe_u2.fastq.gz" | wc -l) / 4 ))

{
  printf "metric\tsingle_end_as_shipped\tnative_paired\n"
  printf "reads_in\t%s\t%s\n"                "$SE_IN" "$PE_IN"
  printf "reads_surviving\t%s\t%s\n"         "$SE_OUT" "$((PE_OUT + PE_ORPHANS))"
  printf "surviving_as_pairs_or_se\t%s\t%s\n" "$SE_OUT" "$PE_OUT"
  printf "surviving_orphans\tNA\t%s\n"       "$PE_ORPHANS"
  printf "adapter_trimmed_reads\t%s\t%s\n"   "$SE_ADAPT_READS" "$PE_ADAPT_READS"
  printf "adapter_trimmed_bases\t%s\t%s\n"   "$SE_ADAPT_BASES" "$PE_ADAPT_BASES"
  printf "survivors_in_both_modes\t%s\t%s\n" "$BOTH" "$BOTH"
  printf "survivors_only_this_mode\t%s\t%s\n" "$ONLY_SE" "$ONLY_PE"
} | tee "$OUT"

echo
echo "wrote $OUT"
