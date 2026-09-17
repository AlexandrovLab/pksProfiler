nextflow.enable.dsl = 2

process krakenPrefilter {
    tag "$sampleID"
    label 'process_high_disk'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/prefilter" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.prefilter_env}"
    input:
    tuple val(sampleID), path(reads)
    output:
    tuple val(sampleID), path("${sampleID}.primary.fastq.gz"), emit: primary
    tuple val(sampleID), path("${sampleID}.non_target.fastq.gz"), emit: non_target
    tuple val(sampleID), path("${sampleID}.krakenuniq.report.txt"), path("${sampleID}.krakenuniq.output.txt"), emit: taxonomy
    tuple val(sampleID), path("${sampleID}.prefilter.qc.tsv"), emit: qc
    script:
    """
    set -euo pipefail
    NODES="${params.kraken_db}/taxonomy/nodes.dmp"
    if [[ ! -s "\$NODES" ]]; then
      NODES="${params.kraken_db}/taxDB"
    fi
    test -s "\$NODES" || { echo "ERROR: Missing taxonomy parents (taxonomy/nodes.dmp or taxDB) in ${params.kraken_db}" >&2; exit 1; }
    gzip -dc "${reads}" > "${sampleID}.non_host.fastq"
    krakenuniq --db "${params.kraken_db}" --threads "${task.cpus}" \
      --report-file "${sampleID}.krakenuniq.report.txt" \
      --output "${sampleID}.krakenuniq.output.txt" "${sampleID}.non_host.fastq"
    python "${params.scripts}/route_kraken_reads.py" \
      --reads "${sampleID}.non_host.fastq" --kraken-output "${sampleID}.krakenuniq.output.txt" \
      --nodes "\$NODES" --target-taxid "${params.primary_taxid}" \
      --primary "${sampleID}.primary.fastq.gz" --non-target "${sampleID}.non_target.fastq.gz" \
      --qc "${sampleID}.prefilter.qc.tsv" --sample "${sampleID}"
    """
}

process buildClbDiamondDb {
    label 'process_low'
    scratch true
    conda "${params.prefilter_env}"
    input:
    path clb_proteins
    output:
    path "clb_reference.dmnd"
    script:
    """
    set -euo pipefail
    diamond makedb --in "${clb_proteins}" --db clb_reference
    """
}

process diamondRescue {
    tag "$sampleID"
    label 'process_low'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/prefilter" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.prefilter_env}"
    input:
    tuple val(sampleID), path(reads)
    path clb_db
    output:
    tuple val(sampleID), path("${sampleID}.rescued.fastq.gz"), emit: reads
    tuple val(sampleID), path("${sampleID}.diamond.tsv"), emit: matches
    tuple val(sampleID), path("${sampleID}.diamond.qc.tsv"), emit: qc
    script:
    """
    set -euo pipefail
    if [[ \$(gzip -dc "${reads}" | wc -l) -eq 0 ]]; then
      : > "${sampleID}.diamond.tsv"
    else
      gzip -dc "${reads}" > "${sampleID}.non_target.fastq"
      diamond blastx --db "${clb_db}" --query "${sampleID}.non_target.fastq" \
        --threads "${task.cpus}" --evalue "${params.diamond_evalue}" \
        --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
        --out "${sampleID}.diamond.tsv"
    fi
    python "${params.scripts}/select_diamond_reads.py" \
      --reads "${reads}" --diamond "${sampleID}.diamond.tsv" \
      --output "${sampleID}.rescued.fastq.gz" --qc "${sampleID}.diamond.qc.tsv" \
      --sample "${sampleID}" --max-evalue "${params.diamond_evalue}" \
      --min-aa "${params.diamond_min_aa}" \
      --min-query-coverage "${params.diamond_min_query_coverage}"
    """
}

process mergePksCandidates {
    tag "$sampleID"
    label 'process_low'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/prefilter" }, mode: 'copy', enabled: params.save_intermediates, saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.fastp_env}"
    input:
    tuple val(sampleID), path(primary), path(rescued)
    output:
    tuple val(sampleID), path("${sampleID}.pks_candidates.fastq.gz"), emit: reads
    tuple val(sampleID), path("${sampleID}.candidate_merge.qc.tsv"), emit: qc
    script:
    """
    set -euo pipefail
    { gzip -dc "${primary}"; gzip -dc "${rescued}"; } | gzip -c > "${sampleID}.pks_candidates.fastq.gz"
    PRIMARY=\$(gzip -dc "${primary}" | awk 'END { print int(NR / 4) }')
    RESCUED=\$(gzip -dc "${rescued}" | awk 'END { print int(NR / 4) }')
    TOTAL=\$(gzip -dc "${sampleID}.pks_candidates.fastq.gz" | awk 'END { print int(NR / 4) }')
    [[ \$((PRIMARY + RESCUED)) -eq "\$TOTAL" ]] || { echo "ERROR: Candidate merge count mismatch" >&2; exit 1; }
    printf "Sample\tMetric\tValue\n" > "${sampleID}.candidate_merge.qc.tsv"
    printf "%s\tcandidate_primary_reads\t%s\n" "${sampleID}" "\$PRIMARY" >> "${sampleID}.candidate_merge.qc.tsv"
    printf "%s\tcandidate_rescued_reads\t%s\n" "${sampleID}" "\$RESCUED" >> "${sampleID}.candidate_merge.qc.tsv"
    printf "%s\tcandidate_total_reads\t%s\n" "${sampleID}" "\$TOTAL" >> "${sampleID}.candidate_merge.qc.tsv"
    """
}

process sampleBracken {
    tag "$sampleID"
    label 'process_high_disk'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/taxonomy" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.krakenuniq_bracken_env}"
    input:
    tuple val(sampleID), path(kraken_report), path(kraken_output)
    output:
    tuple val(sampleID), path("${sampleID}.bracken.G.report.txt"), path("${sampleID}.bracken.S.report.txt"), emit: reports
    tuple val(sampleID), path("${sampleID}.bracken.G.krakenreport.txt"), path("${sampleID}.bracken.S.krakenreport.txt"), emit: kraken_reports
    script:
    """
    set -euo pipefail
    for LEVEL in G S; do
      bracken -d "${params.kraken_db}" -i "${kraken_report}" \
        -o "${sampleID}.bracken.\${LEVEL}.report.txt" \
        -w "${sampleID}.bracken.\${LEVEL}.krakenreport.txt" \
        -r "${params.bracken_read_length}" -l "\$LEVEL" -t 2
    done
    """
}
