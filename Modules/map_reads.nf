process mapReads {

    scratch true
    label 'map_reads'
    publishDir(
        "${params.mapped_reads_dir}",
        mode: 'copy',
        enabled: params.save_intermediates
    )
    conda "${params.minimap2_env}"

    input:
    tuple val(sampleID), path(reads_fastq)

    output:
    tuple val(sampleID), path("${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"), emit: reads
    tuple val(sampleID), path("${sampleID}.depletion.qc.tsv"), emit: qc

    script:
    """
    set -euo pipefail

    if ! gzip -t "${reads_fastq}" >/dev/null 2>&1; then
        echo "ERROR: Corrupt gzip input for ${sampleID}: ${reads_fastq}" >&2
        exit 1
    fi

    echo "Running host depletion"

    QC="${sampleID}.depletion.qc.tsv"

    minimap2 -2 -ax sr -t "${task.cpus}" \
        "${params.hg38_db}" \
        "${reads_fastq}" |
    samtools fastq -@ "${task.cpus}" -f 4 -F 2304 |
    gzip -c > "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz"

    AFTER_HG38=\$(gzip -dc "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz" | awk 'END { print int(NR / 4) }')

    minimap2 -2 -ax sr -t "${task.cpus}" \
        "${params.t2t_phix_db}" \
        "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz" |
    samtools fastq -@ "${task.cpus}" -f 4 -F 2304 |
    gzip -c > "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"

    AFTER_T2T_PHIX=\$(gzip -dc "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz" | awk 'END { print int(NR / 4) }')
    printf "Sample\tMetric\tValue\n" > "\$QC"
    printf "%s\treads_after_hg38\t%s\n" "${sampleID}" "\$AFTER_HG38" >> "\$QC"
    printf "%s\treads_after_t2t_phix\t%s\n" "${sampleID}" "\$AFTER_T2T_PHIX" >> "\$QC"
    """
}
