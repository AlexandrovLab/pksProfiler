process mapReads {

    scratch true
    label 'map_reads'
    publishDir("${params.mapped_reads_dir}", mode: 'copy')
    conda "${params.minimap2_env}"

    input:
    tuple val(sampleID), path(reads_fastq)

    output:
    tuple val(sampleID), path("${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz")

    script:
    """
    set -euo pipefail

    if ! gzip -t "${reads_fastq}" >/dev/null 2>&1; then
        echo "ERROR: Corrupt gzip input for ${sampleID}: ${reads_fastq}" >&2
        exit 1
    fi

    echo "Running host depletion"

    minimap2 -2 -ax sr -t "${task.cpus}" \
        "${params.hg38_db}" \
        "${reads_fastq}" |
    samtools fastq -@ "${task.cpus}" -f 4 -F 2304 |
    gzip -c > "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz"

    minimap2 -2 -ax sr -t "${task.cpus}" \
        "${params.t2t_phix_db}" \
        "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz" |
    samtools fastq -@ "${task.cpus}" -f 4 -F 2304 |
    gzip -c > "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"
    """
}

