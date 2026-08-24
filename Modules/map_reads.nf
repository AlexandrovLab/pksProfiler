process mapReads {

    scratch true
    label 'mapBothReads'
    publishDir("${params.mapped_reads_dir}", mode: 'copy')
    conda "${params.minimap2_env}"

    input:
    tuple val(sampleID), path(r1_fastq), path(r2_fastq)

    output:
    tuple val(sampleID),
        path("${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"),
        path("${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz")

    script:
    """
    set -euo pipefail

    echo "Running host depletion on R1"

    minimap2 -2 -ax sr -t 16 \
        "${params.hg38_db}" \
        "${r1_fastq}" |
    samtools fastq -@ 16 -f 4 -F 256 |
    gzip > "${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz"

    minimap2 -2 -ax sr -t 16 \
        "${params.t2t_phix_db}" \
        "${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz" |
    samtools fastq -@ 16 -f 4 -F 256 |
    gzip > "${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"

    echo "Running host depletion on R2"

    minimap2 -2 -ax sr -t 16 \
        "${params.hg38_db}" \
        "${r2_fastq}" |
    samtools fastq -@ 16 -f 4 -F 256 |
    gzip > "${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz"

    minimap2 -2 -ax sr -t 16 \
        "${params.t2t_phix_db}" \
        "${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz" |
    samtools fastq -@ 16 -f 4 -F 256 |
    gzip > "${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"
    """
}


