nextflow.enable.dsl=2

process extractReads {
    scratch true
    label 'extract_reads'
    publishDir("${params.unmapped_bam_dir}", mode: 'copy')
    conda "${params.samtools_env}"

    input:
	tuple val(sampleID), path(bam)

    output:
	tuple val(sampleID), path("${sampleID}.UNMAPPED.fastq.gz")
    
	script:
	"""
	set -euo pipefail

	READS="${sampleID}.UNMAPPED.fastq.gz"

    # Retain every primary unmapped alignment record, regardless of
    # whether its mate is mapped, unmapped, or absent from the BAM.
    samtools view \
        -@ "${task.cpus}" \
        -f 4 \
        -F 2304 \
        -u \
        "${bam}" |
    samtools fastq \
        -@ "${task.cpus}" \
        -N \
        - |
    gzip -c > "\$READS"
	"""
}
