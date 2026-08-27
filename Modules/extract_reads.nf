nextflow.enable.dsl=2

process extractReads {
    scratch true
    label 'extract_reads'
    publishDir(
        "${params.unmapped_bam_dir}",
        mode: 'copy',
        enabled: params.save_intermediates
    )
    conda "${params.samtools_env}"

    input:
	tuple val(sampleID), path(bam)

    output:
	tuple val(sampleID), path("${sampleID}.UNMAPPED.fastq.gz"), emit: reads
	tuple val(sampleID), path("${sampleID}.extract.qc.tsv"), emit: qc
    
	script:
	"""
	set -euo pipefail

	READS="${sampleID}.UNMAPPED.fastq.gz"
	QC="${sampleID}.extract.qc.tsv"

	INPUT_READS=\$(samtools view -c -F 2304 "${bam}")
	UNMAPPED_READS=\$(samtools view -c -f 4 -F 2304 "${bam}")

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

	printf "Sample\tMetric\tValue\n" > "\$QC"
	printf "%s\tbam_input_primary_records\t%s\n" "${sampleID}" "\$INPUT_READS" >> "\$QC"
	printf "%s\textracted_unmapped_reads\t%s\n" "${sampleID}" "\$UNMAPPED_READS" >> "\$QC"
	"""
}
