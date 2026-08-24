nextflow.enable.dsl=2

process extractReads {
    scratch true
    label 'extract_reads'
    publishDir("${params.unmapped_bam_dir}", mode: 'copy')
    conda "${params.samtools_env}"

    input:
	tuple val(sampleID), path(bam)

    output:
	tuple val(sampleID),
      path("${sampleID}.R1.UNMAPPED.fastq.gz"),
      path("${sampleID}.R2.UNMAPPED.fastq.gz")
    
	script:
	"""
	set -euo pipefail

	R1="${sampleID}.R1.UNMAPPED.fastq.gz"
	R2="${sampleID}.R2.UNMAPPED.fastq.gz"

	samtools view -f 4 -O BAM "${bam}" | samtools bam2fq \
	-1 "\$R1" \
	-2 "\$R2"
	"""
}