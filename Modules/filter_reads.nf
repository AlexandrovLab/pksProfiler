nextflow.enable.dsl=2

process filterReads {
    scratch true
    label 'filter_reads'
    publishDir("${params.unmapped_bam_dir}", mode: 'copy')
    conda "${params.fastp_env}"
	maxRetries 2


	input:
	tuple val(sampleID), path(fastq1), path(fastq2)

	output:
	tuple val(sampleID),
      path("${sampleID}.R1.UNMAPPED.FASTP.FILTERED.fastq.gz"),
      path("${sampleID}.R2.UNMAPPED.FASTP.FILTERED.fastq.gz")

	script:
	"""
	set -euo pipefail

	R1="${sampleID}.R1.UNMAPPED.FASTP.FILTERED.fastq.gz"
	R2="${sampleID}.R2.UNMAPPED.FASTP.FILTERED.fastq.gz"

    if [[ -f "${params.unmapped_bam_dir}/\$R1" && -f "${params.unmapped_bam_dir}/\$R2" ]]; then
        echo "Skipping filterReads: Found existing \$R1, \$R2"
        if [[ ! -f "\$R1" ]]; then
            ln -s "${params.unmapped_bam_dir}/\$R1" . 2>/dev/null || cp "${params.unmapped_bam_dir}/\$R1" .
        fi
        if [[ ! -f "\$R2" ]]; then
            ln -s "${params.unmapped_bam_dir}/\$R2" . 2>/dev/null || cp "${params.unmapped_bam_dir}/\$R2" .
        fi
        exit 0
    fi

    fastp -l 45 --adapter_fasta ${params.adapters} --cut_tail \
        -i ${fastq1} -w 4 -o "\$R1"

    fastp -l 45 --adapter_fasta ${params.adapters} --cut_tail \
        -i ${fastq2} -w 4 -o "\$R2"
    """
}

