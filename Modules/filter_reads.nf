nextflow.enable.dsl=2

process filterReads {
    scratch true
    label 'filter_reads'
    publishDir("${params.unmapped_bam_dir}", mode: 'copy')
    conda "${params.fastp_env}"
	maxRetries 2


	input:
	tuple val(sampleID), path(fastq_files)

	output:
	tuple val(sampleID), path("${sampleID}.UNMAPPED.FASTP.FILTERED.fastq.gz")

	script:
	def input_list = fastq_files instanceof Collection ? fastq_files : [fastq_files]
	def inputs = input_list.collect { input_fastq -> "\"${input_fastq}\"" }.join(' ')
	def paired_merge_commands = input_list.withIndex().collect { input_fastq, index ->
	        def mate = index + 1
	        """gzip -dc \"${input_fastq}\" | awk -v mate=\"${mate}\" '
	            NR % 4 == 1 {
	                split(\$0, header_parts, " ")
	                read_id=header_parts[1]

	                if (read_id !~ /\\/[12]\$/) {
	                    read_id=read_id "/" mate
	                }

	                \$0=read_id substr(\$0, length(header_parts[1]) + 1)
	            }

	            { print }
	        '"""
	    }.join('\n')
	def merge_reads = input_list.size() == 1 ?
	    "gzip -dc ${inputs} | gzip -c > \"\$MERGED\"" :
	    "{\n${paired_merge_commands}\n} | gzip -c > \"\$MERGED\""

	"""
	set -euo pipefail

	MERGED="${sampleID}.UNMAPPED.merged.fastq.gz"
	FILTERED="${sampleID}.UNMAPPED.FASTP.FILTERED.fastq.gz"

	for input_fastq in ${inputs}; do
	    if ! gzip -t "\$input_fastq" >/dev/null 2>&1; then
	        echo "ERROR: Corrupt gzip input for ${sampleID}: \$input_fastq" >&2
	        exit 1
	    fi
	done

	# BAM input contributes one stream. For paired FASTQ input, append /1 and
	# /2 when needed so downstream HMM queries retain distinct mate IDs.
	${merge_reads}

    fastp \
        -l 45 \
        --adapter_fasta "${params.adapters}" \
        --cut_tail \
        -i "\$MERGED" \
        -w "${task.cpus}" \
        -o "\$FILTERED"
    """
}
