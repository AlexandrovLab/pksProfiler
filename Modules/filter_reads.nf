nextflow.enable.dsl=2

process filterReads {
    scratch true
    label 'filter_reads'
	publishDir(
	    "${params.unmapped_bam_dir}",
	    mode: 'copy',
	    enabled: params.save_intermediates
	)
    conda "${params.fastp_env}"
	maxRetries 2


	input:
	tuple val(sampleID), path(fastq_files)

	output:
	tuple val(sampleID), path("${sampleID}.UNMAPPED.FASTP.FILTERED.fastq.gz"), emit: reads
	tuple val(sampleID), path("${sampleID}.filter.qc.tsv"), emit: qc

	script:
	def input_list = fastq_files instanceof Collection ? fastq_files : [fastq_files]
	def inputs = input_list.collect { input_fastq -> "\"${input_fastq}\"" }.join(' ')
	def tagged_input_commands = input_list.withIndex().collect { input_fastq, index ->
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
	def fastp_input = input_list.size() == 1 ?
	    "-i ${inputs}" :
	    "--stdin"
	def fastp_prefix = input_list.size() == 1 ?
	    "" :
	    "{\n${tagged_input_commands}\n} |"

	"""
	set -euo pipefail

	FILTERED="${sampleID}.UNMAPPED.FASTP.FILTERED.fastq.gz"
	FASTP_JSON="${sampleID}.fastp.json"
	FASTP_HTML="${sampleID}.fastp.html"
	QC="${sampleID}.filter.qc.tsv"

	for input_fastq in ${inputs}; do
	    if ! gzip -t "\$input_fastq" >/dev/null 2>&1; then
	        echo "ERROR: Corrupt gzip input for ${sampleID}: \$input_fastq" >&2
	        exit 1
	    fi
	done

	# Single streams go directly to fastp. Paired FASTQs are streamed as one
	# unpaired input after mate suffixes are added; no merged file is written.
	${fastp_prefix} fastp \
	        -l 45 \
	        --adapter_fasta "${params.adapters}" \
	        --cut_tail \
	        ${fastp_input} \
	        -w "${task.cpus}" \
	        --compression 4 \
	        --json "\$FASTP_JSON" \
	        --html "\$FASTP_HTML" \
	        -o "\$FILTERED"

	if ! gzip -t "\$FILTERED" >/dev/null 2>&1; then
	    echo "ERROR: Corrupt fastp output for ${sampleID}: \$FILTERED" >&2
	    exit 1
	fi

	FILTER_INPUT_READS=\$(awk '
	    /"before_filtering"[[:space:]]*:/ { section=1; next }
	    section == 1 && /"total_reads"[[:space:]]*:/ {
	        value=\$0
	        gsub(/[^0-9]/, "", value)
	        print value
	        exit
	    }
	' "\$FASTP_JSON")
	FILTERED_READS=\$(awk '
	    /"after_filtering"[[:space:]]*:/ { section=1; next }
	    section == 1 && /"total_reads"[[:space:]]*:/ {
	        value=\$0
	        gsub(/[^0-9]/, "", value)
	        print value
	        exit
	    }
	' "\$FASTP_JSON")

	if [[ ! "\$FILTER_INPUT_READS" =~ ^[0-9]+\$ || ! "\$FILTERED_READS" =~ ^[0-9]+\$ ]]; then
	    echo "ERROR: Missing read counts in fastp JSON for ${sampleID}" >&2
	    exit 1
	fi

	printf "Sample\tMetric\tValue\n" > "\$QC"
	printf "%s\tfilter_input_reads\t%s\n" "${sampleID}" "\$FILTER_INPUT_READS" >> "\$QC"
	printf "%s\treads_after_fastp\t%s\n" "${sampleID}" "\$FILTERED_READS" >> "\$QC"
    """
}
