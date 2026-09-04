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
		tuple val(sampleID), path(alignment)

    output:
	tuple val(sampleID), path("${sampleID}.UNMAPPED.fastq.gz"), emit: reads
	tuple val(sampleID), path("${sampleID}.extract.qc.tsv"), emit: qc
    
	script:
	"""
	set -euo pipefail

		READS="${sampleID}.UNMAPPED.fastq.gz"
		QC="${sampleID}.extract.qc.tsv"

		FORMAT="\$(htsfile "${alignment}" 2>/dev/null || true)"
		case "\$FORMAT" in
		    *CRAM*) DETECTED=cram ;;
		    *BAM*)  DETECTED=bam ;;
		    *)
		        echo "ERROR: ${sampleID}: input is neither BAM nor CRAM according to htsfile: \$FORMAT" >&2
		        exit 1
		        ;;
		esac

		if [[ "${params.input_data_type}" == "bam" && "\$DETECTED" != "bam" ]]; then
		    echo "ERROR: ${sampleID}: BAM was requested but the input is CRAM" >&2
		    exit 1
		fi
		if [[ "${params.input_data_type}" == "cram" && "\$DETECTED" != "cram" ]]; then
		    echo "ERROR: ${sampleID}: CRAM was requested but the input is BAM" >&2
		    exit 1
		fi

		REFERENCE_ARGS=()
		if [[ "\$DETECTED" == "cram" && -n "${params.cram_reference ?: ''}" ]]; then
		    if [[ ! -f "${params.cram_reference ?: ''}" ]]; then
		        echo "ERROR: CRAM reference does not exist: ${params.cram_reference ?: ''}" >&2
		        exit 1
		    fi
		    python "${params.scripts}/validate_cram_reference.py" \
		        --alignment "${alignment}" \
		        --reference "${params.cram_reference}"
		    REFERENCE_ARGS=(-T "${params.cram_reference}")
		fi

		samtools quickcheck -v "${alignment}"

    # Retain every primary unmapped alignment record, regardless of
    # whether its mate is mapped, unmapped, or absent. The alignment is
    # decoded only once; the extracted-read count comes from the FASTQ.
    if ! samtools view \
        -@ "${task.cpus}" \
	        -f 4 \
	        -F 2304 \
	        -u \
	        "\${REFERENCE_ARGS[@]}" \
	        "${alignment}" |
    samtools fastq \
        -@ "${task.cpus}" \
        -N \
	        -o "\$READS" \
	        -; then
	        echo "ERROR: Could not decode ${alignment}. For CRAM, provide the matching --cram_reference if it is not embedded or cached." >&2
	        exit 1
	    fi

		gzip -t "\$READS"
		UNMAPPED_READS=\$(gzip -dc "\$READS" | awk 'END { print int(NR / 4) }')

	printf "Sample\tMetric\tValue\n" > "\$QC"
	printf "%s\textracted_unmapped_reads\t%s\n" "${sampleID}" "\$UNMAPPED_READS" >> "\$QC"
	"""
}
