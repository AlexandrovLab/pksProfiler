nextflow.enable.dsl=2

process extractReads {
    scratch true
    label 'extract_reads'
    publishDir(
        { "${params.sample_dir}/${sampleID}/intermediates" },
        mode: 'copy',
        enabled: params.save_intermediates.toString().toBoolean()
    )
    conda "${params.samtools_env}"

    input:
		tuple val(sampleID), path(alignment)

    output:
	tuple val(sampleID), path("${sampleID}.UNMAPPED.fastq.gz"), emit: reads
	tuple val(sampleID), path("${sampleID}.extract.qc.tsv"), emit: qc
    
	script:
	"""
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # cram_reference ${params.dep_digest?.cram_reference}  scripts ${params.dep_digest?.scripts}
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
		if [[ "\$DETECTED" == "cram" ]]; then
		    # Ludmil, revised report finding 10: this used to fall through to a
		    # referenceless decode whenever --cram_reference was unset, which
		    # could succeed on an embedded/cached reference or fail with the
		    # generic decode error below -- unpredictable, and stricter than
		    # preflight's own sample-sheet check, which already requires
		    # --cram_reference for anything that looks like CRAM input. Both now
		    # agree: CRAM always requires it, checked here before decode is even
		    # attempted rather than left to samtools to fail on.
		    if [[ -z "${params.cram_reference ?: ''}" ]]; then
		        echo "ERROR: ${sampleID}: CRAM input requires --cram_reference" >&2
		        exit 1
		    fi
		    if [[ ! -f "${params.cram_reference ?: ''}" ]]; then
		        echo "ERROR: CRAM reference does not exist: ${params.cram_reference ?: ''}" >&2
		        exit 1
		    fi
		    python "${params.scripts}/validate_cram_reference.py" \
		        --alignment "${alignment}" \
		        --reference "${params.cram_reference}"
		    REFERENCE_ARGS=(--reference "${params.cram_reference}")
		fi

		samtools quickcheck -v "${alignment}"

    # Retain every primary unmapped alignment record, regardless of
    # whether its mate is mapped, unmapped, or absent. The alignment is
    # decoded only once; the extracted-read count comes from the FASTQ.
    # samtools fastq applies -f/-F itself, so the intermediate uncompressed-BAM
    # pipe through samtools view is unnecessary; dropping it removes ~60 GB of
    # pipe traffic per sample. The read count is tee'd off the live stream
    # instead of decompressing the finished file a second time.
    mkfifo extract.count.fifo
    awk 'END { print int(NR / 4) }' < extract.count.fifo > extract.count &
    EXTRACT_COUNTER=\$!

    # No -o / -0 / -s: every category -- READ1, READ2, READ_OTHER and singletons
    # -- goes to stdout, which is captured below. -o writes only READ1/READ2 and
    # -0 only READ_OTHER, so naming either one silently drops the rest. Records
    # with neither mate bit set (unpaired input) or with both set are READ_OTHER;
    # they were discarded until 2026-09-20. See Ludmil's audit, F01.
    if ! samtools fastq \
        -@ "${task.cpus}" \
	        -N \
	        -f 4 \
	        -F 2304 \
	        "\${REFERENCE_ARGS[@]}" \
	        "${alignment}" |
    tee extract.count.fifo |
    bgzip -@ "${task.cpus}" -c > "\$READS"; then
	        echo "ERROR: Could not decode ${alignment}. For CRAM, provide the matching --cram_reference if it is not embedded or cached." >&2
	        exit 1
	    fi

		wait "\$EXTRACT_COUNTER"
		rm -f extract.count.fifo

		bgzip -t "\$READS"
		UNMAPPED_READS=\$(cat extract.count)

	# F09: a library denominator that does not come from the extraction itself. The QC
	# summary used to fall back to the post-extraction count, so the table could not show
	# how much of the library was dropped -- which is why F01 stayed invisible in QC.
	#
	# idxstats reads the index only, so this costs milliseconds rather than a second
	# decode of a 66 GiB CRAM. It counts alignment *records*, secondary and supplementary
	# included, and it needs an index; both are why the metric is named for what it is,
	# and why it is simply absent when it cannot be had rather than being stood in for.
	INPUT_RECORDS=""
	if samtools idxstats "\${REFERENCE_ARGS[@]}" "${alignment}" > input.idxstats 2>/dev/null; then
	    INPUT_RECORDS=\$(awk '{ total += \$3 + \$4 } END { print total + 0 }' input.idxstats)
	fi

	TOTAL_PRIMARY=""
	if [[ "${params.exact_input_counts}" == "true" ]]; then
	    TOTAL_PRIMARY=\$(samtools view -c -@ "${task.cpus}" -F 2304 "\${REFERENCE_ARGS[@]}" "${alignment}")
	fi

	printf "Sample\tMetric\tValue\n" > "\$QC"
	if [[ -n "\$INPUT_RECORDS" ]]; then
	    printf "%s\tinput_alignment_records\t%s\n" "${sampleID}" "\$INPUT_RECORDS" >> "\$QC"
	fi
	if [[ -n "\$TOTAL_PRIMARY" ]]; then
	    printf "%s\ttotal_primary_reads\t%s\n" "${sampleID}" "\$TOTAL_PRIMARY" >> "\$QC"
	fi
	printf "%s\textracted_unmapped_reads\t%s\n" "${sampleID}" "\$UNMAPPED_READS" >> "\$QC"
	"""
}
