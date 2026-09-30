nextflow.enable.dsl=2

process extractReads {
    scratch true
    label 'extract_reads'
    publishDir(
        { "${params.sample_dir}/${sampleID}/intermediates" },
        mode: 'copy',
        enabled: params.save_intermediates.toString().toBoolean()
    )
    // E3 (Ludmil, 25 Sep): the input's reference and upstream-aligner provenance, always
    // published: by_sample/<sample>/provenance/input_header.tsv.
    publishDir(
        { "${params.sample_dir}/${sampleID}/provenance" },
        mode: 'copy',
        pattern: "*.header_provenance.tsv",
        saveAs: { "input_header.tsv" }
    )
    conda "${params.samtools_env}"

    input:
		tuple val(sampleID), path(alignment)

    output:
	tuple val(sampleID), path("${sampleID}.UNMAPPED.fastq.gz"), emit: reads
	tuple val(sampleID), path("${sampleID}.extract.qc.tsv"), emit: qc
	tuple val(sampleID), path("${sampleID}.header_provenance.tsv"), emit: provenance
    
	script:
	"""
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # cram_reference ${params.dep_digest?.cram_reference}  scripts ${params.dep_digest?.scripts}
	set -euo pipefail
	# F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
	${ToolVersions.record(task, params, 'samtools bgzip:htslib python')}

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

		# An unaligned BAM/CRAM (samtools import, a basecaller's uBAM) has no @SQ lines,
		# and the quickcheck below refuses it with "had no targets in header", which
		# does not say what is wrong. Extraction keeps reads an aligner flagged unmapped
		# against a reference; with no reference there is nothing to extract from. Say so,
		# before CRAM validation or quickcheck. The header is read once and kept.
		if ! samtools view -H "${alignment}" > input.header.sam 2> input.header.err; then
		    echo "ERROR: ${sampleID}: could not read the header of ${alignment}: \$(head -c 500 input.header.err)" >&2
		    exit 1
		fi
		if ! grep -q '^@SQ' input.header.sam; then
		    echo "ERROR: ${sampleID}: ${alignment} has no @SQ lines in its header: it is an unaligned \$DETECTED (e.g. from samtools import or a basecaller), which the input quickcheck rejects. pksProfiler extracts the reads an aligner flagged unmapped against a reference. Convert it to FASTQ and run it with --input_data_type fastq (docs/running/bam.md)." >&2
		    exit 1
		fi
		# E3: @SQ fingerprint (md5 of ordered SN:LN), M5-set fingerprint, contig count,
		# build guess from chr1's length, and the @PG chain (ID/PN/VN along PP, CL never recorded).
		# E1 (Ludmil, 25 Sep): the same step decides whether low-score mapped records are
		# extracted too (--extract_low_score_mapped auto|on|off; auto = on for bwa mem).
		python "${params.scripts}/header_provenance.py" sample \
		    --header input.header.sam --sample "${sampleID}" \
		    --low-score-mode "${params.extract_low_score_mapped}" \
		    --min-as "${params.extract_min_as}" \
		    --decision extract_low_score.decision \
		    --output "${sampleID}.header_provenance.tsv"
		IFS=\$'\\t' read -r LOW_SCORE LOW_SCORE_STATE LOW_SCORE_REASON < extract_low_score.decision
		echo "extractReads: ${sampleID}: low-score mapped extraction \$LOW_SCORE_REASON" >&2

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
    # With low-score extraction off, samtools fastq applies -f/-F itself, so the
    # intermediate uncompressed-BAM pipe through samtools view is unnecessary;
    # dropping it removes ~60 GB of pipe traffic per sample. (On, samtools view
    # filters during the decode and only the selected records cross the pipe.)
    # The read count is tee'd off the live stream instead of decompressing the
    # finished file a second time.
    mkfifo extract.count.fifo
    awk 'END { print int(NR / 4) }' < extract.count.fifo > extract.count &
    EXTRACT_COUNTER=\$!

    if [[ "\$LOW_SCORE" == "on" ]]; then
        # E1 (Ludmil, 25 Sep): the same single decode, plus primary MAPPED records whose
        # AS:i is below --extract_min_as. bwa mem calls a read unmapped only when its
        # best local score is below -T; GDC aligned TCGA with -T 0, so about a third of
        # genuine pks reads come out "mapped" to a short spurious hit and FLAG 4 alone
        # never sees them. AS < 30 (bwa's default -T, Mutographs' setting) is exactly the
        # set -T 30 would have left unmapped (validation/E1: synthetic 99.9-100% of the
        # hidden reads recovered; real pilot, matched normals gained no clb read).
        #
        # samtools view does the decode and passes on only FLAG 4 records, records with
        # AS below the threshold, and mapped records with no AS at all -- the last only to
        # be counted; the awk drops them. What reaches samtools fastq is FLAG 4 plus
        # low-score records, so read names, /1 /2 and reverse-strand orientation are
        # samtools fastq -N's, exactly as in the off branch. The awk counts each class.
        if ! samtools view -h \
            -@ "${task.cpus}" \
            -F 2304 \
            -e '(flag & 4) || !exists([AS]) || [AS] < ${params.extract_min_as}' \
            "\${REFERENCE_ARGS[@]}" \
            "${alignment}" |
        awk -F '\\t' -v min_as="${params.extract_min_as}" -v counts=extract.classes '
            /^@/ { print; next }
            {
                if (int(\$2 / 4) % 2 == 1) { unmapped++; print; next }
                score = ""
                for (i = 12; i <= NF; i++) if (\$i ~ /^AS:[if]:/) { score = substr(\$i, 6); break }
                if (score == "") { no_as++; next }
                if (score + 0 < min_as + 0) { low++; print; next }
                other++
            }
            END { printf "%d\\t%d\\t%d\\t%d\\n", unmapped, low, no_as, other > counts }' |
        samtools fastq \
            -@ "${task.cpus}" \
            -N \
            -F 2304 \
            - |
        tee extract.count.fifo |
        bgzip -@ "${task.cpus}" -c > "\$READS"; then
            echo "ERROR: Could not decode ${alignment}. For CRAM, provide the matching --cram_reference if it is not embedded or cached." >&2
            exit 1
        fi
    else
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
    fi

		wait "\$EXTRACT_COUNTER"
		rm -f extract.count.fifo

		bgzip -t "\$READS"
		EXTRACTED_READS=\$(cat extract.count)
		# extracted_reads is what extraction wrote, and what enters fastp
		# (filter_input_reads). extracted_unmapped_reads keeps its meaning, FLAG 4 records
		# only; extracted_low_score_mapped_reads is E1's addition. Off: all of it is FLAG 4.
		if [[ "\$LOW_SCORE" == "on" ]]; then
		    read -r UNMAPPED_READS LOW_SCORE_READS NO_AS_RECORDS OTHER_RECORDS < extract.classes
		    if [[ "\$OTHER_RECORDS" != 0 || \$((UNMAPPED_READS + LOW_SCORE_READS)) != "\$EXTRACTED_READS" ]]; then
		        echo "ERROR: ${sampleID}: extraction classes do not add up: FLAG 4 \$UNMAPPED_READS + low-score \$LOW_SCORE_READS != \$EXTRACTED_READS FASTQ reads (unclassified: \$OTHER_RECORDS)" >&2
		        exit 1
		    fi
		    echo "extractReads: ${sampleID}: \$UNMAPPED_READS FLAG 4 + \$LOW_SCORE_READS low-score mapped reads (AS < ${params.extract_min_as}); \$NO_AS_RECORDS mapped primary records had no AS tag and were not added" >&2
		else
		    UNMAPPED_READS=\$EXTRACTED_READS
		    LOW_SCORE_READS=0
		    NO_AS_RECORDS=""
		fi

	# F09: a library denominator that does not come from the extraction itself. The QC
	# summary used to fall back to the post-extraction count, so the table could not show
	# how much of the library was dropped -- which is why F01 stayed invisible in QC.
	#
	# idxstats reads the index only, so this costs milliseconds rather than a second
	# decode of a 66 GiB CRAM. It counts alignment *records*, secondary and supplementary
	# included, and it needs an index; both are why the metric is named for what it is,
	# and why it is simply absent when it cannot be had rather than being stood in for.
	# E2 (25 Sep): the per-contig counts are kept as a few reference categories
	# (canonical / mito / alt / unplaced / decoy / viral / other, plus unplaced
	# unmapped) instead of being collapsed straight to the total.
	HAVE_IDXSTATS=false
	if samtools idxstats "\${REFERENCE_ARGS[@]}" "${alignment}" > input.idxstats 2>/dev/null; then
	    python "${params.scripts}/summarize_idxstats.py" --idxstats input.idxstats \
	        --sample "${sampleID}" > input.idxstats.metrics
	    HAVE_IDXSTATS=true
	fi

	TOTAL_PRIMARY=""
	# E5: Groovy boolean coercion, so "True", "1", "yes" etc. are not silently false.
	if [[ "${params.exact_input_counts.toString().toBoolean()}" == "true" ]]; then
	    TOTAL_PRIMARY=\$(samtools view -c -@ "${task.cpus}" -F 2304 "\${REFERENCE_ARGS[@]}" "${alignment}")
	fi

	printf "Sample\tMetric\tValue\n" > "\$QC"
	if [[ "\$HAVE_IDXSTATS" == "true" ]]; then
	    cat input.idxstats.metrics >> "\$QC"
	fi
	if [[ -n "\$TOTAL_PRIMARY" ]]; then
	    printf "%s\ttotal_primary_reads\t%s\n" "${sampleID}" "\$TOTAL_PRIMARY" >> "\$QC"
	fi
	printf "%s\textracted_unmapped_reads\t%s\n" "${sampleID}" "\$UNMAPPED_READS" >> "\$QC"
	printf "%s\textracted_low_score_mapped_reads\t%s\n" "${sampleID}" "\$LOW_SCORE_READS" >> "\$QC"
	printf "%s\textracted_reads\t%s\n" "${sampleID}" "\$EXTRACTED_READS" >> "\$QC"
	if [[ -n "\$NO_AS_RECORDS" ]]; then
	    printf "%s\textract_low_score_no_as_records\t%s\n" "${sampleID}" "\$NO_AS_RECORDS" >> "\$QC"
	fi
	printf "%s\textract_low_score_mode\t%s\n" "${sampleID}" "\$LOW_SCORE_STATE" >> "\$QC"
	"""
}


// E3 (Ludmil, 25 Sep): one cohort table of every input's header provenance, plus a
// counts-only summary (distinct reference fingerprints, distinct aligner chains) that
// RUN_REPORT.txt renders as INPUT ALIGNMENTS. Both under cohort/provenance/.
process inputHeaderSummary {
    label 'process_low'
    publishDir "${params.cohort_dir}/provenance", mode: 'copy'
    conda "${params.samtools_env}"

    input:
    path(provenance_files)
    path(file_list)

    output:
    path "input_headers.tsv", emit: table
    path "input_header_summary.tsv", emit: summary

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python')}
    python "${params.scripts}/header_provenance.py" summarize \
        --inputs-file "${file_list}" \
        --table input_headers.tsv \
        --summary input_header_summary.tsv
    """
}
