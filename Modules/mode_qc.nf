nextflow.enable.dsl = 2

// Mode-specific QC (Ludmil, 26 Sep mode memo sections 3 and 4). Each process writes one
// Sample/Metric/Value fragment that build_qc_summary.py folds into pks.qc.summary.tsv
// (its WES / RNA columns); rnaResidualHuman's also reaches rna_evidence.tsv. Definitions:
// docs/running/wes_rna.md.

// WES: duplicate rate and capture on-target counts of the INPUT alignment. One decode of
// the primary mapped records (-F 2308), teed into flagstat and, with --wes_capture_bed,
// into a -L count. Alignment input only: a FASTQ sheet has no alignment to count.
// Separate from extractReads (which decodes the unmapped reads only) so the WGS path is
// untouched; for WES inputs (a few GB) the second decode is cheap.
process wesInputQC {
    label 'sample_stage'
    tag "$sampleID"
    label 'process_low'
    scratch true
    conda "${params.samtools_env}"

    input:
    tuple val(sampleID), path(alignment)

    output:
    tuple val(sampleID), path("${sampleID}.wes_input.qc.tsv"), emit: qc

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # cram_reference ${params.dep_digest?.cram_reference}  wes_capture_bed ${params.dep_digest?.wes_capture_bed}  scripts ${params.dep_digest?.scripts}
    set -euo pipefail

    FORMAT="\$(htsfile "${alignment}" 2>/dev/null || true)"
    REFERENCE_ARGS=()
    case "\$FORMAT" in
        *CRAM*) REFERENCE_ARGS=(--reference "${params.cram_reference ?: ''}") ;;
    esac

    CAPTURE_BED="${params.wes_capture_bed ?: ''}"
    ON_TARGET_ARGS=()
    if [[ -n "\$CAPTURE_BED" ]]; then
        mkfifo on_target.fifo
        samtools view -c -L "\$CAPTURE_BED" on_target.fifo > on_target.count &
        ON_TARGET_COUNTER=\$!
        samtools view -@ "${task.cpus}" -u -F 2308 "\${REFERENCE_ARGS[@]}" "${alignment}" |
            tee on_target.fifo |
            samtools flagstat -O tsv - > input.flagstat.tsv
        wait "\$ON_TARGET_COUNTER"
        rm -f on_target.fifo
        ON_TARGET_ARGS=(--on-target on_target.count)
    else
        samtools view -@ "${task.cpus}" -u -F 2308 "\${REFERENCE_ARGS[@]}" "${alignment}" |
            samtools flagstat -O tsv - > input.flagstat.tsv
    fi

    python "${params.scripts}/wes_input_qc.py" --sample "${sampleID}" \
        --flagstat input.flagstat.tsv "\${ON_TARGET_ARGS[@]}" \
        --output "${sampleID}.wes_input.qc.tsv"
    """
}

// WES: duplicate rate among pks-island reads of the profiling BAM, by position collapse
// (the input's 0x400 flags do not survive extraction to FASTQ). scripts/pks_duplicate_rate.py.
process wesPksDuplicates {
    label 'sample_stage'
    tag "$sampleID"
    label 'process_low'
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(bam), path(bai)

    output:
    tuple val(sampleID), path("${sampleID}.wes_pks_duplicates.qc.tsv"), emit: qc

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    samtools view -F 2308 "${bam}" \
        "${params.pks_contig}:${params.pks_start_1based}-${params.pks_end_1based}" |
        python "${params.scripts}/pks_duplicate_rate.py" --sample "${sampleID}" \
            --output "${sampleID}.wes_pks_duplicates.qc.tsv"
    """
}

// RNA: residual human after every depletion pass. A capped systematic sample of the final
// host-depleted reads, realigned to the splice index with the human pass's own preset,
// index and junctions plus --rna_residual_human_args (short reads only). A measurement:
// the reads that go on to profiling are not changed. scripts/rna_residual_human.py.
process rnaResidualHuman {
    label 'sample_stage'
    tag "$sampleID"
    // loads the same GRCh38 splice index as mapReads
    label 'map_reads'
    scratch true
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${sampleID}.residual_human.qc.tsv"), emit: qc

    script:
    def human = HostDepletion.humanPass('rna', params.read_type, params.minimap2_frag_mode,
                                        null, null,
                                        params.rna_splice_db, params.rna_splice_db_long,
                                        params.rna_junctions_bed)
    def relaxed = HostDepletion.isLong(params.read_type) ? '' : (params.rna_residual_human_args ?: '')
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # rna_splice_db ${params.dep_digest?.rna_splice_db}  rna_junctions_bed ${params.dep_digest?.rna_junctions_bed}
    # rna_splice_db_long ${params.dep_digest?.rna_splice_db_long}  scripts ${params.dep_digest?.scripts}
    set -euo pipefail

    python "${params.scripts}/rna_residual_human.py" subsample --input "${reads}" \
        --cap ${params.rna_residual_human_sample_reads} \
        --output residual.sample.fastq --stats residual.stats

    echo "Residual-human check: minimap2 -ax ${human.preset} ${human.args} ${relaxed} against ${human.index}"
    minimap2 -ax ${human.preset} ${human.args} ${relaxed} --secondary=no -t "${task.cpus}" \
        "${human.index}" residual.sample.fastq |
        python "${params.scripts}/rna_residual_human.py" count --sam - --sample "${sampleID}" \
            --stats residual.stats --output "${sampleID}.residual_human.qc.tsv"
    """
}
