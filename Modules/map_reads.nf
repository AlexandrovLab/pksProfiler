process mapReads {

    tag "$sampleID"
    scratch true
    label 'map_reads'
    publishDir(
        { "${params.sample_dir}/${sampleID}/intermediates" },
        mode: 'copy',
        enabled: params.save_intermediates.toString().toBoolean()
    )
    conda "${params.minimap2_env}"

    input:
    tuple val(sampleID), path(reads_fastq)

    output:
    tuple val(sampleID), path("${sampleID}.host_depleted.fastq.gz"), emit: reads
    tuple val(sampleID), path("${sampleID}.depletion.qc.tsv"), emit: qc

    script:
    // E0 (Ludmil, 25 Sep): with one FASTQ, -ax sr treats consecutive same-named reads
    // as a pair (fragment mode), so what is depleted depends on read order and
    // therefore on the input route. --frag=no makes every read independent.
    // Every pass's preset, index and flags come from lib/HostDepletion.groovy, which
    // carries frag_arg into each short-read call.
    //
    // Ludmil, 26 Sep mode memo section 4: for --sample_type rna the human pass is
    // splice-aware (-ax splice:sr against --rna_splice_db, plus -j --rna_junctions_bed when
    // given), because -ax sr leaves exon-junction human reads undepleted. The T2T+phiX and
    // pangenome passes stay -ax sr: there is no junction annotation for those references.
    // The QC metric names (reads_after_hg38, ...) are unchanged in every mode.
    // Long reads (--read_type ont|hifi): map-ont/map-hifi for DNA passes, splice/splice:hq
    // for rna's human pass, against the *_long indexes; --frag is omitted (it only changes
    // how minimap2 treats short paired reads).
    def human = HostDepletion.humanPass(params.sample_type, params.read_type, params.minimap2_frag_mode,
                                        params.hg38_db, params.hg38_db_long,
                                        params.rna_splice_db, params.rna_splice_db_long,
                                        params.rna_junctions_bed)
    def t2t = HostDepletion.t2tPass(params.read_type, params.minimap2_frag_mode,
                                    params.t2t_phix_db, params.t2t_phix_db_long)
    def pangenome = HostDepletion.pangenomePass(params.read_type, params.minimap2_frag_mode)
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # hg38_db ${params.dep_digest?.hg38_db}  pangenome_db ${params.dep_digest?.pangenome_db}  t2t_phix_db ${params.dep_digest?.t2t_phix_db}
    # rna_splice_db ${params.dep_digest?.rna_splice_db}  rna_junctions_bed ${params.dep_digest?.rna_junctions_bed}
    # hg38_db_long ${params.dep_digest?.hg38_db_long}  t2t_phix_db_long ${params.dep_digest?.t2t_phix_db_long}  rna_splice_db_long ${params.dep_digest?.rna_splice_db_long}
    set -euo pipefail

    if ! gzip -t "${reads_fastq}" >/dev/null 2>&1; then
        echo "ERROR: Corrupt gzip input for ${sampleID}: ${reads_fastq}" >&2
        exit 1
    fi

    echo "Running host depletion"

    QC="${sampleID}.depletion.qc.tsv"

    # The read count is taken from a tee off the live stream rather than by
    # decompressing the finished file: gzip -dc is single-threaded and was
    # stalling a 16-thread minimap2 behind it. mkfifo + an explicit wait keeps
    # the count deterministic (process substitution is not waited on by bash).
    mkfifo hg38.count.fifo
    awk 'END { print int(NR / 4) }' < hg38.count.fifo > hg38.count &
    HG38_COUNTER=\$!

    echo "Human pass: minimap2 -ax ${human.preset} against ${human.index}"
    minimap2 -2 -ax ${human.preset} ${human.args} --secondary=no -t "${task.cpus}" \
        "${human.index}" \
        "${reads_fastq}" |
    samtools fastq -@ "${task.cpus}" -f 4 -F 2304 |
    tee hg38.count.fifo |
    bgzip -@ "${task.cpus}" -c > "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz"

    wait "\$HG38_COUNTER"
    AFTER_HG38=\$(cat hg38.count)
    rm -f hg38.count.fifo

    mkfifo t2t.count.fifo
    awk 'END { print int(NR / 4) }' < t2t.count.fifo > t2t.count &
    T2T_COUNTER=\$!

    minimap2 -2 -ax ${t2t.preset} ${t2t.args} --secondary=no -t "${task.cpus}" \
        "${t2t.index}" \
        "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz" |
    samtools fastq -@ "${task.cpus}" -f 4 -F 2304 |
    tee t2t.count.fifo |
    bgzip -@ "${task.cpus}" -c > "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"

    wait "\$T2T_COUNTER"
    AFTER_T2T_PHIX=\$(cat t2t.count)
    rm -f t2t.count.fifo

    if [[ -n "${params.pangenome_db ?: ''}" ]]; then
        echo "Running optional human-pangenome depletion"

        if [[ -d "${params.pangenome_db}" ]]; then
            mapfile -t PANGENOME_INDEXES < <(
                find "${params.pangenome_db}" -type f -name '*.mmi' -print | sort
            )
        elif [[ -f "${params.pangenome_db}" ]]; then
            PANGENOME_INDEXES=("${params.pangenome_db}")
        else
            echo "ERROR: Pangenome index path does not exist: ${params.pangenome_db}" >&2
            exit 1
        fi

        if [[ "\${#PANGENOME_INDEXES[@]}" -eq 0 ]]; then
            echo "ERROR: No .mmi files found under ${params.pangenome_db}" >&2
            exit 1
        fi

        cp "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz" \
            "${sampleID}.pangenome.current.fastq.gz"

        for mmi in "\${PANGENOME_INDEXES[@]}"; do
            echo "Running minimap2 on \$mmi"

            minimap2 -2 -ax ${pangenome.preset} ${pangenome.args} --secondary=no -t "${task.cpus}" \
                "\$mmi" \
                "${sampleID}.pangenome.current.fastq.gz" |
            samtools fastq -@ "${task.cpus}" -f 4 -F 2304 |
            bgzip -@ "${task.cpus}" -c > "${sampleID}.pangenome.next.fastq.gz"

            mv "${sampleID}.pangenome.next.fastq.gz" \
                "${sampleID}.pangenome.current.fastq.gz"
        done

        mv "${sampleID}.pangenome.current.fastq.gz" \
            "${sampleID}.host_depleted.fastq.gz"

        AFTER_PANGENOME=\$(bgzip -@ "${task.cpus}" -dc "${sampleID}.host_depleted.fastq.gz" | awk 'END { print int(NR / 4) }')
    else
        mv "${sampleID}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz" \
            "${sampleID}.host_depleted.fastq.gz"
    fi

    printf "Sample\tMetric\tValue\n" > "\$QC"
    printf "%s\treads_after_hg38\t%s\n" "${sampleID}" "\$AFTER_HG38" >> "\$QC"
    printf "%s\treads_after_t2t_phix\t%s\n" "${sampleID}" "\$AFTER_T2T_PHIX" >> "\$QC"

    if [[ -n "${params.pangenome_db ?: ''}" ]]; then
        printf "%s\treads_after_pangenome\t%s\n" "${sampleID}" "\$AFTER_PANGENOME" >> "\$QC"
    fi
    """
}
