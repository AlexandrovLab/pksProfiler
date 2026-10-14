nextflow.enable.dsl = 2

process pksProfiler_align {
    label 'pks_align'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}" }, mode: 'copy',
        // counts and coverage are the two files people open; the rest is alignment detail
        saveAs: { fn -> (fn.endsWith('.counts.txt') || fn.endsWith('.coverage.txt'))
                        ? fn - "${sampleID}." : "alignment/" + (fn - "${sampleID}.") }
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID),
          path("${sampleID}.coverage.txt"),
          path("${sampleID}.coverage_50bp_bins.bedGraph"),
          path("${sampleID}.counts.txt"),
	      path("${sampleID}.sorted.bam"),
	      path("${sampleID}.sorted.bam.bai"),
          emit: profile
    tuple val(sampleID), path("${sampleID}.alignment.qc.tsv"), emit: qc
    // Ludmil, 25 Sep A5: divergence from IHE3034 of the clb-supporting reads. QC only.
    tuple val(sampleID), path("${sampleID}.ihe3034_divergence.qc.tsv"),
          path("${sampleID}.ihe3034_divergence_by_gene.tsv"), emit: divergence

    script:
    def bedtools_cov = "${sampleID}.coverage.txt"
    // A9 (Ludmil, 25 Sep): named for what it is -- 50 bp bins, not per-base depth, which is
    // coverage.txt below.
    def coverage     = "${sampleID}.coverage_50bp_bins.bedGraph"
    def counts       = "${sampleID}.counts.txt"
    def bam          = "${sampleID}.sorted.bam"
	def bai          = "${sampleID}.sorted.bam.bai"
    def qc           = "${sampleID}.alignment.qc.tsv"
    def divergence_qc    = "${sampleID}.ihe3034_divergence.qc.tsv"
    def divergence_genes = "${sampleID}.ihe3034_divergence_by_gene.tsv"
    // Owner decision, 29 Sep: long reads (--read_type ont|hifi) are profiled with
    // minimap2 -ax map-ont/map-hifi against the IHE3034 genome FASTA (the sequence the
    // Bowtie2 index was built from), with the same MAPQ filter; featureCounts gets -L.
    // Every output keeps its name, so nothing downstream changes.
    def long_read          = HostDepletion.isLong(params.read_type)
    def long_preset        = HostDepletion.profilingPreset(params.read_type)
    def featurecounts_long = long_read ? "-L" : ""

    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_annotation ${params.dep_digest?.pks_annotation}  pks_genome ${params.dep_digest?.pks_genome}  scripts ${params.dep_digest?.scripts}
    # pks_reference_fasta ${params.dep_digest?.pks_reference_fasta}
	    set -euo pipefail
	    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
	    ${ToolVersions.record(task, params, (long_read ? 'minimap2' : 'bowtie2') + ' samtools featureCounts:subread bamCoverage:deeptools python')}

	# Build and validate the exact manuscript-facing clbA-clbS annotation.
    awk -F '\t' '
        BEGIN { OFS="\t" }
        \$0 !~ /^#/ && \$3 == "gene" && \$9 ~ /(^|;)Name=clb[A-S](;|\$)/ { print }
    ' "${params.pks_genome_annotation}" > clb_genes.gff

    python "${params.scripts}/validate_featurecounts.py" --annotation clb_genes.gff

	if ! gzip -t "${reads}" >/dev/null 2>&1; then
        echo "ERROR: Corrupt gzip input for ${sampleID}: ${reads}" >&2
        exit 1
    fi

    if [[ "${long_read}" == "true" ]]; then
    echo "minimap2 -ax ${long_preset} alignment (sample: ${sampleID})"
    minimap2 -ax ${long_preset} --secondary=no -t "${task.cpus}" \
        "${params.pks_reference_fasta}" "${reads}" \
      | samtools view -@ "${task.cpus}" -b -F 4 -q ${params.pks_align_min_mapq} - \
      | samtools sort -@ "${task.cpus}" -o "${bam}" -
    else
    echo "Bowtie2 Alignment (sample: ${sampleID})"
    # F18: bowtie2 streams straight into the sorter. The SAM was written to disk, kept
    # as a process output, published beside the BAM that duplicates it, and read by
    # nothing -- main.nf destructured it and dropped it. At cohort scale that is one
    # redundant alignment file per sample.
    # A4 (25 Sep): mode, preset, score-min and MAPQ are named params, recorded in
    # RUN_REPORT; their defaults are the values this command always ran with.
    bowtie2 -x "${params.pks_genome}" -q -U "${reads}" \
        --seed 42 --threads "${task.cpus}" \
        --${params.pks_align_bowtie2_mode} --${params.pks_align_bowtie2_preset} \
        --score-min ${params.pks_align_bowtie2_score_min} --no-unal \
      | samtools view -@ "${task.cpus}" -bS -q ${params.pks_align_min_mapq} - \
      | samtools sort -@ "${task.cpus}" -o "${bam}" -
    fi
    samtools index -@ "${task.cpus}" "${bam}" "${bai}"

    MAPPED_READS=\$(samtools view -c -F 4 "${bam}")

    # Always create a valid 19-gene count table. A successfully
    # processed sample with no aligned reads must remain in summaries
    # as an explicit zero rather than disappearing as an empty file.
    featureCounts \
        -T "${task.cpus}" \
        -a clb_genes.gff \
        -o "${counts}" \
        -t gene \
        -F GFF \
        -g Name \
        --largestOverlap ${featurecounts_long} \
        "${bam}"

    python "${params.scripts}/validate_featurecounts.py" \
        --counts "${counts}" \
        --summary "${counts}.summary"

    # F14: one rule, and it is featureCounts'. --largestOverlap already assigns each
    # read to exactly one clb gene; the QC count is the sum of that assignment, so it
    # cannot disagree with the matrix it sits beside -- it *is* the matrix column
    # total. This used to be recomputed with `bedtools intersect -u`, which counted
    # every read touching any clb gene: a different quantity, higher than the counts
    # table, with nothing saying so. num_clb_genes_align was already derived from this
    # same file.
    CLB_READS=\$(awk -F '\t' '
        \$1 ~ /^clb[A-S]\$/ { total += \$NF + 0 }
        END { print total + 0 }
    ' "${counts}")

    CLB_GENES_DETECTED=\$(awk -F '\t' '
        \$1 ~ /^clb[A-S]\$/ && (\$NF + 0) > 0 { count++ }
        END { print count + 0 }
    ' "${counts}")

    printf "Sample\tMetric\tValue\n" > "${qc}"
    # F09: computed since v0.0.1 and used only for a zero check. It is the reads that
    # mapped to the reference at all -- the denominator the clb counts sit inside.
    # A6 (25 Sep): named for the MAPQ >= 40 filter (samtools view -q 40 above) it carries.
    printf "%s\treads_mapped_ihe3034_mapq40\t%s\n" "${sampleID}" "\$MAPPED_READS" >> "${qc}"
    printf "%s\treads_clb_genes_align\t%s\n" "${sampleID}" "\$CLB_READS" >> "${qc}"
    printf "%s\tnum_clb_genes_align\t%s\n" "${sampleID}" "\$CLB_GENES_DETECTED" >> "${qc}"

    # Ludmil, 25 Sep A5: NM and NM / alignment columns of the reads featureCounts assigned
    # to a clb gene -- divergence FROM IHE3034, and nothing more. MAPQ >= 40 on this one
    # reference is not uniqueness across microbes, and a low divergence is not organism
    # identity; specificity is the competitive-mapping lane (pks_specificity.nf). This
    # reads the BAM only: no count, tier or call below depends on it.
    samtools view "${bam}" \
      | python "${params.scripts}/ihe3034_divergence_qc.py" \
            --sample "${sampleID}" --annotation clb_genes.gff \
            --min-mapq ${params.pks_align_min_mapq} \
            --qc-out "${divergence_qc}" --genes-out "${divergence_genes}"

    if [[ "\$MAPPED_READS" -eq 0 ]]; then
        echo "No confidently mapped reads for ${sampleID}; recording zero clb counts."
        : > "${coverage}"
        : > "${bedtools_cov}"
    else
        # F13: raw depth, not RPKM. RPKM's "per million mapped reads" denominator here
        # is reads that mapped to the island reference at MAPQ >= 40 -- not the library
        # -- so the number was the share of island-mapped reads in a bin while the unit
        # name claimed library-normalised abundance. Two samples were not comparable on
        # it, and it read as though they were. Raw depth says what it is.
        #
        # This track feeds the circos figure only. Breadth and the evidence tiers come
        # from samtools depth in pks_targeted.nf, so no called result changes.
        # A9 (Ludmil, 25 Sep): --binSize 50 is deepTools' default, now stated, and in the
        # file name: each line is the mean depth over a 50 bp bin of the whole IHE3034
        # genome, where coverage.txt is per-base depth over the island only.
        bamCoverage \
			--numberOfProcessors "${task.cpus}" \
            -b "${bam}" \
            -o "${coverage}" \
            --binSize 50 \
            --normalizeUsing None \
            --outFileFormat bedgraph

        python "${params.scripts}/validate_bedgraph.py" "${coverage}"

        # F18: per-base depth over the island, not the whole 5.1 Mb reference.
        # `bedtools genomecov -d` wrote one line per base of the genome -- about 5.1
        # million lines per sample, of which the 50,768 bp island is 1%. Nothing reads
        # the file: classifyPksReadEvidence takes it as an input, passes it straight to
        # its own output, and recomputes depth from the BAM. This is the same tool and
        # the same filters that recomputation uses, so the two now agree by construction.
        # Ludmil, revised report finding 6: pks_shift+1 dropped the true first base
        # of the island (2193827). pks_start_1based/pks_end_1based are already the
        # correct 1-based inclusive bounds -- a samtools region needs no +/-1 at all.
        samtools depth -aa -s -Q 40 \
            -r "${params.pks_contig}:${params.pks_start_1based}-${params.pks_end_1based}" \
            "${bam}" > "${bedtools_cov}"
    fi
    """
}
