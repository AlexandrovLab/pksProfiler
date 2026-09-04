nextflow.enable.dsl = 2

process pksProfiler_align {
    label 'pks_align'
    scratch true
    publishDir "${params.pks_dir}", mode: 'copy'
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID),
          path("${sampleID}.coverage.txt"),
          path("${sampleID}.coverage.bedgraph"),
          path("${sampleID}.counts.txt"),
	      path("${sampleID}.sorted.bam"),
	      path("${sampleID}.sorted.bam.bai"),
	      path("${sampleID}.sam"),
          emit: profile
    tuple val(sampleID), path("${sampleID}.alignment.qc.tsv"), emit: qc

    script:
    def bedtools_cov = "${sampleID}.coverage.txt"
    def coverage     = "${sampleID}.coverage.bedgraph"
    def counts       = "${sampleID}.counts.txt"
    def bam          = "${sampleID}.sorted.bam"
	def bai          = "${sampleID}.sorted.bam.bai"
    def sam          = "${sampleID}.sam"
    def qc           = "${sampleID}.alignment.qc.tsv"

    """
	    set -euo pipefail

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

    echo "Bowtie2 Alignment (sample: ${sampleID})"
    bowtie2 -x "${params.pks_genome}" -q -U "${reads}" \
        --seed 42 --threads "${task.cpus}" --very-sensitive --no-unal -S "${sam}"

    samtools view -@ "${task.cpus}" -bS -q 40 "${sam}" | samtools sort -@ "${task.cpus}" -o "${bam}" -
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
        --largestOverlap \
        "${bam}"

    python "${params.scripts}/validate_featurecounts.py" \
        --counts "${counts}" \
        --summary "${counts}.summary"

    awk -F '\t' '
        BEGIN { OFS="\t" }
        \$0 !~ /^#/ && \$3 == "gene" && \$9 ~ /(^|;)Name=clb[A-S](;|\$)/ {
            gene=""
            n=split(\$9, attributes, ";")
            for (i=1; i<=n; i++) {
                if (attributes[i] ~ /^Name=/) {
                    sub(/^Name=/, "", attributes[i])
                    gene=attributes[i]
                }
            }
            if (gene != "") print \$1, \$4-1, \$5, gene
        }
    ' "${params.pks_genome_annotation}" > clb_genes.qc.bed

    CLB_READS=\$(
        bedtools bamtobed -i "${bam}" |
        bedtools intersect -a - -b clb_genes.qc.bed -u |
        cut -f4 |
        sort -u |
        wc -l
    )

    CLB_GENES_DETECTED=\$(awk -F '\t' '
        \$1 ~ /^clb[A-S]\$/ && (\$NF + 0) > 0 { count++ }
        END { print count + 0 }
    ' "${counts}")

    printf "Sample\tMetric\tValue\n" > "${qc}"
    printf "%s\treads_clb_genes_align\t%s\n" "${sampleID}" "\$CLB_READS" >> "${qc}"
    printf "%s\tnum_clb_genes_align\t%s\n" "${sampleID}" "\$CLB_GENES_DETECTED" >> "${qc}"

    if [[ "\$MAPPED_READS" -eq 0 ]]; then
        echo "No confidently mapped reads for ${sampleID}; recording zero clb counts."
        : > "${coverage}"
        : > "${bedtools_cov}"
    else
        bamCoverage \
			--numberOfProcessors "${task.cpus}" \
            -b "${bam}" \
            -o "${coverage}" \
            --normalizeUsing RPKM \
            --outFileFormat bedgraph

        python "${params.scripts}/validate_bedgraph.py" "${coverage}"

        bedtools genomecov \
            -ibam "${bam}" \
            -d \
            > "${bedtools_cov}"
    fi
    """
}
