nextflow.enable.dsl = 2

process pksProfiler_align {
    label 'process_medium'
    scratch true
    publishDir "${params.pks_dir}", mode: 'copy'
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(r1), path(r2)

    output:
    tuple val(sampleID),
          path("${sampleID}.coverage.txt"),
          path("${sampleID}.coverage.bedgraph"),
          path("${sampleID}.counts.txt"),
	      path("${sampleID}.sorted.bam"),
	      path("${sampleID}.sorted.bam.bai"),
	      path("${sampleID}.sam")

    script:
    def bedtools_cov = "${sampleID}.coverage.txt"
    def coverage     = "${sampleID}.coverage.bedgraph"
    def counts       = "${sampleID}.counts.txt"
    def bam          = "${sampleID}.sorted.bam"
	def bai          = "${sampleID}.sorted.bam.bai"
    def sam          = "${sampleID}.sam"

    """
    set -euo pipefail

	# Validate the manuscript-facing clb annotation before processing samples.
    CLB_GENE_COUNT=\$(awk -F '\t' '
        \$0 !~ /^#/ &&
        \$3 == "gene" &&
        \$9 ~ /(^|;)Name=clb[A-S](;|$)/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${params.pks_genome_annotation}")

    if [[ "\$CLB_GENE_COUNT" -ne 19 ]]; then
        echo "ERROR: Expected exactly 19 clb genes in ${params.pks_genome_annotation}; found \$CLB_GENE_COUNT." >&2
        exit 1
    fi

	if ! gzip -t "${r1}" >/dev/null 2>&1; then
        echo "ERROR: Corrupt gzip input for ${sampleID}: ${r1}" >&2
        exit 1
    fi

    if ! gzip -t "${r2}" >/dev/null 2>&1; then
        echo "ERROR: Corrupt gzip input for ${sampleID}: ${r2}" >&2
        exit 1
    fi


    # Merge gzipped FASTQs safely
    zcat "${r1}" "${r2}" | gzip -c > "${sampleID}.trimmed.fastq.gz"

    echo "Bowtie2 Alignment (sample: ${sampleID})"
    bowtie2 -x "${params.pks_genome}" -q "${sampleID}.trimmed.fastq.gz" \\
        --seed 42 --threads 1 --very-sensitive --no-unal -S "${sam}"

    samtools view -bS -q 40 "${sam}" | samtools sort -o "${bam}" -
    samtools index "${bam}" -o "${bai}"

    MAPPED_READS=\$(samtools view -c -F 4 "${bam}")

    # Always create a valid 19-gene count table. A successfully
    # processed sample with no aligned reads must remain in summaries
    # as an explicit zero rather than disappearing as an empty file.
    featureCounts \
        -a "${params.pks_genome_annotation}" \
        -o "${counts}" \
        -t gene \
        -F GFF \
        -g Name \
        "${bam}"

    if [[ "\$MAPPED_READS" -eq 0 ]]; then
        echo "No confidently mapped reads for ${sampleID}; recording zero clb counts."
        : > "${coverage}"
        : > "${bedtools_cov}"
    else
        bamCoverage \
            -b "${bam}" \
            -o "${coverage}" \
            --normalizeUsing RPKM \
            --outFileFormat bedgraph

        bedtools genomecov \
            -ibam "${bam}" \
            -d \
            > "${bedtools_cov}"
    fi
    """
}
