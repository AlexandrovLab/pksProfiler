nextflow.enable.dsl=2

process plotPKS {

    label 'process_low'
    scratch true
    publishDir "${params.pks_coverage_plots_dir}", mode: 'copy'
    conda "${params.pks_align_env}"

    input:
    path coverage_file

    output:
    path "*.pks.circos.pdf", optional: true

    script:
    """
	coverage_basename="\$(basename "${coverage_file}")"
	sample_name="\${coverage_basename%.coverage.bedgraph}"
	output_pdf="\${sample_name}.pks.circos.pdf"

    Rscript "${params.scripts}/plotPKS.R" ${coverage_file} "${params.pks_cytoband}" "\${output_pdf}"
    """
}

process masterTableAlign {
    label 'process_low'
    scratch true
    publishDir "${params.pks_counts_dir}", mode: 'copy'
    conda "${params.pks_align_env}"

    input:
    path(count_files)

    output:
    path "pks.gene.counts.align.txt"

	script:
    def inputs = count_files.collect { file -> file.toString() }.join(' ')


    """
	echo "Building Gene-by-Sample alignment summary"
    python "${params.scripts}/mergeGeneCounts.py" ${inputs} pks.gene.counts.align.txt
    """
}



process masterTableHMM {
    label 'process_medium'
    scratch true
    publishDir "${params.pks_counts_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    path(count_files)

    output:
    path "pks.gene.counts.hmm.txt"

    script:
    def inputs = count_files.collect { file -> file.toString() }.join(' ')

    """
    python3 "${params.scripts}/build_hmm_matrix.py" \
      --inputs ${inputs} \
      --out pks.gene.counts.hmm.txt \
      --strip-suffix ".hmm_counts.tsv"
    """
}

process masterQCSummary {
    label 'process_low'
    scratch true
    publishDir "${params.pks_qc_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    path(qc_files)
    path(qc_script)

    output:
    path "pks.qc.summary.tsv"

    script:
    def inputs = qc_files.collect { file -> "\"${file}\"" }.join(' ')

    """
    set -euo pipefail

    python3 "${qc_script}" \
        --inputs ${inputs} \
        --output pks.qc.summary.tsv
    """
}
