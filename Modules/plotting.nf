nextflow.enable.dsl=2

process plotPKS {
    label 'sample_stage'
    label 'process_low'
    scratch true
    publishDir "${params.sample_dir}", mode: 'copy',
        // no sampleID in scope here: the circos plot is named <sample>.pks.circos.pdf
        saveAs: { fn -> (fn - '.pks.circos.pdf') + '/figures/circos.pdf' }
    conda "${params.pks_align_env}"

    input:
    path coverage_file

    output:
    path "*.pks.circos.pdf", optional: true

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_cytoband ${params.dep_digest?.pks_cytoband}  scripts ${params.dep_digest?.scripts}
	coverage_basename="\$(basename "${coverage_file}")"
	sample_name="\${coverage_basename%.coverage.bedgraph}"
	output_pdf="\${sample_name}.pks.circos.pdf"

    Rscript "${params.scripts}/plotPKS.R" ${coverage_file} "${params.pks_cytoband}" "\${output_pdf}" \
        "${params.pks_contig}" "${params.pks_shift}" \
        "${params.pks_shift.toString().toInteger() + params.pks_island_len.toString().toInteger()}"
    """
}

process masterTableAlign {
    label 'process_low'
    scratch true
    publishDir "${params.pks_counts_dir}", mode: 'copy'
    conda "${params.pks_align_env}"

    input:
    path(count_files)
    path(manifest)

    output:
    path "pks.gene.counts.align.txt"

	script:
    // The manifest pairs each staged counts file with the sample identifier that came
    // from the sample sheet, so the merge never reconstructs a name (F03). It is also
    // one filename on the command line instead of one per sample (F02).
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
	echo "Building Gene-by-Sample alignment summary"
    python "${params.scripts}/mergeGeneCounts.py" \
        --manifest "${manifest}" \
        --output pks.gene.counts.align.txt
    """
}



process masterTableHMM {
    label 'process_medium'
    scratch true
    publishDir "${params.pks_counts_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    path(count_files)
    path(manifest)

    output:
    path "pks.gene.counts.hmm.txt"

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    python3 "${params.scripts}/build_hmm_matrix.py" \
      --manifest "${manifest}" \
      --out pks.gene.counts.hmm.txt
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
    path(expected_samples)
    path(fragment_list)

    output:
    path "pks.qc.summary.tsv"

    script:
    // F02: one filename per fragment on the command line overflows the OS argument
    // limit on a cohort of this size. The fragments arrive as one list file.
    """
    set -euo pipefail

    python3 "${qc_script}" \
        --inputs-file "${fragment_list}" \
        --expected-samples "${expected_samples}" \
        --output pks.qc.summary.tsv
    """
}


// ─── One page for the whole cohort ────────────────────────────────────────────

process cohortReport {
    label 'process_low'
    publishDir "${params.cohort_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    path(results_marker)
    path(report_script)
    path(expected_samples)

    output:
    path "pks_cohort_report.html", emit: report
    path "pks_cohort_report.tsv",  emit: table

    script:
    // Runs against the published tree rather than staged channels: the report is a
    // view of what the run produced, and reading it back is what makes every value
    // traceable to a file beside it. It derives nothing.
    // The task writes to its own working directory and publishDir moves the files to
    // params.cohort_dir. An earlier version also wrote them into a `cohort/` subdirectory,
    // which put them one level below where the output declarations look: the script
    // exited 0, the report existed, and Nextflow failed the task for a missing output.
    // Ludmil, revised report finding 5: `results_marker` makes the task wait for the
    // channels mixed into it (see main.nf's cohort_report_gate); `expected_samples` is
    // that same wait turned into a filter, so a directory left over from an earlier
    // run at this --outdir can be waited past but never reported on.
    """
    set -euo pipefail
    python3 "${report_script}" \
        --results "${params.outdir}" \
        --output pks_cohort_report.html \
        --table  pks_cohort_report.tsv \
        --expected-samples "${expected_samples}"
    """
}


// ─── One table with every stage joined ────────────────────────────────────────

process masterSummary {
    label 'process_low'
    publishDir "${params.cohort_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    path(results_marker)
    path(report_script)
    path(expected_samples)

    output:
    path "pks.master_summary.tsv", emit: table

    script:
    // Dead-or-unwired-code item d-2: build_master_summary.py was a correct, working
    // script but only ever documented as something you run yourself after a pipeline
    // run. Automated here on the same pattern as cohortReport, one process up:
    // results_marker gates on every optional lane the script actually reads from
    // (see main.nf's master_summary_gate), and expected_samples filters out anything
    // left over from an older run at this --outdir. The script itself is unchanged
    // in what it reads or how -- still safe to re-run by hand against a finished
    // results tree, which is why --expected-samples stayed optional there.
    """
    set -euo pipefail
    python3 "${report_script}" \
        --results "${params.outdir}" \
        --output pks.master_summary.tsv \
        --expected-samples "${expected_samples}"
    """
}
