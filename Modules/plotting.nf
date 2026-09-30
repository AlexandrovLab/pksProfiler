nextflow.enable.dsl=2

process plotPKS {
    label 'sample_stage'
    label 'process_low'
    scratch true
    // Ludmil, 2026-09-25 report P0-5b: the sample identifier travels with the plot so
    // cohortReport can be handed the figure as a declared input. It used to arrive as a
    // bare path and the sample was recovered from the filename, here and in saveAs.
    publishDir { "${params.sample_dir}/${sampleID}/figures" }, mode: 'copy',
        saveAs: { "circos.pdf" }
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(coverage_file)

    output:
    tuple val(sampleID), path("*.pks.circos.pdf"), optional: true

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_cytoband ${params.dep_digest?.pks_cytoband}  scripts ${params.dep_digest?.scripts}
	# F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
	${ToolVersions.record(task, params, 'Rscript:r-base')}
	coverage_basename="\$(basename "${coverage_file}")"
	sample_name="\${coverage_basename%.coverage.bedgraph}"
	output_pdf="\${sample_name}.pks.circos.pdf"

    # P1 (25 Sep): the island interval is the canonical one, converted to the
    # bedgraph's 0-based half-open convention -- not pks_shift/pks_island_len.

    Rscript "${params.scripts}/plotPKS.R" ${coverage_file} "${params.pks_cytoband}" "\${output_pdf}" \
        "${params.pks_contig}" \
        "\$(( ${params.pks_start_1based} - 1 ))" "${params.pks_end_1based}"
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
	# F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
	${ToolVersions.record(task, params, 'python')}
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python3:python')}
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
    val(expected_stages)

    output:
    path "pks.qc.summary.tsv"

    script:
    // F02: one filename per fragment on the command line overflows the OS argument
    // limit on a cohort of this size. The fragments arrive as one list file.
    """
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python3:python')}

    python3 "${qc_script}" \
        --inputs-file "${fragment_list}" \
        --expected-samples "${expected_samples}" \
        --expected-stages "${expected_stages}" \
        --output pks.qc.summary.tsv
    """
}


// ─── One page for the whole cohort ────────────────────────────────────────────

process cohortReport {
    label 'process_low'
    publishDir "${params.cohort_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    path(report_inputs, stageAs: 'inputs/*/*')
    path(inputs_manifest)
    path(report_script)
    path(expected_samples)

    output:
    path "pks_cohort_report.html", emit: report
    path "pks_cohort_report.tsv",  emit: table

    script:
    // Ludmil, 2026-09-25 report P0-5b: the report reads its declared inputs, not the
    // published tree. It used to run `--results ${params.outdir}` and scan by_sample/;
    // the gate made the task wait for every producer, but Nextflow publishes
    // asynchronously, so a published copy could still be missing when it ran, and a
    // file an earlier run left at the same --outdir could be read in its place.
    // `report_inputs` is every file the page displays, each staged into its own
    // inputs/<n>/ directory -- per-sample files share basenames across samples once
    // published, and cohort-level ones could collide with a sample's -- and
    // `inputs_manifest` names each by sample and kind (main.nf's reportManifestText).
    // Figures are staged only so their existence comes from this run; the page links
    // their published copies and never opens them.
    // The task writes to its own working directory and publishDir moves the files to
    // params.cohort_dir. An earlier version also wrote them into a `cohort/` subdirectory,
    // which put them one level below where the output declarations look: the script
    // exited 0, the report existed, and Nextflow failed the task for a missing output.
    // Ludmil, revised report finding 5: `expected_samples` is the run's sample list; the
    // report has exactly those rows. 26 Sep mode memo section 6: --sample-type gives every
    // row its mode, including one whose evidence file never arrived.
    """
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python3:python')}
    python3 "${report_script}" \
        --inputs-manifest "${inputs_manifest}" \
        --output pks_cohort_report.html \
        --table  pks_cohort_report.tsv \
        --sample-type "${params.sample_type}" \
        --expected-samples "${expected_samples}"
    """
}


// ─── One table with every stage joined ────────────────────────────────────────

process masterSummary {
    label 'process_low'
    publishDir "${params.cohort_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    path(summary_inputs, stageAs: 'inputs/*/*')
    path(inputs_manifest)
    path(report_script)
    path(expected_samples)

    output:
    path "pks.master_summary.tsv", emit: table

    script:
    // Dead-or-unwired-code item d-2: build_master_summary.py was a correct, working
    // script but only ever documented as something you run yourself after a pipeline
    // run. Automated here on the same pattern as cohortReport, one process up.
    // Ludmil, 2026-09-25 report P0-5b, mirrored from cohortReport: every file the
    // script reads is staged (`summary_inputs`, one inputs/<n>/ directory each) and
    // named in `inputs_manifest`, so it never reads the published tree. --results is
    // still there for hand runs against a finished results tree; the pipeline never
    // passes it. expected_samples filters to the run's own samples.
    """
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python3:python')}
    python3 "${report_script}" \
        --inputs-manifest "${inputs_manifest}" \
        --output pks.master_summary.tsv \
        --sample-type "${params.sample_type}" \
        --expected-samples "${expected_samples}"
    """
}
