nextflow.enable.dsl = 2

// ─── MLST on an assembled unit (a MAG bin, or a tumour contig set) ──────────────
process mlstTypeAssembly {
    tag "${sampleID}:${label}"
    label 'process_low'
    publishDir { "${params.sample_dir}/${sampleID}/strain/${label}" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.mlst_env}"

    input:
    tuple val(sampleID), val(label), path(assembly)

    output:
    tuple val(sampleID), val(label), path("${label}.mlst.tsv"), emit: mlst

    script:
    """
    set -euo pipefail
    # An assembly below the floor cannot carry all seven housekeeping loci; emit an
    # empty result so downstream reports "insufficient" rather than skipping silently.
    bases=\$(grep -v '^>' "${assembly}" | tr -d '\\n' | wc -c)
    if [[ "\$bases" -lt "${params.typing_min_assembly_bp}" ]]; then
        : > "${label}.mlst.tsv"
    else
        mlst --scheme "${params.mlst_scheme}" --quiet "${assembly}" > "${label}.mlst.tsv"
    fi
    """
}

// ─── Join the ST to clonal complex and phylogroup from the reference panel ──────
process assignStrainType {
    tag "${sampleID}:${label}"
    label 'process_low'
    publishDir { "${params.sample_dir}/${sampleID}/strain/${label}" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), val(label), path(mlst_tsv)

    output:
    tuple val(sampleID), val(label), path("${label}.strain_type.tsv"), emit: typing

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    python "${params.scripts}/assign_strain_type.py" \\
        --mlst "${mlst_tsv}" \\
        --lookup "${params.st_phylogroup_lookup}" \\
        --sample "${sampleID}" \\
        --label "${label}" \\
        --output "${label}.strain_type.tsv"
    """
}

// ─── Collect every typed unit for a sample into one table ───────────────────────
process strainTypeSummary {
    tag "$sampleID"
    label 'process_low'
    publishDir { "${params.sample_dir}/${sampleID}/strain" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(typing_tables)

    output:
    tuple val(sampleID), path("${sampleID}.strain_types.tsv"), emit: summary

    script:
    """
    set -euo pipefail
    head -n 1 \$(ls ${typing_tables} | head -n 1) > "${sampleID}.strain_types.tsv"
    for f in ${typing_tables}; do
        tail -n +2 "\$f" >> "${sampleID}.strain_types.tsv"
    done
    """
}

// ─── Entry point: type a channel of assembled units ─────────────────────────────
workflow strainTyping {
    take:
    units   // tuple val(sampleID), val(label), path(assembly)

    main:
    mlstTypeAssembly(units)
    assignStrainType(mlstTypeAssembly.out.mlst)

    per_sample_ch = assignStrainType.out.typing
        .map { sampleID, _label, table -> tuple(sampleID, table) }
        .groupTuple(by: 0)
    strainTypeSummary(per_sample_ch)

    emit:
    typing  = assignStrainType.out.typing
    summary = strainTypeSummary.out.summary
}
