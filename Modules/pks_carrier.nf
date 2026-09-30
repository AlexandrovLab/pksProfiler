nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------------
// V8b -- CARRIER ATTRIBUTION.
//
// Ludmil, 25 Sep report, section 4:
//   "V8b. Carrier attribution. Separately compare genuine pks reads against alternative
//    pks-positive carriers/strains (pks+ E. coli, Klebsiella, Citrobacter, Enterobacter,
//    etc.). This asks which organism is the likely carrier and should not be conflated
//    with locus-detection specificity."
//
// So this lane is deliberately NOT the V8a lane with a different panel. Its panel is
// pks-POSITIVE throughout (ref/specificity_panels/v8b_carrier_panel.tsv, screened by
// scripts/build_specificity_panels.py --panel v8b, which refuses a member that does not
// carry the island). Its outputs are named carrier_attribution, its verdicts are
// CARRIER_ATTRIBUTED / CARRIER_GENUS_ONLY / CARRIER_AMBIGUOUS / CARRIER_UNDETERMINED, and
// none of V8a's words -- decoy, confirmed_true, SURVIVES, specificity -- appear in them.
// A read that fits another carrier better than IHE3034 is evidence about the ORGANISM; it
// is not evidence against the locus, and must never be reported as if it were.
//
// Ludmil, same report, section 5: "Do not present NM as proving uniqueness or carrier
// identity. It is a divergence metric relative to IHE3034." The same caution applies to
// the alignment-score margins here, which is why scripts/carrier_attribution.py reports
// strain-level and genus-level attribution side by side and never prints one alone.
//
// Optional and off by default: main.nf runs this lane only under
// --pks_carrier_attribution true.
// ---------------------------------------------------------------------------------

params.pks_carrier_panel_index = null   // minimap2 index/FASTA: v8b_panel.mmi, the pks+ carrier panel
params.pks_carrier_map         = null   // v8b_carrier_map.tsv: contig -> carrier key/organism/accession

// Thresholds. Placeholders pending calibration, configurable here, never hardcoded in
// the script; defaults mirror scripts/carrier_attribution.py's own argparse defaults.
params.pks_carrier_min_delta_as            = 1.0
params.pks_carrier_tie_margin              = 0.0
params.pks_carrier_min_attributed_reads    = 10
params.pks_carrier_min_top_carrier_fraction = 0.8

params.scripts = params.scripts ?: "${moduleDir}/../scripts"
params.pks_carrier_minimap2_env = params.minimap2_env ?: "${moduleDir}/../conda_envs/minimap2_env.yml"
params.pks_align_env = params.pks_align_env ?: "${moduleDir}/../conda_envs/pks_align_env.yml"

def pksCarrierAttributionValidateParams() {
    def missing = []
    if (!params.pks_carrier_panel_index) missing << "--pks_carrier_panel_index (v8b_panel.mmi: the pks-positive carrier panel index)"
    if (!params.pks_carrier_map)         missing << "--pks_carrier_map (v8b_carrier_map.tsv: which carrier each panel contig belongs to)"
    if (missing) {
        log.error "V8b carrier-attribution lane: missing required parameter(s):"
        missing.each { log.error "  ${it}" }
        exit 1, "The V8b carrier-attribution lane cannot run without its carrier panel. Build it with " +
                 "scripts/build_specificity_panels.py --panel v8b (see " +
                 "ref/specificity_panels/v8b_carrier_panel.tsv for the candidate carriers)."
    }
}

process alignCarrierPanel {
    label 'sample_stage'
    label 'process_medium'
    tag "$sampleID"
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/carrier_attribution" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_carrier_minimap2_env}"

    input:
    tuple val(sampleID), path(reads_fastq)

    output:
    tuple val(sampleID), path("${sampleID}.v8b_carrier.sam"), emit: sam

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # v8b_panel ${params.dep_digest?.pks_carrier_panel}
    set -euo pipefail
    # Every panel member carries the island, so a read's alignments to different carriers
    # are near-identical by construction and minimap2's default -p 0.8 secondary filter
    # would drop most of them -- leaving one arbitrary winner and no runner-up to measure
    # a margin against. -N 20 -p 0.1 keeps one best alignment per carrier visible, which
    # is exactly what carrier attribution needs and what locus specificity does not.
    minimap2 -a -x sr -N 20 -p 0.1 --secondary=yes -t "${task.cpus}" \
        "${params.pks_carrier_panel_index}" "${reads_fastq}" \
        > "${sampleID}.v8b_carrier.sam"
    """
}

process attributeCarrier {
    label 'sample_stage'
    tag "$sampleID"
    publishDir { "${params.sample_dir}/${sampleID}/carrier_attribution" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(carrier_sam)

    output:
    tuple val(sampleID),
          path("${sampleID}.carrier_attribution.tsv"),
          path("${sampleID}.carrier_per_read.tsv"), emit: attribution

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    python3 "${params.scripts}/carrier_attribution.py" \
        --sample "${sampleID}" \
        --sam "${carrier_sam}" \
        --carrier-map "${params.pks_carrier_map}" \
        --min-delta-as "${params.pks_carrier_min_delta_as}" \
        --tie-margin "${params.pks_carrier_tie_margin}" \
        --min-attributed-reads "${params.pks_carrier_min_attributed_reads}" \
        --min-top-carrier-fraction "${params.pks_carrier_min_top_carrier_fraction}" \
        --output-summary "${sampleID}.carrier_attribution.tsv" \
        --output-per-read "${sampleID}.carrier_per_read.tsv"
    """
}

workflow pksCarrierAttribution {
    take:
    supporting_reads   // tuple(sampleID, pks.fastq.gz) -- extractPksIslandReads' emit in pks_taxa.nf

    main:
    pksCarrierAttributionValidateParams()
    alignCarrierPanel(supporting_reads)
    attributeCarrier(alignCarrierPanel.out.sam)

    emit:
    attribution = attributeCarrier.out.attribution
}

// ---------------------------------------------------------------------------------
// Standalone entry point, for attributing an existing cohort's retained island reads
// without re-running the pipeline:
//   nextflow run Modules/pks_carrier.nf --sample sheet.csv --sample_dir out/ \
//       --pks_carrier_panel_index v8b_panel.mmi --pks_carrier_map v8b_carrier_map.tsv
// sheet.csv columns: sample,reads_fastq
// ---------------------------------------------------------------------------------
workflow {
    if (!params.sample) {
        exit 1, "Missing required parameter: --sample <sheet.csv> (columns: sample,reads_fastq)"
    }
    if (!params.sample_dir) {
        exit 1, "Missing required parameter: --sample_dir <output directory>"
    }

    supporting_reads = Channel.fromPath(params.sample, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.reads_fastq, checkIfExists: true)) }

    pksCarrierAttribution(supporting_reads)
}
