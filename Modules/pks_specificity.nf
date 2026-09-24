nextflow.enable.dsl = 2

// Optional, on-demand secondary check on top of the pipeline's own MAPQ>=40-vs-IHE3034
// read-level pks call (classifyPksReadEvidence in pks_targeted.nf). Not part of the
// default per-sample pipeline and not included from main.nf: the reference panels this
// module depends on (--pks_specificity_true_index/decoy_index/concatenated_index/
// true_contigs/decoy_contigs, and a combined-panel KrakenUniq DB for Arm A) are Phase 2
// deliverables that do not exist yet (see scripts/specificity_decision.py's module
// docstring and this repo's CLAUDE.md for the phase plan). Wiring an unconditional
// import into main.nf for indices that cannot exist yet would only add new failure
// surface to the default pipeline entry point for a feature nobody can run; this module
// is invoked standalone instead: `nextflow run Modules/pks_specificity.nf --sample ...`.
//
// Ludmil's concern this answers: MAPQ >= 40 against IHE3034 alone shows confident
// mapping to *that* reference, not uniqueness against the broader microbial universe.
// clbB/clbC/clbH/clbI/clbJ/clbK/clbN/clbO carry NRPS/PKS megasynthase domains shared
// with unrelated bacterial secondary-metabolite clusters -- the same domain-homology
// failure M1 already caught once at the assembled-MAG-bin level (see
// tests/test_pks_positive_definition.py). This checks for it at the read level, for
// tumour samples too sparse to ever assemble a bin, on exactly the read set that
// already passed the pipeline's own filter -- never a re-run of primary detection.
//
// Two independent arms decide, and full method/threshold documentation lives in
// scripts/competitive_realign_score.py (Arm B) and scripts/specificity_decision.py
// (Arm A + the decision rule), not duplicated here.

// ---------------------------------------------------------------------------------
// Phase 2 inputs. All null until the reference panel exists; each is validated below
// with an error naming exactly what Phase 2 must produce.
// ---------------------------------------------------------------------------------
params.pks_specificity_true_index          = null  // minimap2 index/FASTA: true pks+ reference(s) alone
params.pks_specificity_decoy_index         = null  // minimap2 index/FASTA: pks-negative decoy panel alone
params.pks_specificity_concatenated_index  = null  // minimap2 index/FASTA: true + decoy panel, one index
params.pks_specificity_true_contigs        = null  // text file, one contig name per line: which contigs
                                                    // in the CONCATENATED index are the true reference
params.pks_specificity_decoy_contigs       = null  // text file, one contig name per line: which contigs
                                                    // in the CONCATENATED index are decoys
params.pks_specificity_arm_a_target_clade  = null  // exact Arm-A assigned_clade value that counts as a
                                                    // specific pks+ resolution (not genus/family-collapsed);
                                                    // depends on how Phase 2's KrakenUniq DB names it

// ---------------------------------------------------------------------------------
// Decision-rule thresholds. Placeholders pending Phase 3 calibration against real
// positive/negative controls -- configurable here, never hardcoded in the scripts.
// Defaults mirror scripts/specificity_decision.py's own argparse defaults exactly.
// ---------------------------------------------------------------------------------
params.pks_specificity_arm_a_min_confidence      = 0.65
params.pks_specificity_arm_a_min_unique_kmers    = 3
params.pks_specificity_retention_threshold       = 0.8
params.pks_specificity_retention_ci_lower_threshold = 0.6
params.pks_specificity_bootstrap_iterations      = 1000
params.pks_specificity_bootstrap_seed            = 0
params.pks_specificity_max_tier_drop             = 1
params.pks_specificity_non_megasynthase_min_genes = 2
params.pks_specificity_absorber_fraction_threshold = 0.2
params.pks_specificity_absorber_near_threshold_margin = 5.0
params.pks_specificity_min_supporting_reads      = 10
params.pks_specificity_min_delta_as              = 0.0

// Tier thresholds: must match classifyPksReadEvidence's own params exactly, or tier
// recomputation silently uses a different definition than the tier it is comparing
// against. Not redeclared with different names/defaults -- these ARE
// params.tumor_*_min_* from main.nf when this module is run alongside the rest of the
// pipeline; the standalone entry point below lets them be overridden the same way.
params.pks_specificity_multi_gene_reads       = 5
params.pks_specificity_multi_gene_genes       = 3
params.pks_specificity_multi_gene_breadth     = .01
params.pks_specificity_broad_island_reads     = 30
params.pks_specificity_broad_island_genes     = 8
params.pks_specificity_broad_island_breadth   = .075
params.pks_specificity_extensive_island_reads = 100
params.pks_specificity_extensive_island_genes = 10
params.pks_specificity_extensive_island_breadth = .15
params.pks_specificity_island_length          = null  // required: must match --pks_island_len_1based

// Defaulted here (rather than assumed from main.nf) so this module also runs when
// invoked standalone, without main.nf's own param declarations in scope.
params.scripts = params.scripts ?: "${moduleDir}/../scripts"
params.pks_specificity_minimap2_env = params.minimap2_env ?: "${moduleDir}/../conda_envs/minimap2_env.yml"
params.pks_align_env = params.pks_align_env ?: "${moduleDir}/../conda_envs/pks_align_env.yml"

def pksSpecificityValidateParams() {
    def missing = []
    if (!params.pks_specificity_true_index)         missing << "--pks_specificity_true_index (Phase 2: minimap2 index of the true pks+ reference(s) alone)"
    if (!params.pks_specificity_decoy_index)        missing << "--pks_specificity_decoy_index (Phase 2: minimap2 index of the pks-negative decoy panel alone)"
    if (!params.pks_specificity_concatenated_index) missing << "--pks_specificity_concatenated_index (Phase 2: minimap2 index of true reference + decoy panel combined)"
    if (!params.pks_specificity_true_contigs)       missing << "--pks_specificity_true_contigs (Phase 2: text file listing which contigs in the concatenated index are the true reference)"
    if (!params.pks_specificity_decoy_contigs)      missing << "--pks_specificity_decoy_contigs (Phase 2: text file listing which contigs in the concatenated index are decoys)"
    if (!params.pks_specificity_arm_a_target_clade) missing << "--pks_specificity_arm_a_target_clade (Phase 2/3: the exact Arm-A assigned_clade value that counts as a specific pks+ resolution)"
    if (!params.pks_specificity_island_length)      missing << "--pks_specificity_island_length (must match --pks_island_len_1based)"
    if (missing) {
        log.error "pks_specificity module: missing required parameter(s):"
        missing.each { log.error "  ${it}" }
        exit 1, "pks_specificity module cannot run without the Phase 2 reference panel " +
                 "(see Modules/pks_specificity.nf and scripts/specificity_decision.py for what each parameter is)."
    }
}

process alignConcatenatedIndex {
    label 'sample_stage'
    label 'process_medium'
    tag "$sampleID"
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/specificity_check" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_specificity_minimap2_env}"

    input:
    tuple val(sampleID), path(reads_fastq)

    output:
    tuple val(sampleID), path("${sampleID}.b1_concatenated.sam"), emit: sam

    script:
    """
    set -euo pipefail
    # -p 0.1 (default is 0.8): a decoy-side hit far weaker than the primary would
    # otherwise be dropped as "not a real secondary candidate", but Arm B needs
    # AS_best(decoy) reported whenever any alignment to the decoy family exists at
    # all, not only when it is competitive. -N 5 is comfortably above a 2-family
    # (true vs decoy panel) best-per-family need.
    minimap2 -a -x sr -N 5 -p 0.1 --secondary=yes -t "${task.cpus}" \
        "${params.pks_specificity_concatenated_index}" "${reads_fastq}" \
        > "${sampleID}.b1_concatenated.sam"
    """
}

process alignSeparateIndices {
    label 'sample_stage'
    label 'process_medium'
    tag "$sampleID"
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/specificity_check" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_specificity_minimap2_env}"

    input:
    tuple val(sampleID), path(reads_fastq)

    output:
    tuple val(sampleID), path("${sampleID}.b2_true.sam"), path("${sampleID}.b2_decoy.sam"), emit: sam

    script:
    """
    set -euo pipefail
    minimap2 -a -x sr -t "${task.cpus}" "${params.pks_specificity_true_index}" "${reads_fastq}" \
        > "${sampleID}.b2_true.sam"
    minimap2 -a -x sr -t "${task.cpus}" "${params.pks_specificity_decoy_index}" "${reads_fastq}" \
        > "${sampleID}.b2_decoy.sam"
    """
}

process scoreCompetitiveRealignment {
    label 'sample_stage'
    tag "$sampleID"
    publishDir { "${params.sample_dir}/${sampleID}/specificity_check" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_specificity_minimap2_env}"

    input:
    tuple val(sampleID), path(b1_sam), path(b2_true_sam), path(b2_decoy_sam)

    output:
    tuple val(sampleID), path("${sampleID}.arm_b.tsv"), emit: arm_b

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    python3 "${params.scripts}/competitive_realign_score.py" \
        --b1-sam "${b1_sam}" \
        --b2-true-sam "${b2_true_sam}" \
        --b2-decoy-sam "${b2_decoy_sam}" \
        --true-contigs "${params.pks_specificity_true_contigs}" \
        --decoy-contigs "${params.pks_specificity_decoy_contigs}" \
        --min-delta-as "${params.pks_specificity_min_delta_as}" \
        --output "${sampleID}.arm_b.tsv"
    """
}

process computeSpecificityVerdict {
    label 'sample_stage'
    tag "$sampleID"
    publishDir { "${params.sample_dir}/${sampleID}/specificity_check" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(read_evidence_tsv), path(read_gene_tsv), path(arm_a_tsv), path(arm_b_tsv)

    output:
    tuple val(sampleID), path("${sampleID}.specificity_summary.tsv"), path("${sampleID}.specificity_per_read.tsv"), emit: verdict

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    python3 "${params.scripts}/specificity_decision.py" \
        --sample "${sampleID}" \
        --read-evidence "${read_evidence_tsv}" \
        --read-gene-assignment "${read_gene_tsv}" \
        --arm-a "${arm_a_tsv}" \
        --arm-b "${arm_b_tsv}" \
        --arm-a-target-clade "${params.pks_specificity_arm_a_target_clade}" \
        --arm-a-min-confidence "${params.pks_specificity_arm_a_min_confidence}" \
        --arm-a-min-unique-kmers "${params.pks_specificity_arm_a_min_unique_kmers}" \
        --island-length "${params.pks_specificity_island_length}" \
        --multi-gene-reads "${params.pks_specificity_multi_gene_reads}" \
        --multi-gene-genes "${params.pks_specificity_multi_gene_genes}" \
        --multi-gene-breadth "${params.pks_specificity_multi_gene_breadth}" \
        --broad-island-reads "${params.pks_specificity_broad_island_reads}" \
        --broad-island-genes "${params.pks_specificity_broad_island_genes}" \
        --broad-island-breadth "${params.pks_specificity_broad_island_breadth}" \
        --extensive-island-reads "${params.pks_specificity_extensive_island_reads}" \
        --extensive-island-genes "${params.pks_specificity_extensive_island_genes}" \
        --extensive-island-breadth "${params.pks_specificity_extensive_island_breadth}" \
        --retention-threshold "${params.pks_specificity_retention_threshold}" \
        --retention-ci-lower-threshold "${params.pks_specificity_retention_ci_lower_threshold}" \
        --bootstrap-iterations "${params.pks_specificity_bootstrap_iterations}" \
        --bootstrap-seed "${params.pks_specificity_bootstrap_seed}" \
        --max-tier-drop "${params.pks_specificity_max_tier_drop}" \
        --non-megasynthase-min-genes "${params.pks_specificity_non_megasynthase_min_genes}" \
        --absorber-fraction-threshold "${params.pks_specificity_absorber_fraction_threshold}" \
        --absorber-near-threshold-margin "${params.pks_specificity_absorber_near_threshold_margin}" \
        --min-supporting-reads "${params.pks_specificity_min_supporting_reads}" \
        --output-summary "${sampleID}.specificity_summary.tsv" \
        --output-per-read "${sampleID}.specificity_per_read.tsv"
    """
}

workflow pksSpecificityCheck {
    take:
    supporting_reads      // tuple(sampleID, pks.fastq.gz) -- extractPksIslandReads' emit in pks_taxa.nf
    read_evidence         // tuple(sampleID, read_evidence.tsv) -- classifyPksReadEvidence's emit
    read_gene_assignment  // tuple(sampleID, read_clb_gene.tsv) -- extractPksIslandReads' emit
    arm_a_results         // tuple(sampleID, arm_a.tsv) -- Phase 2/3: a real KrakenUniq run against the
                          // combined panel, parsed into read_id/assigned_clade/confidence/
                          // unique_kmer_count. Not run by this module: deriving a real per-read
                          // confidence and unique-kmer floor from KrakenUniq's raw per-read output
                          // is exactly the kind of thing Phase 3 calibrates, not something to freeze
                          // into pipeline code before a single real panel exists. This workflow
                          // consumes the TSV directly (see scripts/specificity_decision.py's
                          // docstring for the exact schema) so Phase 2/3 can produce it however it
                          // is ultimately built.

    main:
    pksSpecificityValidateParams()

    alignConcatenatedIndex(supporting_reads)
    alignSeparateIndices(supporting_reads)

    b1_and_b2 = alignConcatenatedIndex.out.sam.join(alignSeparateIndices.out.sam, by: 0)
        .map { sampleID, b1_sam, b2_true_sam, b2_decoy_sam -> tuple(sampleID, b1_sam, b2_true_sam, b2_decoy_sam) }
    scoreCompetitiveRealignment(b1_and_b2)

    decision_input = read_evidence
        .join(read_gene_assignment, by: 0)
        .join(arm_a_results, by: 0)
        .join(scoreCompetitiveRealignment.out.arm_b, by: 0)
        .map { sampleID, evidence, geneAssign, armA, armB -> tuple(sampleID, evidence, geneAssign, armA, armB) }
    computeSpecificityVerdict(decision_input)

    emit:
    arm_b   = scoreCompetitiveRealignment.out.arm_b
    verdict = computeSpecificityVerdict.out.verdict
}

// ---------------------------------------------------------------------------------
// Standalone entry point: `nextflow run Modules/pks_specificity.nf --sample sheet.csv
// --sample_dir out/ [threshold overrides...]`. Not reachable from main.nf. sheet.csv
// columns: sample,reads_fastq,read_evidence,read_gene_assignment,arm_a_tsv
// ---------------------------------------------------------------------------------
workflow {
    if (!params.sample) {
        exit 1, "Missing required parameter: --sample <sheet.csv> " +
                 "(columns: sample,reads_fastq,read_evidence,read_gene_assignment,arm_a_tsv)"
    }
    if (!params.sample_dir) {
        exit 1, "Missing required parameter: --sample_dir <output directory>"
    }

    sheet = Channel.fromPath(params.sample, checkIfExists: true)
        .splitCsv(header: true)

    supporting_reads = sheet.map { row -> tuple(row.sample, file(row.reads_fastq, checkIfExists: true)) }
    read_evidence = sheet.map { row -> tuple(row.sample, file(row.read_evidence, checkIfExists: true)) }
    read_gene_assignment = sheet.map { row -> tuple(row.sample, file(row.read_gene_assignment, checkIfExists: true)) }
    arm_a_results = sheet.map { row -> tuple(row.sample, file(row.arm_a_tsv, checkIfExists: true)) }

    pksSpecificityCheck(supporting_reads, read_evidence, read_gene_assignment, arm_a_results)
}
