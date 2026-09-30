nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------------
// V8a -- LOCUS-DETECTION SPECIFICITY.
//
// Ludmil, 25 Sep report, section 4:
//   "V8a. Locus-detection specificity. Challenge pks-supporting reads against
//    pks-negative Enterobacterales, close homologous NRPS/PKS sequences, broad gut/oral
//    backgrounds, GRCh38, and UniVec. This asks whether the sequence really belongs to
//    the pks/clb locus."
//
// That is the ONLY question this lane answers. The separate question -- which organism
// carries the island -- is V8b, Modules/pks_carrier.nf, against a panel of pks-POSITIVE
// carriers, with its own vocabulary and its own outputs. The report explicitly says the
// two "should not be conflated", so nothing here says "carrier" and nothing there says
// "specificity", and neither one's verdict is ever folded into the other's column.
//
// Optional and off by default: main.nf runs this lane only under
// --pks_locus_specificity true, and only when the V8a panel paths are given. Build the
// panel with scripts/build_specificity_panels.py --panel v8a (its member list and
// rationale are ref/specificity_panels/v8a_locus_decoy_panel.tsv); every member of that
// panel is screened against the canonical island and refused if it carries clb sequence.
//
// The concern this answers: MAPQ >= 40 against IHE3034 alone shows confident mapping to
// *that* reference, not uniqueness against the broader microbial universe.
// clbB/clbC/clbH/clbI/clbJ/clbK/clbN/clbO carry NRPS/PKS megasynthase domains shared
// with unrelated bacterial secondary-metabolite clusters -- the same domain-homology
// failure M1 already caught once at the assembled-MAG-bin level (see
// tests/test_pks_positive_definition.py). This checks for it at the read level, for
// tumour samples too sparse to ever assemble a bin, on exactly the read set that already
// passed the pipeline's own filter -- never a re-run of primary detection.
//
// Two independent arms decide, and full method/threshold documentation lives in
// scripts/competitive_realign_score.py (Arm B) and scripts/specificity_decision.py
// (Arm A + the decision rule), not duplicated here.
// ---------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------
// V8a panel inputs. All null until the panel is built; each is validated below with an
// error naming exactly what scripts/build_specificity_panels.py --panel v8a writes.
// ---------------------------------------------------------------------------------
params.pks_specificity_true_index          = null  // minimap2 index/FASTA: true pks+ reference(s) alone
params.pks_specificity_decoy_index         = null  // minimap2 index/FASTA: V8a decoy panel alone
params.pks_specificity_concatenated_index  = null  // minimap2 index/FASTA: true + decoy panel, one index
params.pks_specificity_true_contigs        = null  // text file, one contig name per line: which contigs
                                                    // in the CONCATENATED index are the true reference
params.pks_specificity_decoy_contigs       = null  // text file, one contig name per line: which contigs
                                                    // in the CONCATENATED index are V8a decoys

// Arm A is OPTIONAL. It needs a KrakenUniq database over the same true+decoy taxon set
// and a per-sample reclassification TSV; when --pks_specificity_arm_a_dir names a
// directory holding <sample>.arm_a.tsv, both arms must confirm a read, which is the
// stricter rule. Without it the lane runs Arm B alone and says so in its output
// (arms_used=arm_b_only), rather than silently reporting a two-arm verdict it did not
// compute, or refusing to run at all.
params.pks_specificity_arm_a_dir           = null
params.pks_specificity_arm_a_target_clade  = null  // required only with --pks_specificity_arm_a_dir:
                                                    // the exact assigned_clade value that counts as a
                                                    // specific pks+ resolution (not genus/family-collapsed)

// ---------------------------------------------------------------------------------
// Decision-rule thresholds. Placeholders pending calibration against real positive/
// negative controls -- configurable here, never hardcoded in the scripts. Defaults
// mirror scripts/specificity_decision.py's own argparse defaults exactly.
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

// Staged in place of a real Arm A table when Arm A was not run. Its name is what
// computeLocusSpecificityVerdict switches on, so an arm-B-only run is a visible, named
// state rather than an absent input Nextflow would have no channel element for.
def pksSpecificityNoArmA() {
    return file("${moduleDir}/../ref/specificity_panels/NO_ARM_A.tsv", checkIfExists: true)
}

def pksLocusSpecificityValidateParams() {
    def missing = []
    if (!params.pks_specificity_true_index)         missing << "--pks_specificity_true_index (v8a_true_index.mmi: the true pks+ reference alone)"
    if (!params.pks_specificity_decoy_index)        missing << "--pks_specificity_decoy_index (v8a_panel.mmi: the pks-negative V8a decoy panel alone)"
    if (!params.pks_specificity_concatenated_index) missing << "--pks_specificity_concatenated_index (v8a_locus_concatenated.mmi: true reference + V8a decoy panel in one index)"
    if (!params.pks_specificity_true_contigs)       missing << "--pks_specificity_true_contigs (v8a_true_contigs.txt: which contigs in the concatenated index are the true reference)"
    if (!params.pks_specificity_decoy_contigs)      missing << "--pks_specificity_decoy_contigs (v8a_decoy_contigs.txt: which contigs in the concatenated index are decoys)"
    if (!params.pks_specificity_island_length)      missing << "--pks_specificity_island_length (must match --pks_island_len_1based)"
    if (params.pks_specificity_arm_a_dir && !params.pks_specificity_arm_a_target_clade) {
        missing << "--pks_specificity_arm_a_target_clade (required with --pks_specificity_arm_a_dir: the exact Arm-A assigned_clade value that counts as a specific pks+ resolution)"
    }
    if (missing) {
        log.error "V8a locus-specificity lane: missing required parameter(s):"
        missing.each { log.error "  ${it}" }
        exit 1, "The V8a locus-specificity lane cannot run without its decoy panel. Build it with " +
                 "scripts/build_specificity_panels.py --panel v8a (see " +
                 "ref/specificity_panels/v8a_locus_decoy_panel.tsv for what goes in it)."
    }
}

process alignConcatenatedIndex {
    label 'sample_stage'
    label 'process_medium'
    tag "$sampleID"
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/locus_specificity" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_specificity_minimap2_env}"

    input:
    tuple val(sampleID), path(reads_fastq)

    output:
    tuple val(sampleID), path("${sampleID}.b1_concatenated.sam"), emit: sam

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # v8a_panel ${params.dep_digest?.pks_specificity_panel}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'minimap2')}
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
    publishDir { "${params.sample_dir}/${sampleID}/locus_specificity" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_specificity_minimap2_env}"

    input:
    tuple val(sampleID), path(reads_fastq)

    output:
    tuple val(sampleID), path("${sampleID}.b2_true.sam"), path("${sampleID}.b2_decoy.sam"), emit: sam

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # v8a_panel ${params.dep_digest?.pks_specificity_panel}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'minimap2')}
    minimap2 -a -x sr -t "${task.cpus}" "${params.pks_specificity_true_index}" "${reads_fastq}" \
        > "${sampleID}.b2_true.sam"
    minimap2 -a -x sr -t "${task.cpus}" "${params.pks_specificity_decoy_index}" "${reads_fastq}" \
        > "${sampleID}.b2_decoy.sam"
    """
}

process scoreCompetitiveRealignment {
    label 'sample_stage'
    tag "$sampleID"
    publishDir { "${params.sample_dir}/${sampleID}/locus_specificity" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    // pks_align_env, not the minimap2 env of the two alignment steps above: this runs
    // python3 (and samtools through it), and minimap2_env.yml has no python, so python3
    // came from whatever the launching shell had on PATH. Found by the F16 tool-version
    // record (30 Sep), which refuses exactly that.
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(b1_sam), path(b2_true_sam), path(b2_decoy_sam)

    output:
    tuple val(sampleID), path("${sampleID}.arm_b.tsv"), emit: arm_b

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python3:python samtools')}
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

process computeLocusSpecificityVerdict {
    label 'sample_stage'
    tag "$sampleID"
    publishDir { "${params.sample_dir}/${sampleID}/locus_specificity" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), path(read_evidence_tsv), path(read_gene_tsv), path(arm_a_tsv), path(arm_b_tsv)

    output:
    tuple val(sampleID),
          path("${sampleID}.locus_specificity_summary.tsv"),
          path("${sampleID}.locus_specificity_per_read.tsv"), emit: verdict

    script:
    // Arm A is optional: without it the decision is Arm B alone, and the summary's
    // arms_used column says so. The staged placeholder's name is the switch, because a
    // Nextflow channel element cannot simply be absent for one sample and present for
    // another once the join has happened.
    def arm_a_args = arm_a_tsv.name == 'NO_ARM_A.tsv' ? "" : (
        "--arm-a ${arm_a_tsv} " +
        "--arm-a-target-clade '${params.pks_specificity_arm_a_target_clade}' " +
        "--arm-a-min-confidence ${params.pks_specificity_arm_a_min_confidence} " +
        "--arm-a-min-unique-kmers ${params.pks_specificity_arm_a_min_unique_kmers}")
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python3:python')}
    python3 "${params.scripts}/specificity_decision.py" \
        --sample "${sampleID}" \
        --read-evidence "${read_evidence_tsv}" \
        --read-gene-assignment "${read_gene_tsv}" \
        --arm-b "${arm_b_tsv}" \
        ${arm_a_args} \
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
        --output-summary "${sampleID}.locus_specificity_summary.tsv" \
        --output-per-read "${sampleID}.locus_specificity_per_read.tsv"
    """
}

workflow pksLocusSpecificity {
    take:
    supporting_reads      // tuple(sampleID, pks.fastq.gz) -- extractPksIslandReads' emit in pks_taxa.nf
    read_evidence         // tuple(sampleID, read_evidence.tsv) -- classifyPksReadEvidence's emit
    read_gene_assignment  // tuple(sampleID, read_clb_gene.tsv) -- extractPksIslandReads' emit
    arm_a_results         // tuple(sampleID, arm_a.tsv | NO_ARM_A.tsv) -- an optional KrakenUniq
                          // reclassification against the combined panel, parsed into
                          // read_id/assigned_clade/confidence/unique_kmer_count (see
                          // scripts/specificity_decision.py's docstring for the exact schema).
                          // Not run by this module: deriving a per-read confidence and
                          // unique-k-mer floor from KrakenUniq's raw output is a calibration
                          // question, not something to freeze into pipeline code.

    main:
    pksLocusSpecificityValidateParams()

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
    computeLocusSpecificityVerdict(decision_input)

    emit:
    arm_b   = scoreCompetitiveRealignment.out.arm_b
    verdict = computeLocusSpecificityVerdict.out.verdict
}

// ---------------------------------------------------------------------------------
// Standalone entry point, for re-running the lane on an existing cohort's retained
// outputs without re-running the pipeline:
//   nextflow run Modules/pks_specificity.nf --sample sheet.csv --sample_dir out/ [overrides]
// sheet.csv columns: sample,reads_fastq,read_evidence,read_gene_assignment[,arm_a_tsv]
// arm_a_tsv is optional per row; an empty cell runs that sample Arm-B-only.
// ---------------------------------------------------------------------------------
workflow {
    if (!params.sample) {
        exit 1, "Missing required parameter: --sample <sheet.csv> " +
                 "(columns: sample,reads_fastq,read_evidence,read_gene_assignment[,arm_a_tsv])"
    }
    if (!params.sample_dir) {
        exit 1, "Missing required parameter: --sample_dir <output directory>"
    }

    sheet = Channel.fromPath(params.sample, checkIfExists: true)
        .splitCsv(header: true)

    supporting_reads = sheet.map { row -> tuple(row.sample, file(row.reads_fastq, checkIfExists: true)) }
    read_evidence = sheet.map { row -> tuple(row.sample, file(row.read_evidence, checkIfExists: true)) }
    read_gene_assignment = sheet.map { row -> tuple(row.sample, file(row.read_gene_assignment, checkIfExists: true)) }
    arm_a_results = sheet.map { row ->
        tuple(row.sample, row.arm_a_tsv ? file(row.arm_a_tsv, checkIfExists: true) : pksSpecificityNoArmA())
    }

    pksLocusSpecificity(supporting_reads, read_evidence, read_gene_assignment, arm_a_results)
}
