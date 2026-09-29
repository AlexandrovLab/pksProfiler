nextflow.enable.dsl = 2

// Invoked twice (main.nf aliases it plotPKSTaxa and plotCommunityTaxa). Both used to
// publish figures/taxa_barplots.pdf, so a run with --pks_taxa and --pks_community_taxa
// kept whichever finished last. `scope` names the plot and its published file:
//   pks_island -> figures/pks_island_taxa_barplots.pdf  (Bracken on island-aligned reads)
//   community  -> figures/community_taxa_barplots.pdf   (Bracken on the whole non-host report)
process plotBrackenTaxa {
  label 'sample_stage'
  label 'process_medium'
  scratch true
  publishDir { "${params.sample_dir}/${sampleID}/figures" }, mode: 'copy', saveAs: { "${scope}_taxa_barplots.pdf" }
  conda "${params.krakenuniq_bracken_env}"   // must include python + pandas + matplotlib

  input:
  tuple val(sampleID), path(genus_report), path(species_report)
  val(scope)

  output:
  tuple val(sampleID), path("${sampleID}.${scope}_taxa_barplots.pdf")

  script:
  """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
  set -euo pipefail
  case "${scope}" in pks_island|community) ;; *) echo "unknown taxa plot scope: ${scope}" >&2; exit 1 ;; esac

  python3 "${params.scripts}/plot_bracken_taxa.py" \
    --sample "${sampleID}" \
    --genus "${genus_report}" \
    --species "${species_report}" \
    --scope "${scope}" \
    --out "${sampleID}.${scope}_taxa_barplots.pdf" \
    --top 20
  """
}
