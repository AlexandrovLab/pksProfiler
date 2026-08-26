nextflow.enable.dsl = 2

process plotBrackenTaxa {
  label 'process_medium'
  scratch true
  publishDir "${params.pks_taxonomy_plots_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"   // must include python + pandas + matplotlib

  input:
  tuple val(sampleID), path(genus_report), path(species_report)

  output:
  tuple val(sampleID), path("${sampleID}.pks_island_taxa_barplots.pdf")

  script:
  """
  set -euo pipefail

  python3 "${params.scripts}/plot_bracken_taxa.py" \
    --sample "${sampleID}" \
    --genus "${genus_report}" \
    --species "${species_report}" \
    --out "${sampleID}.pks_island_taxa_barplots.pdf" \
    --top 20
  """
}
