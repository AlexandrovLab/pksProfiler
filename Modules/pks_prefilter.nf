nextflow.enable.dsl = 2

process krakenPrefilter {
    label 'sample_stage'
    tag "$sampleID"
    label 'prefilter_kraken'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/prefilter" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.prefilter_env}"
    input:
    tuple val(sampleID), path(reads)
    output:
    tuple val(sampleID), path("${sampleID}.primary.fastq.gz"), emit: primary
    tuple val(sampleID), path("${sampleID}.non_target.fastq.gz"), emit: non_target
    tuple val(sampleID), path("${sampleID}.krakenuniq.report.txt"), path("${sampleID}.krakenuniq.output.txt"), emit: taxonomy
    tuple val(sampleID), path("${sampleID}.prefilter.qc.tsv"), emit: qc
    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # kraken_db ${params.dep_digest?.kraken_db}  scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    NODES="${params.kraken_db}/taxonomy/nodes.dmp"
    if [[ ! -s "\$NODES" ]]; then
      NODES="${params.kraken_db}/taxDB"
    fi
    test -s "\$NODES" || { echo "ERROR: Missing taxonomy parents (taxonomy/nodes.dmp or taxDB) in ${params.kraken_db}" >&2; exit 1; }
    gzip -dc "${reads}" > "${sampleID}.non_host.fastq"
    krakenuniq --db "${params.kraken_db}" --threads "${task.cpus}" \
      --report-file "${sampleID}.krakenuniq.report.txt" \
      --output "${sampleID}.krakenuniq.output.txt" "${sampleID}.non_host.fastq"
    python "${params.scripts}/route_kraken_reads.py" \
      --reads "${sampleID}.non_host.fastq" --kraken-output "${sampleID}.krakenuniq.output.txt" \
      --nodes "\$NODES" --target-taxid "${params.primary_taxid}" \
      --primary "${sampleID}.primary.fastq.gz" --non-target "${sampleID}.non_target.fastq.gz" \
      --qc "${sampleID}.prefilter.qc.tsv" --sample "${sampleID}"
    """
}

process buildClbDiamondDb {
    label 'process_low'
    scratch true
    conda "${params.prefilter_env}"
    input:
    path clb_proteins
    output:
    path "clb_reference.dmnd"
    script:
    """
    set -euo pipefail
    diamond makedb --in "${clb_proteins}" --db clb_reference
    """
}

process diamondRescue {
    label 'sample_stage'
    tag "$sampleID"
    label 'process_low'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/prefilter" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.prefilter_env}"
    input:
    tuple val(sampleID), path(reads)
    path clb_db
    output:
    tuple val(sampleID), path("${sampleID}.rescued.fastq.gz"), emit: reads
    tuple val(sampleID), path("${sampleID}.diamond.tsv"), emit: matches
    tuple val(sampleID), path("${sampleID}.diamond.qc.tsv"), emit: qc
    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    if [[ \$(gzip -dc "${reads}" | wc -l) -eq 0 ]]; then
      : > "${sampleID}.diamond.tsv"
    else
      gzip -dc "${reads}" > "${sampleID}.non_target.fastq"
      diamond blastx --db "${clb_db}" --query "${sampleID}.non_target.fastq" \
        --threads "${task.cpus}" --evalue "${params.diamond_evalue}" \
        --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
        --out "${sampleID}.diamond.tsv"
    fi
    python "${params.scripts}/select_diamond_reads.py" \
      --reads "${reads}" --diamond "${sampleID}.diamond.tsv" \
      --output "${sampleID}.rescued.fastq.gz" --qc "${sampleID}.diamond.qc.tsv" \
      --sample "${sampleID}" --max-evalue "${params.diamond_evalue}" \
      --min-aa "${params.diamond_min_aa}" \
      --min-query-coverage "${params.diamond_min_query_coverage}"
    """
}

// Not yet invoked from the workflow. Joining diamondRescue.out.matches with
// krakenPrefilter.out.taxonomy by sampleID, and passing the kraken database path
// through, is wiring in main.nf's workflow {} block -- out of scope here. Once that
// join is added, this process turns diamond.tsv's discarded sseqid and
// krakenPrefilter's per-read taxid into a table build_cohort_report.py already knows
// how to read (it looks for this exact published path and renders nothing if it is
// absent).
process diamondRescueTaxonomy {
    label 'sample_stage'
    tag "$sampleID"
    label 'process_low'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/prefilter" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.prefilter_env}"

    input:
    tuple val(sampleID), path(diamond_matches), path(kraken_output)

    output:
    path "${sampleID}.diamond_rescue_taxonomy.tsv"

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # kraken_db ${params.dep_digest?.kraken_db}  scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    python "${params.scripts}/summarize_diamond_rescue_taxonomy.py" \
      --sample "${sampleID}" --diamond "${diamond_matches}" \
      --kraken-output "${kraken_output}" --kraken-db "${params.kraken_db}" \
      --output "${sampleID}.diamond_rescue_taxonomy.tsv"
    """
}

process mergePksCandidates {
    label 'sample_stage'
    tag "$sampleID"
    label 'process_low'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/prefilter" }, mode: 'copy', enabled: params.save_intermediates.toString().toBoolean(), saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.fastp_env}"
    input:
    tuple val(sampleID), path(primary), path(rescued)
    output:
    tuple val(sampleID), path("${sampleID}.pks_candidates.fastq.gz"), emit: reads
    tuple val(sampleID), path("${sampleID}.candidate_merge.qc.tsv"), emit: qc
    script:
    """
    set -euo pipefail
    { gzip -dc "${primary}"; gzip -dc "${rescued}"; } | gzip -c > "${sampleID}.pks_candidates.fastq.gz"
    PRIMARY=\$(gzip -dc "${primary}" | awk 'END { print int(NR / 4) }')
    RESCUED=\$(gzip -dc "${rescued}" | awk 'END { print int(NR / 4) }')
    TOTAL=\$(gzip -dc "${sampleID}.pks_candidates.fastq.gz" | awk 'END { print int(NR / 4) }')
    [[ \$((PRIMARY + RESCUED)) -eq "\$TOTAL" ]] || { echo "ERROR: Candidate merge count mismatch" >&2; exit 1; }
    printf "Sample\tMetric\tValue\n" > "${sampleID}.candidate_merge.qc.tsv"
    printf "%s\tcandidate_primary_reads\t%s\n" "${sampleID}" "\$PRIMARY" >> "${sampleID}.candidate_merge.qc.tsv"
    printf "%s\tcandidate_rescued_reads\t%s\n" "${sampleID}" "\$RESCUED" >> "${sampleID}.candidate_merge.qc.tsv"
    printf "%s\tcandidate_total_reads\t%s\n" "${sampleID}" "\$TOTAL" >> "${sampleID}.candidate_merge.qc.tsv"
    """
}

process sampleBracken {
    label 'sample_stage'
    tag "$sampleID"
    label 'process_high_disk'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/taxonomy" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.krakenuniq_bracken_env}"
    input:
    tuple val(sampleID), path(kraken_report), path(kraken_output)
    output:
    tuple val(sampleID), path("${sampleID}.bracken.G.report.txt"), path("${sampleID}.bracken.S.report.txt"), emit: reports
    tuple val(sampleID), path("${sampleID}.bracken.G.krakenreport.txt"), path("${sampleID}.bracken.S.krakenreport.txt"), emit: kraken_reports
    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # kraken_db ${params.dep_digest?.kraken_db}
    set -euo pipefail

    # Ludmil, revised report finding 4: this called bracken -t 2 unconditionally.
    # Bracken itself fails when no taxon at the requested rank reaches its read
    # threshold, so a thin community-taxonomy sample turned into a task failure
    # instead of an explicit below-threshold/no-call state. pks_taxa.nf's Bracken
    # process already guards exactly this way -- same pre-check, same threshold --
    # mirrored here rather than invented fresh.
    GENUS_READS=\$(awk -F '\\t' '
      \$8 == "genus" && \$2 ~ /^[0-9]+\$/ && \$2+0 > best {
        best = \$2+0
      }
      END {
        print best+0
      }
    ' "${kraken_report}")

    SPECIES_READS=\$(awk -F '\\t' '
      \$8 == "species" && \$2 ~ /^[0-9]+\$/ && \$2+0 > best {
        best = \$2+0
      }
      END {
        print best+0
      }
    ' "${kraken_report}")

    for LEVEL in G S; do
      bracken_output="${sampleID}.bracken.\${LEVEL}.report.txt"
      bracken_kraken_report="${sampleID}.bracken.\${LEVEL}.krakenreport.txt"

      if [[ "\$LEVEL" == "G" ]]; then
        LVL_READS="\$GENUS_READS"
      else
        LVL_READS="\$SPECIES_READS"
      fi

      if [[ "\$LVL_READS" -lt 2 ]]; then
        echo "Skipping Bracken level \$LEVEL for ${sampleID}: no taxon reaches 2 reads (best=\$LVL_READS)"
        : > "\$bracken_output"
        : > "\$bracken_kraken_report"
        continue
      fi

      bracken -d "${params.kraken_db}" -i "${kraken_report}" \
        -o "\$bracken_output" \
        -w "\$bracken_kraken_report" \
        -r "${params.bracken_read_length}" -l "\$LEVEL" -t 2
    done
    """
}
