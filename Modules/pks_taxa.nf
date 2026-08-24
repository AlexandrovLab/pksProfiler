process extractPksIslandReads {
  label 'process_medium'
  scratch true
  publishDir "${params.pks_dir}", mode: 'copy'
  conda "${params.pks_align_env}"

  input:
  tuple val(sampleID), path(bam), path(bai)

  output:
  tuple val(sampleID), path("${sampleID}.pks.fastq.gz")

  script:
  """
  set -euo pipefail

  CHR="${params.pks_contig}"
  START=\$(( ${params.pks_shift} - 1 ))
  END=\$(( ${params.pks_shift} + ${params.pks_island_len} ))

  printf "%s\\t%s\\t%s\\n" "\$CHR" "\$START" "\$END" > pks_island.bed

  # Pull alignments overlapping the island
  samtools view -b -L pks_island.bed "${bam}" > "${sampleID}.pks.bam"

  # Convert to a single FASTQ stream (mates/singletons all included)
  samtools fastq "${sampleID}.pks.bam" | gzip -c > "${sampleID}.pks.fastq.gz"

  # Ensure file exists even if no reads
  [[ -s "${sampleID}.pks.fastq.gz" ]] || : > "${sampleID}.pks.fastq.gz"
  """
}

process Bracken {
  scratch true
  label 'process_high_disk'
  publishDir "${params.pks_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  tuple val(sampleID), path(fastq_gz)

  output:
  tuple val(sampleID),
    path("${sampleID}.krakenuniq.report.txt"),
    path("${sampleID}.classified.fasta"),
    path("${sampleID}.unclassified.fasta"),
    path("${sampleID}.bracken.G.report.txt"),
    path("${sampleID}.bracken.S.report.txt"),
    path("${sampleID}.bracken.G.krakenreport.txt"),
    path("${sampleID}.bracken.S.krakenreport.txt"),
    path("${sampleID}.bracken.G.mpa.krakenreport.txt"),
    path("${sampleID}.bracken.S.mpa.krakenreport.txt")

  script:
  """
  set -euo pipefail

  REPORT="${sampleID}.krakenuniq.report.txt"
  OUTPUT="${sampleID}.krakenuniq.output.txt"
  CLASSIFIED="${sampleID}.classified.fasta"
  UNCLASSIFIED="${sampleID}.unclassified.fasta"

  # Decompress to a plain FASTQ for your krakenuniq usage
  zcat "${fastq_gz}" > "${sampleID}.pks.fastq"

  # If no reads, write empty outputs and exit successfully
  if [[ ! -s "${sampleID}.pks.fastq" ]]; then
    echo "No PKS reads for ${sampleID}; writing empty outputs."
    : > "\$REPORT"
    : > "\$OUTPUT"
    : > "\$CLASSIFIED"
    : > "\$UNCLASSIFIED"
    for lvl in G S; do
      : > "${sampleID}.bracken.\${lvl}.report.txt"
      : > "${sampleID}.bracken.\${lvl}.krakenreport.txt"
      : > "${sampleID}.bracken.\${lvl}.mpa.krakenreport.txt"
    done
    exit 0
  fi

  krakenuniq --db "${params.kraken_db}" --threads "${task.cpus}" \
    --report-file "\$REPORT" --output "\$OUTPUT" \
    --classified-out "\$CLASSIFIED" --unclassified-out "\$UNCLASSIFIED" \
    "${sampleID}.pks.fastq"

  cat "\$REPORT"

  # Count reads per taxonomic level from Kraken report
  GENUS_READS=\$(awk '\$4=="G" && \$2>0 {sum+=\$2} END {print sum+0}' "\$REPORT")
  SPECIES_READS=\$(awk '\$4=="S" && \$2>0 {sum+=\$2} END {print sum+0}' "\$REPORT")

  for lvl in G S; do
    bracken_output="${sampleID}.bracken.\${lvl}.report.txt"
    bracken_kraken_report="${sampleID}.bracken.\${lvl}.krakenreport.txt"
    bracken_kraken_mpa_report="${sampleID}.bracken.\${lvl}.mpa.krakenreport.txt"

    # Select reads for this level
    if [[ "\${lvl}" == "G" ]]; then
      LVL_READS="\${GENUS_READS}"
    else
      LVL_READS="\${SPECIES_READS}"
    fi

    # If no reads for this level, write empty outputs and skip
    if [[ "\${LVL_READS}" -lt 2 ]]; then
      echo "Skipping Bracken level \${lvl}: reads=\${LVL_READS}"
      : > "\${bracken_output}"
      : > "\${bracken_kraken_report}"
      : > "\${bracken_kraken_mpa_report}"
      continue
    fi

	bracken \
    -d "${params.kraken_db}" \
    -i "\$REPORT" \
    -o "\${bracken_output}" \
    -w "\${bracken_kraken_report}" \
    -r ${params.bracken_read_length} \
    -l "\${lvl}" \
    -t 2

    kreport2mpa.py -r "\${bracken_kraken_report}" \
    -o "\${bracken_kraken_mpa_report}" --display-header
  done


  """
}

process process_bracken {

  scratch true
  publishDir "${params.pks_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  path bracken_files

  output:
  tuple path("bracken.genus.mpa.report.txt"),
        path("bracken.species.mpa.report.txt")

  script:
  def genus_files   = bracken_files.findAll { it.name.endsWith('.G.mpa.krakenreport.txt') }
  def species_files = bracken_files.findAll { it.name.endsWith('.S.mpa.krakenreport.txt') }

  def genus_str   = genus_files.collect { "\"${it}\"" }.join(' ')
  def species_str = species_files.collect { "\"${it}\"" }.join(' ')

  """
  set -euo pipefail

  if [ -n "${genus_str}" ]; then
    combine_mpa.py --input ${genus_str} --output bracken.genus.mpa.report.txt
  else
    echo "No genus files found." > bracken.genus.mpa.report.txt
  fi

  if [ -n "${species_str}" ]; then
    combine_mpa.py --input ${species_str} --output bracken.species.mpa.report.txt
  else
    echo "No species files found." > bracken.species.mpa.report.txt
  fi
  """
}




