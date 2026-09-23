process extractPksIslandReads {
  label 'process_medium'
  scratch true
  publishDir { "${params.sample_dir}/${sampleID}/taxonomy" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
  conda "${params.pks_align_env}"

  input:
  tuple val(sampleID), path(bam), path(bai)

  output:
  tuple val(sampleID),
        path("${sampleID}.pks.fastq.gz"),
        path("${sampleID}.read_clb_gene.tsv")

  script:
  """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_annotation ${params.dep_digest?.pks_annotation}  scripts ${params.dep_digest?.scripts}
  set -euo pipefail

  # Ludmil, revised report finding 6: this was already the correct 0-based BED
  # conversion of the 1-based annotation (pks_shift-1 .. pks_shift+len) -- the
  # only one of the affected call sites that got it right on its own. Migrated
  # to the canonical params directly so there is one source of truth, not
  # independently re-derived +/-1 arithmetic in every module.
  CHR="${params.pks_contig}"
  START=\$(( ${params.pks_start_1based} - 1 ))
  END=${params.pks_end_1based}

  printf "%s\\t%s\\t%s\\n" \
    "\$CHR" "\$START" "\$END" \
    > pks_island.bed

  # Extract alignments overlapping the complete pks island
  samtools view \
    -b \
    -L pks_island.bed \
    "${bam}" \
    > "${sampleID}.pks.bam"

  # F14: one definition of the clb intervals and of read-to-gene assignment, shared
  # with the alignment lane. The awk this replaces took every feature with a Name,
  # not only clbA-clbS, so a read over a neighbouring gene was assigned to it here and
  # ignored there.
  python3 "${params.scripts}/clb_gene_bed.py" \
    --annotation "${params.pks_genome_annotation}" \
    --output clb_genes.bed

  bedtools bamtobed -i "${sampleID}.pks.bam" |
    bedtools intersect -a - -b clb_genes.bed -wo |
    python3 "${params.scripts}/assign_reads_to_genes.py" \
      --output "${sampleID}.read_clb_gene.tsv" \
      --min-overlap "${params.min_gene_overlap_bp}" \
      > /dev/null

  # Convert overlapping alignments to a single FASTQ stream
  samtools fastq "${sampleID}.pks.bam" |
    gzip -c > "${sampleID}.pks.fastq.gz"
  """
}


process Bracken {
  scratch true
  label 'process_high_disk'
  publishDir { "${params.sample_dir}/${sampleID}/taxonomy" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
  conda "${params.krakenuniq_bracken_env}"

  input:
  tuple val(sampleID), path(fastq_gz), path(read_gene_tsv)

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
        path("${sampleID}.bracken.S.mpa.krakenreport.txt"),
        path("${sampleID}.clb_species_support.tsv"), emit: reports
  tuple val(sampleID), path("${sampleID}.taxonomy.qc.tsv"), emit: qc

  script:
  """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # kraken_db ${params.dep_digest?.kraken_db}  scripts ${params.dep_digest?.scripts}
  set -euo pipefail

  REPORT="${sampleID}.krakenuniq.report.txt"
  OUTPUT="${sampleID}.krakenuniq.output.txt"
  CLASSIFIED="${sampleID}.classified.fasta"
  UNCLASSIFIED="${sampleID}.unclassified.fasta"
  SPECIES_MATRIX="${sampleID}.clb_species_support.tsv"
  TAXONOMY_QC="${sampleID}.taxonomy.qc.tsv"

  # Decompress PKS reads for KrakenUniq
  zcat "${fastq_gz}" > "${sampleID}.pks.fastq"

  # A valid sample can contain no PKS reads
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

    # Write a valid header-only species-by-clb matrix
    {
      printf "Species\\tTaxID"

      for gene in {A..S}; do
        printf "\\tclb%s" "\$gene"
      done

      printf "\\tTotal\\n"
    } > "\$SPECIES_MATRIX"

    # F10: a sample with no island reads and a sample that was never classified used
    # to be the same absence. pks.qc.summary.tsv has a row for every sample, so the
    # distinction belongs there.
    printf "Sample\\tMetric\\tValue\\n" > "\$TAXONOMY_QC"
    printf "%s\\ttaxonomy_status\\tno_pks_reads\\n" "${sampleID}" >> "\$TAXONOMY_QC"
    printf "%s\\tclb_species_reported\\t0\\n" "${sampleID}" >> "\$TAXONOMY_QC"

    exit 0
  fi

  krakenuniq \
    --db "${params.kraken_db}" \
    --threads "${task.cpus}" \
    --report-file "\$REPORT" \
    --output "\$OUTPUT" \
    --classified-out "\$CLASSIFIED" \
    --unclassified-out "\$UNCLASSIFIED" \
    "${sampleID}.pks.fastq"

  # Direct per-read KrakenUniq assignments retain the connection
  # between taxon and clb gene. Bracken outputs remain separate because
  # Bracken estimates aggregate abundance and has no per-read identity.
  python "${params.scripts}/build_clb_species_matrix.py" \
    --read-gene "${read_gene_tsv}" \
    --kraken-output "\$OUTPUT" \
    --kraken-db "${params.kraken_db}" \
    --output "\$SPECIES_MATRIX"

  # Apply the same conservative support threshold used by Bracken.
  # Reads classified only above the requested rank do not satisfy this
  # requirement until at least two reads support a target-rank node.
  # F10: the largest read count on a single taxon, not the sum across all of them.
  # Bracken's -t 2 admits a taxon with two reads; summing the rank let two species
  # with one read each through, and Bracken then reported neither. The pre-check now
  # asks the question Bracken will ask.
  GENUS_READS=\$(awk -F '\\t' '
    \$8 == "genus" && \$2 ~ /^[0-9]+\$/ && \$2+0 > best {
      best = \$2+0
    }

    END {
      print best+0
    }
  ' "\$REPORT")

  SPECIES_READS=\$(awk -F '\\t' '
    \$8 == "species" && \$2 ~ /^[0-9]+\$/ && \$2+0 > best {
      best = \$2+0
    }

    END {
      print best+0
    }
  ' "\$REPORT")

  for lvl in G S; do
    bracken_output="${sampleID}.bracken.\${lvl}.report.txt"
    bracken_kraken_report="${sampleID}.bracken.\${lvl}.krakenreport.txt"
    bracken_kraken_mpa_report="${sampleID}.bracken.\${lvl}.mpa.krakenreport.txt"

    if [[ "\$lvl" == "G" ]]; then
      LVL_READS="\$GENUS_READS"
    else
      LVL_READS="\$SPECIES_READS"
    fi

    if [[ "\$LVL_READS" -lt 2 ]]; then
      echo "Skipping Bracken level \$lvl: no taxon reaches 2 reads (best=\$LVL_READS)"

      : > "\$bracken_output"
      : > "\$bracken_kraken_report"
      : > "\$bracken_kraken_mpa_report"

      continue
    fi

    bracken \
      -d "${params.kraken_db}" \
      -i "\$REPORT" \
      -o "\$bracken_output" \
      -w "\$bracken_kraken_report" \
      -r ${params.bracken_read_length} \
      -l "\$lvl" \
      -t 2

    kreport2mpa.py \
      -r "\$bracken_kraken_report" \
      -o "\$bracken_kraken_mpa_report" \
      --display-header
  done

  # F10: how this sample ended up, for the one table that has every sample in it.
  SPECIES_ROWS=\$(awk 'END { print (NR > 1 ? NR - 1 : 0) }' "\$SPECIES_MATRIX")
  if [[ "\$SPECIES_ROWS" -gt 0 ]]; then
    STATUS="species_identified"
  elif [[ "\$SPECIES_READS" -lt 2 ]]; then
    STATUS="below_rank_threshold"
  else
    STATUS="no_species_identified"
  fi

  printf "Sample\\tMetric\\tValue\\n" > "\$TAXONOMY_QC"
  printf "%s\\ttaxonomy_status\\t%s\\n" "${sampleID}" "\$STATUS" >> "\$TAXONOMY_QC"
  printf "%s\\tclb_species_reported\\t%s\\n" "${sampleID}" "\$SPECIES_ROWS" >> "\$TAXONOMY_QC"
  """
}


process process_bracken {
  scratch true
  publishDir "${params.pks_taxonomy_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  path bracken_files

  output:
  tuple path("bracken.genus.mpa.report.txt"),
        path("bracken.species.mpa.report.txt")

  script:
  // F02. combine_mpa.py is KrakenTools' and takes its inputs as arguments; it has no
  // manifest option and we do not fork it. What we can stop doing is building that
  // argument list in Groovy, which put one filename per sample into .command.sh with
  // no idea of the limit. The lists are built here with find, and the size is checked
  // against ARG_MAX before the call, so an oversized cohort fails with a message that
  // says what happened instead of "Argument list too long" from the shell.
  """
  set -euo pipefail

  find . -maxdepth 1 -name '*.G.mpa.krakenreport.txt' -printf '%f\\n' | sort > genus.list
  find . -maxdepth 1 -name '*.S.mpa.krakenreport.txt' -printf '%f\\n' | sort > species.list

  arg_limit=\$(( \$(getconf ARG_MAX) / 2 ))

  combine_rank() {
      local list=\$1 out=\$2 label=\$3
      if [[ ! -s "\$list" ]]; then
          echo "No \$label files found." > "\$out"
          return 0
      fi
      local bytes
      bytes=\$(wc -c < "\$list")
      if (( bytes >= arg_limit )); then
          echo "ERROR: \$(wc -l < "\$list") \$label files (\$bytes bytes of paths) exceed" >&2
          echo "       half of ARG_MAX (\$arg_limit). combine_mpa.py takes files as arguments" >&2
          echo "       and cannot read a manifest; this cohort needs it run in batches." >&2
          return 1
      fi
      # shellcheck disable=SC2046
      combine_mpa.py --output "\$out" --input \$(cat "\$list")
  }

  combine_rank genus.list   bracken.genus.mpa.report.txt   genus
  combine_rank species.list bracken.species.mpa.report.txt species
  """
}


process combineClbTaxonomySupport {
  scratch true
  publishDir "${params.pks_taxonomy_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  path species_support_files
  path combine_script
  path species_file_list

  output:
  path("pks.clb_species_support.tsv")

  script:
  // F02: one list file rather than one argument per sample.
  """
  set -euo pipefail

  python "${combine_script}" \
    --species-files-from "${species_file_list}" \
    --species-output pks.clb_species_support.tsv
  """
}
