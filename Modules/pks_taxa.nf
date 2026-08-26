process extractPksIslandReads {
  label 'process_medium'
  scratch true
  publishDir "${params.pks_dir}", mode: 'copy'
  conda "${params.pks_align_env}"

  input:
  tuple val(sampleID), path(bam), path(bai)

  output:
  tuple val(sampleID),
        path("${sampleID}.pks.fastq.gz"),
        path("${sampleID}.read_clb_gene.tsv"),
        path("${sampleID}.clb_gene_support.tsv")

  script:
  """
  set -euo pipefail

  CHR="${params.pks_contig}"
  START=\$(( ${params.pks_shift} - 1 ))
  END=\$(( ${params.pks_shift} + ${params.pks_island_len} ))

  printf "%s\\t%s\\t%s\\n" \
    "\$CHR" "\$START" "\$END" \
    > pks_island.bed

  # Extract alignments overlapping the complete pks island
  samtools view \
    -b \
    -L pks_island.bed \
    "${bam}" \
    > "${sampleID}.pks.bam"

  # Convert clb gene coordinates from GFF to BED
  awk -F '\\t' '
    BEGIN {
      OFS="\\t"
    }

    \$0 !~ /^#/ && \$3 == "gene" {
      gene=""
      n=split(\$9, attributes, ";")

      for (i=1; i<=n; i++) {
        if (attributes[i] ~ /^Name=/) {
          sub(/^Name=/, "", attributes[i])
          gene=attributes[i]
        }
      }

      if (gene != "") {
        print \$1, \$4-1, \$5, gene
      }
    }
  ' "${params.pks_genome_annotation}" > clb_genes.bed

  # Assign each read/template to the clb gene with the greatest
  # total aligned-base overlap
  bedtools bamtobed \
    -i "${sampleID}.pks.bam" |
    bedtools intersect \
      -a - \
      -b clb_genes.bed \
      -wo |
    awk '
      BEGIN {
        OFS="\\t"
      }

      {
        read_id=\$4
        gene=\$10
        overlap=\$11+0
        total[read_id SUBSEP gene] += overlap
      }

      END {
        for (key in total) {
          split(key, fields, SUBSEP)

          read_id=fields[1]
          gene=fields[2]
          overlap=total[key]

          if (!(read_id in best_overlap) ||
              overlap > best_overlap[read_id] ||
              (overlap == best_overlap[read_id] &&
               gene < best_gene[read_id])) {
            best_overlap[read_id]=overlap
            best_gene[read_id]=gene
          }
        }

        for (read_id in best_gene) {
          print read_id, best_gene[read_id], best_overlap[read_id]
        }
      }
    ' > "${sampleID}.read_clb_gene.tmp.tsv"

  printf "read_id\\tGene\\toverlap_bp\\n" \
    > "${sampleID}.read_clb_gene.tsv"

  sort -k1,1 "${sampleID}.read_clb_gene.tmp.tsv" \
    >> "${sampleID}.read_clb_gene.tsv"

  rm -f "${sampleID}.read_clb_gene.tmp.tsv"

  # Overall clb support irrespective of taxonomy. These values are
  # derived from the same read-to-gene assignments used to construct
  # the species-by-gene support matrix, so the column totals reconcile.
  awk -F '\t' -v sample="${sampleID}" '
    BEGIN {
      OFS="\t"
    }

    NR > 1 &&
    length(\$2) == 4 &&
    substr(\$2, 1, 3) == "clb" &&
    index("ABCDEFGHIJKLMNOPQRS", substr(\$2, 4, 1)) > 0 {
      count[\$2]++
      total++
    }

    END {
      printf "Sample"

      for (i=65; i<=83; i++) {
        printf "%sclb%c", OFS, i
      }

      printf "%sTotal\\n%s", OFS, sample

      for (i=65; i<=83; i++) {
        gene=sprintf("clb%c", i)
        printf "%s%d", OFS, count[gene]+0
      }

      printf "%s%d\\n", OFS, total+0
    }
  ' "${sampleID}.read_clb_gene.tsv" \
    > "${sampleID}.clb_gene_support.tsv"

  # Convert overlapping alignments to a single FASTQ stream
  samtools fastq "${sampleID}.pks.bam" |
    gzip -c > "${sampleID}.pks.fastq.gz"
  """
}


process Bracken {
  scratch true
  label 'process_high_disk'
  publishDir "${params.pks_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  tuple val(sampleID), path(fastq_gz), path(read_gene_tsv), path(gene_support_tsv)

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
        path("${sampleID}.clb_species_support.tsv")

  script:
  """
  set -euo pipefail

  REPORT="${sampleID}.krakenuniq.report.txt"
  OUTPUT="${sampleID}.krakenuniq.output.txt"
  CLASSIFIED="${sampleID}.classified.fasta"
  UNCLASSIFIED="${sampleID}.unclassified.fasta"
  SPECIES_MATRIX="${sampleID}.clb_species_support.tsv"

  test -s "${gene_support_tsv}"

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

  # Count reads reported at genus and species levels
  GENUS_READS=\$(awk -F '\\t' '
    \$8 == "genus" && \$2 ~ /^[0-9]+\$/ {
      sum += \$2
    }

    END {
      print sum+0
    }
  ' "\$REPORT")

  SPECIES_READS=\$(awk -F '\\t' '
    \$8 == "species" && \$2 ~ /^[0-9]+\$/ {
      sum += \$2
    }

    END {
      print sum+0
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

    # Bracken requires at least two reads at the requested level
    if [[ "\$LVL_READS" -lt 2 ]]; then
      echo "Skipping Bracken level \$lvl: reads=\$LVL_READS"

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
  """
}


process process_bracken {
  scratch true
  publishDir "${params.pks_summary_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  path bracken_files

  output:
  tuple path("bracken.genus.mpa.report.txt"),
        path("bracken.species.mpa.report.txt")

  script:
  def genus_files = bracken_files.findAll { file ->
    file.name.endsWith('.G.mpa.krakenreport.txt')
  }

  def species_files = bracken_files.findAll { file ->
    file.name.endsWith('.S.mpa.krakenreport.txt')
  }

  def genus_str = genus_files
    .collect { file -> "\"${file}\"" }
    .join(' ')

  def species_str = species_files
    .collect { file -> "\"${file}\"" }
    .join(' ')

  """
  set -euo pipefail

  if [[ -n "${genus_str}" ]]; then
    combine_mpa.py \
      --input ${genus_str} \
      --output bracken.genus.mpa.report.txt
  else
    echo "No genus files found." \
      > bracken.genus.mpa.report.txt
  fi

  if [[ -n "${species_str}" ]]; then
    combine_mpa.py \
      --input ${species_str} \
      --output bracken.species.mpa.report.txt
  else
    echo "No species files found." \
      > bracken.species.mpa.report.txt
  fi
  """
}


process combineClbTaxonomySupport {
  scratch true
  publishDir "${params.pks_summary_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  path gene_support_files
  path species_support_files
  path combine_script

  output:
  tuple path("pks.clb_gene_support.tsv"),
        path("pks.clb_species_support.tsv")

  script:
  def gene_inputs = gene_support_files
    .collect { file -> "\"${file}\"" }
    .join(' ')

  def species_inputs = species_support_files
    .collect { file -> "\"${file}\"" }
    .join(' ')

  """
  set -euo pipefail

  python "${combine_script}" \
    --gene-files ${gene_inputs} \
    --species-files ${species_inputs} \
    --gene-output pks.clb_gene_support.tsv \
    --species-output pks.clb_species_support.tsv
  """
}
