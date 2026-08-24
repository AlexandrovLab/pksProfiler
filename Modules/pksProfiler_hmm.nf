nextflow.enable.dsl = 2

process pksProfiler_hmm {
    label 'pks_hmm'
    scratch true
    publishDir "${params.pks_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), path(r1), path(r2)

    output:
    tuple val(sampleID),
          path("${sampleID}.nhmmscan.tblout"),
          path("${sampleID}.nhmmscan.log"),
          path("${sampleID}.hmm_counts.tsv"),
          path("${sampleID}.hits.read_ids.txt"),
          path("${sampleID}.hmm.filtered.fa")

    script:
    def tblout      = "${sampleID}.nhmmscan.tblout"
    def logfile     = "${sampleID}.nhmmscan.log"
    def counts_tsv  = "${sampleID}.hmm_counts.tsv"
    def read_ids    = "${sampleID}.hits.read_ids.txt"
    def filtered_fa = "${sampleID}.hmm.filtered.fa"

    """
    set -euo pipefail

	if ! gzip -t "${r1}" >/dev/null 2>&1; then
        echo "ERROR: Corrupt gzip input for ${sampleID}: ${r1}" >&2
        exit 1
    fi

    if ! gzip -t "${r2}" >/dev/null 2>&1; then
        echo "ERROR: Corrupt gzip input for ${sampleID}: ${r2}" >&2
        exit 1
    fi

    # Merge gzipped FASTQs safely (concatenated gzip can be finicky depending on tools)
    zcat "${r1}" "${r2}" | gzip -c > "${sampleID}.merged.fastq.gz"

	ls -lh "${sampleID}.merged.fastq.gz"

    # FASTQ -> FASTA
    seqtk seq -a "${sampleID}.merged.fastq.gz" > "${sampleID}.merged.fa"

	ls -lh "${sampleID}.merged.fa"
    echo "nhmmscan --cpu "${task.cpus}" --tblout "${tblout}" "${params.hmm_model}" "${sampleID}.merged.fa" > "${logfile}" 2>&1"

    echo "Running nhmmscan on sample: ${sampleID}"
    nhmmscan --cpu "${task.cpus}" \\
        --tblout "${tblout}" \\
        "${params.hmm_model}" \\
        "${sampleID}.merged.fa" \\
        > "${logfile}" 2>&1

    # Fail early with a helpful message if tblout didn't get created
    if [[ ! -s "${tblout}" ]]; then
        echo "ERROR: nhmmscan did not produce a non-empty ${tblout}. See ${logfile}." >&2
        exit 1
    fi

    E="${params.hmm_evalue}"

    # Counts: best hit per query under E-value threshold
    awk -v e="\$E" '
      \$0 ~ /^#/ { next }
      NF < 14 { next }
      {
        t=\$1; q=\$3; eval=\$13+0; score=\$14+0;
        if (eval <= e) {
          if (!(q in bestE) || eval < bestE[q] || (eval == bestE[q] && score > bestScore[q])) {
            bestE[q]=eval; bestScore[q]=score; bestT[q]=t;
          }
        }
      }
      END {
        for (q in bestT) cnt[bestT[q]]++;
        for (g in cnt) printf "%s\\t%d\\n", g, cnt[g];
      }
    ' "${tblout}" | sort -k1,1 > "${counts_tsv}.tmp"

    printf "Gene\\tCount\\n" > "${counts_tsv}"
    cat "${counts_tsv}.tmp" >> "${counts_tsv}"
    rm -f "${counts_tsv}.tmp"

    # Read IDs that pass threshold
    awk -v e="\$E" '
      \$0 ~ /^#/ { next }
      NF >= 14 { q=\$3; eval=\$13+0; if (eval <= e) print q }
    ' "${tblout}" | sort -u > "${read_ids}"

    # Filter FASTA to reads with hits
    if [[ -s "${read_ids}" ]]; then
        seqkit grep -f "${read_ids}" "${sampleID}.merged.fa" > "${filtered_fa}"
    else
        : > "${filtered_fa}"
    fi
    """
}
