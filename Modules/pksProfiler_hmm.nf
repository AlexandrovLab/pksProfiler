nextflow.enable.dsl = 2

process pksProfiler_hmm {
    label 'pks_hmm'
    scratch true
    publishDir "${params.pks_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID),
          path("${sampleID}.nhmmscan.tblout"),
          path("${sampleID}.nhmmscan.log"),
          path("${sampleID}.hmm_counts.tsv"),
          path("${sampleID}.hits.read_ids.txt"),
          path("${sampleID}.hmm.filtered.fa"),
          emit: profile
    tuple val(sampleID), path("${sampleID}.hmm.qc.tsv"), emit: qc

    script:
    def tblout      = "${sampleID}.nhmmscan.tblout"
    def logfile     = "${sampleID}.nhmmscan.log"
    def counts_tsv  = "${sampleID}.hmm_counts.tsv"
    def read_ids    = "${sampleID}.hits.read_ids.txt"
    def filtered_fa = "${sampleID}.hmm.filtered.fa"
    def qc           = "${sampleID}.hmm.qc.tsv"

    """
    set -euo pipefail

    if ! gzip -t "${reads}" >/dev/null 2>&1; then
        echo "ERROR: Corrupt gzip input for ${sampleID}: ${reads}" >&2
        exit 1
    fi

    # Convert FASTQ to FASTA
    seqtk seq -a "${reads}" > "${sampleID}.merged.fa"

    # A valid sample may contain no reads after upstream filtering
    if [[ ! -s "${sampleID}.merged.fa" ]]; then
        echo "No reads available for ${sampleID}; recording zero HMM counts." > "${logfile}"

        : > "${tblout}"
        : > "${read_ids}"
        : > "${filtered_fa}"

        printf "Gene\\tCount\\n" > "${counts_tsv}"
        for gene in {A..S}; do
            printf "clb%s\\t0\\n" "\$gene" >> "${counts_tsv}"
        done

        printf "Sample\\tMetric\\tValue\\n" > "${qc}"
        printf "%s\\treads_clb_genes_hmm\\t0\\n" "${sampleID}" >> "${qc}"
        printf "%s\\tnum_clb_genes_hmm\\t0\\n" "${sampleID}" >> "${qc}"

        exit 0
    fi

    echo "Running nhmmscan on sample: ${sampleID}"

    if [[ "${params.hmm_chunking}" == "true" && "${task.cpus}" -gt 1 ]]; then
        mkdir -p hmm_chunks hmm_chunk_results
        : > hmm_chunk_results/pids.txt

        # Split into at most one chunk per allocated CPU
        seqkit split2 \\
            --by-part "${task.cpus}" \\
            --out-dir hmm_chunks \\
            "${sampleID}.merged.fa"

        # Run one single-threaded nhmmscan process per chunk
        for chunk in hmm_chunks/*; do
            chunk_name=\$(basename "\$chunk")

            nhmmscan --cpu 1 \\
                --tblout "hmm_chunk_results/\${chunk_name}.tblout" \\
                "${params.hmm_model}" \\
                "\$chunk" \\
                > "hmm_chunk_results/\${chunk_name}.log" 2>&1 &

            echo \$! >> hmm_chunk_results/pids.txt
        done

        # Wait for all chunks and record whether any failed
        chunk_status=0
        while read -r pid; do
            if ! wait "\$pid"; then
                chunk_status=1
            fi
        done < hmm_chunk_results/pids.txt

        if [[ "\$chunk_status" -ne 0 ]]; then
            echo "ERROR: At least one nhmmscan chunk failed for ${sampleID}." >&2
            exit 1
        fi

        # Merge tblout files, retaining the header from the first file only
        first_chunk=true
        : > "${tblout}"
        : > "${logfile}"

        for part in hmm_chunk_results/*.tblout; do
            if [[ "\$first_chunk" == "true" ]]; then
                cp "\$part" "${tblout}"
                first_chunk=false
            else
                { grep -v '^#' "\$part" || true; } >> "${tblout}"
            fi
        done

        # Combine the individual chunk logs
        for part_log in hmm_chunk_results/*.log; do
            printf '\\n===== %s =====\\n' "\$part_log" >> "${logfile}"
            cat "\$part_log" >> "${logfile}"
        done
    else
        # Original non-chunked execution
        nhmmscan --cpu "${task.cpus}" \\
            --tblout "${tblout}" \\
            "${params.hmm_model}" \\
            "${sampleID}.merged.fa" \\
            > "${logfile}" 2>&1
    fi

    # A failed search must not be mistaken for a valid zero-hit sample
    if [[ ! -s "${tblout}" ]]; then
        echo "ERROR: nhmmscan did not produce a non-empty ${tblout}. See ${logfile}." >&2
        exit 1
    fi

    E="${params.hmm_evalue}"

    # Select the best qualifying HMM hit for each query and count by gene
    awk -v e="\$E" '
      \$0 ~ /^#/ { next }
      NF < 14 { next }
      {
        t=\$1;
        sub(/[.]cds[.]aln/, "", t);
        q=\$3;
        eval=\$13+0;
        score=\$14+0;

        if (eval <= e) {
          if (!(q in bestE) || eval < bestE[q] ||
              (eval == bestE[q] && score > bestScore[q])) {
            bestE[q]=eval;
            bestScore[q]=score;
            bestT[q]=t;
          }
        }
      }
      END {
        # Initialize clbA through clbS to zero
        for (i=0; i<19; i++) {
          gene=sprintf("clb%c", 65+i);
          cnt[gene]=0;
        }

        # Count the best qualifying hit for each query
        for (q in bestT) {
          cnt[bestT[q]]++;
        }

        for (g in cnt) {
          printf "%s\\t%d\\n", g, cnt[g];
        }
      }
    ' "${tblout}" | sort -k1,1 > "${counts_tsv}.tmp"

    printf "Gene\\tCount\\n" > "${counts_tsv}"
    cat "${counts_tsv}.tmp" >> "${counts_tsv}"
    rm -f "${counts_tsv}.tmp"

    # Record read identifiers with at least one qualifying hit
    awk -v e="\$E" '
      \$0 ~ /^#/ { next }
      NF >= 14 {
        q=\$3;
        eval=\$13+0;
        if (eval <= e) print q;
      }
    ' "${tblout}" | sort -u > "${read_ids}"

    # Retain the FASTA sequences corresponding to qualifying reads
    if [[ -s "${read_ids}" ]]; then
        seqkit grep \
            -f "${read_ids}" \
            "${sampleID}.merged.fa" \
            > "${filtered_fa}"
    else
        : > "${filtered_fa}"
    fi

    HMM_READS=\$(wc -l < "${read_ids}")
    HMM_GENES_DETECTED=\$(awk -F '\t' '
        \$1 ~ /^clb[A-S]\$/ && (\$2 + 0) > 0 { count++ }
        END { print count + 0 }
    ' "${counts_tsv}")

    printf "Sample\\tMetric\\tValue\\n" > "${qc}"
    printf "%s\\treads_clb_genes_hmm\\t%s\\n" "${sampleID}" "\$HMM_READS" >> "${qc}"
    printf "%s\\tnum_clb_genes_hmm\\t%s\\n" "${sampleID}" "\$HMM_GENES_DETECTED" >> "${qc}"
    """
}
