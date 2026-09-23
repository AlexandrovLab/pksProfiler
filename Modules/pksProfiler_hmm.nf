nextflow.enable.dsl = 2

process pksProfiler_hmm {
    label 'pks_hmm'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/hmm" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
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
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # hmm_model ${params.dep_digest?.hmm_model}  scripts ${params.dep_digest?.scripts}
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
        # Ludmil, revised report finding 2: the zero-read branch omitted this metric
        # entirely, so a true zero sample showed NA in cohort QC instead of 0.
        printf "%s\\thmm_ambiguous_reads\\t0\\n" "${sampleID}" >> "${qc}"
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

    # F07: a typed parser, not awk coercion. `\$13+0` turned a non-numeric E-value
    # into 0, which passes any stringency threshold and became the best hit in the
    # file; short rows were skipped in silence; and ties were settled by whichever
    # line came first, so the answer depended on row order and therefore on chunking.
    #
    # The rule now: lowest E-value, then highest score, then longest aligned length on
    # the read -- the same idea as featureCounts' --largestOverlap -- and a read still
    # tied across two genes after all three is ambiguous and counts for neither, which
    # is what featureCounts does with a perfect tie.
    # Written on one line deliberately. With backslash continuations this rendered with
    # every argument shifted by one position -- the script path became the tblout name --
    # and the first the run knew of it was python3 reporting a path that does not exist,
    # on a task whose failure the per-sample error strategy ignores. One line cannot shift.
    BEST_HIT="${params.scripts}/hmm_best_hit.py"
    if [[ ! -f "\$BEST_HIT" ]]; then
        echo "ERROR: HMM best-hit parser not found at \$BEST_HIT" >&2
        exit 1
    fi
    python3 "\$BEST_HIT" --tblout "${tblout}" --evalue "\$E" --counts-out "${counts_tsv}" --read-ids-out "${read_ids}" --ambiguous-out "${sampleID}.hmm.ambiguous_reads.txt"

    # Retain the FASTA sequences corresponding to qualifying reads
    if [[ -s "${read_ids}" ]]; then
        seqkit grep \
            -f "${read_ids}" \
            "${sampleID}.merged.fa" \
            > "${filtered_fa}"
    else
        : > "${filtered_fa}"
    fi

    # Ludmil, revised report finding 1: this used to be `wc -l < read_ids`, but
    # read_ids (and filtered_fa, above) hold every *qualifying* read -- assigned
    # plus ambiguous -- while the counts matrix only tallies assigned reads. The
    # manuscript-facing metric must equal the matrix total, the same invariant
    # already held for the alignment lane; ambiguity is reported separately below,
    # not folded into this count. The qualifying-read FASTA keeps its current name
    # and content -- still useful downstream -- this only changes what gets counted.
    HMM_READS=\$(awk -F '\t' 'NR>1 { sum += \$2 } END { print sum+0 }' "${counts_tsv}")
    HMM_GENES_DETECTED=\$(awk -F '\t' '
        \$1 ~ /^clb[A-S]\$/ && (\$2 + 0) > 0 { count++ }
        END { print count + 0 }
    ' "${counts_tsv}")

    printf "Sample\\tMetric\\tValue\\n" > "${qc}"
    printf "%s\\treads_clb_genes_hmm\\t%s\\n" "${sampleID}" "\$HMM_READS" >> "${qc}"
    AMBIGUOUS=\$(wc -l < "${sampleID}.hmm.ambiguous_reads.txt")
    printf "%s\\thmm_ambiguous_reads\\t%s\\n" "${sampleID}" "\$AMBIGUOUS" >> "${qc}"
    printf "%s\\tnum_clb_genes_hmm\\t%s\\n" "${sampleID}" "\$HMM_GENES_DETECTED" >> "${qc}"
    """
}
