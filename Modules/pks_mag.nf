nextflow.enable.dsl=2

include { strainTyping as magStrainTyping } from './pks_typing.nf'
include { strainTyping as tumorStrainTyping } from './pks_typing.nf'

// ─── Assembly ─────────────────────────────────────────────────────────────────

process megahitAssemble {
    label 'mag_assembly'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/assembly" },
        mode: 'copy', enabled: params.save_intermediates
    conda "${projectDir}/conda_envs/megahit_env.yml"

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${sampleID}.contigs.fa"), emit: contigs
    tuple val(sampleID), path("${sampleID}.contig_count.txt"), emit: contig_count

    script:
    def reads_list = reads instanceof List ? reads : [reads]
    if (!(reads_list.size() in [1, 2])) {
        error "MAG assembly for ${sampleID} requires one single/interleaved FASTQ or an R1/R2 pair; received ${reads_list.size()} files"
    }
    def r_flag = reads_list.size() == 2
        ? "-1 ${reads_list[0]} -2 ${reads_list[1]}"
        : "-r ${reads_list[0]}"
    """
    set -euo pipefail
    megahit ${r_flag} -t ${task.cpus} -o megahit_out
    cp megahit_out/final.contigs.fa ${sampleID}.contigs.fa
    grep -c '^>' ${sampleID}.contigs.fa > ${sampleID}.contig_count.txt || \
        printf '0\n' > ${sampleID}.contig_count.txt
    """
}

// ─── Read-to-contig alignment ─────────────────────────────────────────────────

process alignToContigs {
    label 'mag_binning'
    scratch true
    conda "${projectDir}/conda_envs/minimap2_env.yml"

    input:
    tuple val(sampleID), path(reads), path(contigs)

    output:
    tuple val(sampleID), path("${sampleID}.contigs.sorted.bam"),
          path("${sampleID}.contigs.sorted.bam.bai"), emit: bam

    script:
    def reads_list = reads instanceof List ? reads : [reads]
    if (!(reads_list.size() in [1, 2])) {
        error "MAG contig alignment for ${sampleID} requires one single/interleaved FASTQ or an R1/R2 pair; received ${reads_list.size()} files"
    }
    def read_args = reads_list.join(' ')
    """
    set -euo pipefail
    minimap2 -ax sr -t ${task.cpus} ${contigs} ${read_args} \
        | samtools view -F 4 -bS \
        | samtools sort -@ ${task.cpus} -o ${sampleID}.contigs.sorted.bam
    samtools index ${sampleID}.contigs.sorted.bam
    """
}

// ─── Contig coverage depth ────────────────────────────────────────────────────

process jgiContigDepths {
    label 'mag_binning'
    scratch true
    conda "${projectDir}/conda_envs/metabat2_env.yml"

    input:
    tuple val(sampleID), path(bam), path(bai)

    output:
    tuple val(sampleID), path("${sampleID}.depth.txt"), emit: depth

    script:
    """
    set -euo pipefail
    jgi_summarize_bam_contig_depths --outputDepth ${sampleID}.depth.txt ${bam}
    """
}

// ─── MetaBAT2 binning ─────────────────────────────────────────────────────────

process metabat2Bin {
    label 'mag_binning'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/bins" },
        mode: 'copy', enabled: params.save_intermediates
    conda "${projectDir}/conda_envs/metabat2_env.yml"

    input:
    tuple val(sampleID), path(contigs), path(depth)

    output:
    // Bin FASTAs are optional, but the status file is mandatory so every
    // successfully assembled sample remains visible even when no MAG is recovered.
    tuple val(sampleID), path("bins/bin.*.fa"), emit: bins, optional: true
    tuple val(sampleID), path("${sampleID}.assembly_binning.tsv"), emit: status

    script:
    """
    set -euo pipefail
    mkdir -p bins

    CONTIG_COUNT=\$(grep -c '^>' "${contigs}" || true)
    if [[ "\$CONTIG_COUNT" -gt 0 ]]; then
        metabat2 -i ${contigs} -a ${depth} -o bins/bin -t ${task.cpus} || {
            rc=\$?
            if [[ \${rc} -eq 1 ]]; then
                echo "MetaBAT2 exited 1 without recoverable bins" >&2
            else
                echo "MetaBAT2 failed with exit code \${rc}" >&2
                exit \${rc}
            fi
        }
    else
        echo "Assembly produced no contigs; skipping MetaBAT2" >&2
    fi

    BIN_COUNT=\$(find bins -maxdepth 1 -type f -name 'bin.*.fa' | wc -l)
    printf "sample\tcontig_count\tbin_count\n" > "${sampleID}.assembly_binning.tsv"
    printf "%s\t%s\t%s\n" "${sampleID}" "\$CONTIG_COUNT" "\$BIN_COUNT" \
        >> "${sampleID}.assembly_binning.tsv"
    """
}

// ─── CheckM2 quality assessment ───────────────────────────────────────────────

process checkm2Predict {
    label 'mag_binning'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/checkm2" }, mode: 'copy'
    conda "${projectDir}/conda_envs/checkm2_env.yml"

    input:
    tuple val(sampleID), path(bins)

    output:
    tuple val(sampleID), path("checkm2_out/quality_report.tsv"), emit: report

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # checkm2_db ${params.dep_digest?.checkm2_db}
    set -euo pipefail
    mkdir -p bin_input
    for b in ${bins}; do ln -s \$(realpath \$b) bin_input/; done
    checkm2 predict --threads ${task.cpus} \
        --input bin_input \
        -x fa \
        --output-directory checkm2_out \
        --database_path ${params.checkm2_db}
    """
}

// ─── GTDB-Tk taxonomy ─────────────────────────────────────────────────────────

process gtdbtkClassify {
    label 'mag_gtdbtk'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/gtdbtk" }, mode: 'copy'
    conda "${projectDir}/conda_envs/gtdbtk_env.yml"

    input:
    tuple val(sampleID), path(bins)

    output:
    tuple val(sampleID), path("gtdbtk_out/gtdbtk.bac120.summary.tsv"), emit: summary

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # gtdbtk_db ${params.dep_digest?.gtdbtk_db}
    set -euo pipefail
    mkdir -p bin_input gtdbtk_out
    for b in ${bins}; do ln -s \$(realpath \$b) bin_input/; done
    GTDBTK_DATA_PATH=${params.gtdbtk_db} gtdbtk classify_wf \
        --genome_dir bin_input \
        --out_dir gtdbtk_out \
        --cpus ${task.cpus} \
        --extension fa
    # Stub if no bacterial bins (only archaea)
    [[ -f gtdbtk_out/gtdbtk.bac120.summary.tsv ]] || \
        printf "user_genome\tclassification\n" > gtdbtk_out/gtdbtk.bac120.summary.tsv
    """
}

// ─── Prokka annotation (ALL bins — must run before hmmsearch so locus_tags match) ──

process prokkaAnnotate {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/annotation" }, mode: 'copy',
        saveAs: { fn -> fn.tokenize('/').last() }
    conda "${projectDir}/conda_envs/prokka_env.yml"

    input:
    tuple val(sampleID), val(binID), path(bin_fa)

    output:
    tuple val(sampleID), val(binID), path("prokka_out/${binID}.gff"), emit: gff
    tuple val(sampleID), val(binID), path("prokka_out/${binID}.faa"), emit: faa_for_hmm

    script:
    """
    set -euo pipefail
    prokka --outdir prokka_out --prefix ${binID} \
        --metagenome --cpus ${task.cpus} --force --quiet \
        ${bin_fa}
    """
}

// ─── hmmsearch vs colibactin protein HMM (per bin, on Prokka proteins) ──────────

process hmmsearchClb {
    label 'mag_hmm'
    scratch true
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), val(binID), path(proteins)

    output:
    tuple val(sampleID), val(binID), path("${binID}.tblout"), emit: tblout
    tuple val(sampleID), val(binID), path("${binID}.clb_gene_count.txt"), emit: clb_gene_count

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # clb_protein_hmm ${params.dep_digest?.clb_protein_hmm}
    set -euo pipefail
    hmmsearch --cpu ${task.cpus} \
        --tblout ${binID}.tblout \
        -E ${params.hmm_protein_evalue} \
        ${params.clb_protein_hmm} \
        ${proteins}
    # M3: distinct clb genes, not hit lines. `grep -vc '^#'` counted every
    # protein-to-model match, so one gene matched by three predicted proteins counted
    # three times and a single line anywhere made the bin pks-positive. This counts the
    # same quantity build_mag_summary.py reports as clb_genes_detected -- field 3 is the
    # query model, field 5 the E-value -- so the two agree by construction.
    awk '!/^#/ && NF>=19 && (\$5+0) <= ${params.hmm_protein_evalue} { genes[\$3]=1 }
         END { print length(genes)+0 }' ${binID}.tblout > ${binID}.clb_gene_count.txt \
      || echo 0 > ${binID}.clb_gene_count.txt
    """
}

// ─── Community-wide prophage prediction (ALL bins) ───────────────────────────

process genomadProphages {
    label 'mag_prophage'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/community/prophages" }, mode: 'copy'
    conda "${projectDir}/conda_envs/genomad_env.yml"

    input:
    tuple val(sampleID), val(binID), path(bin_fa)

    output:
    tuple val(sampleID), val(binID), path("${binID}.virus_summary.tsv"), emit: summary

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # genomad_db ${params.dep_digest?.genomad_db}
    set -euo pipefail
    genomad end-to-end --cleanup --threads ${task.cpus} --splits ${task.cpus} \
        ${bin_fa} genomad_out ${params.genomad_db}
    summary=\$(find genomad_out -type f -name '*_virus_summary.tsv' -print -quit)
    if [[ -z "\$summary" ]]; then
        echo "geNomad did not produce a virus summary for ${binID}" >&2
        exit 1
    fi
    cp "\$summary" ${binID}.virus_summary.tsv
    """
}

// ─── Sample-level producer/lysogen hypothesis tables ─────────────────────────

process communityProphageSummary {
    label 'mag_prophage'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/community" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), path(gtdbtk_summary), path(gffs), path(tblouts), path(virus_summaries)

    output:
    path "${sampleID}.community_prophage_inventory.tsv", emit: inventory
    path "${sampleID}.pks_community_interactions.tsv", emit: interactions
    path "${sampleID}.pks_island_mobility.tsv", emit: mobility

    script:
    """
    set -euo pipefail
    mkdir -p gff_dir tblout_dir genomad_dir
    for f in ${gffs}; do ln -s \$(realpath \$f) gff_dir/; done
    for f in ${tblouts}; do ln -s \$(realpath \$f) tblout_dir/; done
    for f in ${virus_summaries}; do ln -s \$(realpath \$f) genomad_dir/; done
    python ${projectDir}/scripts/build_community_prophage.py \
        --sample ${sampleID} \
        --gtdbtk ${gtdbtk_summary} \
        --gff-dir gff_dir \
        --tblout-dir tblout_dir \
        --genomad-dir genomad_dir \
        --evalue ${params.hmm_protein_evalue} \
        --min-clb-genes ${params.mag_min_clb_genes} \
        --min-provirus-length ${params.provirus_min_length_bp} \
        --min-provirus-hallmarks ${params.provirus_min_hallmarks} \
        --prophage-out ${sampleID}.community_prophage_inventory.tsv \
        --interaction-out ${sampleID}.pks_community_interactions.tsv \
        --mobility-out ${sampleID}.pks_island_mobility.tsv
    """
}
// ─── Genomic context extraction (pks+ bins only) ──────────────────────────────

process extractGenomicContext {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/community/genomic_context" }, mode: 'copy',
        saveAs: { "${binID}.context.tsv" }
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), val(binID), path(gff), path(tblout)

    output:
    tuple val(sampleID), val(binID), path("${binID}.context.tsv"), emit: context

    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/extract_genomic_context.py \
        --gff ${gff} \
        --tblout ${tblout} \
        --evalue ${params.hmm_protein_evalue} \
        --window ${params.context_window_bp} \
        --out ${binID}.context.tsv
    """
}

// ─── Per-sample MAG outcome status ─────────────────────────────────────────────

process magSampleStatus {
    label 'mag_hmm'
    publishDir { "${params.sample_dir}/${sampleID}/genomes" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }

    input:
    tuple val(sampleID), val(contig_count), val(bin_count), val(pks_positive_bins)

    output:
    path "${sampleID}.mag_status.tsv", emit: status

    script:
    def outcome = contig_count == 0
        ? "no_contigs"
        : bin_count == 0
            ? "contigs_no_bins"
            : pks_positive_bins == 0
                ? "bins_no_pks"
                : "pks_positive_bins"
    """
    set -euo pipefail
    printf "sample\\tstatus\\tcontig_count\\tbin_count\\tpks_positive_bin_count\\n" \
        > "${sampleID}.mag_status.tsv"
    printf "%s\\t%s\\t%s\\t%s\\t%s\\n" \
        "${sampleID}" "${outcome}" "${contig_count}" "${bin_count}" "${pks_positive_bins}" \
        >> "${sampleID}.mag_status.tsv"
    """
}


// ─── Per-sample MAG summary table ─────────────────────────────────────────────

process magSummaryTable {
    label 'mag_hmm'
    publishDir { "${params.sample_dir}/${sampleID}/genomes" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_hmm_env}"
    input:
    tuple val(sampleID), path(checkm2_report), path(gtdbtk_summary),
          path(tblouts), path(contexts)

    output:
    path "${sampleID}.pks_mag_summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    mkdir -p tblout_dir context_dir
    for f in ${tblouts}; do ln -s "\$(realpath "\$f")" "tblout_dir/\$(basename "\$f")"; done
    for f in ${contexts}; do ln -s "\$(realpath "\$f")" "context_dir/\$(basename "\$f")"; done
    python ${projectDir}/scripts/build_mag_summary.py \
        --checkm2 ${checkm2_report} \
        --gtdbtk ${gtdbtk_summary} \
        --tblout_dir tblout_dir \
        --context_dir context_dir \
        --sample ${sampleID} \
        --evalue ${params.hmm_protein_evalue} \
        --out ${sampleID}.pks_mag_summary.tsv
    """
}

// ─── pksMAG workflow ──────────────────────────────────────────────────────────

workflow pksMAG {
    take:
    reads  // tuple val(sampleID), path(reads)

    main:

    // 1. Assembly
    megahitAssemble(reads)

    assembly_contig_counts_ch = megahitAssemble.out.contig_count
        .map { sampleID, count_file ->
            def count_text = count_file.text.trim()
            if (!(count_text ==~ /^[0-9]+$/)) {
                error "Invalid contig count for ${sampleID}: ${count_text}"
            }
            tuple(sampleID, count_text.toInteger())
        }

    no_contig_counts_ch = assembly_contig_counts_ch
        .filter { _sampleID, contig_count -> contig_count == 0 }
        .map { sampleID, _contig_count -> tuple(sampleID, 0, 0) }

    nonempty_contigs_ch = megahitAssemble.out.contigs
        .join(assembly_contig_counts_ch, by: 0)
        .filter { _sampleID, _contigs, contig_count -> contig_count > 0 }
        .map { sampleID, contigs, _contig_count -> tuple(sampleID, contigs) }

    // 2. Align reads only when assembly produced at least one contig.
    reads_contigs_ch = reads.join(nonempty_contigs_ch, by: 0)
    alignToContigs(reads_contigs_ch)

    // 3. Contig coverage depth
    jgiContigDepths(alignToContigs.out.bam)

    // 4. Bin contigs
    contigs_depth_ch = nonempty_contigs_ch
        .join(jgiContigDepths.out.depth, by: 0)
    metabat2Bin(contigs_depth_ch)

    // Fan-out: one entry per bin → tuple(sampleID, binID, bin_fa)
    bins_flat_ch = metabat2Bin.out.bins
        .transpose()
        .map { sampleID, bin_fa ->
            def binID = bin_fa.baseName
            tuple(sampleID, binID, bin_fa)
        }

    // 5. CheckM2 quality (all bins together per sample)
    checkm2Predict(metabat2Bin.out.bins)

    // 6. GTDB-Tk taxonomy (all bins together per sample)
    gtdbtkClassify(metabat2Bin.out.bins)

    // 7. Prokka annotation (ALL bins — before hmmsearch so locus_tags are consistent)
    prokkaAnnotate(bins_flat_ch)

    // Explicit provirus calls for every recovered community MAG.
    genomadProphages(bins_flat_ch)

    // 7b. Strain typing per bin: MLST -> ST -> clonal complex / phylogroup.
    if (params.enable_strain_typing.toString().toBoolean()) {
        magStrainTyping(bins_flat_ch.map { sampleID, binID, bin_fa -> tuple(sampleID, binID, bin_fa) })
    }

    // 8. hmmsearch vs clb protein HMM (per bin, on prokka proteins)
    hmmsearchClb(prokkaAnnotate.out.faa_for_hmm)

    // 9. Filter to pks+ bins on the distinct-gene count (avoids reading tblouts in the
    // driver JVM). M3: one definition of pks-positive, shared with the status count
    // below and with community producer selection -- params.mag_min_clb_genes.
    pks_pos_tblout_ch = hmmsearchClb.out.tblout
        .join(hmmsearchClb.out.clb_gene_count, by: [0, 1])
        .filter { sampleID, binID, tblout, gene_count ->
            (gene_count.text.trim() as Integer) >= (params.mag_min_clb_genes as Integer)
        }
        .map { sampleID, binID, tblout, gene_count -> tuple(sampleID, binID, tblout) }

    // 10. Genomic context (prokka GFF + tblout joined per pks+ bin; locus_tags now match)
    extractGenomicContext(prokkaAnnotate.out.gff.join(pks_pos_tblout_ch, by: [0, 1]))

    // 11. Aggregate per sample and build summary
    tblouts_per_sample_ch = hmmsearchClb.out.tblout
        .map { sampleID, binID, tblout -> tuple(sampleID, tblout) }
        .groupTuple(by: 0)

    contexts_per_sample_ch = extractGenomicContext.out.context
        .map { sampleID, binID, ctx -> tuple(sampleID, ctx) }
        .groupTuple(by: 0)

    gffs_per_sample_ch = prokkaAnnotate.out.gff
        .map { sampleID, binID, gff -> tuple(sampleID, gff) }
        .groupTuple(by: 0)

    virus_summaries_per_sample_ch = genomadProphages.out.summary
        .map { sampleID, binID, summary -> tuple(sampleID, summary) }
        .groupTuple(by: 0)

    binned_assembly_counts_ch = metabat2Bin.out.status
        .map { sampleID, status_file ->
            def lines = status_file.text.readLines().findAll { it.trim() }
            if (lines.size() != 2) {
                error "Malformed assembly/binning status for ${sampleID}: ${status_file}"
            }
            def fields = lines[1].split("\\t", -1)
            if (fields.size() != 3 ||
                !(fields[1] ==~ /^[0-9]+$/) ||
                !(fields[2] ==~ /^[0-9]+$/)) {
                error "Invalid assembly/binning counts for ${sampleID}: ${lines[1]}"
            }
            tuple(sampleID, fields[1].toInteger(), fields[2].toInteger())
        }

    assembly_binning_counts_ch = no_contig_counts_ch
        .mix(binned_assembly_counts_ch)

    pks_positive_bin_counts_ch = hmmsearchClb.out.clb_gene_count
        .map { sampleID, _binID, gene_count ->
            tuple(sampleID, gene_count.text.trim().toInteger() >= (params.mag_min_clb_genes as Integer) ? 1 : 0)
        }
        .groupTuple(by: 0)
        .map { sampleID, flags -> tuple(sampleID, flags.sum() as Integer) }

    mag_status_input_ch = assembly_binning_counts_ch
        .join(pks_positive_bin_counts_ch, by: 0, remainder: true)
        .map { sampleID, contig_count, bin_count, pks_positive_bins ->
            tuple(sampleID, contig_count, bin_count, pks_positive_bins ?: 0)
        }

    magSampleStatus(mag_status_input_ch)


    // remainder: true keeps all samples even when contexts_per_sample_ch has no
    // entry (i.e. zero pks+ bins); null → [] so the script receives an empty file list.
    summary_input_ch = checkm2Predict.out.report
        .join(gtdbtkClassify.out.summary, by: 0)
        .join(tblouts_per_sample_ch, by: 0)
        .join(contexts_per_sample_ch, by: 0, remainder: true)
        .map { sampleID, checkm2, gtdbtk, tblouts, contexts ->
            tuple(sampleID, checkm2, gtdbtk, tblouts, contexts ?: [])
        }

    magSummaryTable(summary_input_ch)

    community_input_ch = gtdbtkClassify.out.summary
        .join(gffs_per_sample_ch, by: 0)
        .join(tblouts_per_sample_ch, by: 0)
        .join(virus_summaries_per_sample_ch, by: 0)

    communityProphageSummary(community_input_ch)
}

// ─── Tumor-WGS eligibility and contig analysis ───────────────────────────────

process tumorEligibilityStatus {
    label 'mag_hmm'
    publishDir "${params.cohort_dir}/tumor_eligibility", mode: 'copy'

    input:
    tuple val(sampleID), val(clbGenes), val(pksReads), val(nonhostReads)

    output:
    path "${sampleID}.tumor_wgs_screen.tsv", emit: status

    script:
    def eligible = clbGenes >= params.tumor_min_clb_genes.toString().toInteger() &&
        pksReads >= params.tumor_min_pks_reads.toString().toInteger() &&
        nonhostReads >= params.tumor_min_nonhost_reads.toString().toInteger()
    def reason = eligible ? "eligible" : "below_threshold"
    """
    set -euo pipefail
    printf "sample\\tstatus\\tclb_genes\\tpks_reads\\tnonhost_reads\\tmin_clb_genes\\tmin_pks_reads\\tmin_nonhost_reads\\n" > "${sampleID}.tumor_wgs_screen.tsv"
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \
        "${sampleID}" "${reason}" "${clbGenes}" "${pksReads}" "${nonhostReads}" \
        "${params.tumor_min_clb_genes}" "${params.tumor_min_pks_reads}" \
        "${params.tumor_min_nonhost_reads}" >> "${sampleID}.tumor_wgs_screen.tsv"
    """
}

process prokkaTumorContigs {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/community/annotation" }, mode: 'copy'
    conda "${projectDir}/conda_envs/prokka_env.yml"

    input:
    tuple val(sampleID), path(contigs)

    output:
    tuple val(sampleID), path("prokka_out/${sampleID}.gff"), emit: gff
    tuple val(sampleID), path("prokka_out/${sampleID}.faa"), emit: faa

    script:
    """
    set -euo pipefail
    prokka --outdir prokka_out --prefix ${sampleID} --metagenome \
        --cpus ${task.cpus} --force --quiet ${contigs}
    """
}

process hmmsearchTumorContigs {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/community/annotation" }, mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), path(proteins)

    output:
    tuple val(sampleID), path("${sampleID}.contigs.clb.tblout"), emit: tblout

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # clb_protein_hmm ${params.dep_digest?.clb_protein_hmm}
    set -euo pipefail
    hmmsearch --cpu ${task.cpus} --tblout ${sampleID}.contigs.clb.tblout \
        -E ${params.hmm_protein_evalue} ${params.clb_protein_hmm} ${proteins}
    """
}

process genomadTumorContigs {
    label 'mag_prophage'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/community/prophages" }, mode: 'copy'
    conda "${projectDir}/conda_envs/genomad_env.yml"

    input:
    tuple val(sampleID), path(contigs)

    output:
    tuple val(sampleID), path("${sampleID}.contigs.virus_summary.tsv"), emit: summary

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # genomad_db ${params.dep_digest?.genomad_db}
    set -euo pipefail
    genomad end-to-end --cleanup --threads ${task.cpus} --splits ${task.cpus} \
        ${contigs} genomad_out ${params.genomad_db}
    summary=\$(find genomad_out -type f -name '*_virus_summary.tsv' -print -quit)
    [[ -n "\$summary" ]] || { echo "geNomad produced no tumor contig summary" >&2; exit 1; }
    cp "\$summary" ${sampleID}.contigs.virus_summary.tsv
    """
}

process tumorContigContext {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/community" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), path(gff), path(tblout), path(virusSummary)

    output:
    path "${sampleID}.contig_pks_context.tsv", emit: context
    path "${sampleID}.contig_pks_summary.tsv", emit: summary
    path "${sampleID}.contig_pks_mobility.tsv", emit: mobility

    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/extract_genomic_context.py \
        --gff ${gff} --tblout ${tblout} --evalue ${params.hmm_protein_evalue} \
        --window ${params.context_window_bp} \
        --genomad ${virusSummary} \
        --min-provirus-length ${params.provirus_min_length_bp} \
        --min-provirus-hallmarks ${params.provirus_min_hallmarks} \
        --out ${sampleID}.contig_pks_context.tsv
    awk -F '\\t' 'BEGIN{OFS="\\t"} NR==1{print "sample","locus_tag","clb_gene","contig","has_integrase","has_transposase","nearby_trna","in_prophage","prophage_id","prophage_virus_score"; next} {print "${sampleID}",\$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' \
        ${sampleID}.contig_pks_context.tsv > ${sampleID}.contig_pks_mobility.tsv
    awk -F '\\t' -v sample="${sampleID}" -v minimum="${params.tumor_min_clb_genes}" '
        BEGIN { OFS="\\t" }
        NR > 1 && !seen[\$4 SUBSEP \$2]++ { genes[\$4]++ }
        END {
            print "sample","contig","distinct_clb_genes","pks_class"
            for (contig in genes) {
                class = genes[contig] == 19 ? "complete" : (genes[contig] >= minimum ? "partial" : "fragment")
                print sample,contig,genes[contig],class
            }
        }' ${sampleID}.contig_pks_context.tsv > ${sampleID}.contig_pks_summary.tsv
    """
}

workflow tumorWGS {
    take:
    reads

    main:
    // Eligibility is calculated from pks candidate evidence, but assembly receives
    // the complete host-depleted read stream joined back by sample ID in main.nf.
    megahitAssemble(reads)
    tumor_nonempty_contigs_ch = megahitAssemble.out.contigs
        .join(megahitAssemble.out.contig_count, by: 0)
        .filter { sampleID, contigs, countFile -> countFile.text.trim().toInteger() > 0 }
        .map { sampleID, contigs, countFile -> tuple(sampleID, contigs) }

    if (params.enable_strain_typing.toString().toBoolean()) {
        tumorStrainTyping(tumor_nonempty_contigs_ch.map { sampleID, contigs -> tuple(sampleID, "tumor_contigs", contigs) })
    }

    prokkaTumorContigs(tumor_nonempty_contigs_ch)
    hmmsearchTumorContigs(prokkaTumorContigs.out.faa)
    genomadTumorContigs(tumor_nonempty_contigs_ch)
    tumor_context_input_ch = prokkaTumorContigs.out.gff
        .join(hmmsearchTumorContigs.out.tblout, by: 0)
        .join(genomadTumorContigs.out.summary, by: 0)
    tumorContigContext(tumor_context_input_ch)
}
