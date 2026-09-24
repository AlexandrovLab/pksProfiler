nextflow.enable.dsl=2

include { strainTyping as magStrainTyping } from './pks_typing.nf'

// ─── Assembly ─────────────────────────────────────────────────────────────────

process megahitAssemble {
    label 'sample_stage'
    label 'mag_assembly'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/assembly" },
        mode: 'copy', enabled: params.save_intermediates.toString().toBoolean()
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
    label 'sample_stage'
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
    label 'sample_stage'
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
    label 'sample_stage'
    label 'mag_binning'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/bins" },
        mode: 'copy', enabled: params.save_intermediates.toString().toBoolean()
    conda "${projectDir}/conda_envs/metabat2_env.yml"

    input:
    tuple val(sampleID), path(contigs), path(depth)

    output:
    // Bin FASTAs are optional, but the status file is mandatory so every
    // successfully assembled sample remains visible even when no MAG is recovered.
    // Numeric suffix only: MetaBAT2 names real bins bin.<N>.fa, so this glob
    // deliberately excludes bin.unbinned.fa below -- pooling that file in here
    // would make magSummaryTable/GTDB-Tk/CheckM2 treat unplaced contigs as a
    // genome, exactly the kind of unsupported positivity claim this codebase
    // has just finished removing elsewhere (--tumor_enable_mags, same commit
    // series: no evidence MAG binning works on a tiny, low-diversity fraction).
    tuple val(sampleID), path("bins/bin.[0-9]*.fa"), emit: bins, optional: true
    // U1: contigs MetaBAT2 could not confidently place in any bin are common when
    // abundance is low or binning is underpowered relative to community complexity
    // (confirmed in real testing, 2026-09-24: metabat2 --unbinned pooled every
    // input contig when total assembly size fell under MetaBAT2's 200 kb
    // --minClsSize floor) and, until now, were never checked for pks-island
    // alignment evidence at all -- a real detection blind spot when the read-level
    // lane still shows signal for the same sample. --unbinned pools them so the
    // same canonical-locus alignment used for a real bin (below) can still be
    // tried against them.
    tuple val(sampleID), path("bins/bin.unbinned.fa"), emit: unbinned, optional: true
    tuple val(sampleID), path("${sampleID}.assembly_binning.tsv"), emit: status

    script:
    """
    set -euo pipefail
    mkdir -p bins

    CONTIG_COUNT=\$(grep -c '^>' "${contigs}" || true)
    if [[ "\$CONTIG_COUNT" -gt 0 ]]; then
        metabat2 -i ${contigs} -a ${depth} -o bins/bin -t ${task.cpus} --unbinned || {
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

    BIN_COUNT=\$(find bins -maxdepth 1 -type f -name 'bin.[0-9]*.fa' | wc -l)
    printf "sample\tcontig_count\tbin_count\n" > "${sampleID}.assembly_binning.tsv"
    printf "%s\t%s\t%s\n" "${sampleID}" "\$CONTIG_COUNT" "\$BIN_COUNT" \
        >> "${sampleID}.assembly_binning.tsv"
    """
}

// ─── CheckM2 quality assessment ───────────────────────────────────────────────

process checkm2Predict {
    label 'sample_stage'
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
    label 'sample_stage'
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
    label 'sample_stage'
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
    label 'sample_stage'
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
    # query model, field 5 the E-value -- so the two agree by construction. This count is
    # a cheap pre-filter only (which bins are worth annotating context for); M1's
    # positivity call is decided by magBinLocusEvidence below, not by this count alone.
    awk '!/^#/ && NF>=19 && (\$5+0) <= ${params.hmm_protein_evalue} { genes[\$3]=1 }
         END { print length(genes)+0 }' ${binID}.tblout > ${binID}.clb_gene_count.txt \
      || echo 0 > ${binID}.clb_gene_count.txt
    """
}

// ─── Canonical-locus alignment: is the island actually IN this unit? ─────────
// Shared by a real MetaBAT2 bin and, since U1, the per-sample pool of contigs
// MetaBAT2 never placed in any bin (unitID "unbinned" -- see metabat2Bin above
// and pksMAG's locus_alignment_units_ch). One process, generalized on unitID
// rather than duplicated, since Nextflow forbids invoking the same process
// twice in one workflow and the alignment/classification logic is identical
// either way -- only what FASTA is being scored differs.

process alignMagBinToCanonicalReference {
    label 'sample_stage'
    tag "${sampleID}:${unitID}"
    label 'targeted_alignment'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/locus_alignment" }, mode: 'copy'
    conda "${params.targeted_alignment_env}"

    input:
    tuple val(sampleID), val(unitID), path(unit_fa)
    path reference

    output:
    tuple val(sampleID), val(unitID), path("${unitID}.vs_IHE3034.paf"), emit: paf

    script:
    """
    set -euo pipefail
    paf="${unitID}.vs_IHE3034.paf"
    if [[ -s "${unit_fa}" ]]; then
        minimap2 -x asm10 -c -t ${task.cpus} "${reference}" "${unit_fa}" > "\$paf"
    else
        : > "\$paf"
    fi
    """
}

process magBinLocusEvidence {
    label 'sample_stage'
    tag "${sampleID}:${unitID}"
    label 'mag_hmm'
    publishDir { "${params.sample_dir}/${sampleID}/genomes/locus_alignment" }, mode: 'copy'
    conda "${params.pks_align_env}"

    input:
    tuple val(sampleID), val(unitID), path(paf)

    output:
    tuple val(sampleID), val(unitID), path("${unitID}.locus_evidence.tsv"), emit: evidence

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_annotation ${params.dep_digest?.pks_annotation}
    set -euo pipefail
    # Ludmil, revised report finding 6: island-start/-end here are 0-based
    # half-open (PAF's own convention, matching summarize_tumor_pks_contigs.py) --
    # the correct conversion from the 1-based inclusive annotation is
    # pks_start_1based-1 .. pks_end_1based, not the bare pks_shift this used to
    # pass, which was a 5th call site sharing the same off-by-one his report
    # found in the other four.
    # summarize_mag_bin_locus_evidence.py's --bin-id is an opaque label, not
    # necessarily a real bin: unitID is "unbinned" for the pooled-contig pseudo-unit
    # (U1), and the script's own bin_id output column carries that value through
    # unchanged -- no script change needed to support it.
    python ${projectDir}/scripts/summarize_mag_bin_locus_evidence.py \
        --sample ${sampleID} --bin-id ${unitID} --paf ${paf} \
        --gff ${params.pks_genome_annotation} --contig ${params.pks_contig} \
        --island-start \$(( ${params.pks_start_1based} - 1 )) \
        --island-end ${params.pks_end_1based} \
        --min-identity ${params.mag_locus_min_identity} --min-mapq ${params.mag_locus_min_mapq} \
        --multi-gene-genes ${params.mag_locus_multi_gene_min_genes} \
        --multi-gene-breadth ${params.mag_locus_multi_gene_min_breadth} \
        --broad-island-genes ${params.mag_locus_broad_island_min_genes} \
        --broad-island-breadth ${params.mag_locus_broad_island_min_breadth} \
        --extensive-island-genes ${params.mag_locus_extensive_island_min_genes} \
        --extensive-island-breadth ${params.mag_locus_extensive_island_min_breadth} \
        --output ${unitID}.locus_evidence.tsv
    """
}

// ─── Community-wide prophage prediction (ALL bins) ───────────────────────────

process genomadProphages {
    label 'sample_stage'
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
    label 'sample_stage'
    label 'mag_prophage'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/community" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), path(gtdbtk_summary), path(gffs), path(tblouts), path(virus_summaries), path(locus_evidence)

    output:
    path "${sampleID}.community_prophage_inventory.tsv", emit: inventory
    path "${sampleID}.pks_community_interactions.tsv", emit: interactions
    path "${sampleID}.pks_island_mobility.tsv", emit: mobility

    script:
    """
    set -euo pipefail
    mkdir -p gff_dir tblout_dir genomad_dir locus_dir
    for f in ${gffs}; do ln -s \$(realpath \$f) gff_dir/; done
    for f in ${tblouts}; do ln -s \$(realpath \$f) tblout_dir/; done
    for f in ${virus_summaries}; do ln -s \$(realpath \$f) genomad_dir/; done
    for f in ${locus_evidence}; do ln -s \$(realpath \$f) locus_dir/; done
    python ${projectDir}/scripts/build_community_prophage.py \
        --sample ${sampleID} \
        --gtdbtk ${gtdbtk_summary} \
        --gff-dir gff_dir \
        --tblout-dir tblout_dir \
        --genomad-dir genomad_dir \
        --locus-dir locus_dir \
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
    label 'sample_stage'
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
    label 'sample_stage'
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
    label 'sample_stage'
    label 'mag_hmm'
    publishDir { "${params.sample_dir}/${sampleID}/genomes" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_hmm_env}"
    input:
    tuple val(sampleID), path(checkm2_report), path(gtdbtk_summary),
          path(tblouts), path(contexts), path(locus_evidence)

    output:
    path "${sampleID}.pks_mag_summary.tsv", emit: summary

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_reference_fasta ${params.dep_digest?.pks_reference_fasta}
    set -euo pipefail
    mkdir -p tblout_dir context_dir locus_dir
    for f in ${tblouts}; do ln -s "\$(realpath "\$f")" "tblout_dir/\$(basename "\$f")"; done
    for f in ${contexts}; do ln -s "\$(realpath "\$f")" "context_dir/\$(basename "\$f")"; done
    for f in ${locus_evidence}; do ln -s "\$(realpath "\$f")" "locus_dir/\$(basename "\$f")"; done
    python ${projectDir}/scripts/build_mag_summary.py \
        --checkm2 ${checkm2_report} \
        --gtdbtk ${gtdbtk_summary} \
        --tblout_dir tblout_dir \
        --context_dir context_dir \
        --locus_dir locus_dir \
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

    // U1: the pooled contigs MetaBAT2 never placed in any bin -- one pseudo-unit per
    // sample, "unbinned", standing in for binID the same way bins_flat_ch's real bin
    // IDs do. Confirmed missing today: a contig can carry real pks-island alignment
    // evidence and never surface anywhere once MetaBAT2 declines to bin it (common
    // when abundance is low or the community is underpowered for composition-based
    // binning), even though the read-level pksProfilerAlign/pksProfilerHMM lane may
    // still show signal for the same sample. optional:true on metabat2Bin.out.unbinned
    // means a sample with nothing unbinned (or no contigs at all) simply contributes no
    // entry here, not an empty-string pseudo-unit.
    unbinned_flat_ch = metabat2Bin.out.unbinned
        .map { sampleID, unbinned_fa -> tuple(sampleID, "unbinned", unbinned_fa) }

    // Fed through the exact same alignment + classification processes as a real bin,
    // mixed into one channel first: Nextflow forbids invoking a process twice in one
    // workflow, and the alignment/classification logic does not differ by unitID, only
    // the FASTA being scored does.
    locus_alignment_units_ch = bins_flat_ch.mix(unbinned_flat_ch)

    // 5. CheckM2 quality (all bins together per sample)
    checkm2Predict(metabat2Bin.out.bins)

    // 6. GTDB-Tk taxonomy (all bins together per sample)
    gtdbtkClassify(metabat2Bin.out.bins)

    // 6b. M1: is the colibactin island actually IN this bin? Align every bin's own
    // assembly to the canonical IHE3034 locus and score it the same way read-level
    // evidence already is -- genes an alignment actually covers, and breadth of the
    // island those alignments span -- rather than trusting HMM domain hits alone.
    // Runs for every bin, not only HMM-flagged ones: it is cheap next to
    // CheckM2/GTDB-Tk/Prokka, and a bin whose true positive genes fell just under
    // params.mag_min_clb_genes should not be invisible to this check. Since U1, also
    // runs once more for the sample's unbinned pool via the same mixed channel.
    pks_reference_ch = Channel.value(file(params.pks_reference_fasta, checkIfExists: true))
    alignMagBinToCanonicalReference(locus_alignment_units_ch, pks_reference_ch)
    magBinLocusEvidence(alignMagBinToCanonicalReference.out.paf)

    // U1: split the combined evidence back into real-bin and unbinned-pool streams.
    // Kept apart from here on -- a positive unbinned call means "this sample carries
    // clb sequence outside any recovered genome," not "this genome carries the
    // island," and folding it into the per-bin counts below would conflate the two
    // (pks_positive_bin_counts_ch, community_input_ch's taxonomy-driven interaction
    // table). It still reaches magSummaryTable/build_mag_summary.py, which reports it
    // as its own row -- see summary_input_ch below.
    bin_locus_evidence_ch = magBinLocusEvidence.out.evidence
        .filter { _sampleID, unitID, _evidence -> unitID != "unbinned" }
        .map { sampleID, _unitID, evidence -> tuple(sampleID, evidence) }
    unbinned_locus_evidence_ch = magBinLocusEvidence.out.evidence
        .filter { _sampleID, unitID, _evidence -> unitID == "unbinned" }
        .map { sampleID, _unitID, evidence -> tuple(sampleID, evidence) }

    // 7. Prokka annotation (ALL bins — before hmmsearch so locus_tags are consistent)
    prokkaAnnotate(bins_flat_ch)

    // Explicit provirus calls for every recovered community MAG.
    genomadProphages(bins_flat_ch)

    // 7b. Strain typing per bin: MLST -> ST -> clonal complex / phylogroup.
    // Declared outside the if so main.nf's masterSummary gate has something to
    // read regardless of whether this run enabled strain typing.
    def strain_typing_summary_ch = channel.empty()
    if (params.enable_strain_typing.toString().toBoolean()) {
        magStrainTyping(bins_flat_ch.map { sampleID, binID, bin_fa -> tuple(sampleID, binID, bin_fa) })
        strain_typing_summary_ch = magStrainTyping.out.summary
    }

    // 8. hmmsearch vs clb protein HMM (per bin, on prokka proteins)
    hmmsearchClb(prokkaAnnotate.out.faa_for_hmm)

    // 9. Filter to pks+ bins on the distinct-gene count (avoids reading tblouts in the
    // driver JVM). M3: one definition of pks-positive, shared with the status count
    // below and with community producer selection -- params.mag_min_clb_genes. This is
    // a cheap HMM pre-filter for which bins are worth genomic-context extraction; the
    // actual positivity call is magBinLocusEvidence's alignment-confirmed tier, below.
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

    // Real bins only -- what communityProphageSummary reads below (a taxonomy-driven,
    // per-genome analysis the unbinned pseudo-unit has no taxonomy or annotation for).
    locus_evidence_per_sample_ch = bin_locus_evidence_ch
        .groupTuple(by: 0)

    // U1: bins plus the unbinned pool, for magSummaryTable/build_mag_summary.py only --
    // it reports the unbinned call as its own distinct row (unit_type=unbinned), never
    // merged into a bin's.
    mag_summary_locus_evidence_ch = bin_locus_evidence_ch
        .mix(unbinned_locus_evidence_ch)
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

    // M1: a bin counts toward pks_positive_bin_count only when magBinLocusEvidence's
    // alignment-confirmed tier clears the same bar read-level evidence uses to call a
    // sample positive (multi_gene or above) -- not on the HMM gene count alone, which
    // domain-homology hits to the megasynthases can clear without island carriage.
    // Real bins only (bin_locus_evidence_ch): U1's unbinned pool is not a genome, so a
    // positive call there must never inflate this count -- see mag_summary_locus_evidence_ch
    // above for where the unbinned call is reported instead.
    pks_positive_bin_counts_ch = bin_locus_evidence_ch
        .map { sampleID, evidence_file ->
            def lines = evidence_file.text.readLines().findAll { it.trim() }
            if (lines.size() != 2) {
                error "Malformed locus evidence for ${sampleID}: ${evidence_file}"
            }
            def tier = lines[1].split("\\t", -1)[2]
            tuple(sampleID, (tier in ["multi_gene", "broad_island", "extensive_island"]) ? 1 : 0)
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
    // locus_evidence here is mag_summary_locus_evidence_ch (bins + the U1 unbinned
    // pool), not locus_evidence_per_sample_ch -- build_mag_summary.py's locus_dir
    // needs the unbinned call too, to report it as its own row.
    summary_input_ch = checkm2Predict.out.report
        .join(gtdbtkClassify.out.summary, by: 0)
        .join(tblouts_per_sample_ch, by: 0)
        .join(contexts_per_sample_ch, by: 0, remainder: true)
        .join(mag_summary_locus_evidence_ch, by: 0, remainder: true)
        .map { sampleID, checkm2, gtdbtk, tblouts, contexts, locus_evidence ->
            tuple(sampleID, checkm2, gtdbtk, tblouts, contexts ?: [], locus_evidence ?: [])
        }

    magSummaryTable(summary_input_ch)

    community_input_ch = gtdbtkClassify.out.summary
        .join(gffs_per_sample_ch, by: 0)
        .join(tblouts_per_sample_ch, by: 0)
        .join(virus_summaries_per_sample_ch, by: 0)
        .join(locus_evidence_per_sample_ch, by: 0, remainder: true)
        .map { sampleID, gtdbtk, gffs, tblouts, virus_summaries, locus_evidence ->
            tuple(sampleID, gtdbtk, gffs, tblouts, virus_summaries, locus_evidence ?: [])
        }

    communityProphageSummary(community_input_ch)

    emit:
    // Consumed by main.nf's masterSummary gate (build_master_summary.py's
    // MAG/community/strain-typing columns) -- not by cohortReport, which never
    // reads this lane's output.
    mag_summary       = magSummaryTable.out.summary
    community_summary = communityProphageSummary.out.inventory
                            .mix(communityProphageSummary.out.interactions)
                            .mix(communityProphageSummary.out.mobility)
    // strainTypeSummary emits tuple(val(sampleID), path(...)) -- strip the sampleID so
    // this matches mag_summary/community_summary's bare-path shape. Left as the tuple,
    // main.nf's optional_lane_gate.mix(...) silently accepts it, and masterSummary's
    // path-typed collected input later fails with "Not a valid path value: '<sampleID>'"
    // the moment strain typing is enabled -- caught live on the 2026-09-24 validate1 run.
    strain_summary    = strain_typing_summary_ch.map { sampleID, table -> table }
    // U1: the unbinned pool's own canonical-locus call, one file per sample that had
    // anything left unbinned. Bare-path, matching mag_summary/community_summary/
    // strain_summary's shape (see the strain_summary comment above for why that shape
    // matters to main.nf's optional_lane_gate). Already folded into mag_summary itself
    // (pks_mag_summary.tsv carries it as a distinct row -- see build_mag_summary.py);
    // emitted again here, named, so it is not only discoverable by reading a table.
    unbinned_locus_evidence = unbinned_locus_evidence_ch.map { sampleID, evidence -> evidence }
}
