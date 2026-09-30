nextflow.enable.dsl = 2

process classifyPksReadEvidence {
    label 'sample_stage'
    tag "$sampleID"
    label 'targeted_alignment'
    publishDir { "${params.sample_dir}/${sampleID}" }, mode: 'copy', pattern: "*_evidence.tsv", saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_align_env}"
    input:
    // depletion_qc: mapReads' fragment, the metagenome reads-per-million denominator.
    // mode_qc: rnaResidualHuman's fragment for --sample_type rna (memo section 4), else [].
    tuple val(sampleID), path(raw_coverage), path(counts), path(bam), path(bai), path(qc), path(depletion_qc), path(mode_qc, stageAs: 'mode_qc/*')
    output:
    // Ludmil, 26 Sep mode memo: one schema per mode, and a file name that says which --
    // read_evidence.tsv (wgs, metagenome), wes_evidence.tsv, rna_evidence.tsv. The prefix
    // is classify_tumor_pks_evidence.py's EVIDENCE_KIND with "_evidence" dropped.
    tuple val(sampleID), path(raw_coverage), path(counts), path("${sampleID}.${params.sample_type.toString() in ['wes', 'rna'] ? params.sample_type : 'read'}_evidence.tsv"), emit: evidence
    // GAMEPLAN rc1 QC terminal state: evidence_outcome (positive|negative|untiered),
    // a QC fragment like every other stage's, so build_qc_summary.py can say how the
    // sample ended. Not published.
    tuple val(sampleID), path("${sampleID}.evidence.qc.tsv"), emit: qc
    script:
    def evidence_file = "${sampleID}.${params.sample_type.toString() in ['wes', 'rna'] ? params.sample_type : 'read'}_evidence.tsv"
    def mode_args = "--mode ${params.sample_type}" +
        (params.sample_type.toString() == 'rna' ? " --rna-library-type ${params.rna_library_type}" +
            " --rna-transcript-min-reads ${params.rna_transcript_min_reads}" : '') +
        (mode_qc ? " --residual-human-qc ${mode_qc}" : '') +
        // Long-read WGS: evidence only. The tiers were calibrated on short-read Bowtie2.
        (params.sample_type.toString() == 'wgs' && HostDepletion.isLong(params.read_type) ? ' --no-tiers' : '')
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'samtools python')}
    # Ludmil, revised report finding 6: pks_shift+1 dropped the true first base of
    # the island. pks_start_1based/pks_end_1based are already correct 1-based
    # inclusive bounds; --island-length must match the row count this now produces
    # (50768, not the old pks_island_len=50767), or classify_tumor_pks_evidence.py's
    # own row-count check raises rather than silently disagreeing.
    samtools depth -aa -s -Q 40 -r "${params.pks_contig}:${params.pks_start_1based}-${params.pks_end_1based}" "${bam}" > "${sampleID}.island.raw_depth.tsv"
    python "${params.scripts}/classify_tumor_pks_evidence.py" --sample "${sampleID}" --qc "${qc}" --counts "${counts}" --depth "${sampleID}.island.raw_depth.tsv" --island-length "${params.pks_island_len_1based}" --multi-gene-reads "${params.tumor_multi_gene_min_pks_reads}" --multi-gene-genes "${params.tumor_multi_gene_min_clb_genes}" --multi-gene-breadth "${params.tumor_multi_gene_min_breadth}" --broad-island-reads "${params.tumor_broad_island_min_pks_reads}" --broad-island-genes "${params.tumor_broad_island_min_clb_genes}" --broad-island-breadth "${params.tumor_broad_island_min_breadth}" --extensive-island-reads "${params.tumor_extensive_island_min_pks_reads}" --extensive-island-genes "${params.tumor_extensive_island_min_clb_genes}" --extensive-island-breadth "${params.tumor_extensive_island_min_breadth}" ${mode_args} --depletion-qc "${depletion_qc}" --output "${evidence_file}" --qc-fragment "${sampleID}.evidence.qc.tsv"
    """
}

process targetedPksRecruit {
    label 'sample_stage'
    tag "$sampleID"
    label 'targeted_recruit'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/contigs/recruitment" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.targeted_recruit_env}"
    input:
    tuple val(sampleID), path(reads)
    output:
    tuple val(sampleID), path("${sampleID}.pks.R1.fastq.gz"), path("${sampleID}.pks.R2.fastq.gz"), path("${sampleID}.pks.single.fastq.gz"), emit: reads
    tuple val(sampleID), path("${sampleID}.pks_recruitment.tsv"), emit: stats
    script:
    // Ludmil, 25 Sep A1/A2: the canonical island's span inside the island+flank recruitment
    // reference (indices/pks_island: flank, island, flank). main.nf checks the index FASTA's
    // header (island=..., flank=...) against these params; the script checks the SAM's @SQ LN.
    def recruit_flank = params.pks_recruit_flank_bp.toString().toInteger()
    def recruit_island_start = recruit_flank + 1
    def recruit_island_end = recruit_flank + params.pks_island_len_1based.toString().toInteger()
    def recruit_reference_length = recruit_island_end + recruit_flank
    def flank_only_arg = params.recruit_include_flank_only.toString().toBoolean() ? '' : ' --exclude-flank-only'
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_recruit_index ${params.dep_digest?.pks_recruit_index}  scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'bowtie2 python')}
    # Ludmil, 25 Sep A1: default best-alignment reporting (no -k 1), so MAPQ means something
    # and recover_recruited_mates.py can hold every read to --min-mapq.
    bowtie2 --very-sensitive-local --no-unal --threads ${task.cpus} -x "${params.pks_recruit_index}" -U "${reads}" -S "${sampleID}.pks_recruitment.sam"
    python "${params.scripts}/recover_recruited_mates.py" --sam "${sampleID}.pks_recruitment.sam" --fastq "${reads}" --r1 "${sampleID}.pks.R1.fastq.gz" --r2 "${sampleID}.pks.R2.fastq.gz" --single "${sampleID}.pks.single.fastq.gz" --stats "${sampleID}.pks_recruitment.tsv" --min-aligned-bases "${params.pks_recruit_min_aligned}" --min-identity "${params.pks_recruit_min_identity}" --min-query-coverage "${params.recruit_min_query_coverage}" --min-mapq "${params.recruit_min_mapq}" --boundary-tolerance-bp "${params.recruit_boundary_tolerance_bp}" --island-start ${recruit_island_start} --island-end ${recruit_island_end} --reference-length ${recruit_reference_length}${flank_only_arg}
    """
}

process targetedPksMegahit {
    label 'sample_stage'
    tag "$sampleID"
    label 'targeted_assembly'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/contigs/megahit" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}.targeted_" }
    conda "${projectDir}/conda_envs/megahit_env.yml"
    input:
    tuple val(sampleID), path(r1), path(r2), path(single)
    output:
    tuple val(sampleID), path("${sampleID}.targeted_megahit.contigs.fa"), emit: contigs
    script:
    """
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'megahit')}
    pairs=\$(gzip -dc "${r1}" | awk 'END {print int(NR/4)}')
    singles=\$(gzip -dc "${single}" | awk 'END {print int(NR/4)}')
    if [[ \$pairs -eq 0 && \$singles -eq 0 ]]; then : > "${sampleID}.targeted_megahit.contigs.fa"; exit 0; fi
    args=()
    if [[ \$pairs -gt 0 ]]; then args+=(-1 "${r1}" -2 "${r2}"); fi
    if [[ \$singles -gt 0 ]]; then args+=(-r "${single}"); fi
    megahit "\${args[@]}" -t ${task.cpus} --min-contig-len ${params.targeted_min_contig_len} -o megahit_out
    cp megahit_out/final.contigs.fa "${sampleID}.targeted_megahit.contigs.fa"
    """
}

process targetedPksSpades {
    label 'sample_stage'
    tag "$sampleID"
    label 'targeted_assembly'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/contigs/metaspades" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}.targeted_" }
    conda "${params.targeted_spades_env}"
    input:
    tuple val(sampleID), path(r1), path(r2), path(single)
    output:
    tuple val(sampleID), path("${sampleID}.targeted_metaspades.contigs.fa"), emit: contigs
    script:
    """
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'metaspades.py:spades')}
    pairs=\$(gzip -dc "${r1}" | awk 'END {print int(NR/4)}')
    singles=\$(gzip -dc "${single}" | awk 'END {print int(NR/4)}')
    if [[ \$pairs -eq 0 && \$singles -eq 0 ]]; then : > "${sampleID}.targeted_metaspades.contigs.fa"; exit 0; fi
    args=()
    if [[ \$pairs -gt 0 ]]; then args+=(-1 "${r1}" -2 "${r2}"); fi
    if [[ \$singles -gt 0 ]]; then args+=(-s "${single}"); fi
    # Every FASTQ here is Phred+33 (samtools fastq / fastp). metaSPAdes otherwise guesses the
    # offset from the quality characters and aborts ("Failed to determine offset!") when a
    # small recruited set carries only characters valid in both encodings, e.g. all 'I' or
    # NovaSeq-binned all 'F'.
    metaspades.py "\${args[@]}" --phred-offset 33 -t ${task.cpus} -m ${task.memory.toGiga().intValue()} -o spades_out
    # metaSPAdes has no --min-contig-len, so apply MEGAHIT's floor here. Without this the
    # two assemblers are compared on different length distributions and assembler_agreement
    # reports discordance that is an artefact of the floor, not of the assemblies.
    awk -v minlen=${params.targeted_min_contig_len} '
        /^>/ { if (name != "" && length(seq) >= minlen) print name ORS seq; name=\$0; seq=""; next }
        { seq = seq \$0 }
        END { if (name != "" && length(seq) >= minlen) print name ORS seq }
    ' spades_out/contigs.fasta > "${sampleID}.targeted_metaspades.contigs.fa"
    """
}

process alignTargetedContigsToCanonicalReference {
    label 'sample_stage'
    tag "${sampleID}:${assembler}"
    label 'targeted_alignment'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/contigs/${assembler}" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.targeted_alignment_env}"
    input:
    tuple val(sampleID), val(assembler), path(contigs)
    path reference
    output:
    tuple val(sampleID), val(assembler), path(contigs), path("${sampleID}.${assembler}.vs_IHE3034.paf"), emit: aligned
    script:
    """
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'minimap2')}
    paf="${sampleID}.${assembler}.vs_IHE3034.paf"
    if [[ -s "${contigs}" ]]; then minimap2 -x asm10 -c -t ${task.cpus} "${reference}" "${contigs}" > "\$paf"; else : > "\$paf"; fi
    """
}

process summarizeTargetedPksEvidence {
    label 'sample_stage'
    tag "$sampleID"
    label 'targeted_alignment'
    publishDir { "${params.sample_dir}/${sampleID}/contigs/final_evidence" }, mode: 'copy', pattern: "*.tsv", saveAs: { fn -> fn - "${sampleID}." }
    // the validation plot belongs with the other figures, not buried in final_evidence/
    publishDir { "${params.sample_dir}/${sampleID}/figures" }, mode: 'copy', pattern: "*.svg", saveAs: { "contig_validation.svg" }
    conda "${params.targeted_alignment_env}"
    input:
    tuple val(sampleID), path(megahit_contigs), path(megahit_paf), path(metaspades_contigs), path(metaspades_paf), path(raw_coverage), path(read_evidence), path(recruitment_stats)
    output:
    tuple val(sampleID), path("${sampleID}.final_pks_evidence.tsv"), emit: evidence
    tuple val(sampleID), path("${sampleID}.raw_depth_megahit_metaspades_IHE3034.svg"), emit: plots
    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_annotation ${params.dep_digest?.pks_annotation}  pks_reference_fasta ${params.dep_digest?.pks_reference_fasta}  scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python')}
    # Ludmil, revised report finding 6: island-start/-end here are 0-based half-open
    # (PAF's own convention -- summarize_tumor_pks_contigs.py takes ts/te straight
    # from PAF fields with no adjustment), so the correct conversion from the
    # 1-based inclusive annotation is pks_start_1based-1 .. pks_end_1based, not the
    # bare pks_shift this used to pass, which silently dropped the first base.
    python "${params.scripts}/summarize_tumor_pks_contigs.py" --sample "${sampleID}" --megahit-contigs "${megahit_contigs}" --metaspades-contigs "${metaspades_contigs}" --megahit-paf "${megahit_paf}" --metaspades-paf "${metaspades_paf}" --raw-depth "${raw_coverage}" --read-evidence "${read_evidence}" --recruitment-stats "${recruitment_stats}" --gff "${params.pks_genome_annotation}" --contig "${params.pks_contig}" --region-start "${params.pks_plot_region_start}" --region-end "${params.pks_plot_region_end}" --island-start "${params.pks_start_1based.toString().toInteger() - 1}" --island-end "${params.pks_end_1based}" --min-aligned-bp "${params.tumor_contig_min_aligned_bp}" --output-tsv "${sampleID}.final_pks_evidence.tsv" --output-svg "${sampleID}.raw_depth_megahit_metaspades_IHE3034.svg"
    """
}

workflow targetedPksAssembly {
    take:
    reads
    profiles
    main:
    targetedPksRecruit(reads)
    targetedPksMegahit(targetedPksRecruit.out.reads)
    targetedPksSpades(targetedPksRecruit.out.reads)
    reference_ch = Channel.value(file(params.pks_reference_fasta, checkIfExists: true))
    megahit_ch = targetedPksMegahit.out.contigs.map { sampleID, contigs -> tuple(sampleID, 'megahit', contigs) }
    metaspades_ch = targetedPksSpades.out.contigs.map { sampleID, contigs -> tuple(sampleID, 'metaspades', contigs) }
    alignTargetedContigsToCanonicalReference(megahit_ch.mix(metaspades_ch), reference_ch)
    megahit_aligned_ch = alignTargetedContigsToCanonicalReference.out.aligned.filter { sampleID, assembler, contigs, paf -> assembler == 'megahit' }.map { sampleID, assembler, contigs, paf -> tuple(sampleID, contigs, paf) }
    metaspades_aligned_ch = alignTargetedContigsToCanonicalReference.out.aligned.filter { sampleID, assembler, contigs, paf -> assembler == 'metaspades' }.map { sampleID, assembler, contigs, paf -> tuple(sampleID, contigs, paf) }
    combined_ch = megahit_aligned_ch.join(metaspades_aligned_ch, by: 0).join(profiles, by: 0).join(targetedPksRecruit.out.stats, by: 0)
    summarizeTargetedPksEvidence(combined_ch)
    emit:
    evidence = summarizeTargetedPksEvidence.out.evidence
    plots = summarizeTargetedPksEvidence.out.plots
    alignments = alignTargetedContigsToCanonicalReference.out.aligned
}
