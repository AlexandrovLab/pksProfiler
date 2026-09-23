nextflow.enable.dsl = 2

process classifyPksReadEvidence {
    tag "$sampleID"
    label 'targeted_alignment'
    publishDir { "${params.sample_dir}/${sampleID}" }, mode: 'copy', saveAs: { fn -> fn - "${sampleID}." }
    conda "${params.pks_align_env}"
    input:
    tuple val(sampleID), path(raw_coverage), path(counts), path(bam), path(bai), path(qc)
    output:
    tuple val(sampleID), path(raw_coverage), path(counts), path("${sampleID}.read_evidence.tsv"), emit: evidence
    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    # Ludmil, revised report finding 6: pks_shift+1 dropped the true first base of
    # the island. pks_start_1based/pks_end_1based are already correct 1-based
    # inclusive bounds; --island-length must match the row count this now produces
    # (50768, not the old pks_island_len=50767), or classify_tumor_pks_evidence.py's
    # own row-count check raises rather than silently disagreeing.
    samtools depth -aa -s -Q 40 -r "${params.pks_contig}:${params.pks_start_1based}-${params.pks_end_1based}" "${bam}" > "${sampleID}.island.raw_depth.tsv"
    python "${params.scripts}/classify_tumor_pks_evidence.py" --sample "${sampleID}" --qc "${qc}" --counts "${counts}" --depth "${sampleID}.island.raw_depth.tsv" --island-length "${params.pks_island_len_1based}" --multi-gene-reads "${params.tumor_multi_gene_min_pks_reads}" --multi-gene-genes "${params.tumor_multi_gene_min_clb_genes}" --multi-gene-breadth "${params.tumor_multi_gene_min_breadth}" --broad-island-reads "${params.tumor_broad_island_min_pks_reads}" --broad-island-genes "${params.tumor_broad_island_min_clb_genes}" --broad-island-breadth "${params.tumor_broad_island_min_breadth}" --extensive-island-reads "${params.tumor_extensive_island_min_pks_reads}" --extensive-island-genes "${params.tumor_extensive_island_min_clb_genes}" --extensive-island-breadth "${params.tumor_extensive_island_min_breadth}" --output "${sampleID}.read_evidence.tsv"
    """
}

process targetedPksRecruit {
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
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_recruit_index ${params.dep_digest?.pks_recruit_index}  scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    bowtie2 --very-sensitive-local -k 1 --no-unal --threads ${task.cpus} -x "${params.pks_recruit_index}" -U "${reads}" -S "${sampleID}.pks_recruitment.sam"
    python "${params.scripts}/recover_recruited_mates.py" --sam "${sampleID}.pks_recruitment.sam" --fastq "${reads}" --r1 "${sampleID}.pks.R1.fastq.gz" --r2 "${sampleID}.pks.R2.fastq.gz" --single "${sampleID}.pks.single.fastq.gz" --stats "${sampleID}.pks_recruitment.tsv" --min-aligned-bases "${params.pks_recruit_min_aligned}" --min-identity "${params.pks_recruit_min_identity}"
    """
}

process targetedPksMegahit {
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
    pairs=\$(gzip -dc "${r1}" | awk 'END {print int(NR/4)}')
    singles=\$(gzip -dc "${single}" | awk 'END {print int(NR/4)}')
    if [[ \$pairs -eq 0 && \$singles -eq 0 ]]; then : > "${sampleID}.targeted_metaspades.contigs.fa"; exit 0; fi
    args=()
    if [[ \$pairs -gt 0 ]]; then args+=(-1 "${r1}" -2 "${r2}"); fi
    if [[ \$singles -gt 0 ]]; then args+=(-s "${single}"); fi
    metaspades.py "\${args[@]}" -t ${task.cpus} -m ${task.memory.toGiga().intValue()} -o spades_out
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
    paf="${sampleID}.${assembler}.vs_IHE3034.paf"
    if [[ -s "${contigs}" ]]; then minimap2 -x asm10 -c -t ${task.cpus} "${reference}" "${contigs}" > "\$paf"; else : > "\$paf"; fi
    """
}

process summarizeTargetedPksEvidence {
    tag "$sampleID"
    label 'targeted_alignment'
    publishDir { "${params.sample_dir}/${sampleID}/contigs/final_evidence" }, mode: 'copy', pattern: "*.tsv", saveAs: { fn -> fn - "${sampleID}." }
    // the validation plot belongs with the other figures, not buried in final_evidence/
    publishDir { "${params.sample_dir}/${sampleID}/figures" }, mode: 'copy', pattern: "*.svg", saveAs: { "contig_validation.svg" }
    conda "${params.targeted_alignment_env}"
    input:
    tuple val(sampleID), path(megahit_contigs), path(megahit_paf), path(metaspades_contigs), path(metaspades_paf), path(raw_coverage), path(read_evidence)
    output:
    tuple val(sampleID), path("${sampleID}.final_pks_evidence.tsv"), emit: evidence
    tuple val(sampleID), path("${sampleID}.raw_depth_megahit_metaspades_IHE3034.svg"), emit: plots
    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_annotation ${params.dep_digest?.pks_annotation}  pks_reference_fasta ${params.dep_digest?.pks_reference_fasta}  scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    # Ludmil, revised report finding 6: island-start/-end here are 0-based half-open
    # (PAF's own convention -- summarize_tumor_pks_contigs.py takes ts/te straight
    # from PAF fields with no adjustment), so the correct conversion from the
    # 1-based inclusive annotation is pks_start_1based-1 .. pks_end_1based, not the
    # bare pks_shift this used to pass, which silently dropped the first base.
    python "${params.scripts}/summarize_tumor_pks_contigs.py" --sample "${sampleID}" --megahit-contigs "${megahit_contigs}" --metaspades-contigs "${metaspades_contigs}" --megahit-paf "${megahit_paf}" --metaspades-paf "${metaspades_paf}" --raw-depth "${raw_coverage}" --read-evidence "${read_evidence}" --gff "${params.pks_genome_annotation}" --contig "${params.pks_contig}" --region-start "${params.pks_plot_region_start}" --region-end "${params.pks_plot_region_end}" --island-start "${params.pks_start_1based.toString().toInteger() - 1}" --island-end "${params.pks_end_1based}" --min-aligned-bp "${params.tumor_contig_min_aligned_bp}" --output-tsv "${sampleID}.final_pks_evidence.tsv" --output-svg "${sampleID}.raw_depth_megahit_metaspades_IHE3034.svg"
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
    combined_ch = megahit_aligned_ch.join(metaspades_aligned_ch, by: 0).join(profiles, by: 0)
    summarizeTargetedPksEvidence(combined_ch)
    emit:
    evidence = summarizeTargetedPksEvidence.out.evidence
    plots = summarizeTargetedPksEvidence.out.plots
    alignments = alignTargetedContigsToCanonicalReference.out.aligned
}
