nextflow.enable.dsl=2

include { strainTyping as magStrainTyping } from './pks_typing.nf'

// ─── Mate pairs, rebuilt after host depletion ────────────────────────────────
// Ludmil, 26 Sep mode memo section 5: the MAG lane already took -1/-2 but received one
// host-depleted FASTQ, so every metagenome assembled as single-end. Depletion stays
// read-independent (compatible with --frag=no); the pairs are reconstructed afterwards,
// from the mate suffixes, into R1/R2 plus singletons. Every read is kept: pairs +
// singletons == the host-depleted stream, so assembly still sees all of it.

process pairDepletedReads {
    label 'sample_stage'
    label 'mag_binning'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/assembly" },
        mode: 'copy', pattern: "*.mate_pairing.tsv", saveAs: { "mate_pairing.tsv" }
    conda "${params.samtools_env}"

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${sampleID}.pairs.R1.fastq.gz"), path("${sampleID}.pairs.R2.fastq.gz"),
          path("${sampleID}.pairs.single.fastq.gz"), emit: reads
    tuple val(sampleID), path("${sampleID}.mate_pairing.tsv"), emit: stats

    script:
    def reads_list = reads instanceof List ? reads : [reads]
    if (reads_list.size() != 1) {
        error "Mate pairing for ${sampleID} expects the one host-depleted FASTQ; received ${reads_list.size()} files"
    }
    // external sort buffer: half the task's memory, at least 1 GB (integer GB; Groovy `/`
    // would give a BigDecimal)
    def sort_gb = Math.max(1, (task.memory ? task.memory.toGiga().intValue() : 2).intdiv(2))
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # scripts ${params.dep_digest?.scripts}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python')}
    python "${params.scripts}/recover_recruited_mates.py" --all-reads \\
        --fastq "${reads_list[0]}" \\
        --r1 "${sampleID}.pairs.R1.fastq.gz" --r2 "${sampleID}.pairs.R2.fastq.gz" \\
        --single "${sampleID}.pairs.single.fastq.gz" \\
        --stats "${sampleID}.mate_pairing.tsv" \\
        --threads ${task.cpus} --sort-memory ${sort_gb}G \\
        --tmpdir .
    """
}

// ─── Assembly ─────────────────────────────────────────────────────────────────

process megahitAssemble {
    label 'sample_stage'
    label 'mag_assembly'
    scratch true
    publishDir { "${params.sample_dir}/${sampleID}/genomes/assembly" },
        mode: 'copy', pattern: "*.contig*", enabled: params.save_intermediates.toString().toBoolean()
    conda "${projectDir}/conda_envs/megahit_env.yml"

    input:
    tuple val(sampleID), path(r1), path(r2), path(single)

    output:
    tuple val(sampleID), path("${sampleID}.contigs.fa"), emit: contigs
    tuple val(sampleID), path("${sampleID}.contig_count.txt"), emit: contig_count

    script:
    // Pairs as -1/-2, singletons (unsuffixed, orphaned or single-end reads) as -r; each
    // only when it has reads, since MEGAHIT refuses an empty input file.
    """
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'megahit')}
    pairs=\$(gzip -dc "${r1}" | awk 'END {print int(NR/4)}')
    singles=\$(gzip -dc "${single}" | awk 'END {print int(NR/4)}')
    args=()
    if [[ \$pairs -gt 0 ]]; then args+=(-1 "${r1}" -2 "${r2}"); fi
    if [[ \$singles -gt 0 ]]; then args+=(-r "${single}"); fi
    if [[ \${#args[@]} -eq 0 ]]; then
        : > ${sampleID}.contigs.fa
    else
        megahit "\${args[@]}" -t ${task.cpus} -o megahit_out
        cp megahit_out/final.contigs.fa ${sampleID}.contigs.fa
    fi
    grep -c '^>' ${sampleID}.contigs.fa > ${sampleID}.contig_count.txt || \\
        printf '0\\n' > ${sampleID}.contig_count.txt
    """
}

// ─── Read-to-contig alignment ─────────────────────────────────────────────────

process alignToContigs {
    label 'sample_stage'
    label 'mag_binning'
    scratch true
    conda "${projectDir}/conda_envs/minimap2_env.yml"

    input:
    tuple val(sampleID), path(r1), path(r2), path(single), path(contigs)

    output:
    tuple val(sampleID), path("${sampleID}.contigs.sorted.bam"),
          path("${sampleID}.contigs.sorted.bam.bai"), emit: bam

    script:
    // Pairs are aligned as pairs (two files: minimap2 -ax sr pairs them) and singletons
    // on their own; the two sorted BAMs are merged, so MetaBAT2's depths count every read.
    """
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'minimap2 samtools')}
    pairs=\$(gzip -dc "${r1}" | awk 'END {print int(NR/4)}')
    singles=\$(gzip -dc "${single}" | awk 'END {print int(NR/4)}')
    parts=()
    if [[ \$pairs -gt 0 ]]; then
        minimap2 -ax sr -t ${task.cpus} ${contigs} "${r1}" "${r2}" \\
            | samtools view -F 4 -b \\
            | samtools sort -@ ${task.cpus} -o paired.bam
        parts+=(paired.bam)
    fi
    if [[ \$singles -gt 0 ]]; then
        minimap2 -ax sr -t ${task.cpus} ${contigs} "${single}" \\
            | samtools view -F 4 -b \\
            | samtools sort -@ ${task.cpus} -o single.bam
        parts+=(single.bam)
    fi
    if [[ \${#parts[@]} -eq 1 ]]; then
        mv "\${parts[0]}" ${sampleID}.contigs.sorted.bam
    else
        samtools merge -@ ${task.cpus} -o ${sampleID}.contigs.sorted.bam "\${parts[@]}"
    fi
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'jgi_summarize_bam_contig_depths:metabat2')}
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'metabat2')}
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'checkm2')}
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'gtdbtk')}
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

// ─── U1 follow-up: stand-in for a sample with zero real bins ──────────────────
// checkm2Predict/gtdbtkClassify both key off metabat2Bin.out.bins (real bins only,
// optional: true), so a sample with zero real bins never gets an entry in either
// channel. Confirmed via a real (non-mocked) metagenome --enable_mags execution: a
// small synthetic assembly under MetaBAT2's 200 kb --minClsSize floor pooled every
// contig into the unbinned pool (see metabat2Bin above) and had ZERO real bins --
// exactly the "abundance is low or binning is underpowered" case U1's own comment
// names as the reason the unbinned pool exists. Before this fix, summary_input_ch's
// strict join on checkm2/gtdbtk silently dropped that sample before its
// mag_summary_locus_evidence_ch entry (a real, positive unbinned call) ever got a
// chance to keep it via that join's own remainder: true -- pks_mag_summary.tsv was
// never written for it at all, the exact blind spot U1 was written to close, just
// one join upstream of where U1 closed it. This supplies an empty but correctly
// shaped pair so magSummaryTable still runs for exactly the samples this affects.
// Deliberately no publish directive here: unlike checkm2Predict/gtdbtkClassify, no
// real tool ran, and publishing a fake checkm2/gtdbtk directory would misrepresent
// what happened.
process stubBinlessMagQuality {
    label 'sample_stage'
    label 'process_low'
    tag "$sampleID"

    input:
    val sampleID

    output:
    tuple val(sampleID), path("stub.checkm2_quality_report.tsv"), emit: checkm2
    tuple val(sampleID), path("stub.gtdbtk.bac120.summary.tsv"), emit: gtdbtk

    script:
    """
    set -euo pipefail
    printf "Name\tCompleteness\tContamination\tGenome_Size\tContig_N50\n" > stub.checkm2_quality_report.tsv
    printf "user_genome\tclassification\n" > stub.gtdbtk.bac120.summary.tsv
    """
}

// ─── rc2: stand-in for a CheckM2/GTDB-Tk task that failed and was ignored ─────
// Mayo metagenome run (v1.0.0-rc1, 2026-10): several gtdbtkClassify tasks exited 1 and,
// under the sample_stage errorStrategy ('ignore' by default), simply never emitted a
// summary. The sample still had bins, tblouts and per-bin locus evidence, so
// summary_input_ch's strict checkm2/gtdbtk join dropped it and the remainder:true joins
// after it, and community_input_ch's, produced a left-orphaned 3-element row
// [sampleID, null, [bin.*.locus_evidence.tsv, ...]] that the next 6-parameter closure
// could not destructure -- Nextflow stopped the whole run ~37 h in. This supplies, for
// exactly the tool(s) that produced nothing, a correctly shaped stand-in, so every bin is
// still reported instead of the sample vanishing or the run crashing:
//   GTDB-Tk: one row per bin, classification "gtdbtk_failed" (taxonomy says so per bin);
//   CheckM2: the report header plus a Status column and NO rows, so every quality field
//            reads NA (unmeasured), never 0. build_mag_summary.py reads that Status header
//            to set mag_qc_status=checkm2_failed.
// Neither stand-in changes how scripts/ parse these tables, so the scripts digest (and with
// it -resume of every earlier task) is untouched by this fix. The failed task itself is
// still listed under RUN_REPORT "WHAT FAILED" (the trace records it as FAILED). Not
// published: no real tool output exists.
process stubFailedMagQuality {
    label 'sample_stage'
    label 'process_low'
    tag "$sampleID"

    input:
    tuple val(sampleID), val(bin_ids), val(checkm2_failed), val(gtdbtk_failed)

    output:
    tuple val(sampleID), path("failed.checkm2_quality_report.tsv"), emit: checkm2, optional: true
    tuple val(sampleID), path("failed.gtdbtk.bac120.summary.tsv"), emit: gtdbtk, optional: true

    script:
    def ids = bin_ids.join(' ')
    """
    set -euo pipefail
    if [[ "${checkm2_failed}" == "true" ]]; then
        # header only: no bin has a measurement. "Status" marks the stand-in.
        printf "Name\\tCompleteness\\tContamination\\tGenome_Size\\tContig_N50\\tStatus\\n" > failed.checkm2_quality_report.tsv
    fi
    if [[ "${gtdbtk_failed}" == "true" ]]; then
        printf "user_genome\\tclassification\\n" > failed.gtdbtk.bac120.summary.tsv
        for b in ${ids}; do
            printf "%s\\tgtdbtk_failed\\n" "\$b" >> failed.gtdbtk.bac120.summary.tsv
        done
    fi
    """
}

// rc2: guard for the row a remainder:true join emits when the LEFT side has no entry for
// a key -- [key, null, <right values>], which is shorter than the row the next closure
// destructures. Such a row is dropped with a warning (the sample's upstream failure is
// already in RUN_REPORT "WHAT FAILED"); anything else passes through unchanged.
def completeJoinRow(row, int width, String where) {
    if (row instanceof List && row.size() == width && row[1] != null) {
        return true
    }
    def key = row instanceof List ? row[0] : row
    log.warn "pksMAG: ${where} skipped for ${key}: an upstream per-sample task produced no output " +
             "(failed and ignored; see RUN_REPORT WHAT FAILED)"
    return false
}

def noRealBinsSentinel() {
    "__stubBinlessMagQuality_no_real_bins_anywhere_sentinel__"
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'prokka')}
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'hmmsearch:hmmer')}
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'minimap2')}
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python')}
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'genomad')}
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
    // P0-5b: with their sampleID, so the cohort report and master summary can stage
    // them as declared inputs rather than find them in the published tree.
    tuple val(sampleID), path("${sampleID}.community_prophage_inventory.tsv"), emit: inventory
    tuple val(sampleID), path("${sampleID}.pks_community_interactions.tsv"), emit: interactions
    tuple val(sampleID), path("${sampleID}.pks_island_mobility.tsv"), emit: mobility

    script:
    """
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python')}
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
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python')}
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
    // P0-5b: with its sampleID -- see communityProphageSummary.
    tuple val(sampleID), path("${sampleID}.pks_mag_summary.tsv"), emit: summary

    script:
    """
    # F15 dependency digests -- a change here must invalidate this task; lib/Provenance.groovy
    # pks_reference_fasta ${params.dep_digest?.pks_reference_fasta}
    set -euo pipefail
    # F16: the version and path of each tool this task runs; lib/ToolVersions.groovy
    ${ToolVersions.record(task, params, 'python')}
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
    reads  // tuple val(sampleID), path(reads): the complete host-depleted stream

    main:

    // 0. Mate pairs, rebuilt from the host-depleted stream (26 Sep memo section 5)
    pairDepletedReads(reads)
    paired_reads = pairDepletedReads.out.reads

    // 1. Assembly
    megahitAssemble(paired_reads)

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
    reads_contigs_ch = paired_reads.join(nonempty_contigs_ch, by: 0)
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


    // U1 follow-up (see stubBinlessMagQuality above): give every sample that has a
    // mag_summary_locus_evidence_ch entry but no real bins (and so no checkm2/gtdbtk
    // output of its own) a stand-in pair, so the strict join two lines down keeps it
    // instead of silently dropping it.
    //
    // The sentinel mix-in below is not decoration: `join(..., remainder: true)` against
    // a right-hand channel that is COMPLETELY empty (every sample in the run has zero
    // real bins -- exactly this fixture) does not null-pad the way it does when the
    // right-hand channel has some entries just not for this key. Confirmed by standalone
    // reproduction and a live run: with a truly empty right-hand side, Nextflow's join
    // silently degenerates and re-emits the left-hand channel's items unchanged (bare
    // sampleID strings, not [sampleID, null] tuples), and the very next 2-parameter
    // closure then throws "Invalid method invocation `call` with arguments: <sampleID>
    // (java.lang.String)" -- a real, live crash, not a hypothetical. A right-hand
    // channel that can never be structurally empty (it always carries at least the
    // sentinel, which no real sampleID can ever equal) keeps the join's normal
    // null-padding behaviour in every case. The same degeneration reaches all the way
    // through summary_input_ch's own join chain below when a run has zero real bins
    // AND zero pks+ bins anywhere (this fixture's exact case: tblouts_per_sample_ch and
    // contexts_per_sample_ch are then ALSO completely empty) -- every channel that
    // chain joins is guarded with the same sentinel key for the same reason, and the
    // sentinel row that survives the joins as a result is filtered back out immediately
    // after, before magSummaryTable (or anything else) ever sees it.
    def NO_REAL_BINS_SENTINEL = noRealBinsSentinel()
    samples_with_bins_ch = checkm2Predict.out.report
        .map { sampleID, _report -> tuple(sampleID, true) }
        .mix(Channel.value(tuple(NO_REAL_BINS_SENTINEL, true)))
    binless_samples_needing_summary_ch = mag_summary_locus_evidence_ch
        .map { sampleID, _evidence -> sampleID }
        .join(samples_with_bins_ch, by: 0, remainder: true)
        .filter { _sampleID, has_bins -> has_bins == null }
        .map { sampleID, _has_bins -> sampleID }
    stubBinlessMagQuality(binless_samples_needing_summary_ch)

    // rc2: samples WITH bins whose checkm2Predict/gtdbtkClassify failed and was ignored
    // get a per-bin "<tool>_failed" stand-in -- see stubFailedMagQuality.
    magQualityFallback(metabat2Bin.out.bins, checkm2Predict.out.report, gtdbtkClassify.out.summary)

    checkm2_report_ch = checkm2Predict.out.report.mix(stubBinlessMagQuality.out.checkm2)
        .mix(magQualityFallback.out.checkm2_stub)
        .mix(Channel.value(tuple(NO_REAL_BINS_SENTINEL, "SENTINEL")))
    gtdbtk_summary_ch = gtdbtkClassify.out.summary.mix(stubBinlessMagQuality.out.gtdbtk)
        .mix(magQualityFallback.out.gtdbtk_stub)
        .mix(Channel.value(tuple(NO_REAL_BINS_SENTINEL, "SENTINEL")))
    tblouts_per_sample_ch = tblouts_per_sample_ch
        .mix(Channel.value(tuple(NO_REAL_BINS_SENTINEL, [])))
    contexts_per_sample_ch = contexts_per_sample_ch
        .mix(Channel.value(tuple(NO_REAL_BINS_SENTINEL, [])))
    mag_summary_locus_evidence_ch = mag_summary_locus_evidence_ch
        .mix(Channel.value(tuple(NO_REAL_BINS_SENTINEL, [])))

    // remainder: true keeps all samples even when contexts_per_sample_ch or
    // tblouts_per_sample_ch have no entry (zero pks+ bins, or -- since the fix just
    // above -- zero bins at all); null → [] so the script receives an empty file list.
    // locus_evidence here is mag_summary_locus_evidence_ch (bins + the U1 unbinned
    // pool), not locus_evidence_per_sample_ch -- build_mag_summary.py's locus_dir
    // needs the unbinned call too, to report it as its own row.
    //
    // rc2: the joins live in magSummaryInput / magCommunityInput (below) so that
    // tests/check_mag_channel_regression.sh can drive the exact same channel logic with
    // a failed GTDB-Tk task, which crashed the rc1 Mayo run.
    magSummaryInput(checkm2_report_ch, gtdbtk_summary_ch, tblouts_per_sample_ch,
                    contexts_per_sample_ch, mag_summary_locus_evidence_ch)
    summary_input_ch = magSummaryInput.out.rows

    magSummaryTable(summary_input_ch)

    // Real GTDB-Tk output plus the rc2 gtdbtk_failed stand-in (not the binless stub:
    // a sample with no real bins has no genomes to build community tables from).
    magCommunityInput(gtdbtkClassify.out.summary.mix(magQualityFallback.out.gtdbtk_stub),
                      gffs_per_sample_ch, tblouts_per_sample_ch,
                      virus_summaries_per_sample_ch, locus_evidence_per_sample_ch)
    community_input_ch = magCommunityInput.out.rows

    communityProphageSummary(community_input_ch)

    emit:
    // Consumed by main.nf's optional_lane_gate, which stages them into masterSummary and
    // cohortReport as declared inputs (Ludmil, 2026-09-25 report P0-5b). Each is
    // tuple(sampleID, file): the scripts no longer find these by globbing the published
    // tree, so the sample identifier has to travel with the file. These used to be
    // bare paths, the shape a pure ordering gate needed; a tuple mixed into that gate
    // failed masterSummary's path-typed input with "Not a valid path value:
    // '<sampleID>'" (caught live on the 2026-09-24 validate1 run). main.nf now tags
    // every one with its kind before the mix, so there is no bare-path gate left to break.
    mag_summary       = magSummaryTable.out.summary
    community_summary = communityProphageSummary.out.inventory
    community_interactions = communityProphageSummary.out.interactions
    island_mobility   = communityProphageSummary.out.mobility
    strain_summary    = strain_typing_summary_ch
    // U1: the unbinned pool's own canonical-locus call, one file per sample that had
    // anything left unbinned. Already folded into mag_summary itself
    // (pks_mag_summary.tsv carries it as a distinct row -- see build_mag_summary.py);
    // emitted again here, named, so it is not only discoverable by reading a table.
    // Neither cohort-level script reads it, so main.nf does not stage it.
    unbinned_locus_evidence = unbinned_locus_evidence_ch
    // 26 Sep memo section 5 / A7: pairs, singletons and mate-tag diagnostics per sample,
    // published as genomes/assembly/mate_pairing.tsv.
    mate_pairing = pairDepletedReads.out.stats
}


// ─── rc2: channel logic factored out of pksMAG so it can be tested on its own ───
// tests/fixtures/mag_channel_regression.nf drives these with fabricated channels; the
// processes above are not involved except stubFailedMagQuality (no conda, no tools).

// Which samples with real bins got no CheckM2 report / no GTDB-Tk summary because the
// task failed (and was ignored)? Emits a stand-in for exactly those.
workflow magQualityFallback {
    take:
    bins      // tuple(sampleID, bin FASTA | [bin FASTAs]) -- metabat2Bin.out.bins
    checkm2   // tuple(sampleID, quality_report.tsv)     -- checkm2Predict.out.report
    gtdbtk    // tuple(sampleID, bac120.summary.tsv)     -- gtdbtkClassify.out.summary

    main:
    // The same never-empty right-hand side trick as pksMAG's NO_REAL_BINS_SENTINEL: a
    // remainder:true join against a completely empty channel does not null-pad.
    def SENTINEL = noRealBinsSentinel()
    bin_ids_ch = bins.map { sampleID, bin_fas ->
        def fas = bin_fas instanceof List ? bin_fas : [bin_fas]
        tuple(sampleID, fas.collect { new File(it.toString()).name.replaceFirst(/\.fa$/, '') }.sort())
    }
    has_checkm2_ch = checkm2.map { sampleID, _report -> tuple(sampleID, true) }
        .mix(Channel.value(tuple(SENTINEL, true)))
    has_gtdbtk_ch = gtdbtk.map { sampleID, _summary -> tuple(sampleID, true) }
        .mix(Channel.value(tuple(SENTINEL, true)))

    // One-parameter closures throughout: a left-orphaned remainder row is shorter than
    // the 4-tuple, and must be filtered, not destructured.
    failed_ch = bin_ids_ch
        .join(has_checkm2_ch, by: 0, remainder: true)
        .join(has_gtdbtk_ch, by: 0, remainder: true)
        .filter { row -> row instanceof List && row.size() == 4 && row[0] != SENTINEL && row[1] != null }
        .filter { row -> row[2] == null || row[3] == null }
        .map { row ->
            def failed = [row[2] == null ? 'checkm2' : null, row[3] == null ? 'gtdbtk' : null].findAll()
            log.warn "pksMAG: ${failed.join(' and ')} produced no output for ${row[0]} (task failed and " +
                     "was ignored); its ${row[1].size()} bin(s) are reported as " +
                     "${failed.collect { it + '_failed' }.join('/')} in pks_mag_summary.tsv"
            tuple(row[0], row[1], row[2] == null, row[3] == null)
        }
    stubFailedMagQuality(failed_ch)

    emit:
    checkm2_stub = stubFailedMagQuality.out.checkm2
    gtdbtk_stub  = stubFailedMagQuality.out.gtdbtk
    failed       = failed_ch
}

// magSummaryTable's input. checkm2/gtdbtk already carry every sample (real, binless stub,
// failed stand-in, sentinel); the remainder joins keep samples with no pks+ bins.
workflow magSummaryInput {
    take:
    checkm2
    gtdbtk
    tblouts
    contexts
    locus_evidence

    main:
    def SENTINEL = noRealBinsSentinel()
    rows = checkm2
        .join(gtdbtk, by: 0)
        .join(tblouts, by: 0, remainder: true)
        .join(contexts, by: 0, remainder: true)
        .join(locus_evidence, by: 0, remainder: true)
        .filter { row -> !(row instanceof List && row[0] == SENTINEL) }
        .filter { row -> completeJoinRow(row, 6, 'magSummaryTable') }
        .map { row -> tuple(row[0], row[1], row[2], row[3] ?: [], row[4] ?: [], row[5] ?: []) }

    emit:
    rows
}

// communityProphageSummary's input: genomes with taxonomy, annotation and provirus calls.
workflow magCommunityInput {
    take:
    gtdbtk
    gffs
    tblouts
    virus_summaries
    locus_evidence

    main:
    // Sentinel on the remainder side so it is never completely empty (see pksMAG).
    def SENTINEL = noRealBinsSentinel()
    rows = gtdbtk
        .join(gffs, by: 0)
        .join(tblouts, by: 0)
        .join(virus_summaries, by: 0)
        .join(locus_evidence.mix(Channel.value(tuple(SENTINEL, []))), by: 0, remainder: true)
        .filter { row -> !(row instanceof List && row[0] == SENTINEL) }
        .filter { row -> completeJoinRow(row, 6, 'communityProphageSummary') }
        .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5] ?: []) }

    emit:
    rows
}
