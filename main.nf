nextflow.enable.dsl = 2

// ---------------- Parameters ----------------
params.help = false
// What to do when one sample fails a per-sample step. 'ignore' lets the remaining
// samples finish and the cohort tables build; the failed sample is recorded in
// pks.qc.summary.tsv as no_qc_produced or incomplete. 'finish' stops the run, which
// also means no cohort table is written at all.
params.sample_failure_strategy = 'ignore'
params.sample = null

params.input_data_type = "auto"        // auto | bam | cram | fastq
params.cram_reference  = null          // optional; used and validated when supplied
params.pks_taxa = false   // set true to run krakenuniq/bracken on pks-island reads
params.save_intermediates = false // publish extracted/filtered/host-depleted FASTQs

params.profiling_method = "bowtie2" // bowtie2 | hmm | both
params.hmm_evalue       = 1e-10
params.hmm_chunking = false
// v0.0.2: population-scale nucleotide profiles (exact_v1 build). The v0.0.1 benchmark
// model, ref/hmm/clb_all_dna.hmm, is retained in the tree and selectable via --hmm_model.
params.hmm_model        = "${projectDir}/ref/hmm/clb_population_dna_exact_v1.hmm"
params.bracken_read_length = null // Must match a read length supported by the selected Bracken database.

// ---------------- v0.0.2: sample-type routing ----------------
params.sample_type = "auto"   // auto | metagenome | tumor_wgs | tumor_wes | tumor_rna

// ---------------- v0.0.2: optional taxonomic prefilter + clb homology rescue ----------------
// off      -> host-depleted reads go straight to profiling (v0.0.1 behaviour, the default)
// balanced -> krakenuniq narrows to primary_taxid, then DIAMOND blastx rescues clb-homologous
//             reads krakenuniq placed outside it, and the two sets are merged.
// auto resolves to balanced for metagenomes when --kraken_db is available, otherwise off.
// Metagenomes carry far more non-target background, so the narrowing pays there; tumour
// WGS keeps v0.0.1 behaviour unless asked. Assembly is never prefiltered in either case.
params.prefilter_mode = "auto"    // auto | off | balanced
params.primary_taxid  = 91347     // Enterobacterales
params.diamond_rescue = true
params.diamond_evalue = 1e-5
params.diamond_min_aa = 25
params.diamond_min_query_coverage = 0.5
// 844 unique exact_v1 training proteins, all 19 genes — replaces the 19 single reference
// proteins so the rescue can recover divergent reads. ref/clb_reference_proteins.faa (19 seqs)
// is retained and selectable via --clb_protein_fasta.
params.clb_protein_fasta = "${projectDir}/ref/clb_reference_proteins_exact_v1.faa"

// Community-level taxonomy: Bracken over the whole non-host krakenuniq report.
// Distinct from --pks_taxa, which is Bracken over reads extracted from the island itself.
params.pks_community_taxa = false

// ---------------- v0.0.2: tumour pks evidence tiers ----------------
// Frozen three-tier system from tier_threshold_development_20260909.
// All three criteria must hold; classification proceeds from the highest tier down.
// Breadth is the fraction of the canonical island covered at >=1x.
params.tumor_multi_gene_min_pks_reads       = 5
params.tumor_multi_gene_min_clb_genes       = 3
params.tumor_multi_gene_min_breadth         = 0.01
params.tumor_broad_island_min_pks_reads     = 30
params.tumor_broad_island_min_clb_genes     = 8
params.tumor_broad_island_min_breadth       = 0.075
params.tumor_extensive_island_min_pks_reads = 100
params.tumor_extensive_island_min_clb_genes = 10
params.tumor_extensive_island_min_breadth   = 0.15

// Tiers that qualify a tumour sample for contig analysis.
// Restricted to the two higher tiers: the depth-control sweep showed targeted assembly
// is underpowered at multi_gene-level depth, so multi_gene classifies but does not assemble.
// Widen with --tumor_contig_tiers multi_gene,broad_island,extensive_island if needed.
params.tumor_contig_tiers = "broad_island,extensive_island"

// ---------------- v0.0.2: contig analysis ----------------
params.tumor_targeted_assembly   = true    // recruit + MEGAHIT/metaSPAdes on eligible tumours
params.tumor_full_contig_context = false   // prophage + neighbouring-gene context; requires --genomad_db
params.tumor_enable_mags         = false   // genome-resolved binning on tumour reads
params.enable_mags               = false   // MAG reconstruction; metagenome sample_type only

params.pks_recruit_min_aligned  = 60
params.pks_recruit_min_identity = 0.90

// Shared contig-length floor for both targeted assemblers. MEGAHIT takes it directly;
// metaSPAdes has no equivalent option and is filtered to the same floor afterwards, so
// assembler_agreement compares like with like.
params.targeted_min_contig_len = 300

// Minimum aligned bases for a contig alignment to count as island evidence.
// Was hardcoded at 500 in summarize_tumor_pks_contigs.py, which discarded 10-77% of
// in-island alignments and 100% of them on both broad_island samples in the
// v0.0.2_functional_test_20260910 tumour arm, reporting those as no_contig_support.
// Set to the assembler floor: below it MEGAHIT cannot contribute by construction, so a
// lower value makes metaSPAdes look artificially divergent.
params.tumor_contig_min_aligned_bp = 200   // T1: total island bases a contig must carry to count as supporting; below clbR (213 bp)

// Legacy read-level eligibility screen inside tumorWGS. Undefined upstream (latent NPE);
// aligned here with the lowest positive tier. Non-host floor disabled by default -- the
// tier classifier is the authoritative gate.
params.tumor_min_clb_genes    = 3
params.tumor_min_pks_reads    = 5
params.tumor_min_nonhost_reads = 0

// ---------------- v0.0.2: MAG / prophage ----------------
params.gtdbtk_db               = null
params.checkm2_db              = null
params.genomad_db              = null
params.clb_protein_hmm         = "${projectDir}/ref/hmm/clb_population_protein_exact_v1.hmm"
params.hmm_protein_evalue      = 1e-5
// M3: one definition of a pks-positive bin, used by the status count, by the choice
// of which bins get genomic context, and by community producer selection. Counted as
// distinct clb genes, which is what build_mag_summary.py reports.
params.mag_min_clb_genes = 3

// M1: clbB/clbC/clbH/clbI/clbJ/clbK/clbN/clbO are multi-domain NRPS/PKS megasynthases
// whose domains are shared with unrelated secondary-metabolite clusters, so an E-value
// cut alone let a non-Enterobacterales bin clear mag_min_clb_genes on megasynthase
// hits by domain homology, not island carriage (v0.0.2_functional_test_20260910,
// ERR525841: 7 of 8 bins called pks-positive this way, including two Bifidobacterium
// bins). Rather than requiring specific genes by name, each bin's own assembly is
// aligned to the canonical IHE3034 locus (params.pks_reference_fasta,
// alignMagBinToCanonicalReference/magBinLocusEvidence in Modules/pks_mag.nf) and
// scored the same way read-level evidence already is: genes an alignment actually
// covers, and breadth of the island those alignments span. Tier names and threshold
// values below are deliberately the same as classify_tumor_pks_evidence.py's read-level
// tiers, so "positive" means the same thing whether the evidence is reads or a bin.
// hmmsearchClb (mag_min_clb_genes, above) stays a cheap HMM pre-filter for which bins
// get genomic-context extraction; it does not decide positivity.
params.mag_locus_min_identity = 0.90
params.mag_locus_min_mapq     = 20
params.mag_locus_multi_gene_min_genes         = 3
params.mag_locus_multi_gene_min_breadth       = .01
params.mag_locus_broad_island_min_genes       = 8
params.mag_locus_broad_island_min_breadth     = .075
params.mag_locus_extensive_island_min_genes   = 10
params.mag_locus_extensive_island_min_breadth = .15

// T3: a geNomad provirus must clear both floors to be reported. geNomad called 134
// "viral contigs" on AA-3850 whose top hits were 369 bp with one gene and one
// hallmark -- a count that tracked assembly fragmentation, not biology.
params.provirus_min_length_bp = 3000
params.provirus_min_hallmarks = 2
// Flanking window for nearby-gene context around each clb hit, both lanes.
params.context_window_bp = 50000

// ---------------- v0.0.2: strain typing ----------------
// MLST on assembled units (MAG bins, tumour contig sets), joined to clonal complex and
// phylogroup via the 45,761-genome reference panel. Reads are never typed: a handful of
// pks reads cannot support a strain call.
params.enable_strain_typing   = true
params.mlst_scheme            = "ecoli_achtman_4"
params.typing_min_assembly_bp = 500000   // below this, MLST cannot recover 7 loci
params.st_phylogroup_lookup   = "${projectDir}/ref/typing/st_phylogroup.tsv"

// profile taxa that map to the pks island:
params.pks_shift = 2193827          // island start in E. coli genome coords
params.pks_island_len = 50767       // island length (0..50767 in your file)
params.pks_contig = 'NC_017628.1'   // contig name in BAM
params.pks_plot_region_start = 2183826   // plotting window around the island
params.pks_plot_region_end   = 2254594

// Ludmil, revised report finding 6: the annotation spans 2,193,827-2,244,594
// inclusive (50,768 bp), but pks_shift/pks_island_len above got a different
// +/-1 adjustment in every module that turned them into a region -- correct in
// pks_taxa.nf (pks_shift-1 for a 0-based BED start), off by one bp at the start
// everywhere else (pksProfiler_align.nf, pks_targeted.nf's two call sites, and
// this session's own pks_mag.nf locus alignment), silently dropping the true
// first base of the island. These two values are the single authoritative
// definition; every correctness-sensitive call site below now derives its own
// coordinate representation (1-based samtools region, or 0-based BED/PAF start)
// from them directly. pks_shift/pks_island_len are UNCHANGED, and still used by
// plotting.nf's coverage-plot window, where a 1 bp offset has no numeric effect.
params.pks_start_1based    = 2193827   // first base of the annotated island (1-based, inclusive)
params.pks_end_1based      = 2244594   // last base of the annotated island (1-based, inclusive)
params.pks_island_len_1based = params.pks_end_1based - params.pks_start_1based + 1   // 50768, not 50767

// Output directories
params.outdir = "${launchDir}/results"

// v0.0.2 output layout. Three top-level trees and no fourth:
//   by_sample/<sample>/   everything one sample produced
//   cohort/               everything that spans samples
//   runs/                 nextflow's own report/timeline/trace
// Replaces pks_per_sample/ + pks_summary/ + mags/ + community_context/ +
// strain_typing/ + prefilter/, which split each sample's artefacts across six trees
// with inconsistent nesting (community_context carried an extra arm level on the
// tumour side only, strain_typing did not, pks_per_sample was flat) and put figures
// in three separate places.
params.sample_dir = "${params.outdir}/by_sample"
params.cohort_dir = "${params.outdir}/cohort"
params.runs_dir   = "${params.outdir}/runs"

params.pks_counts_dir = "${params.cohort_dir}/gene_counts"
params.pks_taxonomy_dir = "${params.cohort_dir}/taxonomy"
params.pks_qc_dir = "${params.cohort_dir}/qc"

// Databases and refs [CHANGE THIS]
params.hg38_db      = null
params.t2t_phix_db  = null
params.pangenome_db = null
// Optional. fastp detects adapters by itself -- by overlap analysis for paired reads
// and a built-in list otherwise -- and matching every read against the 234 sequences
// in ref/known_adapters.fna was a measurable share of the Hartwig extraction runtime
// for no change in what survived. Pass --adapters ref/known_adapters.fna to restore it.
params.adapters     = null
params.kraken_db= null


// PKS + E. coli annotation
params.pks_genome            = "${projectDir}/indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic"
// v0.0.2: targeted recruitment reuses the canonical IHE3034 index and its FASTA.
// Island + 5 kb flank, NOT the whole 5.1 Mb chromosome. Recruiting against the full
// genome made targeted assembly genome-proportional: ~1% of recruited fragments were
// island reads across the v0.0.2 tumour arm, matching the island's 0.99% share of the
// reference. The 5 kb flank carries the mobility architecture the context analysis
// looks for -- integrase at -280 bp, tRNA-Asn at -1713 bp, IS3 at +146 and +1014 --
// which a zero-flank index would truncate away. Rebuild with
// reporting/tumor_edits_2026-09-16/04-build-island-recruit-index.sh.
params.pks_recruit_index     = "${projectDir}/indices/pks_island/pks_island"
params.pks_reference_fasta   = "${params.pks_genome}.fna"
params.pks_genome_annotation = "${projectDir}/ref/annotations/IHE3034.clbA-clbS.gff"
params.pks_cytoband          = "${projectDir}/indices/GCF_000025745.1/genomic_pks.txt"

// Envs
params.samtools_env = "${projectDir}/conda_envs/samtools_env.yml"
params.fastp_env = "${projectDir}/conda_envs/fastp_env.yml"
params.minimap2_env = "${projectDir}/conda_envs/minimap2_env.yml"
params.pks_align_env = "${projectDir}/conda_envs/pks_align_env.yml"
params.pks_hmm_env = "${projectDir}/conda_envs/pks_hmm_env.yml"
params.krakenuniq_bracken_env = "${projectDir}/conda_envs/krakenUniq_bracken_env.yml"
params.prefilter_env          = "${projectDir}/conda_envs/prefilter_env.yml"
params.mlst_env               = "${projectDir}/conda_envs/mlst_env.yml"
params.targeted_recruit_env   = "${projectDir}/conda_envs/targeted_recruit_env.yml"
params.targeted_spades_env    = "${projectDir}/conda_envs/targeted_spades_env.yml"
params.targeted_alignment_env = "${projectDir}/conda_envs/targeted_alignment_env.yml"
params.megahit_env            = "${projectDir}/conda_envs/megahit_env.yml"
params.metabat2_env           = "${projectDir}/conda_envs/metabat2_env.yml"
params.checkm2_env            = "${projectDir}/conda_envs/checkm2_env.yml"
params.gtdbtk_env             = "${projectDir}/conda_envs/gtdbtk_env.yml"
params.genomad_env            = "${projectDir}/conda_envs/genomad_env.yml"
params.prokka_env             = "${projectDir}/conda_envs/prokka_env.yml"
// F09: the cheap denominator comes from the alignment index and counts records,
// including secondary and supplementary. Set this to count primary records exactly,
// which costs a second full decode of every input.
params.exact_input_counts = false

// F14: aligned bases a read must share with a clb gene to be assigned to it. One
// value, used by the alignment QC and the taxonomy lane alike.
params.min_gene_overlap_bp = 1
params.conda_cache_dir = null   // pin conda envs outside the work dir; also read by the run report
params.scripts = "${projectDir}/scripts"

// ---------------- Dependency digests (F15) ----------------
// Nextflow hashes a task from its declared inputs and its rendered script. A reference,
// model or helper that appears only as a path string inside the script is neither, so
// rebuilding an index in place or editing a helper leaves every task hash unchanged and
// `-resume` reuses work that no longer corresponds to the current code or references.
// These digests are interpolated into the affected task scripts, which are hashed, so a
// change to any of them invalidates exactly the tasks that depend on it.
//
// Provenance.code() reads content, for the things small enough to read on every launch.
// Provenance.data() reads name, size and mtime of the top level, for indexes and
// databases that are far too large -- see lib/Provenance.groovy for what that misses.
params.dep_digest = [
    scripts             : Provenance.code(params.scripts),
    hmm_model           : Provenance.code(params.hmm_model),
    clb_protein_hmm     : Provenance.code(params.clb_protein_hmm),
    pks_annotation      : Provenance.code(params.pks_genome_annotation),
    pks_cytoband        : Provenance.code(params.pks_cytoband),
    adapters            : Provenance.code(params.adapters),
    hg38_db             : Provenance.data(params.hg38_db),
    t2t_phix_db         : Provenance.data(params.t2t_phix_db),
    pangenome_db        : Provenance.data(params.pangenome_db),
    pks_genome          : Provenance.data(params.pks_genome),
    pks_recruit_index   : Provenance.data(params.pks_recruit_index),
    pks_reference_fasta : Provenance.data(params.pks_reference_fasta),
    cram_reference      : Provenance.data(params.cram_reference),
    kraken_db           : Provenance.data(params.kraken_db),
    genomad_db          : Provenance.data(params.genomad_db),
    gtdbtk_db           : Provenance.data(params.gtdbtk_db),
    checkm2_db          : Provenance.data(params.checkm2_db),
]

// ---------------- Modules ----------------
include { extractReads } from './Modules/extract_reads.nf'
include { filterReads } from './Modules/filter_reads.nf'
include { mapReads } from './Modules/map_reads.nf'
include { pksProfiler_align as pksProfilerAlign } from './Modules/pksProfiler_align.nf'
include { pksProfiler_hmm as pksProfilerHMM } from './Modules/pksProfiler_hmm.nf'
include { plotPKS; masterTableAlign; masterTableHMM; masterQCSummary; cohortReport } from './Modules/plotting.nf'
include { extractPksIslandReads; Bracken; process_bracken as combinePKSTaxa; combineClbTaxonomySupport } from './Modules/pks_taxa.nf'
include { plotBrackenTaxa as plotPKSTaxa } from './Modules/plot_bracken_taxa.nf'
include { plotBrackenTaxa as plotCommunityTaxa } from './Modules/plot_bracken_taxa.nf'
include { classifyPksReadEvidence; targetedPksAssembly } from './Modules/pks_targeted.nf'
include { pksMAG; tumorWGS } from './Modules/pks_mag.nf'
// Separate alias so the tumour path cannot collide with the metagenome one:
// Nextflow forbids invoking a component twice, and validation alone is a weaker
// guarantee than making the second invocation a different component.
include { pksMAG as tumorPksMAG } from './Modules/pks_mag.nf'
include { krakenPrefilter; buildClbDiamondDb; diamondRescue; mergePksCandidates; sampleBracken } from './Modules/pks_prefilter.nf'

// An empty channel, named for what it gates. A local `def` inside the branch that
// fills it would not be visible where cohortReport is invoked, and Nextflow's parser
// allows no statements at the top level.
def classifyPksReadEvidence_gate() { return channel.empty() }

// ---------------- Workflow ----------------
workflow {

    // ---------- Help ----------
    // First thing the workflow does: `--help` must work without --sample, without
    // databases, and without emitting the prefilter resolution line.
    if (params.help.toString().toBoolean()) {
        log.info """
        pksProfiler ${workflow.manifest.version} -- ${workflow.manifest.description}

        Usage:
          nextflow run ${workflow.manifest.name} --sample SHEET.csv \\
              --hg38_db INDEX --t2t_phix_db INDEX [options]

        Required:
          --sample            sample sheet; needs a unique `patient` column plus your
                              files (bam | cram + cram_reference | fastq1/fastq2)
          --hg38_db           minimap2 index or FASTA for human depletion
          --t2t_phix_db       minimap2 index or FASTA for T2T + PhiX depletion

        Common:
          --input_data_type   auto | bam | cram | fastq            (default: auto)
          --sample_type       auto | metagenome | tumor_wgs |
                              tumor_wes | tumor_rna                (default: auto)
          --profiling_method  bowtie2 | hmm | both                 (default: bowtie2)
          --outdir            results directory                    (default: ./results)
          --save_intermediates  publish extracted/filtered/host-depleted FASTQs

        Optional stages (each off unless stated):
          --enable_mags            recover draft genomes from metagenomes
          --tumor_targeted_assembly  reassemble the island from tumour reads
          --enable_strain_typing   MLST on assembled output        (default: on)
          --prefilter_mode         auto | off | balanced           (default: auto)
          --pks_taxa               taxonomy of island reads; needs --kraken_db

        Evidence tiers are descriptive strata, not a validated presence/absence test.
        See README.md and docs/ for the full parameter list and interpretation limits.
        """.stripIndent()
        exit 0
    }

    // ---------- v0.0.2: boolean flag coercion ----------
    // `--flag false` on the command line arrives as the String "false", which is truthy
    // in Groovy. Every v0.0.2 boolean is resolved once here; only these locals are tested.
    // toString().toBoolean() maps Boolean true/false and the Strings "true"/"false" alike.
    def enable_mags_b               = params.enable_mags.toString().toBoolean()
    def tumor_enable_mags_b         = params.tumor_enable_mags.toString().toBoolean()
    def tumor_full_contig_context_b = params.tumor_full_contig_context.toString().toBoolean()
    def tumor_targeted_assembly_b   = params.tumor_targeted_assembly.toString().toBoolean()
    def diamond_rescue_b            = params.diamond_rescue.toString().toBoolean()
    def pks_community_taxa_b        = params.pks_community_taxa.toString().toBoolean()
    def enable_strain_typing_b      = params.enable_strain_typing.toString().toBoolean()

    // ---------- Equivalence record (F16) ----------
    // Written now, before any task runs, so a run that dies still leaves a record of
    // what it was attempting; finished off in onComplete with the outcome and the
    // output checksums. One file per session: with -resume a results tree is usually
    // the product of several sessions, and a record describing only the last one would
    // misdescribe the cohort.
    def run_record_path = RunRecord.path(params.runs_dir, workflow.sessionId)
    RunRecord.write(run_record_path, RunRecord.start([
        session_id  : workflow.sessionId.toString(),
        code        : RunRecord.gitIdentity(workflow.projectDir, workflow.commitId, workflow.revision) + [
            repository      : workflow.manifest.name,
            manifest_version: workflow.manifest.version,
            script_id       : workflow.scriptId?.toString(),
            project_dir     : workflow.projectDir.toString(),
        ],
        invocation  : [
            command_line    : workflow.commandLine,
            launch_dir      : workflow.launchDir.toString(),
            work_dir        : workflow.workDir.toString(),
            profile         : workflow.profile,
            outdir          : params.outdir.toString(),
            sample_sheet    : [path: params.sample?.toString(), digest: Provenance.data(params.sample)],
            stages_requested: RunRecord.stages(params),
            params          : params.findAll { key, _value -> key != 'dep_digest' }
                                    .collectEntries { key, value -> [(key): value?.toString()] },
        ],
        dependencies: params.dep_digest,
        environment : [
            nextflow_version : workflow.nextflow.version?.toString(),
            nextflow_build   : workflow.nextflow.build?.toString(),
            java_version     : System.getProperty('java.version'),
            os               : "${System.getProperty('os.name')} ${System.getProperty('os.version')}",
            user_name        : workflow.userName,
            // The launch host. Tasks run wherever the executor put them, but on a
            // cluster the head node is often what explains a failure.
            hostname         : java.net.InetAddress.getLocalHost().getHostName(),
            conda_envs_digest: Provenance.code("${projectDir}/conda_envs"),
        ],
        notes       : [
            "stages_requested is what the parameters asked for; runs/trace.txt lists the tasks that actually executed.",
            "Solved tool versions are not captured: conda_envs_digest fingerprints the pinned specifications, not the environments conda resolved from them.",
            "Reference and database digests follow lib/Provenance.groovy -- content for code and models, name/size/mtime for large indexes and databases.",
        ],
    ]))
    log.info "Provenance record: ${run_record_path}"

    // ---------- v0.0.2: resolve the prefilter mode ----------
    def prefilter_mode = params.prefilter_mode.toString()
    if (prefilter_mode == "auto") {
        if (params.sample_type == "metagenome" && params.kraken_db) {
            prefilter_mode = "balanced"
        }
        else {
            prefilter_mode = "off"
        }
        log.info "prefilter_mode auto -> ${prefilter_mode} (sample_type=${params.sample_type})"
    }

    // ---------- Preflight contract (F20) ----------
    // One pass over everything the run needs, before a single task is submitted, and
    // every problem at once. Stage prerequisites used to check only that a flag was
    // set: `--hg38_db /typo/human.mmi` satisfied every check and failed in the first
    // mapReads task, after the whole cohort had been extracted.
    def preflight_status = Preflight.run("${params.scripts}/preflight.py", [
        sample                   : params.sample,
        input_data_type          : params.input_data_type,
        sample_type              : params.sample_type,
        cram_reference           : params.cram_reference,
        profiling_method         : params.profiling_method,
        prefilter_mode           : prefilter_mode,
        hg38_db                  : params.hg38_db,
        t2t_phix_db              : params.t2t_phix_db,
        pangenome_db             : params.pangenome_db,
        pks_genome               : params.pks_genome,
        pks_recruit_index        : params.pks_recruit_index,
        adapters                 : params.adapters,
        pks_reference_fasta      : params.pks_reference_fasta,
        pks_genome_annotation    : params.pks_genome_annotation,
        hmm_model                : params.hmm_model,
        clb_protein_hmm          : params.clb_protein_hmm,
        clb_protein_fasta        : params.clb_protein_fasta,
        st_phylogroup_lookup     : params.st_phylogroup_lookup,
        kraken_db                : params.kraken_db,
        bracken_read_length      : params.bracken_read_length,
        genomad_db               : params.genomad_db,
        gtdbtk_db                : params.gtdbtk_db,
        checkm2_db               : params.checkm2_db,
        pks_taxa                 : params.pks_taxa.toString().toBoolean(),
        pks_community_taxa       : pks_community_taxa_b,
        enable_mags              : enable_mags_b,
        tumor_enable_mags        : tumor_enable_mags_b,
        tumor_full_contig_context: tumor_full_contig_context_b,
        tumor_targeted_assembly  : tumor_targeted_assembly_b,
        enable_strain_typing     : enable_strain_typing_b,
    ], log)

    if (preflight_status != 0) {
        error "Preflight checks failed. Nothing has run; fix the problems above and relaunch."
    }


	// ---------- Required input validation ----------
    if (!params.sample) {
        exit 1, "Missing required parameter: --sample"
    }

    if (!params.hg38_db) {
        exit 1, "Missing required parameter: --hg38_db"
    }

    if (!params.t2t_phix_db) {
        exit 1, "Missing required parameter: --t2t_phix_db"
    }

    if (params.pks_taxa && !params.kraken_db) {
        exit 1, "Taxonomic profiling requires: --kraken_db"
    }
	if (params.pks_taxa && !params.bracken_read_length) {
	    exit 1, "Taxonomic profiling requires: --bracken_read_length"
	}
	if (params.pks_taxa && (
        !(params.bracken_read_length.toString() ==~ /^[0-9]+$/) ||
        params.bracken_read_length.toString().toInteger() <= 0
    )) {
	    exit 1, "--bracken_read_length must be a positive integer"
	}

    def hmm_evalue_text = params.hmm_evalue.toString()

    if (!(
        hmm_evalue_text ==~
        /^(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$/
    )) {
        exit 1, "--hmm_evalue must be a positive number, for example 1e-10"
    }

    if (hmm_evalue_text.toDouble() <= 0) {
        exit 1, "--hmm_evalue must be greater than zero"
    }

    // ---------- STEP 1: Inputs + filtering ----------
    if (!(params.input_data_type in ["auto", "bam", "cram", "fastq"])) {
        exit 1, "Unknown --input_data_type: ${params.input_data_type}. Supported: auto, bam, cram, fastq"
    }

    def alignment_input = params.input_data_type in ["auto", "bam", "cram"]
    def required_sample_columns = alignment_input ? ["patient"] : ["patient", "fastq1"]

    def sample_sheet = channel
        .fromPath(params.sample, checkIfExists: true)
        .splitCsv(header: true)
        .collect()
        .flatMap { rows ->
            if (!rows) {
                error "Sample sheet contains no samples: ${params.sample}"
            }

            def observed_columns = rows[0].keySet()
            if (alignment_input && !("alignment" in observed_columns) && !("bam" in observed_columns)) {
                error "Alignment sample sheets require an alignment column (or legacy bam column)"
            }
            def missing_columns = required_sample_columns.findAll { column ->
                !(column in observed_columns)
            }

            if (missing_columns) {
                error "Sample sheet is missing required column(s) for ${params.input_data_type} input: ${missing_columns.join(', ')}"
            }

            def observed_ids = new HashSet()
            def duplicate_ids = new TreeSet()

            rows.eachWithIndex { row, index ->
                def row_number = index + 2
                def sample_id = row.patient?.toString()?.trim()

                if (!sample_id) {
                    error "Sample sheet row ${row_number} has an empty patient value"
                }

                if (!(sample_id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/)) {
                    error "Invalid patient value '${sample_id}' on row ${row_number}. Use only letters, numbers, periods, underscores, and hyphens; the first character must be alphanumeric."
                }

                if (!observed_ids.add(sample_id)) {
                    duplicate_ids.add(sample_id)
                }

                required_sample_columns
                    .findAll { column -> column != "patient" }
                    .each { column ->
                        if (!row[column]?.toString()?.trim()) {
                            error "Sample sheet row ${row_number} has an empty ${column} value"
                        }
                    }

                if (alignment_input) {
                    def alignment = row.alignment?.toString()?.trim() ?: row.bam?.toString()?.trim()
                    if (!alignment) {
                        error "Sample sheet row ${row_number} has an empty alignment value"
                    }
                }
            }

            if (duplicate_ids) {
                error "Sample identifiers must be unique. Duplicate patient value(s): ${duplicate_ids.join(', ')}"
            }

            rows
        }

    // Every sample that entered the run, captured before the sheet is remapped per
    // input type. masterQCSummary uses this to give a row to samples that produced no
    // QC at all, so a failure is visible in the table instead of silently absent.
    def EXPECTED_SAMPLE_IDS = sample_sheet
        .map { row -> row.patient.toString().trim() }
        .collectFile(name: 'expected_samples.txt', newLine: true, sort: true)

    def QC_FRAGMENTS = channel.empty()

    if (alignment_input) {

			// Prefer patient,alignment; retain patient,bam for compatibility.
			sample_sheet = sample_sheet.map { row ->
			    def alignment = row.alignment?.toString()?.trim() ?: row.bam?.toString()?.trim()
			    // Trimmed, exactly as validation and EXPECTED_SAMPLE_IDS trim it. Passing
			    // the raw value made a sheet with stray whitespace validate as one
			    // identifier and run as another, so the sample read as missing in QC.
			    tuple(row.patient.toString().trim(), file(alignment, checkIfExists: true))
			}

        EXTRACT_OUT = extractReads(sample_sheet)

        EXTRACT_OUT.reads
            .map { sampleID, reads -> tuple(sampleID, [reads]) }
            .set { READS_TO_FILTER }

        QC_FRAGMENTS = QC_FRAGMENTS.mix(
            EXTRACT_OUT.qc.map { _sampleID, qc_file -> qc_file }
        )

    } else if (params.input_data_type == "fastq") {

        def sample_sheet_fastq = sample_sheet
            .map { row ->
                def fastq_files = [
                    file(row.fastq1, checkIfExists: true)
                ]

                def fastq2 = row.fastq2?.toString()?.trim()
                if (fastq2) {
                    fastq_files << file(fastq2, checkIfExists: true)
                }

                tuple(
                    row.patient.toString().trim(),
                    fastq_files
                )
            }

        sample_sheet_fastq.set { READS_TO_FILTER }

    }

    FILTER_OUT = filterReads(READS_TO_FILTER)

    FILTER_OUT.reads
        .set { FILTERED_UNMAPPED_READS }

    QC_FRAGMENTS = QC_FRAGMENTS.mix(
        FILTER_OUT.qc.map { _sampleID, qc_file -> qc_file }
    )

	// ---------- STEP 1b: Host read depletion ----------
	MAP_OUT = mapReads(FILTERED_UNMAPPED_READS)

	MAP_OUT.reads
	    .set { MAPPED_READS }

    QC_FRAGMENTS = QC_FRAGMENTS.mix(
        MAP_OUT.qc.map { _sampleID, qc_file -> qc_file }
    )

    // ---------- v0.0.2: preserve the complete host-depleted stream for MAG assembly ----------
    // Branched before profiling so genome-resolved analysis sees every host-depleted read.
    MAPPED_READS
        .filter { _sampleID, _reads -> params.sample_type == "metagenome" && enable_mags_b }
        .set { MAG_ASSEMBLY_READS }

    if (enable_mags_b) {
        pksMAG(MAG_ASSEMBLY_READS)
    }

	// ---------- v0.0.2: prefilter validation ----------
    def valid_prefilter_modes = ["auto", "off", "balanced"]
    if (!(params.prefilter_mode in valid_prefilter_modes)) {
        exit 1, "Unknown --prefilter_mode: ${params.prefilter_mode}. Supported: auto, off, balanced"
    }
    if (prefilter_mode == "balanced") {
        if (!params.kraken_db) { exit 1, "Balanced prefiltering requires: --kraken_db" }
        if (!diamond_rescue_b) {
            exit 1, "--prefilter_mode balanced requires --diamond_rescue true; without the rescue, clb reads outside taxid ${params.primary_taxid} are discarded"
        }
        if (!file(params.clb_protein_fasta).exists()) {
            exit 1, "clb protein FASTA not found: ${params.clb_protein_fasta}"
        }
    }
    if (!(params.primary_taxid.toString() ==~ /^[1-9][0-9]*$/)) {
        exit 1, "--primary_taxid must be a positive integer (got ${params.primary_taxid})"
    }
    if (params.diamond_evalue.toString().toDouble() <= 0 ||
        params.diamond_min_aa.toString().toInteger() <= 0 ||
        params.diamond_min_query_coverage.toString().toDouble() <= 0 ||
        params.diamond_min_query_coverage.toString().toDouble() > 1) {
        exit 1, "--diamond_evalue and --diamond_min_aa must be positive, and --diamond_min_query_coverage must be in (0,1]"
    }
    if (pks_community_taxa_b && prefilter_mode != "balanced") {
        exit 1, "--pks_community_taxa requires --prefilter_mode balanced so Bracken receives the complete non-host krakenuniq report"
    }
    if (pks_community_taxa_b && !params.bracken_read_length) {
        exit 1, "--pks_community_taxa requires: --bracken_read_length"
    }

	// ---------- v0.0.2: sample-type and contig-analysis validation ----------
    def valid_sample_types = ["auto", "metagenome", "tumor_wgs", "tumor_wes", "tumor_rna"]
    if (!(params.sample_type in valid_sample_types)) {
        exit 1, "Unknown --sample_type: ${params.sample_type}. Supported: ${valid_sample_types.join(', ')}"
    }
    if (enable_mags_b && params.sample_type != "metagenome") {
        exit 1, "--enable_mags requires --sample_type metagenome; tumour assembly is selected with --sample_type tumor_wgs"
    }
    // Both flags invoke pksMAG, and Nextflow forbids invoking a process twice. Tying each
    // to its own sample_type makes them mutually exclusive by construction.
    if (tumor_enable_mags_b && params.sample_type != "tumor_wgs") {
        exit 1, "--tumor_enable_mags requires --sample_type tumor_wgs; metagenome binning is selected with --enable_mags"
    }
    if (params.sample_type == "tumor_wgs" && params.profiling_method == "hmm") {
        exit 1, "--sample_type tumor_wgs requires --profiling_method bowtie2 or both for read-level clb screening"
    }
    if (params.sample_type == "tumor_wgs" && tumor_targeted_assembly_b) {
        if (!file(params.pks_reference_fasta).exists()) {
            exit 1, "Canonical PKS reference not found: ${params.pks_reference_fasta}"
        }
        if (!file("${params.pks_recruit_index}.1.bt2").exists()) {
            exit 1, "PKS recruitment index not found: ${params.pks_recruit_index}.*.bt2"
        }
    }
    if (tumor_full_contig_context_b && !params.genomad_db) {
        exit 1, "--tumor_full_contig_context requires --genomad_db"
    }
    // pksMAG runs checkm2, gtdbtk and genomad unconditionally, so all three are required.
    [enable_mags: enable_mags_b, tumor_enable_mags: tumor_enable_mags_b].each { flag, on ->
        if (on) {
            if (!params.gtdbtk_db)  { exit 1, "--${flag} requires --gtdbtk_db" }
            if (!params.checkm2_db) { exit 1, "--${flag} requires --checkm2_db" }
            if (!params.genomad_db) { exit 1, "--${flag} requires --genomad_db for prophage detection" }
        }
    }
    if ((enable_mags_b || tumor_enable_mags_b || tumor_full_contig_context_b) && !file(params.clb_protein_hmm).exists()) {
        exit 1, "Protein HMM not found: ${params.clb_protein_hmm}"
    }
    // nhmmscan reads the pressed index, not the plain .hmm; hmmsearch needs only the .hmm.
    if (params.profiling_method in ["hmm", "both"]) {
        if (!file(params.hmm_model).exists()) {
            exit 1, "Nucleotide HMM not found: ${params.hmm_model}"
        }
        if (!file("${params.hmm_model}.h3i").exists()) {
            exit 1, "Nucleotide HMM is not pressed: ${params.hmm_model}.h3i missing. Run: hmmpress ${params.hmm_model}"
        }
    }

    if (enable_strain_typing_b) {
        if (!file(params.st_phylogroup_lookup).exists()) {
            exit 1, "ST/phylogroup lookup not found: ${params.st_phylogroup_lookup}"
        }
        if (!(params.typing_min_assembly_bp.toString() ==~ /^[0-9]+$/)) {
            exit 1, "--typing_min_assembly_bp must be a non-negative integer (got ${params.typing_min_assembly_bp})"
        }
    }

    if (!(params.context_window_bp.toString() ==~ /^[0-9]+$/) || params.context_window_bp.toString().toInteger() <= 0) {
        exit 1, "--context_window_bp must be a positive integer (got ${params.context_window_bp})"
    }

    def tier_order = ["multi_gene", "broad_island", "extensive_island"]
    def contig_tiers = params.tumor_contig_tiers.toString().split(",").collect { it.trim() }.findAll { it }
    if (contig_tiers.isEmpty()) {
        exit 1, "--tumor_contig_tiers must name at least one tier of: ${tier_order.join(', ')}"
    }
    contig_tiers.each { tier ->
        if (!(tier in tier_order)) {
            exit 1, "Unknown tier in --tumor_contig_tiers: ${tier}. Supported: ${tier_order.join(', ')}"
        }
    }

    [tumor_multi_gene_min_breadth: params.tumor_multi_gene_min_breadth,
     tumor_broad_island_min_breadth: params.tumor_broad_island_min_breadth,
     tumor_extensive_island_min_breadth: params.tumor_extensive_island_min_breadth].each { name, value ->
        def parsed = value.toString().toDouble()
        if (parsed < 0 || parsed > 1) {
            exit 1, "--${name} must be a fraction between 0 and 1 (got ${value})"
        }
    }
    [tumor_multi_gene_min_clb_genes: params.tumor_multi_gene_min_clb_genes,
     tumor_multi_gene_min_pks_reads: params.tumor_multi_gene_min_pks_reads,
     tumor_broad_island_min_clb_genes: params.tumor_broad_island_min_clb_genes,
     tumor_broad_island_min_pks_reads: params.tumor_broad_island_min_pks_reads,
     tumor_extensive_island_min_clb_genes: params.tumor_extensive_island_min_clb_genes,
     tumor_extensive_island_min_pks_reads: params.tumor_extensive_island_min_pks_reads].each { name, value ->
        if (!(value.toString() ==~ /^[0-9]+$/) || value.toString().toInteger() < 0) {
            exit 1, "--${name} must be a non-negative integer (got ${value})"
        }
    }
    if (params.tumor_multi_gene_min_clb_genes.toString().toInteger() > 19 ||
        params.tumor_broad_island_min_clb_genes.toString().toInteger() > 19 ||
        params.tumor_extensive_island_min_clb_genes.toString().toInteger() > 19) {
        exit 1, "clb gene thresholds cannot exceed 19 (clbA-clbS)"
    }
    // Tiers must be monotonically non-decreasing, or the highest-tier-down cascade is unreachable.
    if (params.tumor_broad_island_min_pks_reads.toString().toInteger() < params.tumor_multi_gene_min_pks_reads.toString().toInteger() ||
        params.tumor_extensive_island_min_pks_reads.toString().toInteger() < params.tumor_broad_island_min_pks_reads.toString().toInteger() ||
        params.tumor_broad_island_min_clb_genes.toString().toInteger() < params.tumor_multi_gene_min_clb_genes.toString().toInteger() ||
        params.tumor_extensive_island_min_clb_genes.toString().toInteger() < params.tumor_broad_island_min_clb_genes.toString().toInteger() ||
        params.tumor_broad_island_min_breadth.toString().toDouble() < params.tumor_multi_gene_min_breadth.toString().toDouble() ||
        params.tumor_extensive_island_min_breadth.toString().toDouble() < params.tumor_broad_island_min_breadth.toString().toDouble()) {
        exit 1, "Tier thresholds must be non-decreasing from multi_gene to broad_island to extensive_island"
    }

	// ---------- STEP 2: Profiling ----------
    def valid_methods = ["bowtie2", "hmm", "both"]

    if (!(params.profiling_method in valid_methods)) {
        exit 1, "Unknown --profiling_method: ${params.profiling_method}. Supported: bowtie2, hmm, both"
    }

    if (params.pks_taxa && params.profiling_method == "hmm") {
        exit 1, "--pks_taxa requires alignment profiling. Use --profiling_method bowtie2 or both."
    }

    def do_align = params.profiling_method in ["bowtie2", "both"]
    def do_hmm   = params.profiling_method in ["hmm", "both"]

    // The cohort report's ordering gate is declared here, not inside `if (do_align)`,
    // because cohortReport reads it unconditionally. Assigned only in the alignment
    // branch, `--profiling_method hmm` left it undefined and the run died while the
    // workflow was still being built: "No such variable: cohort_report_gate". The
    // alignment branch overwrites it with the real evidence channel; read-level evidence
    // comes from the alignment lane, so an HMM-only run genuinely has no gate to wait on
    // and the empty channel is the correct value rather than a placeholder.
    cohort_report_gate = classifyPksReadEvidence_gate()

    // ---------- v0.0.2: optional prefilter ahead of profiling ----------
    // MAG assembly and the tumour re-join both draw on MAPPED_READS, so neither is
    // affected by this narrowing -- only what gets profiled changes.
    if (prefilter_mode == "balanced") {
        KRAKEN_PREFILTER_OUT = krakenPrefilter(MAPPED_READS)
        QC_FRAGMENTS = QC_FRAGMENTS.mix(KRAKEN_PREFILTER_OUT.qc.map { _sampleID, qc_file -> qc_file })

        CLB_DIAMOND_DB = buildClbDiamondDb(file(params.clb_protein_fasta, checkIfExists: true))
        DIAMOND_RESCUE_OUT = diamondRescue(KRAKEN_PREFILTER_OUT.non_target, CLB_DIAMOND_DB)
        QC_FRAGMENTS = QC_FRAGMENTS.mix(DIAMOND_RESCUE_OUT.qc.map { _sampleID, qc_file -> qc_file })

        CANDIDATE_OUT = mergePksCandidates(
            KRAKEN_PREFILTER_OUT.primary
                .join(DIAMOND_RESCUE_OUT.reads)
                .map { sampleID, primary, rescued -> tuple(sampleID, primary, rescued) }
        )
        CANDIDATE_OUT.reads.set { PROFILING_READS }
        QC_FRAGMENTS = QC_FRAGMENTS.mix(CANDIDATE_OUT.qc.map { _sampleID, qc_file -> qc_file })

        if (pks_community_taxa_b) {
            plotCommunityTaxa(sampleBracken(KRAKEN_PREFILTER_OUT.taxonomy).reports)
        }
        // Profiling now runs on a narrowed stream while assembly still uses the complete
        // one, so read-level and assembly-level clb evidence for the same sample are
        // computed on different inputs and may legitimately disagree.
        if (enable_mags_b || tumor_targeted_assembly_b || tumor_full_contig_context_b) {
            log.warn "prefilter_mode=balanced narrows profiling only; assembly uses the complete host-depleted stream, so read-level and assembly-level clb evidence are not computed on the same reads"
        }
    } else {
        MAPPED_READS.set { PROFILING_READS }
    }

    if (do_align) {
        ALIGN_OUT = pksProfilerAlign(PROFILING_READS)

        ALIGN_OUT.profile
            .set { PKS_ALIGN_OUT }

        QC_FRAGMENTS = QC_FRAGMENTS.mix(
            ALIGN_OUT.qc.map { _sampleID, qc_file -> qc_file }
        )

        // ---------- v0.0.2: read-level pks evidence tier ----------
        // Computed for tumour and metagenome alike. For metagenomes it is an annotation
        // only: the tiers were derived from TCGA tumour breadth distributions and are not
        // calibrated for metagenome depth, and gating genome-resolved analysis on a
        // reference-based nucleotide tier would filter out the divergent carriers that
        // lane exists to find. The contig-count gate remains the metagenome gate.
        // The cohort report reads the published tree, so it must not start before the
        // per-sample evidence is there. Depending only on the QC summary, it could run
        // while read_evidence.tsv was still unpublished and report every sample as
        // not_classified. Empty when the lane that writes it did not run.
        if (params.sample_type in ["tumor_wgs", "metagenome"]) {
            read_evidence_input_ch = ALIGN_OUT.profile
                .join(ALIGN_OUT.qc, by: 0)
                .map { sampleID, rawCoverage, bedgraph, counts, bam, bai, qc ->
                    tuple(sampleID, rawCoverage, counts, bam, bai, qc)
                }
            classifyPksReadEvidence(read_evidence_input_ch)
            cohort_report_gate = classifyPksReadEvidence.out.evidence
                .map { _sampleID, _rawCoverage, _counts, evidenceFile -> evidenceFile }
        }

        // ---------- v0.0.2: tumour contig analysis on formally positive samples ----------
        if (params.sample_type == "tumor_wgs") {
            tumor_eligible_profiles_ch = classifyPksReadEvidence.out.evidence
                .filter { sampleID, rawCoverage, counts, evidenceFile ->
                    def lines = evidenceFile.text.readLines().findAll { it?.trim() }
                    if (lines.size() < 2) {
                        exit 1, "Empty evidence file for ${sampleID}: ${evidenceFile}"
                    }
                    def header = lines[0].split("\\t", -1) as List
                    def tier_index = header.indexOf("read_evidence")
                    if (tier_index < 0) {
                        exit 1, "Evidence file for ${sampleID} has no read_evidence column: ${evidenceFile}"
                    }
                    def tier = lines[1].split("\\t", -1)[tier_index]
                    tier in contig_tiers
                }

            tumor_eligible_ids_ch = tumor_eligible_profiles_ch
                .map { sampleID, rawCoverage, counts, evidenceFile -> tuple(sampleID) }
            tumor_assembly_reads_ch = MAPPED_READS
                .join(tumor_eligible_ids_ch, by: 0)
                .map { sampleID, reads -> tuple(sampleID, reads) }

            if (tumor_targeted_assembly_b) {
                targeted_profiles_ch = tumor_eligible_profiles_ch
                    .map { sampleID, rawCoverage, counts, evidenceFile -> tuple(sampleID, rawCoverage, evidenceFile) }
                targetedPksAssembly(tumor_assembly_reads_ch, targeted_profiles_ch)
                // Ludmil, revised report finding 5: the cohort report reads
                // contigs/final_evidence/final_pks_evidence.tsv directly from the
                // published tree. Assembly is the slowest lane in the pipeline, so
                // without this the report could start, and publish, before it lands.
                cohort_report_gate = cohort_report_gate.mix(
                    targetedPksAssembly.out.evidence
                        .map { _sampleID, evidenceFile -> evidenceFile }
                )
            }
            if (tumor_full_contig_context_b) {
                tumorWGS(tumor_assembly_reads_ch)
            }
            if (tumor_enable_mags_b) {
                tumorPksMAG(tumor_assembly_reads_ch)
            }
        }
    }
    if (do_hmm) {
        HMM_OUT = pksProfilerHMM(PROFILING_READS)

        HMM_OUT.profile
            .set { PKS_HMM_OUT }

        QC_FRAGMENTS = QC_FRAGMENTS.mix(
            HMM_OUT.qc.map { _sampleID, qc_file -> qc_file }
        )
    }

    // ---------- STEP 3: Plotting (align only) ----------
    if (do_align) {
        PKS_ALIGN_OUT
            .map { output -> output[2]  }   // bedgraph
            .set { COVERAGE_BEDGRAPH }

        plotPKS(COVERAGE_BEDGRAPH)
    }

	// ---------- STEP 3b: Optional PKS-island taxa profiling (align only) ----------
    if (do_align && params.pks_taxa) {
        PKS_ALIGN_OUT
			.map { sampleID, _covtxt, _bedgraph, _counts, bam, bai ->
			    tuple(sampleID, bam, bai)
			}
            .set { PKS_BAM_FOR_TAXA }

        extractPksIslandReads(PKS_BAM_FOR_TAXA)
            .set { PKS_ISLAND_FASTQ }

		// Bracken emits two channels since F10 added the QC fragment, so the reports
		// channel is addressed by name. `Bracken(...).set { ... }` captured the
		// multi-channel object instead, and the next `.map` on it failed at parse time --
		// which is every --pks_taxa run, and nothing without a live Nextflow caught it.
		Bracken(PKS_ISLAND_FASTQ)
		Bracken.out.reports.set { BRACKEN_PER_SAMPLE }

		// F10: the taxonomy outcome reaches pks.qc.summary.tsv, which already carries a
		// row for every sample, so "no species" stops being indistinguishable from
		// "never classified".
		QC_FRAGMENTS = QC_FRAGMENTS.mix(Bracken.out.qc.map { _sampleID, qc_file -> qc_file })

		// F02: the files, plus one list naming them, instead of one argument each.
		BRACKEN_PER_SAMPLE
		.map { _sampleID, _report, _classified, _unclassified, _brG, _brS, _gk, _sk, _gmpa, _smpa, speciesSupport -> speciesSupport }
		.multiMap { support ->
		    files: support
		    rows:  support.name
		}
		.set { CLB_SUPPORT_MERGE }

		def combine_clb_support_script = file(
		    "${params.scripts}/combine_clb_species_support.py",
		    checkIfExists: true
		)

		combineClbTaxonomySupport(
		    CLB_SUPPORT_MERGE.files.collect(),
		    combine_clb_support_script,
		    CLB_SUPPORT_MERGE.rows.collectFile(name: 'pks.clb_species_support.list',
		                                       newLine: true, sort: true)
		)

		BRACKEN_PER_SAMPLE
		.map { sampleID, _kreport, _classified, _unclassified, brG, brS, _gk, _sk, _gmpa, _smpa, _speciesMatrix ->
			    tuple(sampleID, brG, brS)
		}
	    .set { BRACKEN_GS_REPORTS }

		plotPKSTaxa(BRACKEN_GS_REPORTS)

	
		BRACKEN_PER_SAMPLE
		.map { _sampleID, _report, _classified, _unclassified, _brG, _brS, _gk, _sk, gmpa, smpa, _speciesMatrix ->
		    [gmpa, smpa]
		}
	   .flatten()
       .collect()
       .set { BRACKEN_MPA_FILES }

		combinePKSTaxa(BRACKEN_MPA_FILES)
		

    }

    // ---------- STEP 4: Master tables ----------
    // The sample identifier travels with its counts file into the merge, as a manifest.
    // It used to be dropped here and guessed back from the filename, which collapsed the
    // distinct sheet identifiers `case` and `case.txt` onto one column (F03). Passing a
    // manifest also means the merge command holds one filename, not one per sample, which
    // is what F02 is about.
    if (do_align) {
        PKS_ALIGN_OUT
            .multiMap { output ->
                files: output[3]                                            // counts.txt
                rows:  "${output[0].toString().trim()}\t${output[3].name}"
            }
            .set { ALIGN_MERGE }

        masterTableAlign(
            ALIGN_MERGE.files.collect(),
            ALIGN_MERGE.rows.collectFile(name: 'pks.align.counts.manifest.tsv',
                                         newLine: true, sort: true)
        )
        // Ludmil, revised report finding 5: pks.gene.counts.align.txt is the first
        // thing build_cohort_report.py reads. Nothing previously made cohortReport
        // wait for it -- only for the per-sample read-evidence tier -- so the cohort
        // report could start, and publish an empty gene matrix, before this landed.
        cohort_report_gate = cohort_report_gate.mix(masterTableAlign.out)
    }

    if (do_hmm) {
        PKS_HMM_OUT
            .multiMap { output ->
                files: output[3]                                            // hmm_counts.tsv
                rows:  "${output[0].toString().trim()}\t${output[3].name}"
            }
            .set { HMM_MERGE }

        masterTableHMM(
            HMM_MERGE.files.collect(),
            HMM_MERGE.rows.collectFile(name: 'pks.hmm.counts.manifest.tsv',
                                       newLine: true, sort: true)
        )
    }

    // ---------- STEP 5: Cohort QC ----------
    // F02: the fragments, plus one list naming them, instead of one argument each.
    QC_FRAGMENTS
        .multiMap { fragment ->
            files: fragment
            rows:  fragment.name
        }
        .set { QC_MERGE }

    def qc_summary_script = file(
        "${params.scripts}/build_qc_summary.py",
        checkIfExists: true
    )

    masterQCSummary(
        QC_MERGE.files.collect(),
        qc_summary_script,
        EXPECTED_SAMPLE_IDS,
        QC_MERGE.rows.collectFile(name: 'pks.qc.fragments.list', newLine: true, sort: true)
    )

    // ---------- STEP 6: One page for the cohort ----------
    // Two thousand sample directories and no view across them is not a result anyone
    // can read. This runs last, after the cohort tables exist, and only reads what the
    // pipeline already wrote -- the tier from read_evidence.tsv, the breadth from the
    // depth that produced it. An earlier hand-built version recomputed both with the
    // v0.0.1 thresholds and could disagree with the run it described.
    def cohort_report_script = file("${params.scripts}/build_cohort_report.py",
                                    checkIfExists: true)
    cohortReport(masterQCSummary.out.mix(cohort_report_gate).collect(),
                 cohort_report_script,
                 EXPECTED_SAMPLE_IDS)
}
