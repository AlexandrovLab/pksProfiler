/**
 * Which minimap2 preset and index each host-depletion pass uses, per sample type.
 *
 * Ludmil, 26 Sep mode memo section 4: `-ax sr` is not splice-aware, so for RNA-seq a read
 * spanning an exon-exon junction aligns to the genome in two pieces or not at all, and
 * survives human depletion as if it were non-human. For --sample_type rna the FIRST
 * (human) pass is therefore `-ax splice:sr` against --rna_splice_db (an index built with
 * `-x splice:sr`, or a FASTA), with annotated junctions from --rna_junctions_bed (`-j`)
 * when given. The T2T+phiX and pangenome passes stay `-ax sr`: there is no junction
 * annotation for T2T-CHM13, phiX or the pangenome haplotypes, and those passes exist to
 * catch the residual human (and spike-in) reads GRCh38 missed, not to re-model splicing.
 *
 * Kept here, not inline in Modules/map_reads.nf, so the choice is one table that the
 * module, the preflight contract and tests/check_depletion_plan.sh all read.
 * --minimap2_frag_mode (E0) applies to every short-read call, spliced or not.
 */
class HostDepletion {

    /** minimap2 -x preset for the human pass, by sample type. */
    static final Map<String, String> HUMAN_PRESET = [
        wgs: 'sr', wes: 'sr', metagenome: 'sr', rna: 'splice:sr',
    ].asImmutable()

    /** minimap2 -x preset for the T2T+phiX and pangenome passes. */
    static final String RESIDUAL_PRESET = 'sr'

    static String fragArg(Object fragMode) {
        return fragMode?.toString() == 'no' ? '--frag=no' : ''
    }

    /**
     * The human pass: [preset, index, args]. `args` is every flag between the preset and
     * --secondary=no, already shell-quoted.
     */
    static Map humanPass(Object sampleType, Object fragMode, Object hg38Db,
                         Object rnaSpliceDb, Object rnaJunctionsBed) {
        String mode = sampleType?.toString()
        boolean rna = mode == 'rna'
        List<String> args = []
        String frag = fragArg(fragMode)
        if (frag) args << frag
        if (rna && rnaJunctionsBed) args << "-j \"${rnaJunctionsBed}\"".toString()
        return [preset: HUMAN_PRESET.get(mode, 'sr'),
                index : (rna ? rnaSpliceDb : hg38Db)?.toString(),
                args  : args.join(' ')]
    }

    /** The T2T+phiX pass: [preset, index, args]. */
    static Map t2tPass(Object fragMode, Object t2tPhixDb) {
        return [preset: RESIDUAL_PRESET, index: t2tPhixDb?.toString(), args: fragArg(fragMode)]
    }

    /** The optional pangenome passes: [preset, args]; the indexes are found at run time. */
    static Map pangenomePass(Object fragMode) {
        return [preset: RESIDUAL_PRESET, args: fragArg(fragMode)]
    }

    /** The parameter naming the human index a sample type needs, for required-param checks. */
    static String humanIndexParam(Object sampleType) {
        return sampleType?.toString() == 'rna' ? 'rna_splice_db' : 'hg38_db'
    }
}
