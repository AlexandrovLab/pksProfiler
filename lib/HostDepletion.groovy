/**
 * Which minimap2 preset and index each host-depletion pass uses, per sample type and
 * read type -- and, for long reads, the profiling aligner's preset too.
 *
 * Ludmil, 26 Sep mode memo section 4: `-ax sr` is not splice-aware, so for RNA-seq a read
 * spanning an exon-exon junction aligns to the genome in two pieces or not at all, and
 * survives human depletion as if it were non-human. For --sample_type rna the FIRST
 * (human) pass is therefore spliced -- `-ax splice:sr` against --rna_splice_db (an index
 * built with `-x splice:sr`, or a FASTA) -- with annotated junctions from
 * --rna_junctions_bed (`-j`) when given. The T2T+phiX and pangenome passes stay unspliced
 * (`-ax sr`): there is no junction annotation for T2T-CHM13, phiX or the pangenome
 * haplotypes, and those passes exist to catch the residual human (and spike-in) reads
 * GRCh38 missed, not to re-model splicing.
 *
 * Long reads (owner decision, 29 Sep; --read_type ont | hifi, wgs and rna only): every
 * pass uses the long-read preset -- map-ont / map-hifi for DNA, splice / splice:hq for the
 * RNA human pass -- against the *_long indexes. A minimap2 index carries the k-mer and
 * window settings of the preset it was built with, so a short-read .mmi is not reused for
 * long reads: --hg38_db_long, --t2t_phix_db_long and --rna_splice_db_long are separate
 * parameters (a FASTA works for any preset). --frag only changes how minimap2 treats
 * short paired reads, so --minimap2_frag_mode (E0) reaches every short-read call and no
 * long-read call.
 *
 * Kept here, not inline in Modules/map_reads.nf, so the choice is one table that the
 * module, main.nf's required-parameter check and tests/check_depletion_plan.sh all read
 * (scripts/preflight.py mirrors humanIndexParam/t2tIndexParam).
 */
class HostDepletion {

    static final List<String> READ_TYPES = ['short', 'ont', 'hifi'].asImmutable()
    /** Sample types that accept long reads; metagenome and wes are short-read only. */
    static final List<String> LONG_READ_SAMPLE_TYPES = ['wgs', 'rna'].asImmutable()

    /** Unspliced DNA preset by read type: every non-RNA pass, and long-read profiling. */
    static final Map<String, String> DNA_PRESET = [
        short: 'sr', ont: 'map-ont', hifi: 'map-hifi',
    ].asImmutable()

    /** Spliced preset by read type: the human pass of --sample_type rna only. */
    static final Map<String, String> RNA_PRESET = [
        short: 'splice:sr', ont: 'splice', hifi: 'splice:hq',
    ].asImmutable()

    static String resolveReadType(Object readType) {
        return readType?.toString() ?: 'short'
    }

    static boolean isLong(Object readType) {
        return resolveReadType(readType) != 'short'
    }

    static String fragArg(Object readType, Object fragMode) {
        return (!isLong(readType) && fragMode?.toString() == 'no') ? '--frag=no' : ''
    }

    /**
     * The human pass: [preset, index, args]. `args` is every flag between the preset and
     * --secondary=no, already shell-quoted.
     */
    static Map humanPass(Object sampleType, Object readType, Object fragMode,
                         Object hg38Db, Object hg38DbLong,
                         Object rnaSpliceDb, Object rnaSpliceDbLong, Object rnaJunctionsBed) {
        String rt = resolveReadType(readType)
        boolean rna = sampleType?.toString() == 'rna'
        boolean lr = isLong(rt)
        List<String> args = []
        String frag = fragArg(rt, fragMode)
        if (frag) args << frag
        if (rna && rnaJunctionsBed) args << "-j \"${rnaJunctionsBed}\"".toString()
        def index = rna ? (lr ? rnaSpliceDbLong : rnaSpliceDb) : (lr ? hg38DbLong : hg38Db)
        return [preset: (rna ? RNA_PRESET : DNA_PRESET)[rt],
                index : index?.toString(),
                args  : args.join(' ')]
    }

    /** The T2T+phiX pass: [preset, index, args]. Never spliced. */
    static Map t2tPass(Object readType, Object fragMode, Object t2tPhixDb, Object t2tPhixDbLong) {
        String rt = resolveReadType(readType)
        return [preset: DNA_PRESET[rt],
                index : (isLong(rt) ? t2tPhixDbLong : t2tPhixDb)?.toString(),
                args  : fragArg(rt, fragMode)]
    }

    /**
     * The optional pangenome passes: [preset, args]; the indexes are found at run time.
     * For long reads, build the pangenome .mmi files with the matching preset.
     */
    static Map pangenomePass(Object readType, Object fragMode) {
        String rt = resolveReadType(readType)
        return [preset: DNA_PRESET[rt], args: fragArg(rt, fragMode)]
    }

    /** The parameter naming the human index a run needs, for required-param checks. */
    static String humanIndexParam(Object sampleType, Object readType) {
        String base = sampleType?.toString() == 'rna' ? 'rna_splice_db' : 'hg38_db'
        return isLong(readType) ? "${base}_long".toString() : base
    }

    /** The parameter naming the T2T+phiX index a run needs. */
    static String t2tIndexParam(Object readType) {
        return isLong(readType) ? 't2t_phix_db_long' : 't2t_phix_db'
    }

    /** minimap2 -x preset for long-read profiling against the IHE3034 genome FASTA. */
    static String profilingPreset(Object readType) {
        return DNA_PRESET[resolveReadType(readType)]
    }
}
