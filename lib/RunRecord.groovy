import groovy.json.JsonOutput
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.time.ZonedDateTime
import java.time.format.DateTimeFormatter

/**
 * The equivalence record, launch half: what code, what parameters, what references.
 *
 * F16. Nothing tied the released source to the code that actually produced a cohort.
 * Every run now writes one JSON record per session at launch -- before any task runs,
 * so a run that dies still says what it was attempting. The outcome and the output
 * checksums are added afterwards by scripts/finalize_provenance.py, called from the
 * workflow.onComplete handler in nextflow.config.
 *
 * The split is not aesthetic. Nextflow's script parser rejects a top-level
 * `workflow.onComplete` block, and classes in lib/ are not visible to the config
 * parser, so the completion half cannot call into here. Keeping it in Python also puts
 * the rendering and checksum logic somewhere the test suite can reach.
 *
 * One record per session, not one per results directory: with `-resume` a results tree
 * is usually the product of several sessions, and F15 is exactly the case where a later
 * session ran different code against the same tree.
 */
class RunRecord {

    static final int VERSION = 1

    /** Everything knowable before the first task runs. */
    static Map start(Map context) {
        return [
            record_version: VERSION,
            status        : 'started',
            code          : context.code,
            invocation    : context.invocation,
            dependencies  : context.dependencies,
            environment   : context.environment,
            run           : [
                session_id: context.session_id,
                started   : DateTimeFormatter.ISO_OFFSET_DATE_TIME.format(ZonedDateTime.now().withNano(0)),
            ],
            outputs       : [:],
            notes         : context.notes ?: [],
        ]
    }

    /**
     * Where this session's record lives. Derived from the session id alone, with no
     * timestamp, so the finaliser can find the file this half wrote without the two
     * sharing any state. The launch time is inside the record.
     */
    static Path path(Object runsDir, Object sessionId) {
        Path dir = Paths.get(runsDir.toString(), 'provenance')
        Files.createDirectories(dir)
        return dir.resolve("${sessionId}.json")
    }

    static void write(Path target, Map record) {
        Files.createDirectories(target.getParent())
        target.text = JsonOutput.prettyPrint(JsonOutput.toJson(record)) + "\n"
    }

    /**
     * Which stages the parameters asked for. runs/trace.txt is what actually executed;
     * this is what was requested, which is the part a reader cannot reconstruct.
     */
    static List<String> stages(Map params) {
        def requested = [
            "profiling:${params.profiling_method}".toString(),
            "sample_type:${params.sample_type}".toString(),
            "input:${params.input_data_type}".toString(),
        ]
        [
            'enable_mags', 'tumor_targeted_assembly',
            'diamond_rescue', 'pks_community_taxa',
            'enable_strain_typing', 'pks_taxa',
        ].each { flag ->
            if (params[flag]?.toString()?.toBoolean()) requested << flag
        }
        return requested
    }

    /**
     * Which commit this is. workflow.commitId is set only when Nextflow pulled the
     * project itself; most runs are launched from a clone, where it is null and git is
     * the only source. A dirty working tree is recorded rather than ignored -- "commit
     * X" is not true of a run launched with uncommitted edits, and F16 is precisely
     * about claims of that kind.
     */
    static Map gitIdentity(Object projectDir, Object commitId, Object revision) {
        def identity = [
            commit        : commitId?.toString(),
            commit_source : commitId ? 'nextflow' : 'unknown',
            revision      : revision?.toString(),
            dirty         : false,
            modified_files: [],
        ]
        try {
            def dir = projectDir.toString()
            def head = capture(['git', '-C', dir, 'rev-parse', 'HEAD'])
            if (head) {
                if (!identity.commit) {
                    identity.commit = head
                    identity.commit_source = 'git'
                }
                else if (identity.commit != head) {
                    identity.commit_source = "nextflow (git HEAD differs: ${head})".toString()
                }
                def status = capture(['git', '-C', dir, 'status', '--porcelain'])
                if (status) {
                    identity.dirty = true
                    identity.modified_files = status.readLines().collect { it.trim() }
                }
            }
        }
        catch (Exception ignored) {
            // not a git work tree, or no git on PATH: keep what we have
        }
        return identity
    }

    private static String capture(List<String> command) {
        def process = new ProcessBuilder(command).start()
        def output = process.inputStream.text.trim()
        return process.waitFor() == 0 ? output : null
    }
}
