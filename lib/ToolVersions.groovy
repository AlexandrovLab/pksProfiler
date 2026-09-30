/**
 * The runtime tool-version record (F16, "tool versions"): one line at the top of every
 * process that runs an external tool.
 *
 *     ${ToolVersions.record(task, params, 'minimap2 samtools bgzip:htslib')}
 *
 * renders a call to scripts/record_tool_versions.sh, which asks each binary on PATH --
 * inside the task, after conda activation -- for its version and location, writes them
 * to the task work dir's .tool_versions.tsv, and in strict mode stops the task (exit 78)
 * when one ran from outside the task's conda environment or differs from its pin in
 * the environment's yml. RUN_REPORT's SOFTWARE section reads conda-meta after the run,
 * which is what conda installed; on TSCC a user minimap2 2.24 ran every host-depletion
 * pass while it said 2.31. scripts/finalize_provenance.py reads the records back through
 * the trace's workdir column and compares them.
 *
 * The env passed is the task's own conda directive (task.conda), so the pins compared
 * are always the ones the task was built from. The record goes next to .command.sh
 * rather than into the current directory: under `scratch true` the current directory is
 * node-local and discarded, while $0 is still the work dir's .command.sh. It is a dot
 * file so no output glob can pick it up and publish it. `|| exit $?` because Nextflow runs
 * .command.sh with -u but not -e, and not every script sets -e before this line.
 */
class ToolVersions {

    static final List<String> MODES = ['strict', 'warn']

    /** The exit status record_tool_versions.sh uses in strict mode (EX_CONFIG). */
    static final int EXIT_STATUS = 78

    static String record(task, params, String tools) {
        def env = (task.conda ?: '').toString()
        return "bash \"${params.scripts}/record_tool_versions.sh\" ${params.tool_version_check} " +
               "'${env}' '${task.process}' \"\$(dirname \"\$0\")/.tool_versions.tsv\" ${tools} || exit \$?"
    }

    /** null when valid, else the message main.nf fails with. */
    static String validate(Object mode) {
        def value = mode?.toString()
        return value in MODES ? null :
            "--tool_version_check must be one of ${MODES.join(', ')} (got ${mode})"
    }
}
