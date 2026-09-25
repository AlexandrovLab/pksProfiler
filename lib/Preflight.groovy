import groovy.json.JsonOutput

/**
 * Runs the preflight contract at launch, before any task is submitted.
 *
 * F20. The checks themselves live in scripts/preflight.py rather than here, for the
 * same reason the run report does: a Groovy class in lib/ cannot be reached by the
 * test suite, and validation that nothing tests is validation nobody should trust.
 * This class only marshals the resolved parameters and reports what the script said.
 */
class Preflight {

    /**
     * Returns the script's exit status; 0 means every check passed.
     *
     * logSink takes (level, message) with level one of 'info'/'warn'/'error' and
     * dispatches to Nextflow's own log -- not passed as the log object itself,
     * because Nextflow 24.10.0 (this manifest's declared floor, and the version CI
     * pins) fails to compile a bare reference to the implicit `log` binding used as
     * a value. See the call site in main.nf for the full explanation.
     */
    static int run(Object scriptPath, Map config, Closure logSink) {
        def script = new File(scriptPath.toString())
        if (!script.isFile()) {
            logSink('warn', "Preflight script not found at ${script}; skipping the launch checks")
            return 0
        }

        File spec = File.createTempFile('pksprofiler-preflight', '.json')
        spec.deleteOnExit()
        spec.text = JsonOutput.toJson(config)

        try {
            def process = new ProcessBuilder(['python3', script.absolutePath,
                                              '--config', spec.absolutePath]).start()
            String out = process.inputStream.text
            String err = process.errorStream.text
            int status = process.waitFor()

            // One call, not one per line: Nextflow decorates every log.error with its
            // own banner, and twenty problems became forty lines of "check .nextflow.log".
            String text = out.readLines().findAll { it.trim() }.collect { "  ${it}" }.join('\n')
            if (text) {
                logSink(status == 0 ? 'info' : 'error', "\n${text}")
            }
            if (err?.trim()) logSink('warn', "preflight: ${err.trim()}")
            return status
        }
        catch (Exception problem) {
            // No python3 on the launch node is itself a problem -- every helper the
            // pipeline calls needs it -- but say so plainly rather than as a stack trace.
            logSink('error', "Could not run the preflight checks: ${problem.message}")
            return 1
        }
        finally {
            spec.delete()
        }
    }
}
