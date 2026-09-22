import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.security.MessageDigest

/**
 * Fingerprints for the things a task depends on but does not declare: helper scripts,
 * reference indexes, annotation files and databases.
 *
 * F15. Nextflow hashes a task from its declared inputs and its rendered script. A
 * reference or helper referenced only as a path string inside the script is neither,
 * so rebuilding an index in place, or editing a helper, leaves every hash unchanged
 * and `-resume` reuses tasks that no longer correspond to the current code. Injecting
 * these digests into the script makes the content part of the hash.
 *
 * Two modes, because the inputs differ by orders of magnitude:
 *
 *   code(...)  content of every file, recursively. For the scripts directory -- 22
 *              files, under half a megabyte -- an edit must invalidate, and nothing
 *              else will catch it.
 *   data(...)  name, size and modification time of the top level only. An 8.6 GB
 *              minimap2 index or a GTDB-Tk database cannot be read on every launch,
 *              and a directory walk of one is slow enough on Lustre to matter. A
 *              rebuilt database changes size or timestamp, which this sees; a
 *              byte-level edit that preserves both, which it does not, is not a
 *              failure mode any cheap check catches.
 *
 * A path that does not exist is treated as an index prefix -- Bowtie2 and minimap2
 * indexes are named `<prefix>.1.bt2`, `<prefix>.fna` and so on -- and the siblings
 * sharing that prefix are fingerprinted instead.
 */
class Provenance {

    /** Content of every file below `target`, recursively. For code, not for data. */
    /**
     * Build byproducts, which must not reach the digest.
     *
     * code() hashed everything under scripts/, including __pycache__ -- which is
     * gitignored, so the task hash depended on files that are not source at all. Running
     * any helper script rewrote a .pyc and invalidated all fifteen processes that carry
     * this digest, so a -resume after executing one script re-ran the entire cohort. At
     * 2,054 samples that is hours of compute for a byproduct. It also made the digest
     * non-reproducible across machines -- different interpreter versions or optimisation
     * levels write different .pyc -- which undermines the equivalence record F16 exists
     * to provide.
     *
     * F15 asks that a changed reference or helper script invalidate the cache. A .pyc is
     * neither. Observed 2026-09-22: one regenerated .pyc cost a full cache miss.
     */
    private static boolean generated(File f) {
        def path = f.absolutePath
        return path.contains('/__pycache__/') ||
               path.contains('/.ipynb_checkpoints/') ||
               path.contains('/.git/') ||
               f.name.endsWith('.pyc') ||
               f.name.endsWith('.pyo') ||
               f.name == '.DS_Store'
    }

    static String code(Object target) {
        Path root = resolve(target)
        if (root == null) return 'unset'
        if (!Files.exists(root)) return missing(root)

        def digest = MessageDigest.getInstance('SHA-256')
        def files = []
        if (Files.isDirectory(root)) {
            root.toFile().eachFileRecurse { f ->
                if (f.isFile() && !generated(f)) files << f
            }
        } else {
            files << root.toFile()
        }
        files.sort { it.absolutePath }
        files.each { f ->
            digest.update(f.name.getBytes('UTF-8'))
            digest.update(f.bytes)
        }
        return "${short_(digest)}/content"
    }

    /** Name, size and mtime of the top level. For references and databases. */
    static String data(Object target) {
        Path root = resolve(target)
        if (root == null) return 'unset'
        if (!Files.exists(root)) return prefix(root)

        def digest = MessageDigest.getInstance('SHA-256')
        def entries = Files.isDirectory(root) ? listing(root) : [describe(root)]
        entries.sort()
        entries.each { digest.update(it.getBytes('UTF-8')) }
        return "${short_(digest)}/meta"
    }

    // ---- internals ----

    private static Path resolve(Object target) {
        if (target == null) return null
        String text = target.toString().trim()
        return text ? Paths.get(text) : null
    }

    /** Siblings sharing a basename prefix: how multi-file indexes are named. */
    private static String prefix(Path root) {
        Path parent = root.getParent()
        String base = root.getFileName().toString()
        if (parent == null || !Files.exists(parent)) return missing(root)

        def siblings = []
        parent.toFile().eachFile { f -> if (f.name.startsWith(base)) siblings << describe(f.toPath()) }
        if (!siblings) return missing(root)

        def digest = MessageDigest.getInstance('SHA-256')
        siblings.sort()
        siblings.each { digest.update(it.getBytes('UTF-8')) }
        return "${short_(digest)}/prefix"
    }

    private static List<String> listing(Path dir) {
        def entries = []
        dir.toFile().eachFile { f -> entries << describe(f.toPath()) }
        return entries
    }

    private static String describe(Path p) {
        def file = p.toFile()
        long size = file.isDirectory() ? -1L : file.length()
        return "${file.name}:${size}:${file.lastModified()}"
    }

    private static String missing(Path root) {
        return "absent(${root.getFileName()})"
    }

    private static String short_(MessageDigest digest) {
        return digest.digest().encodeHex().toString().substring(0, 12)
    }
}
