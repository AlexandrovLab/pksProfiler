#!/usr/bin/env bash
# T7 (Ludmil, 28 Sep): "freeze and record versions/checksums of KrakenUniq/Bracken, GTDB,
# CheckM2, geNomad, MLST and any other DB".
#
# RUN_REPORT's DATABASES section listed the PATH each database was read from and nothing
# else. A path is not a version: the GTDB-Tk tree at .../gtdbtk_r232 can be replaced with
# a newer release in place, and every run before and after would report the same line.
# lib/Provenance.data() does fingerprint these, but only the top level (for a directory,
# the names/sizes/mtimes of its immediate children), it is truncated to 12 hex characters,
# and it is a cache key, not something a reader can compare against a release.
#
#   record_database_versions.sh OUT [--manifest FILE] NAME=PATH ...
#
#   OUT        the record, written as TSV
#   --manifest conf/databases.tsv; when given, each database's fingerprint is compared
#              with the frozen one and the verdict goes in the `known` column
#   NAME=PATH  the pipeline parameter and the path it was set to; an unset or empty
#              PATH is skipped, so a caller can pass every database unconditionally
#
# Columns: name, path, kind, status, version, size_bytes, file_count,
#          fingerprint_method, fingerprint, known
#
# FINGERPRINT METHOD, and why there are two.
#
#   sha256:content   sha256 of the file itself. Used for a regular file at or below
#                    PKS_DB_HASH_MAX_BYTES (default 2 GiB). This is the real checksum.
#
#   sha256:listing   sha256 over a sorted listing of `<relative path>\t<size in bytes>`
#                    for every regular file in the tree (or, for one oversized file, that
#                    single line). Used for directories and for files above the limit.
#
# The listing method is used because the databases are too large to checksum on every
# run: KrakenUniq's database.kdb alone is 574 GB, and the GTDB-Tk r232 tree and the
# host-depletion minimap2 indexes are tens of gigabytes more. Reading all of that at
# roughly 1 GB/s is tens of minutes of pure I/O per run, repeated for every run, against
# a shared Lustre filesystem -- the cost is not the hashing, it is dragging the bytes
# across the network for no new information. The listing is O(file count), completes in
# well under a second, and changes whenever a database is rebuilt, re-downloaded or
# extended, because a rebuild changes file sizes or the set of files. What it does NOT
# catch is an edit that rewrites bytes in place while preserving every file's size. No
# cheap check catches that; the frozen `sha256:content` of each database's own small key
# index file, recorded in conf/databases.tsv, is the backstop. Sizes, not mtimes, are
# hashed: an rsync or a restore from backup rewrites mtimes without changing the data,
# and that is not a different database.
#
# VERSION is read from the release marker the database itself ships, never guessed from
# the directory name -- a directory can be renamed or a newer release dropped into the
# old name, which is the failure this exists to catch:
#
#   GTDB-Tk    metadata/metadata.txt, VERSION_DATA=      -> r232
#   geNomad    version.txt                               -> 1.9
#   CheckM2    CONTENTS.json, "version":                 -> 3
#   KrakenUniq no marker file in the release             -> unknown (see conf/databases.tsv)
#
# No `set -e`: a database that cannot be read is recorded with status `unreadable` and
# the others are still recorded. This script never fails a run -- it is a record, not a
# gate -- so it always exits 0 unless its arguments are wrong.
set -uo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: record_database_versions.sh OUT [--manifest FILE] NAME=PATH ..." >&2
    exit 2
fi

OUT=$1
shift

MANIFEST=""
if [[ "${1:-}" == "--manifest" ]]; then
    MANIFEST=${2:-}
    shift 2
fi

# Above this, a regular file gets the listing fingerprint instead of a real checksum.
HASH_MAX_BYTES=${PKS_DB_HASH_MAX_BYTES:-2147483648}
[[ "$HASH_MAX_BYTES" =~ ^[0-9]+$ ]] || HASH_MAX_BYTES=2147483648

sha256() { sha256sum | cut -d' ' -f1; }

# The release the database declares about itself, looking in `root` and then in its
# parent directories, up to `climb` levels. The climb is needed because a parameter can
# point at a file inside the release rather than at its root: --checkm2_db is the
# .dmnd itself, two levels below the CONTENTS.json that names the release.
declared_version() {
    local root=$1 climb=${2:-0} value="" level=0
    while :; do
        if [[ -f "$root/metadata/metadata.txt" ]]; then
            value=$(sed -n 's/^VERSION_DATA=[[:space:]]*//p' "$root/metadata/metadata.txt" | sed -n 1p)
        elif [[ -f "$root/version.txt" ]]; then
            value=$(sed -n 1p "$root/version.txt")
        elif [[ -f "$root/CONTENTS.json" ]]; then
            value=$(sed -n 's/.*"version"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$root/CONTENTS.json" | sed -n 1p)
        fi
        [[ -n "$value" ]] && break
        (( level++ >= climb )) && break
        local parent
        parent=$(dirname "$root")
        [[ "$parent" == "$root" ]] && break
        root=$parent
    done
    printf '%s' "${value//[[:space:]]/}"
}

printf 'name\tpath\tkind\tstatus\tversion\tsize_bytes\tfile_count\tfingerprint_method\tfingerprint\tknown\n' > "$OUT"

for entry in "$@"; do
    name=${entry%%=*}
    path=${entry#*=}

    # Unset parameters arrive as an empty string or Nextflow's "null".
    if [[ -z "$path" || "$path" == "null" ]]; then
        continue
    fi

    kind=missing status=missing version="" size=0 count=0 method="" fingerprint=""

    if [[ -d "$path" ]]; then
        kind=directory
        # One pass, NUL-delimited: names with spaces or newlines cannot split a field.
        listing=$(cd "$path" 2>/dev/null && find . -type f -printf '%P\t%s\n' 2>/dev/null | LC_ALL=C sort)
        if [[ -z "$listing" ]]; then
            status=unreadable
        else
            status=ok
            count=$(printf '%s\n' "$listing" | wc -l)
            size=$(printf '%s\n' "$listing" | awk -F '\t' '{ total += $2 } END { printf "%.0f", total }')
            method=sha256:listing
            fingerprint=$(printf '%s\n' "$listing" | sha256)
            version=$(declared_version "$path")
        fi
    elif [[ -f "$path" ]]; then
        kind=file
        count=1
        size=$(stat -c '%s' "$path" 2>/dev/null || echo 0)
        version=$(declared_version "$(dirname "$path")" 2)
        version=$(declared_version "$(dirname "$path")" 2)
        if [[ ! -r "$path" ]]; then
            status=unreadable
        elif (( size <= HASH_MAX_BYTES )); then
            status=ok
            method=sha256:content
            fingerprint=$(sha256 < "$path")
        else
            status=ok
            method=sha256:listing
            fingerprint=$(printf '%s\t%s\n' "$(basename "$path")" "$size" | sha256)
        fi
    else
        # A bowtie2/minimap2 index is given as a prefix, not an existing file.
        siblings=$(find "$(dirname "$path")" -maxdepth 1 -type f -name "$(basename "$path")*" \
                       -printf '%f\t%s\n' 2>/dev/null | LC_ALL=C sort)
        if [[ -n "$siblings" ]]; then
            kind=prefix status=ok
            count=$(printf '%s\n' "$siblings" | wc -l)
            size=$(printf '%s\n' "$siblings" | awk -F '\t' '{ total += $2 } END { printf "%.0f", total }')
            method=sha256:listing
            fingerprint=$(printf '%s\n' "$siblings" | sha256)
        fi
    fi

    # Compare with the frozen manifest: does this run read the database we froze?
    known=no_manifest
    if [[ -n "$MANIFEST" && -f "$MANIFEST" ]]; then
        known=$(awk -F '\t' -v n="$name" -v f="$fingerprint" '
            !/^#/ && $1 == n {
                found = 1
                if ($8 == f) { print "frozen"; exit }
                print "CHANGED(frozen " substr($8, 1, 16) "...)"
                exit
            }
            END { if (!found) print "not_in_manifest" }' "$MANIFEST")
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$path" "$kind" "$status" "$version" "$size" "$count" \
        "$method" "$fingerprint" "$known" >> "$OUT"
done

exit 0
