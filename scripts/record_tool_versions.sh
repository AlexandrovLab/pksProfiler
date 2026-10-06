#!/usr/bin/env bash
# Record the version and path of every tool a task is about to run, and refuse to run
# one that is not the pinned build from the task's own conda environment.
#
# F16 (Ludmil) asks for tool versions. RUN_REPORT's SOFTWARE section lists what conda
# INSTALLED, read from conda-meta after the run -- not what ran. On TSCC a user
# minimap2 2.24 earlier on PATH than the conda base ran every host-depletion pass
# while the report said 2.31 (29 Sep; see nextflow.config, process.beforeScript). This
# asks the binary itself, from inside the task, after conda activation.
#
#   record_tool_versions.sh MODE ENV PROCESS OUT TOOL[:PACKAGE] ...
#
#   MODE     strict | warn (--tool_version_check)
#   ENV      the task's conda directive: a .yml, whose `- package=version` pins are
#            compared, or an existing environment (path check only), or ''
#   PROCESS  task.process
#   OUT      the record, written as TSV: the task work dir's .tool_versions.tsv, which
#            scripts/finalize_provenance.py reads through the trace's workdir column
#   TOOL     a command the task runs; PACKAGE is the conda package that provides it
#            when the names differ (bgzip:htslib, featureCounts:subread)
#
# Per tool, one of these statuses:
#   ok                inside $CONDA_PREFIX and, when pinned, the pinned version
#   unpinned          inside $CONDA_PREFIX; the yml pins no version for the package
#   version_unread    inside $CONDA_PREFIX; its version output held no version number
#   OUTSIDE_ENV       resolved somewhere other than the task's conda environment
#   VERSION_MISMATCH  reports a version other than its pin
#   NOT_FOUND         not on PATH at all
#   NO_ENV            no conda environment is active in this task
#
# The upper-case statuses stop the task in strict mode with exit 78 (EX_CONFIG), which
# conf/base.config turns into 'terminate': the environment is run-wide, so every other
# sample would hit the same problem. In warn mode they are printed and recorded, and
# RUN_REPORT lists them under TOOL VERSION MISMATCH.
#
# No `set -e`: a tool whose version flag exits non-zero is recorded, not fatal. Every
# pipe ends in `sed -n 1p` or `grep` guarded by `|| true`, never `head`, which under
# pipefail can make the writer die of SIGPIPE.
#
# Reading a version is retried (rc1.1). On the rc1 TCGA WES run (1,304 samples) a few
# tasks per tool recorded version_unread (bamCoverage 11, fastp 8, samtools 4, bowtie2 3,
# bgzip 1; RNA python3 1 of 3) while every other task of the same tool, from the same
# env, read fine. The parse is deterministic, so in those tasks the single call itself
# produced no version line: killed by the fixed 60 s `timeout` at a stalled cold start
# (exec and shared libraries -- or, for bamCoverage, ~30 Python modules, numpy, pysam and
# matplotlib -- loaded from a network-filesystem conda env as a whole cohort starts at
# once), or a transient exec failure. Which one could not be told afterwards: the exit
# status was discarded and `raw` kept only lines with a version number in them, so it was
# empty. Now each read:
#   - captures stdout and stderr together (fastp, bracken and others print to stderr);
#   - runs under `timeout -k 5`, PKS_TOOL_VERSION_TIMEOUT seconds (default 60) on the
#     first attempt and that times the attempt number after;
#   - is retried, up to PKS_TOOL_VERSION_ATTEMPTS (default 3) in all, after a backoff of
#     PKS_TOOL_VERSION_BACKOFF (default 2) times the attempt number seconds, when it
#     exited non-zero or printed nothing and no version was found. A call that exited 0
#     with output but no version number is not retried: that is deterministic;
#   - records `attempts` and `read_exit` (124 = timed out), and for version_unread puts
#     the exit status and the first output line in `raw`, so the next one is diagnosable.
# A fast tool costs one call, as before. Version reading never fails the task by itself:
# version_unread is not a blocking status.
set -uo pipefail

if [[ $# -lt 4 ]]; then
    echo "usage: record_tool_versions.sh MODE ENV PROCESS OUT TOOL[:PACKAGE] ..." >&2
    exit 2
fi

MODE=$1
ENV_SPEC=$2
PROCESS=$3
OUT=$4
shift 4

if [[ "$MODE" != "strict" && "$MODE" != "warn" ]]; then
    echo "ERROR: record_tool_versions.sh: MODE must be strict or warn (got '$MODE')" >&2
    exit 2
fi

PREFIX=${CONDA_PREFIX:-}

# A non-negative integer from the environment, else the default.
setting() {
    local value=${!1:-}
    [[ "$value" =~ ^[0-9]+$ ]] && printf '%s' "$value" || printf '%s' "$2"
}
READ_ATTEMPTS=$(setting PKS_TOOL_VERSION_ATTEMPTS 3)
(( READ_ATTEMPTS >= 1 )) || READ_ATTEMPTS=1
READ_TIMEOUT=$(setting PKS_TOOL_VERSION_TIMEOUT 60)
(( READ_TIMEOUT >= 1 )) || READ_TIMEOUT=1
READ_BACKOFF=$(setting PKS_TOOL_VERSION_BACKOFF 2)
HAVE_TIMEOUT=0
command -v timeout >/dev/null 2>&1 && HAVE_TIMEOUT=1
TIMEOUT=()

# The pin for PACKAGE in the env yml, as written after the package name: "2.31", or ""
# when the yml does not pin it exactly (no version, or a range such as >=1.0).
pin_for() {
    local package=$1
    [[ -f "$ENV_SPEC" && "$ENV_SPEC" == *.y*ml ]] || return 0
    local spec
    spec=$(awk -v p="$package" '
        {
            line = $0
            sub(/#.*/, "", line)
            if (!match(line, /^[[:space:]]*-[[:space:]]*/)) next
            spec = substr(line, RLENGTH + 1)
            gsub(/[[:space:]]/, "", spec)
            sub(/^.*::/, "", spec)
            name = spec
            sub(/[=<>!~].*$/, "", name)
            if (name == p) { print substr(spec, length(name) + 1); exit }
        }' "$ENV_SPEC")
    # "=2.31" and "==2.31" are exact; "=2.31=h1234_0" carries a build after the version.
    case "$spec" in
        ==*) spec=${spec#==} ;;
        =*)  spec=${spec#=} ;;
        *)   return 0 ;;
    esac
    spec=${spec%%=*}
    spec=${spec%.\*}
    spec=${spec%\*}
    [[ "$spec" =~ ^[0-9A-Za-z._+-]+$ ]] && printf '%s' "$spec"
    return 0
}

# Tool-bump (Oct 2026): a few conda builds print a version other than their package
# version -- bracken 3.1p1 says "Bracken v3.0.1", because the upstream v3.1 tag never
# updated its VERSION string. scripts/tool_version_aliases.tsv lists each such build;
# finalize_provenance.py reads the same file, so the task and RUN_REPORT agree.
ALIASES="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tool_version_aliases.tsv"

# The version PACKAGE's binary is known to print at conda version PIN, or "".
alias_for() {
    [[ -n "$2" && -f "$ALIASES" ]] || return 0
    awk -F '\t' -v p="$1" -v v="$2" '!/^#/ && $1 == p && $2 == v { print $3; exit }' "$ALIASES"
}

# True when VERSION satisfies PIN: the pin, a more specific release of it, or its alias.
satisfies() {
    local version=$1 pin=$2 package=$3 alias
    [[ "$version" == "$pin" || "$version" == "$pin".* ]] && return 0
    alias=$(alias_for "$package" "$pin")
    [[ -n "$alias" && ( "$version" == "$alias" || "$version" == "$alias".* ) ]]
}

# What the tool says about its own version. Most take --version; these do not.
version_output() {
    local tool=$1
    case "$tool" in
        featureCounts)                  "${TIMEOUT[@]}" featureCounts -v ;;
        seqkit|diamond)                 "${TIMEOUT[@]}" "$tool" version ;;
        seqtk)                          "${TIMEOUT[@]}" seqtk ;;
        bracken)                        "${TIMEOUT[@]}" bracken -v ;;
        # MetaBAT2 (2.17 and 2.18) aborts on --version; -h prints "(version 2:2.18; ...)"
        # and the depth tool prints "jgi_summarize_bam_contig_depths 2.18 ..." with no args.
        metabat2)                       "${TIMEOUT[@]}" metabat2 -h ;;
        jgi_summarize_bam_contig_depths) "${TIMEOUT[@]}" jgi_summarize_bam_contig_depths ;;
        hmmsearch|hmmscan|nhmmscan|nhmmer|hmmpress|hmmbuild)
                                        "${TIMEOUT[@]}" "$tool" -h ;;
        *)                              "${TIMEOUT[@]}" "$tool" --version ;;
    esac
}

clean() { tr '\t\r\n' '   ' | cut -c1-200; }

# Paths first: bowtie2 prints its own absolute path, and an env named python3.10 would
# otherwise read as a version.
strip_paths() { sed -E 's#[^[:space:]]*/[^[:space:]]*##g'; }

# Ask TOOL for its version, retrying a failed call. Sets VERSION, RAW, ATTEMPTS, READ_EXIT.
read_version() {
    local tool=$1 attempt out rc
    VERSION="" RAW="" ATTEMPTS=0 READ_EXIT=""
    for (( attempt = 1; attempt <= READ_ATTEMPTS; attempt++ )); do
        TIMEOUT=()
        (( HAVE_TIMEOUT )) && TIMEOUT=(timeout -k 5 $(( READ_TIMEOUT * attempt )))
        # Both streams; sed reads to the end, so the tool never sees SIGPIPE.
        out=$(version_output "$tool" </dev/null 2>&1 | sed -n '1,200p'; exit "${PIPESTATUS[0]}")
        rc=$?
        ATTEMPTS=$attempt
        READ_EXIT=$rc
        out=$(printf '%s\n' "$out" | strip_paths)
        RAW=$(printf '%s\n' "$out" | grep -E '[0-9]+\.[0-9]+' | sed -n 1p || true)
        VERSION=$(printf '%s\n' "$RAW" | grep -oE '[0-9]+(\.[0-9]+)+' | sed -n 1p || true)
        [[ -n "$VERSION" ]] && return 0
        # Exit 0 with output but no version number: asking again gives the same answer.
        if (( rc == 0 )) && [[ -n "${out//[[:space:]]/}" ]]; then
            break
        fi
        (( attempt < READ_ATTEMPTS )) && sleep $(( READ_BACKOFF * attempt ))
    done
    local first
    first=$(printf '%s\n' "$out" | grep -vE '^[[:space:]]*$' | sed -n 1p || true)
    RAW="exit $READ_EXIT after $ATTEMPTS attempt(s)$( (( READ_EXIT == 124 )) && printf ' (timed out)')${first:+: $first}"
    return 0
}

printf 'process\ttool\tpackage\tstatus\tversion\tpinned\tpath\tconda_prefix\tenv\traw\tattempts\tread_exit\n' > "$OUT"

problems=()
for entry in "$@"; do
    tool=${entry%%:*}
    package=$tool
    [[ "$entry" == *:* ]] && package=${entry#*:}
    pinned=$(pin_for "$package")

    path=$(command -v "$tool" 2>/dev/null || true)
    version=""
    raw=""
    attempts=0
    read_exit=""
    if [[ -z "$path" ]]; then
        status=NOT_FOUND
    else
        read_version "$tool"
        version=$VERSION raw=$RAW attempts=$ATTEMPTS read_exit=$READ_EXIT
        if [[ -z "$PREFIX" ]]; then
            status=NO_ENV
        elif [[ "$path" != "$PREFIX"/* ]]; then
            status=OUTSIDE_ENV
        elif [[ -n "$pinned" && -n "$version" ]] && ! satisfies "$version" "$pinned" "$package"; then
            status=VERSION_MISMATCH
        elif [[ -z "$version" ]]; then
            status=version_unread
        elif [[ -z "$pinned" ]]; then
            status=unpinned
        else
            status=ok
        fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$PROCESS" "$tool" "$package" "$status" "$version" "$pinned" "$path" "$PREFIX" \
        "$ENV_SPEC" "$(printf '%s' "$raw" | clean)" "$attempts" "$read_exit" >> "$OUT"

    case "$status" in
        NOT_FOUND)        problems+=("$tool: not found on PATH") ;;
        NO_ENV)           problems+=("$tool: ran from $path (${version:-version unread}) with no conda environment active") ;;
        OUTSIDE_ENV)      problems+=("$tool: ran from $path (${version:-version unread}), outside this task's conda environment $PREFIX${pinned:+; pinned $package=$pinned}") ;;
        VERSION_MISMATCH) problems+=("$tool: reports $version from $path; $ENV_SPEC pins $package=$pinned") ;;
    esac
done

if (( ${#problems[@]} )); then
    label=ERROR
    [[ "$MODE" == "warn" ]] && label=WARNING
    {
        echo "$label: tool-version check, $PROCESS (--tool_version_check $MODE):"
        printf '  %s\n' "${problems[@]}"
        echo "  Something earlier on PATH shadows the pinned environment, or the environment is not"
        echo "  the one its yml describes. Launch Nextflow from a clean PATH (docs/hpc.md);"
        echo "  --tool_version_check warn records this and carries on."
    } >&2
    [[ "$MODE" == "strict" ]] && exit 78
fi
exit 0
