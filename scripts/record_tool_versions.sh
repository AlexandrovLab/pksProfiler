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
TIMEOUT=()
if command -v timeout >/dev/null 2>&1; then
    TIMEOUT=(timeout 60)
fi

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

# What the tool says about its own version. Most take --version; these do not.
version_output() {
    local tool=$1
    case "$tool" in
        featureCounts)                  "${TIMEOUT[@]}" featureCounts -v ;;
        seqkit|diamond)                 "${TIMEOUT[@]}" "$tool" version ;;
        seqtk)                          "${TIMEOUT[@]}" seqtk ;;
        bracken)                        "${TIMEOUT[@]}" bracken -v ;;
        hmmsearch|hmmscan|nhmmscan|nhmmer|hmmpress|hmmbuild)
                                        "${TIMEOUT[@]}" "$tool" -h ;;
        *)                              "${TIMEOUT[@]}" "$tool" --version ;;
    esac
}

clean() { tr '\t\r\n' '   ' | cut -c1-200; }

printf 'process\ttool\tpackage\tstatus\tversion\tpinned\tpath\tconda_prefix\tenv\traw\n' > "$OUT"

problems=()
for entry in "$@"; do
    tool=${entry%%:*}
    package=$tool
    [[ "$entry" == *:* ]] && package=${entry#*:}
    pinned=$(pin_for "$package")

    path=$(command -v "$tool" 2>/dev/null || true)
    version=""
    raw=""
    if [[ -z "$path" ]]; then
        status=NOT_FOUND
    else
        # Paths first: bowtie2 prints its own absolute path, and an env named
        # python3.10 would otherwise read as a version.
        raw=$(version_output "$tool" </dev/null 2>&1 | sed -E 's#[^[:space:]]*/[^[:space:]]*##g' |
              grep -E '[0-9]+\.[0-9]+' | sed -n 1p || true)
        version=$(printf '%s\n' "$raw" | grep -oE '[0-9]+(\.[0-9]+)+' | sed -n 1p || true)
        if [[ -z "$PREFIX" ]]; then
            status=NO_ENV
        elif [[ "$path" != "$PREFIX"/* ]]; then
            status=OUTSIDE_ENV
        elif [[ -n "$pinned" && -n "$version" && "$version" != "$pinned" && "$version" != "$pinned".* ]]; then
            status=VERSION_MISMATCH
        elif [[ -z "$version" ]]; then
            status=version_unread
        elif [[ -z "$pinned" ]]; then
            status=unpinned
        else
            status=ok
        fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$PROCESS" "$tool" "$package" "$status" "$version" "$pinned" "$path" "$PREFIX" \
        "$ENV_SPEC" "$(printf '%s' "$raw" | clean)" >> "$OUT"

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
