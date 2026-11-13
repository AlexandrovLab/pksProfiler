#!/usr/bin/env bash
# The Bracken -r read length, decided the same way everywhere.
#
# Bracken's abundance model is built per read length: a database ships one
# database<N>mers.kmer_distrib per N it was built for, and -r selects among them. One
# global value cannot fit a multi-cohort run -- in the 36-cohort infant MAG set only 11
# cohorts record a read length in the archive at all, and the implied single-read lengths
# span ~96-148 bp. So rather than guess, measure the median length of the very records
# the task is about to classify and snap it to whichever N the database actually ships.
# That adapts to any Bracken database instead of a hard-coded ladder.
#
# v3 review item 13: "Centralize Bracken read-length determination across island and
# community taxonomy paths." Modules/pks_taxa.nf (process Bracken, the island path) did
# all of this inline. Modules/pks_prefilter.nf (process sampleBracken, the community
# path) passed `-r "${params.bracken_read_length}"` straight through, and that parameter
# defaults to null -- so enabling --pks_community_taxa without setting it handed Bracken
# an empty -r. The failure mode is a dead task rather than a wrong abundance, which is
# why it explains no shipped result, but the two paths should not disagree about
# something this basic. This is that logic, once, for both.
#
# It lives in a shell script rather than a lib/*.groovy helper on purpose. A Groovy
# helper would have to emit shell through a GString, and a `\$` that is correct inside a
# module's script block arrives in the task script as a literal backslash-dollar, which
# bash reads as an escaped dollar rather than a command substitution. lib/ToolVersions
# avoids the same trap the same way: return a command line, keep the shell in a file.
#
# The two callers differ only in where the lengths come from:
#   fastq   the reads themselves (island path): sequence lines of a FASTQ
#   kraken  the KrakenUniq output (community path), whose column 4 is the sequence
#           length -- that process is given no FASTQ, only the report and output. A
#           paired record reads "a|b"; the first mate's length is the single-read length
#           Bracken's model is indexed by.
#
# Fallback order, identical for both:
#   1. the measured median, snapped to the nearest N the database ships
#   2. --bracken_read_length, when the record set is empty or unmeasurable
#   3. the smallest N the database ships, so a thin sample still yields a call rather
#      than dying on an empty -r
#
# Bounded cost: at most 20k records are inspected, ample for a median and trivial beside
# the classification it accompanies.
#
# Prints the chosen length on stdout; progress and the reasoning go to stderr so the
# caller can capture the value with $(...).
#
#     BRACKEN_RL=$(bash bracken_read_length.sh <db> <fallback|""> <fastq|kraken> <file>)

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

[[ $# -eq 4 ]] || die "usage: bracken_read_length.sh <db> <fallback|\"\"> <fastq|kraken> <file>"

DB=$1
FALLBACK=$2
SOURCE=$3
FILE=$4

[[ -d "$DB" ]] || die "not a database directory: $DB"
case "$SOURCE" in
    fastq|kraken) ;;
    *) die "unknown source '$SOURCE'; expected fastq or kraken" ;;
esac

# The read lengths this database was built for, ascending.
#
# `|| true`: with set -e and pipefail, a database shipping no kmer_distrib makes ls exit
# 2 and kills the script right here -- before the explanatory error at the bottom ever
# runs, so the user gets a bare exit 2 and no reason. An empty AVAIL is a state this
# script is built to report on, not a reason to die silently.
AVAIL=$(ls "$DB"/database*mers.kmer_distrib 2>/dev/null \
        | sed -E 's|.*/database([0-9]+)mers\.kmer_distrib|\1|' \
        | sort -n | tr '\n' ' ' || true)

median_of() {
    sort -n | awk '{v[NR]=$1}
                   END {if (NR) print (NR%2 ? v[(NR+1)/2] : int((v[NR/2]+v[NR/2+1])/2))}'
}

measure() {
    case "$SOURCE" in
        fastq)
            awk 'NR%4==2 {print length($0); n++} n>=20000 {exit}' "$FILE"
            ;;
        kraken)
            # Column 4 is the sequence length; "a|b" for a pair. Skip anything
            # non-numeric rather than guessing at it.
            awk -F '\t' '{ split($4, L, "|")
                           if (L[1] ~ /^[0-9]+$/ && L[1]+0 > 0) { print L[1]+0; n++ } }
                         n>=20000 {exit}' "$FILE"
            ;;
    esac
}

BRACKEN_RL=""
if [[ -s "$FILE" && -n "$AVAIL" ]]; then
    MEASURED=$(measure | median_of || true)
    if [[ -n "${MEASURED:-}" && "$MEASURED" -gt 0 ]]; then
        # Snap to the nearest length the database actually ships.
        BRACKEN_RL=$(awk -v m="$MEASURED" \
            'BEGIN {best=""; bd=1e9}
             {for (i=1; i<=NF; i++) {
                d = ($i > m) ? $i - m : m - $i
                if (d < bd) { bd = d; best = $i }
              }}
             END {print best}' <<< "$AVAIL")
        echo "bracken read length: measured ${MEASURED} bp -> using ${BRACKEN_RL} (available: $AVAIL)" >&2
    fi
fi

if [[ -z "$BRACKEN_RL" ]]; then
    BRACKEN_RL="$FALLBACK"
    if [[ -z "$BRACKEN_RL" ]]; then
        BRACKEN_RL=$(awk '{print $1}' <<< "$AVAIL")
        echo "bracken read length: not measurable and no --bracken_read_length; using ${BRACKEN_RL:-<none>}" >&2
    else
        echo "bracken read length: not measurable; falling back to --bracken_read_length=$BRACKEN_RL" >&2
    fi
fi

[[ -n "$BRACKEN_RL" ]] || die "no Bracken read length could be determined: $DB ships no
       database<N>mers.kmer_distrib and --bracken_read_length was not given"

printf '%s\n' "$BRACKEN_RL"
