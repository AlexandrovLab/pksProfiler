#!/usr/bin/env bash
# Execute lib/HostDepletion.groovy -- the table Modules/map_reads.nf builds every
# host-depletion minimap2 call from -- and compare what it returns with the expected
# preset, index and flags for each sample type. Runs no task and reads no data.
#
# Ludmil, 26 Sep mode memo section 4: RNA's human pass must be splice-aware
# (-ax splice:sr, optional -j junctions) while the T2T+phiX and pangenome passes stay
# -ax sr; E0's --frag=no must reach every short-read call. The unit tests read the
# Groovy source; this runs it, through the same Nextflow the pipeline uses.
#
#     bash tests/check_depletion_plan.sh

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
work=$(mktemp -d "${TMPDIR:-/tmp}/pks-depletion-plan.XXXXXX")
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/lib"
cp "$repo_dir/lib/HostDepletion.groovy" "$work/lib/"

cat > "$work/main.nf" <<'NF'
def show(label, m) {
    println "${label}|${m.preset}|${m.index ?: ''}|${m.args}"
}

workflow {
    // sample_type, frag_mode, hg38_db, rna_splice_db, rna_junctions_bed
    def cases = [
        ['wgs',        'default', 'HG38', null,     null],
        ['wgs',        'no',      'HG38', null,     null],
        ['wes',        'default', 'HG38', null,     null],
        ['metagenome', 'no',      'HG38', 'SPLICE', 'JUNC'],
        ['rna',        'default', 'HG38', 'SPLICE', null],
        ['rna',        'no',      null,   'SPLICE', 'JUNC'],
    ]
    cases.each { c ->
        show("human:${c[0]}:${c[1]}", HostDepletion.humanPass(c[0], c[1], c[2], c[3], c[4]))
    }
    show('t2t:default', HostDepletion.t2tPass('default', 'T2T'))
    show('t2t:no', HostDepletion.t2tPass('no', 'T2T'))
    show('pangenome:no', HostDepletion.pangenomePass('no'))
    println "param:wgs|${HostDepletion.humanIndexParam('wgs')}"
    println "param:rna|${HostDepletion.humanIndexParam('rna')}"
}
NF

cat > "$work/expected.txt" <<'TXT'
human:wgs:default|sr|HG38|
human:wgs:no|sr|HG38|--frag=no
human:wes:default|sr|HG38|
human:metagenome:no|sr|HG38|--frag=no
human:rna:default|splice:sr|SPLICE|
human:rna:no|splice:sr|SPLICE|--frag=no -j "JUNC"
t2t:default|sr|T2T|
t2t:no|sr|T2T|--frag=no
pangenome:no|sr||--frag=no
param:wgs|hg38_db
param:rna|rna_splice_db
TXT

(cd "$work" && nextflow -q run main.nf -work-dir "$work/w" > "$work/actual.raw" 2>&1) || {
    cat "$work/actual.raw" >&2; echo "ERROR: could not run lib/HostDepletion.groovy" >&2; exit 1; }
grep -E '^(human|t2t|pangenome|param):' "$work/actual.raw" > "$work/actual.txt" || true

if ! diff -u "$work/expected.txt" "$work/actual.txt"; then
    echo "ERROR: lib/HostDepletion.groovy returned a different depletion plan (diff above)" >&2
    exit 1
fi
echo "depletion plan: $(wc -l < "$work/expected.txt") cases as expected"
