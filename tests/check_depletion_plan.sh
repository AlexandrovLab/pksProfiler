#!/usr/bin/env bash
# Execute lib/HostDepletion.groovy -- the table Modules/map_reads.nf builds every
# host-depletion minimap2 call from -- and compare what it returns with the expected
# preset, index and flags for each sample type. Runs no task and reads no data.
#
# Ludmil, 26 Sep mode memo section 4: RNA's human pass must be splice-aware
# (-ax splice:sr, optional -j junctions) while the T2T+phiX and pangenome passes stay
# unspliced; E0's --frag=no must reach every short-read call. Owner decision, 29 Sep:
# ont/hifi use map-ont/map-hifi (splice/splice:hq for rna's human pass) against the
# *_long indexes, and never --frag. The unit tests read the
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
    // sample_type, read_type, frag_mode
    def cases = [
        ['wgs', 'short', 'default'], ['wgs', 'short', 'no'], ['wes', 'short', 'default'],
        ['metagenome', 'short', 'no'], ['rna', 'short', 'default'], ['rna', 'short', 'no'],
        ['wgs', 'ont', 'no'], ['wgs', 'hifi', 'default'], ['rna', 'ont', 'no'], ['rna', 'hifi', 'default'],
    ]
    cases.each { c ->
        def label = "${c[0]}:${c[1]}:${c[2]}"
        // hg38, hg38_long, splice, splice_long, junctions (junctions only reach rna)
        show("human:${label}", HostDepletion.humanPass(c[0], c[1], c[2],
            'HG38', 'HG38_LONG', 'SPLICE', 'SPLICE_LONG', c[2] == 'no' ? 'JUNC' : null))
        show("t2t:${label}", HostDepletion.t2tPass(c[1], c[2], 'T2T', 'T2T_LONG'))
        show("pangenome:${label}", HostDepletion.pangenomePass(c[1], c[2]))
        println "param:${label}|${HostDepletion.humanIndexParam(c[0], c[1])}|${HostDepletion.t2tIndexParam(c[1])}|${HostDepletion.profilingPreset(c[1])}"
    }
    println "param:default|${HostDepletion.humanIndexParam('wgs', null)}|${HostDepletion.t2tIndexParam(null)}"
}
NF

cat > "$work/expected.txt" <<'TXT'
human:wgs:short:default|sr|HG38|
t2t:wgs:short:default|sr|T2T|
pangenome:wgs:short:default|sr||
param:wgs:short:default|hg38_db|t2t_phix_db|sr
human:wgs:short:no|sr|HG38|--frag=no
t2t:wgs:short:no|sr|T2T|--frag=no
pangenome:wgs:short:no|sr||--frag=no
param:wgs:short:no|hg38_db|t2t_phix_db|sr
human:wes:short:default|sr|HG38|
t2t:wes:short:default|sr|T2T|
pangenome:wes:short:default|sr||
param:wes:short:default|hg38_db|t2t_phix_db|sr
human:metagenome:short:no|sr|HG38|--frag=no
t2t:metagenome:short:no|sr|T2T|--frag=no
pangenome:metagenome:short:no|sr||--frag=no
param:metagenome:short:no|hg38_db|t2t_phix_db|sr
human:rna:short:default|splice:sr|SPLICE|
t2t:rna:short:default|sr|T2T|
pangenome:rna:short:default|sr||
param:rna:short:default|rna_splice_db|t2t_phix_db|sr
human:rna:short:no|splice:sr|SPLICE|--frag=no -j "JUNC"
t2t:rna:short:no|sr|T2T|--frag=no
pangenome:rna:short:no|sr||--frag=no
param:rna:short:no|rna_splice_db|t2t_phix_db|sr
human:wgs:ont:no|map-ont|HG38_LONG|
t2t:wgs:ont:no|map-ont|T2T_LONG|
pangenome:wgs:ont:no|map-ont||
param:wgs:ont:no|hg38_db_long|t2t_phix_db_long|map-ont
human:wgs:hifi:default|map-hifi|HG38_LONG|
t2t:wgs:hifi:default|map-hifi|T2T_LONG|
pangenome:wgs:hifi:default|map-hifi||
param:wgs:hifi:default|hg38_db_long|t2t_phix_db_long|map-hifi
human:rna:ont:no|splice|SPLICE_LONG|-j "JUNC"
t2t:rna:ont:no|map-ont|T2T_LONG|
pangenome:rna:ont:no|map-ont||
param:rna:ont:no|rna_splice_db_long|t2t_phix_db_long|map-ont
human:rna:hifi:default|splice:hq|SPLICE_LONG|
t2t:rna:hifi:default|map-hifi|T2T_LONG|
pangenome:rna:hifi:default|map-hifi||
param:rna:hifi:default|rna_splice_db_long|t2t_phix_db_long|map-hifi
param:default|hg38_db|t2t_phix_db
TXT

(cd "$work" && nextflow -q run main.nf -work-dir "$work/w" > "$work/actual.raw" 2>&1) || {
    cat "$work/actual.raw" >&2; echo "ERROR: could not run lib/HostDepletion.groovy" >&2; exit 1; }
grep -E '^(human|t2t|pangenome|param):' "$work/actual.raw" > "$work/actual.txt" || true

if ! diff -u "$work/expected.txt" "$work/actual.txt"; then
    echo "ERROR: lib/HostDepletion.groovy returned a different depletion plan (diff above)" >&2
    exit 1
fi
echo "depletion plan: $(wc -l < "$work/expected.txt") cases as expected"
