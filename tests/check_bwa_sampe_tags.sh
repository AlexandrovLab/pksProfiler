#!/usr/bin/env bash
# Does `bwa aln` + `sampe` emit an AS tag? scripts/header_provenance.py's decision for a
# non-`mem` bwa @PG rests on the answer being no.
#
# Ludmil's v3 review, V2: "The present fallback partly relies on the premise that
# unconfirmable legacy `bwa aln` / `samse` / `sampe` inputs do not provide an AS tag
# compatible with the BWA-MEM rescue threshold. Confirm the actual SAM tags emitted by
# those subcommands before deciding whether `unconfirmed` should fail or resolve to
# rescue."
#
# This confirms it by running the aligner rather than citing the manual, on synthetic
# sequence with no reference to any sample. If a future bwa ever starts writing AS from
# sampe, off_auto_no_as_tag stops being the only correct answer and this test fails,
# which is the point of keeping it.
#
#     bash tests/check_bwa_sampe_tags.sh
set -uo pipefail

BWA=${BWA_BIN:-$(command -v bwa || true)}
if [[ -z "$BWA" || ! -x "$BWA" ]]; then
    echo "SKIP: no bwa on PATH and BWA_BIN unset (set BWA_BIN=/path/to/bwa to run this)"
    exit 0
fi

work=$(mktemp -d "${TMPDIR:-/tmp}/bwa-sampe-tags.XXXXXX")
trap 'rm -rf "$work"' EXIT
cd "$work"

python3 - <<'PY'
import random
random.seed(42)
seq = "".join(random.choice("ACGT") for _ in range(20000))
with open("ref.fa", "w") as fh:
    fh.write(">synthetic_ref\n")
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i+60] + "\n")
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
with open("r1.fq", "w") as f1, open("r2.fq", "w") as f2:
    # the last two pairs carry mismatches, so a scoring scheme would have something to say
    for n, (start, mm) in enumerate([(1000, 0), (5000, 0), (9000, 4), (13000, 8)]):
        a = list(seq[start:start+100])
        for k in range(mm):
            a[10 + k*7] = {"A": "C", "C": "A", "G": "T", "T": "G"}[a[10 + k*7]]
        f1.write(f"@pair{n}/1\n{''.join(a)}\n+\n{'I'*100}\n")
        f2.write(f"@pair{n}/2\n{rc(seq[start+200:start+300])}\n+\n{'I'*100}\n")
PY

"$BWA" index ref.fa >/dev/null 2>&1
"$BWA" aln ref.fa r1.fq > r1.sai 2>/dev/null
"$BWA" aln ref.fa r2.fq > r2.sai 2>/dev/null
"$BWA" sampe ref.fa r1.sai r2.sai r1.fq r2.fq > sampe.sam 2>/dev/null
"$BWA" mem ref.fa r1.fq r2.fq > mem.sam 2>/dev/null

sampe_as=$(grep -v '^@' sampe.sam | grep -c 'AS:i:' || true)
sampe_n=$(grep -vc '^@' sampe.sam || true)
mem_as=$(grep -v '^@' mem.sam | grep -c 'AS:i:' || true)
mem_n=$(grep -vc '^@' mem.sam || true)
echo "  bwa: $("$BWA" 2>&1 | awk '/^Version/{print $2}')"
echo "  sampe records with AS:i:  $sampe_as / $sampe_n"
echo "  mem   records with AS:i:  $mem_as / $mem_n"

fail=0
if [[ "$sampe_as" -ne 0 ]]; then
    echo "FAIL: sampe emitted an AS tag on $sampe_as records." >&2
    echo "      header_provenance.py returns off_auto_no_as_tag for a non-mem bwa @PG on" >&2
    echo "      the premise that there is no AS to rescue on. That premise no longer holds;" >&2
    echo "      revisit low_score_decision before trusting the decision for those inputs." >&2
    fail=1
else
    echo "ok: sampe emits no AS tag, so off_auto_no_as_tag is the only possible answer"
fi
if [[ "$mem_as" -eq 0 || "$mem_n" -eq 0 ]]; then
    echo "FAIL: mem emitted no AS tag either -- the test fixture is not exercising bwa properly." >&2
    fail=1
else
    echo "ok: mem emits AS on every record, so the two paths are genuinely different"
fi
exit $fail
