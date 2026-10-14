nextflow.enable.dsl=2
// rc2 regression for the v1.0.0-rc1 Mayo metagenome crash (see tests/test_mag_tool_failure.py):
// drives pksMAG's own channel logic (magQualityFallback, magSummaryInput, magCommunityInput,
// from Modules/pks_mag.nf) with the shapes a failed-and-ignored gtdbtkClassify / checkm2Predict
// leaves behind. Everything is fabricated: no reads, no conda, no real tool.
//
//   s_ok          2 bins, CheckM2 + GTDB-Tk both produced output
//   s_gtdb_fail   2 bins (bin.18, bin.19), GTDB-Tk produced nothing   <- the rc1 crash
//   s_checkm_fail 1 bin, CheckM2 produced nothing
//   s_binless     no real bins, only an unbinned-pool call (binless stub path, unchanged)

include { magQualityFallback; magSummaryInput; magCommunityInput; noRealBinsSentinel } from '../../Modules/pks_mag.nf'
include { magCommunityInput as magCommunityInputRaw } from '../../Modules/pks_mag.nf'

params.outdir = "results"

def mk(String name, String text) {
    def f = file("${workflow.workDir}/fixture_inputs/${name}")
    f.parent.mkdirs()
    f.text = text
    return f
}

def perBin(Map bins, String kind) {
    bins.collect { s, ids -> tuple(s, ids.collect { mk("${s}/${it}.${kind}", "x\n") }) }
}

def baseName(p) {
    p instanceof Path ? p.name : p.toString()
}

workflow {
    def S = noRealBinsSentinel()
    def bins = [
        s_ok: ['bin.1', 'bin.2'], s_gtdb_fail: ['bin.18', 'bin.19'], s_checkm_fail: ['bin.3'],
    ]
    bins_ch = Channel.fromList(bins.collect { s, ids ->
        tuple(s, ids.collect { mk("${s}/${it}.fa", ">c\nACGT\n") }) })
    checkm2_ch = Channel.fromList(['s_ok', 's_gtdb_fail'].collect { s ->
        tuple(s, mk("${s}/quality_report.tsv",
            "Name\tCompleteness\tContamination\tGenome_Size\tContig_N50\n" +
            bins[s].collect { "${it}\t90\t1\t4600000\t50000\n" }.join(''))) })
    gtdbtk_ch = Channel.fromList(['s_ok', 's_checkm_fail'].collect { s ->
        tuple(s, mk("${s}/gtdbtk.bac120.summary.tsv",
            "user_genome\tclassification\n" + bins[s].collect { "${it}\td__Bacteria\n" }.join(''))) })

    tblouts_ch = Channel.fromList(perBin(bins, 'tblout'))
    gffs_ch = Channel.fromList(perBin(bins, 'gff'))
    virus_ch = Channel.fromList(perBin(bins, 'virus_summary.tsv'))
    bin_locus_ch = Channel.fromList(perBin(bins, 'locus_evidence.tsv'))
    contexts_ch = Channel.fromList([tuple('s_ok', [mk('s_ok/bin.1.context.tsv', "x\n")])])
    binless_locus = tuple('s_binless', [mk('s_binless/unbinned.locus_evidence.tsv', "x\n")])
    binless_checkm2 = tuple('s_binless', mk('s_binless/stub.checkm2.tsv', "Name\n"))
    binless_gtdbtk = tuple('s_binless', mk('s_binless/stub.gtdbtk.tsv', "user_genome\tclassification\n"))

    magQualityFallback(bins_ch, checkm2_ch, gtdbtk_ch)

    // Wired exactly as pksMAG wires them, sentinels included.
    checkm2_all = checkm2_ch.mix(Channel.of(binless_checkm2)).mix(magQualityFallback.out.checkm2_stub)
        .mix(Channel.value(tuple(S, "SENTINEL")))
    gtdbtk_all = gtdbtk_ch.mix(Channel.of(binless_gtdbtk)).mix(magQualityFallback.out.gtdbtk_stub)
        .mix(Channel.value(tuple(S, "SENTINEL")))
    tblouts_s = tblouts_ch.mix(Channel.value(tuple(S, [])))
    contexts_s = contexts_ch.mix(Channel.value(tuple(S, [])))
    locus_s = bin_locus_ch.mix(Channel.of(binless_locus)).mix(Channel.value(tuple(S, [])))

    magSummaryInput(checkm2_all, gtdbtk_all, tblouts_s, contexts_s, locus_s)
    magCommunityInput(gtdbtk_ch.mix(magQualityFallback.out.gtdbtk_stub), gffs_ch, tblouts_s, virus_ch, bin_locus_ch)

    // rc1's exact crash shape: community join WITHOUT the fallback (raw GTDB-Tk output only).
    // s_gtdb_fail is left-orphaned here; it must be skipped with a warning, not crash.
    magCommunityInputRaw(gtdbtk_ch, gffs_ch, tblouts_s, virus_ch, bin_locus_ch)

    magSummaryInput.out.rows
        .map { s, c, g, t, x, l -> "summary\t${s}\t${baseName(c)}\t${baseName(g)}\t${t.size()}\t${x.size()}\t${l.size()}" }
        .mix(magCommunityInput.out.rows
            .map { s, g, gf, t, v, l -> "community\t${s}\t${baseName(g)}\t${gf.size()}\t${t.size()}\t${v.size()}\t${l.size()}" })
        .mix(magCommunityInputRaw.out.rows.map { row -> "community_raw\t${row[0]}\t${row.size()}" })
        .mix(magQualityFallback.out.gtdbtk_stub.map { s, f -> "gtdbtk_stub\t${s}\t" + f.text.readLines().drop(1).join(';') })
        .mix(magQualityFallback.out.checkm2_stub.map { s, f -> "checkm2_stub\t${s}\t" + f.text.readLines().drop(1).join(';') })
        .collectFile(name: 'channel_rows.tsv', storeDir: params.outdir, newLine: true, sort: true)
}
