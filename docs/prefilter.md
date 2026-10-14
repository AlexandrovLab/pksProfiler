# Optional read prefilter

`--prefilter_mode off | balanced | auto` (default `auto`)

In a metagenome most reads have nothing to do with this question. `balanced` first keeps only reads
from the relevant bacterial group, then goes back over the rejected ones and **rescues** anything
that still looks like a *clb* gene — so *clb* sequence in a read the classifier placed outside that
group is not thrown away by the first pass. The taxon each rescued read was assigned is recorded
in `prefilter/diamond_rescue_taxonomy.tsv` as read-level taxonomic support, not carrier attribution.
`auto` (the default) resolves to `balanced` for `--sample_type metagenome` whenever `--kraken_db` is
given, and to `off` otherwise (every tumour run, and a metagenome with no Kraken database). So a
metagenome run that passes `--kraken_db` for `--pks_taxa` / `--pks_community_taxa` is prefiltered
unless you add `--prefilter_mode off`. The resolved mode is logged at launch
(`prefilter_mode auto -> balanced`) and recorded in RUN_REPORT.txt.

It only changes what gets **counted**. Assembly still sees every host-depleted read.

---

Back to the [README](../README.md).
