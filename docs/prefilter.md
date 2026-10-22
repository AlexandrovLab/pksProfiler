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
(`prefilter_mode auto -> balanced`).

It only changes what gets **counted**. Assembly still sees every host-depleted read.

## Which mode actually ran

`auto` is a request, not a mode anything runs in, so every output records the **resolved**
value as well — `off` or `balanced`, never `auto`. Three places, and they cannot disagree,
because all three read the one value `main.nf` resolves at launch
(`params.prefilter_mode_resolved`), which is also what the channel graph is built from:

| Where | What you see |
|---|---|
| `RUN_REPORT.txt`, PROFILING section | `prefilter_mode` (requested, e.g. `auto`) and `prefilter_mode_resolved` (`off` / `balanced`) on the next line |
| `runs/provenance/<session>.json` | `invocation.params.prefilter_mode_resolved` |
| `cohort/cohort_qc_summary.tsv` | a `prefilter_mode` column, per sample, written by whichever profiling lane ran |

Per sample as well as per run, because a resumed results tree can be the product of
several sessions: the cohort table is then the only place that says, sample by sample,
whether the *clb* counts beside it came from the complete host-depleted stream or from a
narrowed subset of it. A sample whose profiling never ran shows `NA`; `auto` in that
column is refused outright, since nothing could be reproduced from it.

Raised as rc2 review item 2.1 — before this, only the request was recorded, so two runs
with byte-identical recorded parameters could have profiled different read sets.

---

Back to the [README](../README.md).
