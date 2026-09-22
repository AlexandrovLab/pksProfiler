# Optional read prefilter

`--prefilter_mode off | balanced | auto` (default `auto`)

In a metagenome most reads have nothing to do with this question. `balanced` first keeps only reads
from the relevant bacterial group, then goes back over the rejected ones and **rescues** anything
that still looks like a *clb* gene — so an unusual carrier is not thrown away by the first pass.
`auto` switches this on for metagenomes when a Kraken database is available, and leaves tumour runs
alone.

It only changes what gets **counted**. Assembly still sees every host-depleted read.

---

Back to the [README](../README.md).
