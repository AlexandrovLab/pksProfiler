# Community context, prophages and island mobility

Two questions this answers. Is the *pks* island itself mobile — able to move between bacteria?
And is the organism carrying it plausibly affecting its neighbours?

Both are addressed from assembled sequence, so this page assumes you have either draft genomes
from [genome-resolved analysis](mags.md) or contigs from [island reassembly](assembly.md).

## Why the neighbours matter

Colibactin damages DNA. Bacteria answer DNA damage with the **SOS response**, governed by *recA*
and *lexA* — and the same response can wake dormant prophages, bursting the cell. So a colibactin
producer sitting beside prophage-carrying neighbours is a testable hypothesis about what it is
doing to its community, rather than a coincidence.

That is why the community tables record, for every neighbouring organism, its prophage load and
whether it carries *recA*, *lexA*, and a *clbS*-like self-protection gene.

Colibactin damages DNA. Bacteria respond to DNA damage with the **SOS response**, controlled by
*recA* and *lexA* — and that same response can wake up dormant prophages, bursting the cell. So a
colibactin producer sitting next to prophage-carrying neighbours is a testable hypothesis about
what it is doing to its community.

This runs automatically with step 5; no extra flags. It looks at **every** draft genome in the
sample, not just the one carrying the island — the neighbours are the point.

**What you get**

```text
results/by_sample/<sample>/community/community_prophage_inventory.tsv   every prophage in every bin
results/by_sample/<sample>/community/pks_community_interactions.tsv     producer × each neighbour
results/by_sample/<sample>/community/pks_island_mobility.tsv            is the island itself mobile
```

`pks_community_interactions.tsv` is the one to read. Each row pairs the *clb*-positive bin with one
neighbour and records that neighbour's prophage count, whether it carries `recA` and `lexA`, and
whether it has a *clbS*-like gene (the self-protection gene — a neighbour with one may be immune).

## Where the results go

Everything on this page is written to its own tree, separate from the draft genomes in `mags/`:

```text
results/by_sample/<sample>/community/
├── <sample>/
│   ├── prophages/                                  prophage calls per draft genome
│   ├── genomic_context/                            genes flanking each clb hit
│   ├── <sample>.community_prophage_inventory.tsv   every prophage in every genome
│   ├── <sample>.pks_community_interactions.tsv     producer x each neighbour
│   └── <sample>.pks_island_mobility.tsv            is the island itself mobile
└── tumor_wgs/<sample>/
    ├── annotation/                                 genes called on the contigs
    ├── prophages/                                  prophage calls on the contigs
    ├── <sample>.contig_pks_context.tsv             what sits within the flanking window
    ├── <sample>.contig_pks_mobility.tsv            mobility flags, condensed
    └── contig_pks_summary.tsv             how much island each contig carries
```

`mags/` keeps the genomes themselves — bins, completeness, taxonomy and the per-sample summary.
The split is deliberate: the draft genomes are one product, what surrounds the island is another.

## Tuning

| Flag | Default | Effect |
|---|---|---|
| `--context_window_bp` | 50000 | how far either side of a *clb* hit to look for flanking genes |
| `--community_min_clb_genes` | 3 | *clb* genes a genome needs before it is treated as a producer |
| `--hmm_protein_evalue` | 1e-5 | threshold for calling a *clb* protein present |
| `--genomad_db` | — | required; without it no prophage calls are made |

---

Back to the [README](../../README.md).
