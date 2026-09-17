# Strain typing

Runs automatically whenever there are contigs or MAGs — no flag needed. It compares 7 housekeeping
genes to get a sequence type, then looks that ST up in a panel of 45,761 *E. coli* genomes to get
the clonal complex and phylogroup.

```text
results/by_sample/<sample>/strain/<unit>/strain_type.tsv
results/by_sample/<sample>/strain/strain_types.tsv     all units for the sample
```

```text
sample    unit    status  ST  clonal_complex  phylogroup  panel_pks_fraction  note
STOOL_01  bin.3   typed   73  ST73_Cplx       B2          0.6884
```

**How to read it.** `status` matters as much as the ST:

| status | Meaning |
|---|---|
| `typed` | all 7 genes found and the combination is known |
| `novel_allele_combination` | all 7 found, but the combination is not in the database |
| `insufficient_loci` | too little sequence to type — **no ST is reported rather than a guessed one** |

92% of colibactin-positive genomes in the reference panel are phylogroup **B2**, so a confident
non-B2 result gets a note asking you to double-check it. Turn typing off with
`--enable_strain_typing false`.

---

Back to the [README](../../README.md).
