# Genome-resolved analysis (metagenomes)

For a metagenome the question is not how many reads carry *clb* genes, but **which organism**
carries them — and what that organism's neighbours look like.

## Recovering draft genomes

For a metagenome the interesting question is not *how many reads* but *which organism*. So the
route is different: assemble everything, sort the contigs into draft genomes, then ask which of
those genomes carries *clb* genes.

Needs three databases.

```bash
nextflow run main.nf -profile local \
  --sample samples.csv \
  --input_data_type fastq \
  --sample_type metagenome \
  --enable_mags true \
  --gtdbtk_db    /db/gtdbtk_r232 \
  --checkm2_db   /db/checkm2/uniref100.KO.1.dmnd \
  --genomad_db   /db/genomad_db_v1.9 \
  --kraken_db    /db/krakenuniq \
  --hg38_db      /refs/human-GRC-db.mmi \
  --t2t_phix_db  /refs/human-GCA-phix-db.mmi \
  --outdir       results
```

What happens, in plain terms: everything in the sample is assembled at once, the resulting contigs
are sorted into groups that came from the same organism, each group is scored for how complete a
genome it represents and assigned a species name, and then every group is tested for *clb* genes.
So you end up with a list of draft genomes and know which of them carries the island.

**What you get**

```text
results/by_sample/<sample>/genomes/
├── bins/                        the draft genomes
├── checkm2/                     how complete and how contaminated each one is
├── gtdbtk/                      which species each one is
└── <sample>.mag_summary.tsv     one row per draft genome: species, quality, clb genes found
```

**How to read it.** Open `mag_summary.tsv` and look for a draft genome with *clb* genes in it. Then
check its completeness score before you trust the species name — a genome that is only 40% complete
is a hint, not an identification. **If a sample yields zero contigs, that is reported as a genuine negative**, not
an error.

Once draft genomes are recovered, [community and prophage context](community_context.md) describes what surrounds the island.

---

Back to the [README](../../README.md).
