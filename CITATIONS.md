# Citations

## pksProfiler

A manuscript is in preparation. Until then, cite this repository by URL and commit.

## The biological premise

- Pleguezuelos-Manzano C, Puschhof J, Rosendahl Huber A, *et al.*
  **Mutational signature in colorectal cancer caused by genotoxic pks+ *E. coli*.**
  *Nature* 580, 269–273 (2020). doi:10.1038/s41586-020-2080-8
- Dziubańska-Kusibab PJ, Berger H, Battistini F, *et al.*
  **Colibactin DNA-damage signature indicates mutational impact in colorectal cancer.**
  *Nature Medicine* 26, 1063–1069 (2020). doi:10.1038/s41591-020-0908-2

## Reference models and panels

The *clb* nucleotide and protein profile models shipped in `ref/hmm/` are built from the
Mäklin *et al.* collection of 45,761 *E. coli* genomes, Zenodo `10.5281/zenodo.13374348`.
`ref/hmm/PROVENANCE.md` records the exact build, the lineage-aware train/validation/test split,
and why exact deduplication was chosen over 99% clustering. The same collection supplies the
sequence-type, clonal-complex and phylogroup panel in `ref/typing/`.

## Tools

pksProfiler orchestrates the following. Please cite the ones your run used.

| Stage | Tool |
|---|---|
| workflow engine | **Nextflow** — Di Tommaso P, *et al.* *Nat Biotechnol* 35, 316–319 (2017). doi:10.1038/nbt.3820 |
| read quality filtering | **fastp** — Chen S, *et al.* *Bioinformatics* 34, i884–i890 (2018). doi:10.1093/bioinformatics/bty560 |
| human read removal | **Minimap2** — Li H. *Bioinformatics* 34, 3094–3100 (2018). doi:10.1093/bioinformatics/bty191 |
| alignment profiling | **Bowtie 2** — Langmead B, Salzberg SL. *Nat Methods* 9, 357–359 (2012). doi:10.1038/nmeth.1923 |
| read counting | **featureCounts** — Liao Y, *et al.* *Bioinformatics* 30, 923–930 (2014). doi:10.1093/bioinformatics/btt656 |
| BAM/CRAM handling | **SAMtools** — Danecek P, *et al.* *GigaScience* 10, giab008 (2021). doi:10.1093/gigascience/giab008 |
| profile-model search | **HMMER 3** — Eddy SR. *PLoS Comput Biol* 7, e1002195 (2011). doi:10.1371/journal.pcbi.1002195 |
| translated homology search | **DIAMOND** — Buchfink B, *et al.* *Nat Methods* 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x |
| assembly | **MEGAHIT** — Li D, *et al.* *Bioinformatics* 31, 1674–1676 (2015). doi:10.1093/bioinformatics/btv033 |
| assembly | **metaSPAdes** — Nurk S, *et al.* *Genome Res* 27, 824–834 (2017). doi:10.1101/gr.213959.116 |
| binning | **MetaBAT 2** — Kang DD, *et al.* *PeerJ* 7, e7359 (2019). doi:10.7717/peerj.7359 |
| genome completeness | **CheckM2** — Chklovski A, *et al.* *Nat Methods* 20, 1203–1212 (2023). doi:10.1038/s41592-023-01940-w |
| genome taxonomy | **GTDB-Tk** — Chaumeil P-A, *et al.* *Bioinformatics* 38, 5315–5316 (2022). doi:10.1093/bioinformatics/btac672 |
| gene annotation | **Prokka** — Seemann T. *Bioinformatics* 30, 2068–2069 (2014). doi:10.1093/bioinformatics/btu153 |
| prophage detection | **geNomad** — Camargo AP, *et al.* *Nat Biotechnol* 42, 1303–1312 (2024). doi:10.1038/s41587-023-01953-y |
| taxonomic classification | **KrakenUniq** — Breitwieser FP, *et al.* *Genome Biol* 19, 198 (2018). doi:10.1186/s13059-018-1568-1 |
| abundance estimation | **Bracken** — Lu J, *et al.* *PeerJ Comput Sci* 3, e104 (2017). doi:10.7717/peerj-cs.104 |
| strain typing | **mlst** — Seemann T. https://github.com/tseemann/mlst, using the PubMLST *E. coli* Achtman scheme (Jolley KA, *et al.* *Wellcome Open Res* 3, 124 (2018). doi:10.12688/wellcomeopenres.14826.1) |
| sequence clustering (model build) | **CD-HIT** — Fu L, *et al.* *Bioinformatics* 28, 3150–3152 (2012). doi:10.1093/bioinformatics/bts565 |
| alignment (model build) | **MAFFT** — Katoh K, Standley DM. *Mol Biol Evol* 30, 772–780 (2013). doi:10.1093/molbev/mst010 |
| gene prediction (model build) | **Prodigal** — Hyatt D, *et al.* *BMC Bioinformatics* 11, 119 (2010). doi:10.1186/1471-2105-11-119 |

DOIs were transcribed from the published records and should be verified before use in a
manuscript.
