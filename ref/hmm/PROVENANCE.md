# clb HMM profiles shipped with pksProfiler

## In use by default (v0.0.2)

| File | Used by | Alphabet | Models |
|---|---|---|---|
| `clb_population_protein_exact_v1.hmm` | `hmmsearchClb`, tumour contig HMM search (`hmmsearch`) | amino | 19 |
| `clb_population_dna_exact_v1.hmm` (+ `.h3f .h3i .h3m .h3p`) | `pksProfilerHMM` (`nhmmscan`) | DNA | 19 |

`hmmsearch` reads a plain `.hmm`, so the protein profiles ship without pressed indexes.
`nhmmscan` reads the pressed index, so the nucleotide profiles ship with all four.

SHA256, also in [`audit/SHA256SUMS`](audit/SHA256SUMS) for `sha256sum -c`:

```
6d68ce28736c40408c3d368fd64567f04262012bca0cf00b27a4bcfca1b0272c  clb_population_protein_exact_v1.hmm
95242c9f4ef86030157c8c0c89817915a2fc5df25420d1f8bb6c1e574132627c  clb_population_dna_exact_v1.hmm
```

## Where they came from

Full build record: [`BUILD_COMMANDS_AND_STEPS.txt`](BUILD_COMMANDS_AND_STEPS.txt) in this
directory. Read the note at its head first: the shipped models are the `exact_v1` rebuild
described in its final section, not the 99%-clustered original described in its main section.
The seed HMM and build scripts are not included in this repository.

Source collection: Lawley/Mäklin 45,761-genome AGC archive, Zenodo `10.5281/zenodo.13374348`
(2024 Lancet Microbe). Assemblies count as pks-positive at >=17/19 island genes. The 4,760
pks-positive genomes were split **without lineage leakage** into 2,868 train / 964 validation /
928 test. Only the training partition fit these models.

Per gene: Prodigal 2.6.3 `-p single`; HMMER 3.3.2 candidates at i-E-value <=1e-10 and >=70%
model coverage; best candidate per clbA–clbS; genomes admitted at >=17/19 recovered; paired
protein and CDS retained; MAFFT 7.526 alignment; `hmmbuild`; concatenate; `hmmpress`.

## Why `exact_v1` and not the original build

`exact_v1` removes only **100%-identical** duplicates, after excluding sequences outside
80–120% of each gene's modal length as likely truncations or misassemblies. The original build
clustered at **99% identity** with CD-HIT/CD-HIT-EST instead.

On an island this conserved, 99% clustering was far too aggressive — for clbA (244 aa) it
tolerates roughly two amino acids, so genuine variants were merged away:

Protein sequences (the DNA build retains more, 1,418 at a median of 41 per gene):

| | raw sequences | retained |
|---|---:|---:|
| exact_v1 (100% dedup) | 52,865 | **844** |
| original (99% clustering) | 52,865 | **65** |

Per-gene training depth (protein `NSEQ`): exact_v1 median 23, minimum 7, **no** single-sequence
models. The 99% build had a median of 2 and built clbG, clbL, clbM, clbO and clbP from a
**single sequence** each; its DNA counterpart had eight such genes. A profile HMM trained on one
sequence cannot model positional variation, which is the reason to use an HMM at all.

Exact deduplication still achieves what the clustering was for — suppressing abundant ST73 and
ST95 clones — because identical copies collapse regardless. Clustering additionally discarded
real diversity.

Per-gene audit tables: [`audit/`](audit) — `raw`, `modal_length`, `retained_exact_unique`,
`rejected_length` and `rejected_exact_duplicate` per gene, for both alphabets.

## Retained, not in use

`clb_all_dna.hmm` (+ indexes) is the **v0.0.1 benchmark model** and is deliberately left
untouched — the reference-build README treats it as an immutable benchmark until sensitivity,
specificity, locus-reconstruction and MAG-fragment tests are complete. Select it with
`--hmm_model ref/hmm/clb_all_dna.hmm` to reproduce v0.0.1 numbers.

## Not yet calibrated

Neither profile set carries GA/TC/NC cutoffs, so detection depends entirely on `--hmm_evalue`
(1e-10, nucleotide) and `--hmm_protein_evalue` (1e-5, protein). The 964 validation and 928 test
genomes are still reserved and unused. Treat hit counts as uncalibrated evidence.

Two open questions from the audit tables, worth resolving before calibration:

- **clbJ** lost 687 sequences to the length filter, far more than any other gene — consistent
  with its length fragmenting in assembly, but worth confirming.
- **clbK** yielded only 1,731 raw sequences where every other gene yielded ~2,840. Roughly 1,100
  training genomes produced no clbK call. Genuine absence, or a gap in extraction?
