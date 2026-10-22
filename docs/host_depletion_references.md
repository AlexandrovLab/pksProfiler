# What is in the host-depletion references

`--hg38_db` and `--t2t_phix_db` are the two indexes every run depletes against before
anything is counted, so **everything pksProfiler reports is what these two did not
remove**. Until now the only thing recorded about them anywhere was a path and a size.

rc2 review item 2.7 (Ludmil's 27 Aug ask, reopened): *"Record exactly which human
assemblies/decoys/phiX are present in `human-GRC-db.mmi` / `human-GCA-phix-db.mmi` and
the build command used."*

Everything below was read out of the index files themselves on **2026-10-22** with

```
python scripts/describe_minimap2_index.py <index.mmi> --tsv <dictionary.tsv>
```

which parses the header and the complete sequence dictionary minimap2 writes in the
first bytes of a `.mmi` (see that script's docstring for the format). It is exact: these
are the sequence names and lengths the index carries, not an inference from its filename.

Both files live under
`/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/`, both have mtime
**2025-01-19**, and both are rows in [`conf/databases.tsv`](../conf/databases.tsv).

---

## `human-GRC-db.mmi` — `--hg38_db`

| | |
|---|---|
| Size | 9,147,467,965 bytes |
| Sequences | **705** |
| Total bases | **3,298,430,636** |
| Naming | RefSeq accessions (`NC_`, `NT_`, `NW_`) |
| minimap2 `k` / `w` / `b` | 21 / 11 / 14 |
| `flag` | 0 — no `-H`, sequences stored, names kept |
| Dictionary sha256 | `a06386dc8b08f9547f6cbb80afab223f92161093e1355b367a59166054a010c9` (first 12,301 bytes) |

**Assembly: GRCh38.p14, the RefSeq full assembly (GCF_000001405.40).** The sequence count
and the total length are the published figures for that assembly exactly — 705 scaffolds,
3,298,430,636 bp — and the accessions are its own.

| Class | Sequences | Notes |
|---|---|---|
| `NC_000001.11` … `NC_000024.10` | 24 | chr1–chr22, X, Y |
| `NC_012920.1` | 1 | mitochondrion (rCRS) |
| `NT_…` | 355 | unlocalized / unplaced scaffolds, patch and ALT scaffolds |
| `NW_…` | 325 | unlocalized / unplaced scaffolds, patch and ALT scaffolds |

What is **not** in it:

* **no hs38d1 decoy contigs** (`*_decoy`, `chrUn_JTFH*`) — this is not an analysis set;
* **no EBV decoy** (`chrEBV` / `NC_007605`);
* **no HLA allele contigs** (`HLA-*`);
* **no phiX**.

ALT and patch scaffolds **are** present, as part of the full RefSeq assembly. That is the
relevant point for Ludmil's decoy-absorption question (rc2 §4): depletion is against a
reference that includes ALT/patch sequence, and the dedicated *pks* simulations found no
high-score island reads absorbed by the decoy/HLA/viral references that were tested.

Full dictionary: `human-GRC-db.sequences.tsv` beside this run's validation outputs
(regenerate with the command above).

---

## `human-GCA-phix-db.mmi` — `--t2t_phix_db`

| | |
|---|---|
| Size | 8,933,629,875 bytes |
| Sequences | **25** |
| Total bases | **3,117,292,070** |
| Naming | GenBank accessions (`CP0682…`, `CP086569.2`) |
| minimap2 `k` / `w` / `b` | 21 / 11 / 14 |
| `flag` | 0 — no `-H`, sequences stored, names kept |
| Dictionary sha256 | `ba0747e8d82f8be009b6494ee952e0497fb4c263e373b2d151449e4271b6a6f3` (first 399 bytes) |

**Assembly: T2T-CHM13v2.0 (GCA_009914755.4), complete and alone.** 25 sequences,
3,117,292,070 bp — the published total for CHM13v2.0 exactly:

| Sequence | What it is |
|---|---|
| `CP068277.2` … `CP068255.2` | chr1–chr22 and chrX, CHM13 |
| `CP086569.2` | chrY, from HG002 (as in CHM13v2.0) |
| `CP068254.1` | mitochondrion, 16,569 bp |

### ⚠ There is no phiX in the phiX index

**The file is named `human-GCA-phix-db.mmi`, the parameter is `--t2t_phix_db`, and the QC
column is `reads_after_t2t_phix` — but the index contains no phiX174 sequence.** phiX
(`NC_001422.1`, 5,386 bp) is not among the 25 sequences, and the total base count is
CHM13v2.0's to the base, which it could not be if a 26th sequence were present.

This is a naming claim the data does not support, and it is recorded here rather than
quietly corrected, because it changes what the pass means:

* **what it does**: a second human-depletion pass against T2T-CHM13v2.0, catching reads
  from sequence that GRCh38 does not represent. That is the pass's main value and it is
  unaffected;
* **what it does not do**: remove phiX spike-in. A library with a phiX spike carries
  those reads through to profiling;
* **effect on the *pks* result**: phiX is a 5.4 kb coliphage genome with no homology to
  the 50,768 bp *pks* island; its reads do not reach `clb` gene counts through the
  Bowtie2 MAPQ ≥ 40 + featureCounts path. So this is a provenance and naming defect, not
  a correctness one for the counts — but it should not be described in a methods section
  as phiX depletion, and `reads_after_t2t_phix` should be read as "reads surviving the
  T2T pass".

Resolution for v1.0.0: **documented, not changed.** Rebuilding the index would change
host depletion for every cohort already run and would split provenance mid-flight, which
is exactly what Ludmil's §2 asks us not to do. The option for a later release is to add
phiX to the index (or a third, tiny pass) and rename the parameter; either is a
behaviour change and belongs in a version where it can be validated.

---

## The build command

**It could not be recovered, and nothing in this repository or in the shared database
tree records it.** Both files were built on 2025-01-19, before pksProfiler existed, as
part of the CMPipeline reference set; there is no build script, log or README beside them
(`dbs/` holds the downloads and extracted trees for KrakenUniq, GTDB-Tk, CheckM2 and
geNomad, and no provenance for these two), and no file in either repository names them
except as a parameter value.

Rather than invent one, here is what **is** verifiable, read from the indexes themselves:

| Verifiable | `human-GRC-db.mmi` | `human-GCA-phix-db.mmi` |
|---|---|---|
| minimizer length `k` | 21 | 21 |
| minimizer window `w` | 11 | 11 |
| bucket bits `b` | 14 | 14 |
| homopolymer compression (`-H`) | not used | not used |
| target sequences stored | yes | yes |
| sequence names kept | yes | yes |

`k=21, w=11` is minimap2's **`-x sr`** preset, and `b=14` is the default. So both indexes
were built as

```
minimap2 -x sr -d <index>.mmi <reference>.fa
```

possibly with `-I` and `-t` as well; neither affects the index layout and neither is
stored in the file, so neither can be recovered. The input FASTAs' paths are likewise not
stored. This matters because it is also what the pipeline *uses*: `lib/HostDepletion.groovy`
runs both passes with `-ax sr`, so the index settings and the alignment preset agree.

To re-create equivalent indexes from scratch:

```
# --hg38_db  : GRCh38.p14 RefSeq full assembly, GCF_000001405.40
minimap2 -x sr -d human-GRC-db.mmi GCF_000001405.40_GRCh38.p14_genomic.fna

# --t2t_phix_db : T2T-CHM13v2.0, GCA_009914755.4  (note: no phiX, as above)
minimap2 -x sr -d human-GCA-phix-db.mmi GCA_009914755.4_T2T-CHM13v2.0_genomic.fna
```

These are reconstructions consistent with every verifiable property of the existing
files. They are **not** claimed to be the original commands, and an index built this way
will not be byte-identical (minimap2's version and `-I` batching both change the file).

## Re-checking either index

Identity, cheap — reads a few kB, and is what the table above reports:

```
python scripts/describe_minimap2_index.py <index.mmi>
```

Integrity, expensive — reads the whole 9 GB file, so run it deliberately, not per run:

```
sha256sum /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi
sha256sum /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi
```

`conf/databases.tsv` carries the per-run fingerprint (`sha256:listing`) for both, and
every run compares against it and reports `frozen` / `CHANGED` in RUN_REPORT.txt.

---

Back to the [README](../README.md).
