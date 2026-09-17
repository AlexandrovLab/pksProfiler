# Glossary

Terms used across pksProfiler's documentation and outputs.

These come up constantly in the outputs. If you already know them, skip ahead.

| Term | What it means here |
|---|---|
| **read** | One short DNA fragment from the sequencer, typically 100–150 letters long. A sample contains hundreds of millions of them. |
| **host depletion** | Throwing away the human reads first, so only non-human reads are searched. Most of a tumour sample is human. |
| **breadth** | What fraction of the island's 50,767 letters has at least one read sitting on it. 20% breadth means a fifth of the island was seen. |
| **assembly** | Stitching overlapping reads back into longer stretches of DNA, like reassembling a shredded page. |
| **contig** | One stretch produced by assembly — short for *contiguous sequence*. A good assembly gives long contigs; a poor one gives many short fragments. |
| **binning** | Sorting the contigs from a mixed sample into groups that probably came from the same organism. |
| **MAG** | *Metagenome-assembled genome* — one of those groups, i.e. a draft genome of a single bacterium recovered from a mixture, without ever culturing it. |
| **tier** | pksProfiler's evidence label for a sample, from `negative` up to `extensive_island`. Not a yes/no call — see [Evidence tiers](#evidence-tiers). |
| **MLST / ST** | *Multi-locus sequence typing*. Compares 7 housekeeping genes to a public database and returns a numbered **sequence type**, e.g. ST73. A standard way to name a bacterial strain. |
| **phylogroup** | A major branch of the *E. coli* family tree (A, B1, B2, D, E, F). Colibactin-producers are overwhelmingly **B2**. |
| **prophage** | A virus genome parked inside a bacterial chromosome. Relevant because colibactin damages DNA, which can wake prophages up. |

**"Assembled input"**, where it appears below, simply means *contigs or MAGs* — output from an
assembly step — as opposed to raw reads. Strain typing needs assembled input because you cannot
type a strain from a handful of reads; you need enough sequence to recover all 7 MLST genes.

---

Back to the [README](../README.md).
