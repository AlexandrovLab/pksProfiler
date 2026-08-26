# Lightweight regression checks

Run the repository checks with:

```bash
bash tests/run_checks.sh
```

The checks do not run bioinformatics workloads or require the large host and
taxonomy databases. They verify:

- Nextflow syntax for the main workflow and included modules;
- the 19-gene GFF and HMM reference contracts;
- presence of the complete bundled Bowtie2 index;
- mate-safe taxonomy joins;
- combined species-support table structure;
- positive and negative example-table structure; and
- rejection of duplicate sample identifiers before task execution.

These checks complement, rather than replace, end-to-end runs on representative
BAM, FASTQ, HMM, and taxonomy inputs.
