# Running pksProfiler on HPC systems

pksProfiler uses Nextflow executors to submit each process to the local batch
scheduler. All compute nodes must be able to access the repository, input data,
reference indexes, Conda environments, work directory, and output directory.

## General recommendations

- Launch Nextflow from a small driver allocation rather than performing compute
  work on a login node.
- Put the Nextflow work directory on a shared high-throughput filesystem unless
  the site profile explicitly manages node-local scratch.
- Keep the same parameters and work directory when using `-resume`.
- Store site accounts, queues, projects, and QOS values in a private local config.
- Start with one representative sample before scaling to a cohort.

## Resource expectations

The bundled resource defaults are intentionally conservative and target
large-memory HPC systems processing cancer BAMs or large metagenomes. They
favor avoiding out-of-memory failures over minimizing queue time.

| Process or label | CPUs | Initial memory | Initial time |
|---|---:|---:|---:|
| BAM extraction | 4 | 64 GB | 8 h |
| FASTQ filtering | 4 | 128 GB | 10 h |
| Host depletion | 16 | 64 GB | 50 h |
| Bowtie2/*clb* alignment | 4 | 128 GB | 30 h |
| HMM profiling | 8 | 128 GB | 30 h |
| Low-resource summaries/plots | 4 | 100 GB | 10 h |
| Medium-resource tasks | 4 | 128 GB | 30 h |
| KrakenUniq/Bracken high-disk tasks | 4 | 256 GB | 50 h |

Some memory and time requests increase on retry, subject to the global limits
in `nextflow.config`. These defaults may wait a long time or exceed the limits
of smaller clusters, but changing them is not required for biological
correctness.

To use smaller requests without modifying the repository, create a local
configuration file. For example:

```groovy
// resources.config
process {
    withLabel:filter_reads {
        memory = 16.GB
    }
    withLabel:pks_align {
        memory = 32.GB
    }
    withLabel:pks_hmm {
        memory = 32.GB
    }
}
```

Add it to the run command after the selected profile:

```bash
nextflow run main.nf \
    -profile slurm \
    -c site.config \
    -c resources.config \
    [pipeline options]
```

Site administrators or experienced users should choose overrides appropriate
for their scheduler and data size.

## Generic Slurm

Create `site.config`:

```groovy
process {
    queue = 'compute'
    clusterOptions = '--account=my_account'
}
```

Run with:

```bash
nextflow run main.nf -profile slurm -c site.config [pipeline options]
```

For a long run, submit the Nextflow driver itself:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=pksProfiler-driver
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00

set -euo pipefail

module load nextflow
nextflow run main.nf -profile slurm -c site.config [pipeline options]
```

## UC San Diego TSCC

The included `tscc` profile selects Slurm, the `platinum` partition, and the
configured project account:

```bash
nextflow run main.nf -profile tscc [pipeline options]
```

Users outside that allocation should copy `conf/tscc.config` to a private site
config and change the queue, account, and QOS settings.

## NIH Biowulf

The included `biowulf` profile uses the `norm` partition, conservative scheduler
polling, and 200 GB of `lscratch` per process:

```bash
module load nextflow

nextflow run main.nf \
    -profile biowulf \
    --biowulf_lscratch_gb 200 \
    [pipeline options]
```

Increase `--biowulf_lscratch_gb` for very large FASTQ or BAM inputs. Submit the
Nextflow driver as a Biowulf batch job and follow current NIH guidance for
scheduler polling and local scratch.

This repository currently uses Conda environments. A fully containerized
Biowulf profile requires tested per-process Apptainer/Singularity images and is
future work; enabling Singularity alone does not replace the Conda directives.

Biowulf guidance: <https://hpc.nih.gov/apps/nextflow.html>

## PBS Pro

```bash
nextflow run main.nf -profile pbspro -c site.config [pipeline options]
```

Example site settings:

```groovy
process {
    queue = 'workq'
    clusterOptions = '-A my_project'
}
```

## IBM Spectrum LSF

```bash
nextflow run main.nf -profile lsf -c site.config [pipeline options]
```

Example site settings:

```groovy
process {
    queue = 'normal'
    clusterOptions = '-P my_project'
}
```

## Sun/Oracle Grid Engine

```bash
nextflow run main.nf -profile sge -c site.config [pipeline options]
```

Example site settings:

```groovy
process {
    queue = 'all.q'
    clusterOptions = '-P my_project'
}
```

## Monitoring resources

Always enable Nextflow reports for a benchmark:

```bash
RUN_TAG=$(date +%Y%m%d_%H%M%S)

nextflow run main.nf \
    [pipeline options] \
    -with-trace "run.${RUN_TAG}.trace.txt" \
    -with-report "run.${RUN_TAG}.report.html" \
    -with-timeline "run.${RUN_TAG}.timeline.html"
```

For Slurm, use the trace `native_id` values with `sacct` or `sstat` to compare
requested resources with elapsed time, CPU utilization, peak RSS, and disk I/O.

Nextflow executor documentation:
<https://www.nextflow.io/docs/latest/executor.html>
