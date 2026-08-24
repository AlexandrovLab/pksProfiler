nextflow.enable.dsl=2

process FASTQC {
    scratch true
    label 'fastqc'
    publishDir("${params.fastqc_dir}", mode: 'copy')
    conda "${params.samtools_env}" // Note: This likely should be a fastqc env, not samtools
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(r1_fastq), path(r2_fastq)

    output:
    path("*.html"), emit: html
    path("*.zip"), emit: zip

    script:

    """
    # Defined report naming convention (FastQC: <basename>_fastqc.html/.zip)
    # Process basename using shell
    r1_base=\$(basename ${r1_fastq.toString()} .fastq.gz)
    r2_base=\$(basename ${r2_fastq.toString()} .fastq.gz)
    out1="\${r1_base}_fastqc.html"
    out2="\${r1_base}_fastqc.zip"
    out3="\${r2_base}_fastqc.html"
    out4="\${r2_base}_fastqc.zip"

    # Actual FastQC execution
    # Running FastQC on the extracted R1 and R2 fastq.gz files
    fastqc -t "${task.cpus}" -o ./  ${r1_fastq.toString()} ${r2_fastq.toString()}
    """
}
