process mapReads {

    scratch true
    label 'mapBothReads'
    publishDir("${params.mapped_reads_dir}", mode: 'copy')
    conda "${params.minimap2_env}"

    input:
    tuple val(sampleID), path(r1_fastq), path(r2_fastq)
    val(mmi_files)  // pass the list of .mmi files

    output:
    tuple val(sampleID),
        path("${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz"),
        path("${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz")

    script:
    """
	set -euo pipefail

    out1="${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz"
    out2="${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz"

    # ========= Actual execution =========
    echo "Running on R1"
    # Run minimap2 on hg38 reference
    minimap2 -2 -ax sr -t 16 ${params.hg38_db} ${r1_fastq} -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz

    # Run minimap2 on t2t_phix reference
    minimap2 -2 -ax sr -t 16 ${params.t2t_phix_db} ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz

    # Create temporary file for processing
    cp ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz

    # Convert the mmi_files list to a space-separated string
    mmi_files_list=\$(echo ${mmi_files.join(' ')})
    # Iterate over pangenome database files (list passed from main workflow)
    for mmi in \${mmi_files_list}; do
        echo "Running minimap2 on \${mmi}"
        minimap2 -2 -ax sr -t 16 \${mmi} ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.new.fastq.gz

        # Move output to temporary file for next iteration
        mv ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.new.fastq.gz ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz
    done

    # Final move to save the processed output
    mv ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz

    echo "Running on R2"
    # Run minimap2 on hg38 reference
    minimap2 -2 -ax sr -t 16 ${params.hg38_db} ${r2_fastq} -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz

    # Run minimap2 on t2t_phix reference
    minimap2 -2 -ax sr -t 16 ${params.t2t_phix_db} ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz

    # Create temporary file for processing
    cp ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz

    # Convert the mmi_files list to a space-separated string
    mmi_files_list=\$(echo ${mmi_files.join(' ')})
    # Iterate over pangenome database files (list passed from main workflow)
    for mmi in \${mmi_files_list}; do
        echo "Running minimap2 on \${mmi}"
        minimap2 -2 -ax sr -t 16 \${mmi} ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.new.fastq.gz

        # Move output to temporary file for next iteration
        mv ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.new.fastq.gz ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz
    done

    # Final move to save the processed output
    mv ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz
    """
}



