nextflow.enable.dsl=2

params.samplesheet = "sample_sheet.csv"
params.genomeDir   = "data/humanSTARindex/"
params.gtf         = "data/humanSTARindex/Homo_sapiens.GRCh38.99.gtf"
params.rscript     = "scripts/pat_down.R"

// Create channel from CSV
workflow {

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .set { samples_ch }

    // Step 1: Align reads
    aligned_bams_ch = align_reads(samples_ch)

    // Step 2: Merge BAM files
    merged_bam_ch = merge_bams(aligned_bams_ch)

    // Step 3: Run EMA
    ema_out_ch = run_ema(merged_bam_ch)

    // Step 4: R analysis
    downstream_analysis(ema_out_ch)
}

process align_reads {
    tag { sample.sample }
    publishDir "results/star", mode: 'copy'

    input:
    val sample

    output:
    path "${sample.sample}.bam"

    script:
    def read1_path   = file(sample.read1).toAbsolutePath()
    def read2_path   = file(sample.read2).toAbsolutePath()
    def genome_dir   = file(params.genomeDir).toAbsolutePath()
    def umi_start    = sample.cb_len.toInteger() + 1

    """
    STAR --runThreadN 64 \\
         --genomeDir ${genome_dir} \\
         --readFilesIn ${read2_path} ${read1_path} \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFileNamePrefix ${sample.sample}_ \\
         --outSAMattributes CR UR CY UY CB UB \\
         --soloType CB_UMI_Simple \\
         --soloCBstart 1 \\
         --soloUMIstart ${umi_start} \\
         --soloCBlen ${sample.cb_len} \\
         --soloUMIlen ${sample.umi_len} \\
         --soloUMIdedup 1MM_CR \\
         --soloCBwhitelist None \\
         --soloCBmatchWLtype 1MM_multi \\
         --soloCellReadStats Standard \\
         --readFilesCommand 'gunzip -c' \\
         --soloBarcodeReadLength 0

    mv ${sample.sample}_Aligned.sortedByCoord.out.bam ${sample.sample}.bam
    """
}

process merge_bams {

    publishDir "results/ema_merge", mode: 'copy'

    input:
    path bams

    output:
    path "Aligned.sortedByCoord.merged.out.bam"

    script:
    """
    source .emaenv/bin/activate
    echo "bamFiles:" > bamfiles.yaml
    for bam in ${bams.join(' ')}; do
      echo "  - \${bam}" >> bamfiles.yaml
    done
    ema_merge --bamFiles bamfiles.yaml --threads 100
    """
}

process run_ema {

    publishDir "results/emaout", mode: 'copy'

    input:
    path bam_file

    output:
    path "emaout"

    script:
    """
    source .emaenv/bin/activate
    ema --bamDir ${bam_file} \\
        --sequenceLen 101 \\
        --CellBarcodeLen 12 \\
        --BarcodeTag CB \\
        --gap 200 \\
        --min_read 2000 \\
        --min_cells 200 \\
        --min_genes 200 \\
        --gtfDir ${params.gtf}
    """
}

process downstream_analysis {

    publishDir "results/final", mode: 'copy'

    input:
    path ema_folder

    script:
    """
    Rscript ${params.rscript} ${ema_folder} Macs_sc
    """
}
