nextflow.enable.dsl=2

params.samplesheet = "sample_sheet.csv"
params.genomeDir   = "data/humanSTARindex/"
params.gtf         = "data/humanSTARindex/Homo_sapiens.GRCh38.99.gtf"
params.chromSizes  = "data/chrom.sizes.nochr.filt"
params.flankLength = 5000
params.rscript     = "scripts/pat_down.R"

workflow {

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .set { samples_ch }

    // Step 1: Align reads
    aligned_bams_ch = align_reads(samples_ch)

    // Step 2: Wait for all BAMs then merge
    merged_bam_ch = aligned_bams_ch.collect().ifEmpty([]) | merge_bams

    // Step 3: Filter BAMs to UTR
    filtered_sampled_bam_ch = merged_bam_ch | filter_and_sample_bam

    // Step 4: Run EMA
    ema_out_ch = run_ema(filtered_sampled_bam_ch)

    // Step 5: R Analysis
    downstream_analysis(ema_out_ch)
}

process align_reads {
    tag { sample.sample }
    publishDir "results/star", mode: 'copy'

     // conda "./envs/star_env.yaml" // This YAML should define the STAR environment
    input:
    val sample

    output:
    path "${sample.sample}.bam"

    script:
    def read1_path = file(sample.read1).toAbsolutePath()
    def read2_path = file(sample.read2).toAbsolutePath()
    def genome_dir = file(params.genomeDir).toAbsolutePath()
    def umi_start  = sample.cb_len.toInteger() + 1

    """
    source /home/biolab/miniconda3/etc/profile.d/conda.sh
    conda activate STAR
    STAR --runThreadN 64 \
         --genomeDir ${genome_dir} \
         --readFilesIn ${read2_path} ${read1_path} \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${sample.sample}_ \
         --outSAMattributes CR UR CY UY CB UB \
         --soloType CB_UMI_Simple \
         --soloCBstart 1 \
         --soloUMIstart ${umi_start} \
         --soloCBlen ${sample.cb_len} \
         --soloUMIlen ${sample.umi_len} \
         --soloUMIdedup 1MM_CR \
         --soloCBwhitelist None \
         --soloCBmatchWLtype 1MM_multi \
         --soloCellReadStats Standard \
         --readFilesCommand 'gunzip -c' \
         --soloBarcodeReadLength 0

    mv ${sample.sample}_Aligned.sortedByCoord.out.bam ${sample.sample}.bam
    conda deactivate
    """
}

process merge_bams {
    publishDir "results/ema_merge", mode: 'copy'

    input:
    path bams

    output:
    path "Aligned.sortedByCoord.merged.out.bam"

    script:
    def env_path = file("tools/PeakATail/.emaenv/bin/activate").toAbsolutePath()
    """
    source ${env_path}
    for bam in ${bams.join(' ')}; do
      echo "  - \${bam}" >> bamfiles.yaml
    done
    ema_merge --bamFiles bamfiles.yaml --threads 100
    deactivate
    """
}

process filter_and_sample_bam {
    tag { bam.getBaseName() }
    publishDir 'results/utr_sampled', mode: 'copy'

    input:
    path bam

    output:
    path "${bam.getBaseName().replace('.bam','')}_3utr_selected.bam"

    script:
    def gtf_file = file(params.gtf).toAbsolutePath()
    def chromosome_size = file(params.chromSizes).toAbsolutePath()
    """
    source /home/biolab/miniconda3/etc/profile.d/conda.sh
    conda activate STAR
    
    python /home/biolab/Projects/PeakATail_wd/scripts/gtf2bed.py ${gtf_file}  > regions_3utr_py.bed    
    bedtools intersect \
      -wa -s \
      -abam ${bam} \
      -b regions_3utr_py.bed  \
      | samtools view -b - \
      > ${bam.getBaseName().replace('.bam','')}_3utr_selected.bam

    samtools index ${bam.getBaseName().replace('.bam','')}_3utr_selected.bam

    conda deactivate
    """
}


process run_ema {
    publishDir "results/emaout", mode: 'copy'

    input:
    path bam_file

    output:
    path "emaout"

    script:
    def env_path = file("tools/PeakATail/.emaenv/bin/activate").toAbsolutePath()
    def gtf_file = file(params.gtf).toAbsolutePath()
    """
    source ${env_path}
    ema --bamDir ${bam_file} \
        --sequenceLen 101 \
        --CellBarcodeLen 12 \
        --BarcodeTag CB \
        --gap 200 \
        --min_read 20 \
        --min_cells 20 \
        --min_genes 20 \
        --gtfDir ${gtf_file}
    deactivate
    """
}

process downstream_analysis {
    publishDir "results/final", mode: 'copy'

    input:
    path ema_folder

    script:
    def downstream_r = file(params.rscript).toAbsolutePath()
    """
    Rscript ${downstream_r} ${ema_folder} Macs_sc
    """
}
