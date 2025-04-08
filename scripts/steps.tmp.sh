
STAR --runThreadN 64 \
     --genomeDir data/humanSTARindex/ \
     --readFilesIn data/Mono_day0_S1_R1_001.fastq.gz data/Mono_day0_S1_R2_001.fastq.gz \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes CR UR CY UY CB UB \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 \
     --soloUMIstart 13 \
     --soloCBlen 12 \
     --soloUMIlen 8 \
     --soloUMIdedup 1MM_CR \
     --soloCBwhitelist None \
     --soloCBmatchWLtype 1MM_multi \
     --soloCellReadStats Standard \
     --readFilesCommand 'gunzip -c' \
     --soloBarcodeReadLength 0

## Install PeakATail
   python -m venv .emaenv
   source .emaenv/bin/activate
   pip install -e .
   ema --help

#Merge Bam files
   ema_merge --bamFiles bamfiles.yaml --threads 100

#Run PeakATail
ema --bamDir Aligned.sortedByCoord.merged.out.bam \
    --sequenceLen 101 \
    --CellBarcodeLen 12 \
    --BarcodeTag CB \
    --gap 200 \
    --min_read 2000 \
    --min_cells 200 \
    --min_genes 200 \
    --gtfDir /home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf

#Run downstream analysis:
Rscript /home/sharedFolder/yk_data/pat_down.R /mnt/second/Macrophage_scRNAseq/Run_peakATail/emaout Macs_sc
