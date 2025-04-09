conda activate STAR

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
   source .emaenv/bin/activate
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





#Align reads to genome using STAR for each sample separately
for sample in A549trd_M2mac_day8 A549trd_M2mac_day9 A549trd_mac_day7 untrd_M2mac_day8 untrd_M2mac_day9 untrd_mac_day7;
   do 
   STAR --runThreadN 64 \
     --genomeDir data/humanSTARindex/ \
     --readFilesIn data/"$sample"_R2_001.fastq.gz data/"$sample"_R1_001.fastq.gz \
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
     --soloBarcodeReadLength 0;
done

#Create a bamfiles.yaml from the outputs of the previous step
source .emaenv/bin/activate

ema_merge --bamFiles bamfiles.yaml --threads 100 #This will output a file Aligned.sortedByCoord.merged.out.bam

#Run ema aka PeakATail 
source .emaenv/bin/activate

ema --bamDir data/Aligned.sortedByCoord.merged.out.bam \
    --sequenceLen 101 \
    --CellBarcodeLen 12 \
    --BarcodeTag CB \
    --gap 200 \
    --min_read 2000 \
    --min_cells 200 \
    --min_genes 200 \
    --gtfDir data/humanSTARindex/Homo_sapiens.GRCh38.99.gtf #This will create an output folder named 'emaout'

#Run downstream analysis using the output folder from the previous step. 'Macs_sc' is just to name the run.
Rscript scripts/pat_down.R ./emaout Macs_sc
