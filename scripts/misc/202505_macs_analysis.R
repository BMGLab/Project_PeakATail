#!/usr/bin/Rscript --vanilla
library(Seurat)
library(patchwork)


setwd("/home/biolab/Projects/MacAPAs_wd/")
macs<- readRDS("/home/biolab/Projects/MacAPAs_wd/data/merged_raw_macs.rds")
macs[["percent.mt"]] <- PercentageFeatureSet(macs, pattern = "^MT-")
macs <- NormalizeData(object = macs)
macs <- FindVariableFeatures(macs)
macs <- ScaleData(macs)
macs <- RunPCA(macs)



Idents(macs)<- "orig.ident"
p <- VizDimLoadings(macs, dims = 1:2, reduction = "pca")
ggsave("2025_macs_pca_loadings.pdf", plot = p, width = 7, height = 5)

p<- DimPlot(macs, reduction = "pca") + NoLegend()
ggsave("2025_macs_pca.pdf", plot = p, width = 7, height = 5)

p <- ElbowPlot(macs)
ggsave("2025_macs_elbow_plot.pdf", plot = p, width = 7, height = 5)

p <- DimHeatmap(macs, dims = 1:10, cells = 500, balanced = TRUE)
ggsave("2025_macs_dim_heatmap.pdf", plot = p, width = 20, height = 15)

macs <- FindNeighbors(macs, dims = 1:10)
macs <- FindClusters(macs, resolution = 0.5)

macs <- RunUMAP(macs, dims = 1:10)

p <- DimPlot(macs, reduction = "umap", group.by = "orig.ident", label = T)
ggsave("2025_macs_umap_orig_ident.pdf", plot = p, width = 7, height = 5)

p <- DimPlot(macs, reduction = "umap", label = T)
ggsave("2025_macs_umap.pdf", plot = p, width = 7, height = 5)


maturation <- c("CD68", "MRC1", "CD209", "CD163", "FCGR1A", 
                "CD14", "CCR7", "CD40","TREM2", "STAB1", 
                "VCAN", "MARCO", "ITGAX", "SMDT1", "FOLR2",
                "IRF4", "IRF5", "ITGB2", "CLEC7A", "CD36",
                "MERTK", "NR3C2","NR1H2", "IRF8", "SIGLEC1")


activation <- c( "CD80", "CD86", "IL12B", "IL10", "TNF", 
                "HLA-DRB1", "HLA-DPB1", "HLA-DQB1", "CD274", "PDCD1LG2", "IL-6", 
                "TGFB1", "CXCL10", "CXCL9", "IL2RA", "IL1B", "FGL2",
                "RETN", "TLR2", "TLR4", "VEGFA", "CX3CR1", "CD209" )


pdf("all_activation_ridgeplots.pdf", width = 8, height = 10)  

for (gene in activation) {
  print(RidgePlot(macs, features = gene) + ggtitle(gene))
}

dev.off()

pdf("all_maturation_ridgeplots.pdf", width = 8, height = 10)  

for (gene in maturation) {
  print(RidgePlot(macs, features = gene) + ggtitle(gene))
}

dev.off()


pdf("all_maturation_vlnplots.pdf", width = 8, height = 10)  

for (gene in maturation) {
  print(VlnPlot(macs, features = gene,group.by = "orig.ident") + ggtitle(gene))
}

dev.off()

pdf("all_activation_vlnplots.pdf", width = 8, height = 10)  

for (gene in activation) {
  print(VlnPlot(macs, features = gene, group.by = "orig.ident") + ggtitle(gene))
}

dev.off()



