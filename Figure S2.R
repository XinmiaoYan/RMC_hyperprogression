library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(patchwork)
library(ggplot2)
library(devtools)
library(CelliD)
library(tidyverse)
library(ggpubr)
library(clustree)
library(tidyverse)
library(readxl)
library(colorspace)
library(RColorBrewer)
library(viridis)
library(scales)
library(ggmap)
library(DoubletFinder)
library(sctransform)
library(glmGamPoi)
rm(list = ls())

setwd('/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/')
dir='16_source_data'
dir.create(dir)
setwd(dir)

load('/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/16_source_data/01_FigureS2/00_scRNA_object.Rdata')

dir='01_FigureS2'
dir.create(dir)
setwd(dir)

obj = scRNA_object
pt = 0.7
# my_color = colorRampPalette(brewer.pal(12, "Paired"))(12)
my_color = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
             "#E31A1C", "#FDBF6F", "#CAB2D6","#FF7F00", "#6A3D9A")

#################
## figure S2a  ##
#################
pS2a = DimPlot(obj, label = T, reduction = 'umap',pt.size = pt)
ggsave("01_umap_pS2a.pdf",plot=pS2a, width = 5.5, height = 4.5)

# Extract UMAP coordinates
umap_coords <- Embeddings(obj, reduction = "umap")

# Extract metadata
metadata <- obj@meta.data
table(metadata$seurat_clusters, metadata$seurat_clusters_new)

identical(rownames(umap_coords), rownames(metadata))
# Combine UMAP and metadata
output_df <- cbind(UMAP_1 = umap_coords[, 1],
                   UMAP_2 = umap_coords[, 2],
                   metadata)

head(output_df)
# Write to a tab-separated file
write.table(output_df, file = "01_umap_with_metadata_23880.tsv", sep = "\t", quote = FALSE, row.names = TRUE)


anno_gene = c("PTPRC","CD19", "MS4A1","CD79A","CD79B","BANK1","JCHAIN","XBP1","MZB1",    
              "CD3D","CD3E","IL7R", "CD8A","CD8B", "GZMA","GZMK","GNLY","NKG7", #"CD4",
              "MKI67","BIRC5","CDK1","TOP2A","TYMS",   
              "FCGR3A","CD14","ITGAX","LYZ","FCN1","S100A8","S100A9","C1QB","C1QC","CD68",
              "TPSAB1","TPSB2","LILRA4","CLEC4C",
              "EPCAM","CDH1", "KRT8","COL1A1","COL1A2","COL3A1","COL6A1",   
              "PECAM1","VWF", "ENG",'CLDN5')

#################
## figure S2c  ##
#################

Idents(obj) = 'seurat_clusters'
pS2c = DotPlot(obj, features = unique(anno_gene), cluster.idents = F) + 
  RotatedAxis() + # scale_color_gradient2(low = "lightgray", mid = "white", high = "#FB8604") + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=7),legend.text = element_text(size=7)) +
  FontSize(x.text = 9, y.text = 10,main = 9) +labs(x = NULL, y = NULL)
ggsave("02_bubblePlot_pS2c.pdf",plot=pS2c, width = 11, height = 5)

#################
## figure S2b  ##
#################

Idents(obj) <- "level1_cell"
levels(obj) <- c("B","Plasma", "T&NK", "T_Proliferating",'Myeloid','Mast', "pDCs","Epithelial", "Fibroblasts", "Endothelial")
pS2b = DimPlot(obj, reduction = 'umap',cols = my_color, pt.size = pt) #same as pS2b
ggsave("03_clustering_pS2b.pdf",plot=p2b, width = 6, height = 4.5)

#################
## figure S2d  ##
#################

pS2d = DotPlot(obj, features = unique(anno_gene), cluster.idents = F) + 
  RotatedAxis() +  scale_color_gradient2(low = "lightgray", mid = "white", high = "#FB8604") + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=7),legend.text = element_text(size=7)) +
  FontSize(x.text = 9, y.text = 10,main = 9) +labs(x = NULL, y = NULL)
ggsave("03_pS2d.pdf",plot=pS2d, width = 11, height = 3)

# output data
# Extract the data used for the plot
dot_data <- pS2d$data
# Preview the data
head(dot_data)
# Write to a tab-separated file
write.table(dot_data, file = "03_pS2d_data_23880.tsv", sep = "\t", quote = FALSE, row.names = TRUE)