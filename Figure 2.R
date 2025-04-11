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

dir='01_Figure2'
dir.create(dir)
setwd(dir)

#################
## figure 2b  ##
#################
Idents(obj) <- "level1_cell"
levels(obj) <- c("B","Plasma", "T&NK", "T_Proliferating",'Myeloid','Mast', "pDCs","Epithelial", "Fibroblasts", "Endothelial")
p2b = DimPlot(obj, reduction = 'umap',cols = my_color, pt.size = pt) #same as pS2b
ggsave("01_clustering_p2b.pdf",plot=p2b, width = 6, height = 4.5)

setwd('../')

#################
## figure 2c  ##
#################
load('/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/06_subset/01_epi/00_obj_v3_3561_malignancy.Rdata')

obj = obj_v3
pt = 0.6
pa_color = colorRampPalette(brewer.pal(9, "Paired"))(9)
p2c1 = DimPlot(obj, group.by = 'Patient', reduction = 'umap',pt.size = pt, cols= pa_color)
ggsave("02_p2c1_patient.pdf",plot=p2c1, width = 7, height = 6)

treat_color = rev(c('#FFD966','#a8adb4'))
p2c2 = DimPlot(obj, group.by = 'Treatment', reduction = 'umap', cols = treat_color, pt.size = pt)
ggsave("03_p2c2_Treatment.pdf",plot=p2c2, width = 7, height = 6)

# Extract UMAP coordinates
umap_coords <- Embeddings(obj, reduction = "umap")
# Extract metadata
metadata <- obj@meta.data
# Combine UMAP and metadata
output_df <- cbind(UMAP_1 = umap_coords[, 1],
                   UMAP_2 = umap_coords[, 2],
                   metadata)
head(output_df)
# Write to a tab-separated file
write.table(output_df, file = "01_umap_with_metadata.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
