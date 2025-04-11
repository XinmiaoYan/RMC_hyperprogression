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
library(sctransform)
library(glmGamPoi)
library(openxlsx)
library(future)
library('msigdbr')
library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)
library(msigdf)
library(fgsea)
library(GSEABase)
library(biomaRt)
library(curl)
library(readxl)
library(msigdbr)
library(ggrepel)
rm(list = ls())

setwd('/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/')
dir='16_source_data'
dir.create(dir)
setwd(dir)

load('/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/16_source_data/01_FigureS2/00_scRNA_object.Rdata')

dir='16_FigureS16'
dir.create(dir)
setwd(dir)

#################
## figure S16  ##
#################
load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/2_subprojects/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/07_whole_clean_data_240422/03_similarity/01_ALL/00_obj_TumorMye_6339.Rdata')

table(obj$CellClass, obj$Treatment)
treat_color = c('PostNI'='#FFD966',
                'Baseline' = '#a8adb4',
                'TumorCell' = '#CAB2D6',
                'MyeloidCell' = '#FB9A99')

# treat_color = rev(c('#FFD966','#a8adb4'))
pt = 0.5
p1 = DimPlot(obj, group.by = 'Treatment', reduction = 'umap', cols = treat_color, pt.size = pt)
ggsave("01_pS16_Treatment.pdf",plot=p1, width = 5, height = 4)

p1 = DimPlot(obj, group.by = 'CellClass', reduction = 'umap', cols = treat_color, pt.size = pt)
ggsave("01_pS16_cell.pdf",plot=p1, width = 5, height = 4)

# p1 = DimPlot(obj, group.by = 'CellClass',split.by = 'Treatment',label = F, reduction = 'umap', cols = treat_color)
# ggsave("02_Split.pdf",plot=p1, width = 8, height = 4)


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
write.table(output_df, file = "01_umap_with_metadata.tsv", sep = "\t", quote = FALSE, row.names = TRUE)