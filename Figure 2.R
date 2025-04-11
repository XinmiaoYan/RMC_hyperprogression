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