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

dir='06_FigureS6'
dir.create(dir)
setwd(dir)

###########
## pS6  ##
###########
load('/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/05_RMC_pre_post/03_result/02_whole_umap/06_subset/01_epi/00_obj_v3_3561_malignancy.Rdata')
obj_v3
obj1 = subset(obj_v3, subset = Malignancy == 'Malignant') 
obj = subset(obj1, subset = BiopsySite == 'Liver')
# table(obj$Treatment, obj$Patient)

pa_color = colorRampPalette(brewer.pal(9, "Paired"))(9)
p2 = DimPlot(obj, group.by = 'Patient', reduction = 'umap', cols= pa_color)
ggsave("01_PS6a_patient.pdf",plot=p2, width = 5, height = 4)

my_color = rev(c('#FFD966','#a8adb4'))
p2 = DimPlot(obj, reduction = 'umap', group.by  = 'Treatment',cols= my_color)
ggsave("02_PS6b_group.pdf",plot=p2, width = 5, height = 4)
