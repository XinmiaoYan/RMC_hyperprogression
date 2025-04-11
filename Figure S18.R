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

dir='18_FigureS18'
dir.create(dir)
setwd(dir)

#################
## figure S18  ##
#################

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/2_subprojects/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/07_whole_clean_data_240422/03_similarity/01_ALL/00_obj_TumorMye_6339.Rdata')
# obj = merge(obj_epi,obj_mye)
obj = obj
cell_base_tumor = colnames(subset(obj, subset = (Treatment == 'Baseline' & CellClass == 'TumorCell')))
cell_Post_tumor = colnames(subset(obj, subset = (Treatment == 'PostNI' & CellClass == 'TumorCell')))
cell_base_mye = colnames(subset(obj, subset = (Treatment == 'Baseline' & CellClass == 'MyeloidCell')))
cell_Post_mye = colnames(subset(obj, subset = (Treatment == 'PostNI' & CellClass == 'MyeloidCell')))
# length(cell_base_tumor)
# length(cell_Post_tumor)
# length(cell_base_mye)
# length(cell_Post_mye)
############################################

table(obj$Treatment, obj$CellClass)
varible_genes_3000 = VariableFeatures(obj)
matrix = obj@assays$RNA@scale.data[varible_genes_3000,]
# matrix = test[1:5,1:5]

# 计算欧氏距离
dir = '01_euclidean'
dir.create(dir)
setwd(dir)

Euclidean_data <- read.xlsx("/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/07_whole_clean_data_240422/03_similarity_240520/02_Liver/01_euclidean/01_Euclidean_distance.xlsx", sheet = 1)
res4 = as.data.frame(Euclidean_data[, -1])
dim(res4)
res4 = mutate(res4,
              Cell = case_when(Cell1 %in% cell_base_tumor ~ 'Baseline_Tumor',
                               Cell1 %in% cell_Post_tumor ~ 'PostNI_Tumor',
                               TRUE ~ NA_character_))
table(res4$Cell)
res4$Cell = factor(res4$Cell, level = c('Baseline_Tumor','PostNI_Tumor'))
levels(res4$Cell)


dt = res4 %>% filter(((Cell1 %in% cell_base_tumor) & (Cell2 %in% cell_base_mye))|((Cell1 %in% cell_Post_tumor)&(Cell2 %in% cell_Post_mye)))
dt$Distance = rep('Euclidean', dim(dt)[1])
p1 = ggplot(dt, aes(Distance, Euclidean_distance, fill=Cell, color = Cell)) + 
      geom_split_violin(alpha = .5,trim = TRUE,draw_quantiles = T) + 
      geom_boxplot(color = '#fcfcfc',width = 0.25, linewidth = 0.35,notch = FALSE, show.legend = FALSE,notchwidth = .2, outlier.shape = NA, coef=0) +
      labs(x=NULL,y="Distance to myeloid cells") +
      theme_classic() +
      stat_compare_means(method = "t.test", #wilcox.test
                        show.legend = F,
                        label = "p.format", 
                        hide.ns = TRUE)+
      theme(text = element_text(size = 10)) +
      scale_color_manual(values = rev(c('#FFD966','#a8adb4')))+
      scale_fill_manual(values = rev(c('#FFD966','#a8adb4'))) 

ggsave(file = 'pS181.pdf',plot = p1,  width = 3, height=3) 

# write.xlsx(dt, "01_Euclidean_distance.xlsx", rowNames = T)
dt=dt[,-5]
write.table(dt, file = "01_Euclidean_distance.tsv", sep = "\t", quote = FALSE, row.names = TRUE)


setwd('../')



dir = '02_manhattan'
dir.create(dir)
setwd(dir)
# use heatmap to show the average distance between mye and tumor cells
manhattan_data <- read.xlsx("/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/07_whole_clean_data_240422/03_similarity_240520/02_Liver/02_manhattan/01_Manhattan_distance.xlsx", sheet = 1)
res4 = manhattan_data[, -1]
dim(res4)
res4 = mutate(res4,
              Cell = case_when(Cell1 %in% cell_base_tumor ~ 'Baseline_Tumor',
                               Cell1 %in% cell_Post_tumor ~ 'PostNI_Tumor',
                               TRUE ~ NA_character_))
table(res4$Cell)
res4$Cell = factor(res4$Cell, level = c('Baseline_Tumor','PostNI_Tumor'))
levels(res4$Cell)


dt = res4 %>% filter(((Cell1 %in% cell_base_tumor) & (Cell2 %in% cell_base_mye))|((Cell1 %in% cell_Post_tumor)&(Cell2 %in% cell_Post_mye)))
dt$Distance = rep('Manhattan', dim(dt)[1])
p1 = ggplot(dt, aes(Distance, Manhattan_distance, fill=Cell, color = Cell)) + 
      geom_split_violin(alpha = .5,trim = TRUE,draw_quantiles = T) + 
      geom_boxplot(color = '#fcfcfc',width = 0.25, linewidth = 0.35,notch = FALSE, show.legend = FALSE,notchwidth = .2, outlier.shape = NA, coef=0) +
      labs(x=NULL,y="Distance to myeloid cells") +
      theme_classic() +
      stat_compare_means(method = "t.test", #wilcox.test
                        show.legend = F,
                        label = "p.format", 
                        hide.ns = TRUE)+
      theme(text = element_text(size = 10)) +
      scale_color_manual(values = rev(c('#FFD966','#a8adb4')))+
      scale_fill_manual(values = rev(c('#FFD966','#a8adb4'))) 

ggsave(file = 'p2k2.pdf',plot = p1,  width = 3, height=3) 

# write.xlsx(dt, "01_Manhattan_distance.xlsx", rowNames = T)
dt=dt[,-5]
write.table(dt, file = "01_Manhattan_distance.tsv", sep = "\t", quote = FALSE, row.names = TRUE)