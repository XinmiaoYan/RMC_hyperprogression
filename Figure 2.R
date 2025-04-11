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

#################
## figure 2d  ##
#################
parent_directory <- "/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/2_subprojects/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/07_whole_clean_data_240422/02_Degs_pathways"

# List subdirectories
subfolders <- list.dirs(path = parent_directory, full.names = TRUE, recursive = FALSE)
work_subfolders = subfolders[c(2,4)]

for(i in work_subfolders){
  setwd(i)
  ###
  dir_name = '04_heatmap'
  dir.create(dir_name)
  setwd(dir_name)
  focus_pathways = NULL
  CELL_CYCLE = c('REACTOME_CELL_CYCLE_MITOTIC',
                 'REACTOME_CELL_CYCLE',
                 'REACTOME_MAPK_FAMILY_SIGNALING_CASCADES',
                 'REACTOME_M_PHASE'
                 )
  immune_related_pathways = c('REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL',
                              'REACTOME_INTERFERON_ALPHA_BETA_SIGNALING',
                              'REACTOME_INTERFERON_SIGNALING',
                              'REACTOME_INTERFERON_GAMMA_SIGNALING'
                              )
  focus_pathways$pathway = c(CELL_CYCLE,immune_related_pathways)
  focus_pathways$Category = c(rep('CelCycle',4),rep('ImmuneSystem',4))
  focus_pathways = as.data.frame(focus_pathways)

  tumor_dt = read_excel(paste(i,'/Epi/02_pathways/03_Reactome/01_Reactome.xlsx', sep=''))
  tumor_focus = (tumor_dt%>% filter(ID%in%focus_pathways$pathway))[,c('ID','NES','pvalue','p.adjust')]
  tumor_focus$Celltype = rep('TumorCell',dim(tumor_focus)[1])

  CD8T = read_excel(paste(i,'/CD8T/02_pathways/03_Reactome/01_Reactome.xlsx', sep=''))
  setdiff(focus_pathways$pathway, CD8T$ID)
  CD8T_focus = (CD8T%>% filter(ID%in%focus_pathways$pathway))[,c('ID','NES','pvalue','p.adjust')]
  CD8T_focus$Celltype = rep('CD8T',dim(CD8T_focus)[1])

  # CD4T
  CD4T = read_excel(paste(i,'/CD4T/02_pathways/03_Reactome/01_Reactome.xlsx', sep=''))
  setdiff(focus_pathways$pathway, CD4T$ID)

  CD4T_focus = (CD4T%>% filter(ID%in%focus_pathways$pathway))[,c('ID','NES','pvalue','p.adjust')]
  CD4T_focus$Celltype = rep('CD4T',dim(CD4T_focus)[1])

  Myeloid = read_excel(paste(i,'/Myeloids/02_pathways/03_Reactome/01_Reactome.xlsx', sep=''))
  setdiff(focus_pathways$pathway, Myeloid$ID)
  Myeloid_focus = (Myeloid%>% filter(ID%in%focus_pathways$pathway))[,c('ID','NES','pvalue','p.adjust')]
  Myeloid_focus$Celltype = rep('Myeloids',dim(Myeloid_focus)[1])

  df = rbind(tumor_focus,CD8T_focus,CD4T_focus,Myeloid_focus)
  rownames(df) = NULL
  df$stars <- ifelse(df$pvalue < 0.001, "***",
                      ifelse(df$pvalue < 0.01, "**",
                          ifelse(df$pvalue < 0.05, "*", "")))

  df = mutate(df, 
              Name = case_when(ID == 'REACTOME_CELL_CYCLE_MITOTIC' ~ 'Cell cycle mitotic',
                               ID == 'REACTOME_CELL_CYCLE' ~ 'Cell cycle',
                               ID == 'REACTOME_MAPK_FAMILY_SIGNALING_CASCADES' ~ 'MAPK signaling cascades',
                               ID == 'REACTOME_M_PHASE' ~ 'M phase',
                               ID == 'REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL' ~ 'Immunoregulatory between lymphoid and non-lymphoid cell',
                               ID == 'REACTOME_INTERFERON_ALPHA_BETA_SIGNALING' ~ 'Interferon alpha & beta signaling',
                               ID == 'REACTOME_INTERFERON_SIGNALING' ~ 'Interferon signaling',
                               ID == 'REACTOME_INTERFERON_GAMMA_SIGNALING' ~ 'Interferon gamma signaling',
                               TRUE ~ NA_character_))

  pathways = c(
  'Cell cycle mitotic',
  'Cell cycle',
  'MAPK signaling cascades',
  'M phase',
  'Immunoregulatory between lymphoid and non-lymphoid cell',
  'Interferon alpha & beta signaling',
  'Interferon signaling',
  'Interferon gamma signaling'
  )
df$Celltype <- factor(df$Celltype, levels = c('TumorCell','CD8T','CD4T','Myeloids'))
df$Name <- factor(df$Name, levels = rev(pathways))

heatmap_plot <- ggplot(df, aes(x = Celltype, y = Name)) +
  geom_tile(aes(fill = NES)) +
  scale_fill_gradient2(low = "#55ACEE",mid = "white",high = "#d12d50",midpoint = 0) +
  geom_text(aes(label = stars), color = "black") +  # 使用stars列添加星号
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot=heatmap_plot,file = 'p2d.pdf', width = 6, height = 3)

}

#################
## figure 2e  ##
#################
load('/rsrch3/scratch/genomic_med/xyan4/Project/1_Jianfeng/2_subprojects/6_RMC_RCC/05_RMC_pre_post/03_result/02_whole_umap/06_subset/01_epi/00_obj_v3_3561_malignancy.Rdata')
obj_v3
obj = subset(obj_v3, subset = Malignancy == 'Malignant') #3196

cutoff = 0
obj_S100A9_MKI67 = subset(x = obj, subset = S100A9 > cutoff & MKI67 > cutoff)
obj_S100A9_MKI67$class2 = rep('S100A9+MKI67+',dim(obj_S100A9_MKI67)[2])
obj_other = subset(x = obj, subset = S100A9 > cutoff & MKI67 > cutoff, invert = T)
obj_other$class2 = rep('Other',dim(obj_other)[2])

obj_0 = merge(obj_S100A9_MKI67, obj_other)
table(obj_0$Treatment, obj_0$class2)

cell_number = as.data.frame(cbind(table(obj_0$Treatment, obj_0$class2)))
cell_number$id = rownames(cell_number)
library(reshape)
cell_number = melt(cell_number, id = 'id')
colnames(cell_number)<-c("Sample","Cluster","number")

my_color = rev(c('#FFD966','#a8adb4'))
hist_plot <- ggplot(cell_number,aes(Group,number,fill=Cluster))+
  scale_fill_manual(values = alpha(my_color, 0.9)) +
  #scale_fill_brewer(values = my_color) +
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab('Proportion')+ theme_bw()+
  theme(axis.ticks.length=unit(0.3,'cm'),axis.text.x = element_text(angle = 90)) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + xlab('')+
  guides(fill=guide_legend(title=NULL)) 
ggsave(plot = hist_plot, file = 'p2e.pdf', width = 3, height = 3)
