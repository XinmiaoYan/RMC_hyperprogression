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

dir='07_FigureS7'
dir.create(dir)
setwd(dir)

###########
## pS7  ##
###########
degs = read_excel('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/2_subprojects/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/07_whole_clean_data_240422/02_Degs_pathways/02_Liver_degs/Epi/01_degs/01_degs.xlsx')
df = degs
cutoff = 0.1
df = mutate(df,
            Condition = case_when((avg_log2FC >= cutoff)&(p_val_adj<=0.05) ~ 'Up', 
                                (avg_log2FC < cutoff)&(p_val_adj<=0.05) ~ 'Down', 
                                TRUE ~ 'Stable'))#NA_character_
df$gene = df$...1
df$avg_log2FC = as.numeric(df$avg_log2FC)
df$p_val_adj = as.numeric(df$p_val_adj)
df$logp = -log10(df$p_val_adj)

# df$label=ifelse(df$gene %in% genes,df$gene,"")
df$label=ifelse(df$p_val_adj < 0.05 & abs((df$avg_log2FC)) >= 0.6,df$gene,"")

colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
p <- ggplot(data = df, 
            aes(x = avg_log2FC, 
                y = logp, 
                colour=avg_log2FC,
                label = Condition,
                size = abs(avg_log2FC))) +
geom_point(alpha=0.8) +
scale_color_gradientn(colors = colors) +
geom_vline(xintercept=c(-0.2,0.2),lty=4,col="gray",lwd=0.6) +
geom_hline(yintercept = -log10(0.01),lty=4,col="gray",lwd=0.6) +
labs(x="avg_log2FC",
    y="-log10 (p_val_adj)",
    title=" ")  +
theme_bw()+
theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p1 = p+geom_text_repel(data = df, aes(x = df$avg_log2FC, 
                                y = logp, 
                                label = label),
                    size = 3,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"), 
                    vadjust= -1,
                    #segment.color = "black", 
                    show.legend = FALSE)

ggsave(plot = p1, file = '01_pS7.pdf', width = 14, height = 13)
