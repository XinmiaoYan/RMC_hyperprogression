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
dir='08_FigureS8'
dir.create(dir)
setwd(dir)

#################
## figure S8  ##
#################
load('/rsrch3/scratch/genomic_med/xyan4/Project/1_Jianfeng/2_subprojects/6_RMC_RCC/05_RMC_pre_post/03_result/02_whole_umap/06_subset/01_epi/00_obj_v3_3561_malignancy.Rdata')
obj = subset(obj_v3, subset = Malignancy == 'Malignant' & BiopsySite %in% c('Liver'))

Idents(obj) <- "Treatment"
Idents(obj) = factor(Idents(obj), levels = rev(c("PostNI", "Baseline")))
test_sign <- list(rev(c("PostNI", "Baseline")))
###
vp_case <- function(gene_signature, name, test_sign,nrow, width, height){
  median.stat <- function(x){
    out <- quantile(x, probs = c(0.5))
    names(out) <- c("ymed")
    return(out) 
  }

  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(obj, features = signature,
            pt.size = 0, 
            cols = rev(c('#FFD966','#a8adb4')), 
            fill.by = 'Treatment',
            log = FALSE,
            y.max = y_max
    ) + NoLegend() + stat_compare_means(comparisons = test_sign, label = "p.format") + labs(x = NULL)+
    stat_summary(fun = mean, geom='point', size = 6, colour = "#4C4C4C", shape = 95)
    #geom_boxplot(color = "White")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]])
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  #cowplot::plot_grid(plotlist = plot_list)
  p = wrap_plots(plots = plot_list, nrow = nrow)
  file_name <- paste0(name, ".pdf")
  ggsave(file_name, plot = p, width = width, height = height)
}

focus_genes = c('S100A9', 'CEBPB','IFNGR1','JAK1','MAP3K1','MAP2K2','MAP3K11','EP300') 
# focus_genes = intersect(focus_genes, rownames(obj))
vp_case(gene_signature = focus_genes, name ='01_pS8.pdf', test_sign = test_sign, nrow = 1, width = 15, height = 3)

# Extract expression data and cell identities for GeneX
plot_data <- FetchData(obj, vars = c(focus_genes, "ident"))
head(plot_data)

write.table(plot_data, file = "01_p2S8_data.tsv", sep = "\t", quote = FALSE, row.names = TRUE)