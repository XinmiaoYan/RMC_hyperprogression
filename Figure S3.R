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

dir='03_FigureS3'
dir.create(dir)
setwd(dir)

load('/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/16_source_data/03_FigureS3/01_obj_TNK.Rdata')
obj = obj_TNK
###########
## pS3a  ##
###########

my_color = colorRampPalette(brewer.pal(12, "Paired"))(12)
pS3a = DimPlot(obj, label = T, reduction = 'umap')
ggsave("01_umap_pS3a.pdf",plot=pS3a, width = 5, height = 4)
# table(obj$seurat_clusters)

# Extract UMAP coordinates
umap_coords <- Embeddings(obj, reduction = "umap")
# Extract metadata
metadata <- obj@meta.data
identical(rownames(umap_coords), rownames(metadata))
# Combine UMAP and metadata
output_df <- cbind(UMAP_1 = umap_coords[, 1],
                   UMAP_2 = umap_coords[, 2],
                   metadata)

head(output_df)
# Write to a tab-separated file
write.table(output_df, file = "01_umap_with_metadata.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

###########
## pS3b  ##
###########
annotation_gene = c('CD3D', 'CD3E', 'CD3G', #T     
                    'CD4', 'CD40LG', 'FOXP3',#CD4
                    'CD8A','CD8B', #CD8
                    'NCAM1', 'NKG7', 'GNLY', # NK
                    'CCR7','TCF7',
                    'IL7R', 'MKI67'#ils
                    )
pS3b = DotPlot(obj, features = unique(annotation_gene), cluster.idents = T) +
  RotatedAxis() + #scale_colour_gradientn(colors=c('royalblue','royalblue','springgreen','yellow',"red")) +
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)
ggsave("02_pS3b.pdf",plot=pS3b, width = 5, height = 4)

dot_data <- pS3b$data
# Preview the data
head(dot_data)
# Write to a tab-separated file
write.table(dot_data, file = "02_pS3b_data.tsv", sep = "\t", quote = FALSE, row.names = TRUE)


###########
## pS3c  ##
###########

plots = list()
for (i in annotation_gene){
  plots[[i]] = FeaturePlot(obj, features = i,order = TRUE) + NoLegend()
}

pS3c = wrap_plots(plots = plots, nrow =3)
ggsave("03_pS3c.pdf",plot=pS3c, width = 10, height =6)

###########
## pS3d  ##
###########
Idents(object = obj) <- "level2_cell"
levels(obj)  = c('CD4+T','CD8+T','NK','ProliferatingT')

my_color = colorRampPalette(brewer.pal(12, "Paired"))(6)
pS3d = DimPlot(obj, label = F, reduction = 'umap', cols = my_color)
ggsave("04_pS3d.pdf",plot=pS3d, width = 5, height = 3.5)

###########
## pS3e  ##
###########

pS3e = DotPlot(obj, features =unique(annotation_gene), cluster.idents = F) +
  RotatedAxis() + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)
ggsave("05_pS3e.pdf",plot=pS3e, width = 5, height = 2)

# output data
# Extract the data used for the plot
dot_data <- pS3e$data
# Preview the data
head(dot_data)
# Write to a tab-separated file
write.table(dot_data, file = "05_pS3e_data.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
