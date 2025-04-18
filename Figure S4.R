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

dir='04_FigureS4'
dir.create(dir)
setwd(dir)

load('/rsrch9/home/genomic_med/lwang22_lab/Xinmiao/6_RMC_RCC/06_RMC_revision_240415/03_result/02_whole_umap/06_subset_240422/01_T/00_obj_TNK_5518.Rdata')

###########
## pS4a  ##
###########

df = read.table('/rsrch3/scratch/genomic_med/xyan4/Project/1_Jianfeng/2_subprojects/6_RMC_RCC/03_result/8_sub_cell/05_Epithelial/04_identify_malignant_cell/01_kmeans_df_s_new_k2.txt',row.names = 1)
obj_file_path = '/rsrch3/scratch/genomic_med/xyan4/Project/1_Jianfeng/2_subprojects/6_RMC_RCC/03_result/8_sub_cell/05_Epithelial/03_InferCNV_15s/01_infercnv'

#start
out_path = './'
infercnv_obj = readRDS(paste(obj_file_path,"/run.final.infercnv_obj", sep=''))

expr <- infercnv_obj@expr.data
normal_loc <- as.numeric(unlist(infercnv_obj@reference_grouped_cell_indices))
test_loc <- as.numeric(unlist(infercnv_obj@observation_grouped_cell_indices))


df = df %>% filter(treatment_condition %in% c('RMC_Baseline', 'RMC_Post'))
df = mutate(df,
            sample_name = case_when(orig.ident == '0650101' ~ 'RMC46_Baseline',
                                    orig.ident == '0650201' ~ 'RMC49_Baseline',
                                    orig.ident == '0650401_S7' ~ 'RMC56_Baseline',
                                    orig.ident == '2010101_S8' ~ 'RMC53_Baseline',
                                    orig.ident == '2010201_S10' ~ 'RMC60_Baseline',
                                    orig.ident == '0650801_S14' ~ 'RMC71_Baseline',
                                    orig.ident == '0650901_S15' ~ 'RMC66_Baseline',
                                    orig.ident == '0650501_S9' ~ 'RMC57_Post',
                                    orig.ident == '2010302_S13' ~ 'RMC61_Post',
                                    TRUE ~ NA_character_))
df = mutate(df,
            sample = case_when(orig.ident == '0650101' ~ '0650101',
                               orig.ident == '0650201' ~ '0650201',
                               orig.ident == '0650401_S7' ~ '0650401',
                               orig.ident == '2010101_S8' ~ '2010101',
                               orig.ident == '2010201_S10' ~ '2010201',
                               orig.ident == '0650801_S14' ~ '0650801',
                               orig.ident == '0650901_S15' ~ '0650901',
                               orig.ident == '0650501_S9' ~ '0650501',
                               orig.ident == '2010302_S13' ~ '2010302',
                               TRUE ~ NA_character_))

df = mutate(df,
            Patient = case_when(sample == '0650101' ~ 'RMC46',
                                sample == '0650201' ~ 'RMC49',
                                sample == '0650401' ~ 'RMC56',
                                sample == '2010101' ~ 'RMC53',
                                sample == '2010201' ~ 'RMC60',
                                sample == '0650801' ~ 'RMC71',
                                sample == '0650901' ~ 'RMC66',
                                sample == '0650501' ~ 'RMC57',
                                sample == '2010302' ~ 'RMC61',
                                TRUE ~ NA_character_))
df = mutate(df,
            BiopsySite = case_when(sample == '0650101' ~ 'LymphNode',
                                   sample == '0650201' ~ 'LymphNode',
                                   sample == '0650401' ~ 'Liver',
                                   sample == '2010101' ~ 'Kidney',
                                   sample == '2010201' ~ 'Lung',
                                   sample == '0650801' ~ 'Kidney',
                                   sample == '0650901' ~ 'Liver',
                                   sample == '0650501' ~ 'Liver',
                                   sample == '2010302' ~ 'Liver',
                                   TRUE ~ NA_character_))
df = mutate(df,
            Treatment = case_when(sample == '0650101' ~ 'Baseline',
                                  sample == '0650201' ~ 'Baseline',
                                  sample == '0650401' ~ 'Baseline',
                                  sample == '2010101' ~ 'Baseline',
                                  sample == '2010201' ~ 'Baseline',
                                  sample == '0650801' ~ 'Baseline',
                                  sample == '0650901' ~ 'Baseline',
                                  sample == '0650501' ~ 'PostNI',
                                  sample == '2010302' ~ 'PostNI',
                                  TRUE ~ NA_character_))

table(df$orig.ident, df$Patient)
table(df$treatment_condition, df$Treatment)
table(df$orig.ident, df$BiopsySite)
table(df$orig.ident, df$Treatment)

df = df %>% filter(sample %in% names(table(obj$orig.ident)))
df1 = df

# expr
# start plot the figure
anno.df=data.frame(
CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
class=c(rep("T",length(normal_loc)),rep("Epithelia",length(test_loc)))
)
# head(anno.df)
gn <- rownames(expr)
geneFile <- read.table("/rsrch3/scratch/genomic_med/xyan4/Project/singlecell_reference/3_chr_final.txt",header = F,sep = "\t",stringsAsFactors = F)
# head(geneFile)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
# head(sub_geneFile,4)
# expr[1:4,1:4]

expr1 = expr[,rownames(df1)]
dt = df1 %>% filter(kmeans_class==1)
table(dt$kmeans_class, dt$class)
# colnames(dt) = c("kmeans_class","class", "Patient","sample", "BiopsySite",'Treatment')

dt = dt[order(dt$BiopsySite),]
dt <- dt[order(dt$Treatment),]
dt <- dt[order(dt$Patient),]

dt = arrange(dt, Treatment, BiopsySite, Patient)
expr = expr[,rownames(dt)]

#define annotation style
pa = colorRampPalette(brewer.pal(9, "Paired"))(9)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
c6 = colorRampPalette(brewer.pal(8, 'Set2'))(8)[1:6]
left_anno <- rowAnnotation(df = dt[,c(11,4,10)],
                           col=list(Treatment = c('Baseline'="#a8adb4", 'PostNI'="#FFD966"),
                                   BiopsySite = c('Kidney' = '#f78e31', 'Liver' = '#f0c808', 'LymphNode'='#90be6d'),
                                   Patient = c('RMC46'=pa[1],'RMC49'=pa[2],'RMC53'=pa[3],'RMC57'=pa[4],'RMC61'=pa[5],'RMC66'=pa[6],'RMC71'=pa[7])))

#plot
pdf(paste(out_path,"/01_heatmap_pS4a.pdf",sep=''),width = 10,height = 5)
ht = Heatmap(t(expr)[rownames(dt),], #绘图数据的CB顺序和注释CB顺序保持一致
            col = colorRamp2(c(0.85,1.0,1.15), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
            cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
            column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
            column_gap = unit(2, "mm"),
            heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
            top_annotation = top_anno,left_annotation = left_anno, #添加注释
            row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

#calculate CNV score
###########
## pS4b  ##
###########
df1
## remove 2 from kmeans_class adn relabel 3 into 2
df2 = df1 %>% filter(kmeans_class %in% c(1,2))
df2

expr2 = expr1[,rownames(df2)]
expr2=expr2-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"

CNV_score$CB=rownames(CNV_score)
df2$CB=rownames(df2)
CNV_score=CNV_score%>%inner_join(df2,by="CB")
CNV_score$kmeans_class = as.character(CNV_score$kmeans_class)

color_v=c('#7DB46CFF','#ABD6DFFF')
CNV_score%>%ggplot(aes(kmeans_class,log(CNV_score)))+geom_violin(aes(fill=kmeans_class),color="NA")+
  scale_fill_manual(values = color_v)+ylim(c(-9, -5))+
  theme_bw()
ggsave(paste(out_path,"/02_CNV_level_new_k2.pdf",sep=''),width = 7,height = 5,units = "cm")

####
obj
epi_dt = (df2 %>% filter(class=='Epithelia'))[colnames(obj),]
library(xlsx)
write.xlsx(epi_dt, file = '03_class_result.xlsx', sheetName = "Sheet1", row.names = T)
