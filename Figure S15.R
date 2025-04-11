library(xlsx)

#################
## figure S15  ##
#################

dt = read.xlsx('/rsrch3/scratch/genomic_med/xyan4/Project/1_Jianfeng/2_subprojects/6_RMC_RCC/05_RMC_pre_post/03_result/03_bulk/01_bulk_gsea_result_REACTOME_MS_edit.xlsx',sheet='01_bulk_gsea_result_REACTOME')
cc_green = c('REACTOME_CELL_CYCLE_CHECKPOINTS','REACTOME_CELLULAR_RESPONSE_TO_STARVATION',
'REACTOME_DNA_REPLICATION','REACTOME_S_PHASE','REACTOME_SYNTHESIS_OF_DNA','REACTOME_G2_M_CHECKPOINTS',
'REACTOME_DNA_REPLICATION_PRE_INITIATION','REACTOME_CELL_CYCLE_MITOTIC','REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS','REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION'
)
isa_blue = c('REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL',
'REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS',
'REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR',
'REACTOME_DISEASES_OF_IMMUNE_SYSTEM',
'REACTOME_PD_1_SIGNALING',
'REACTOME_INTERFERON_GAMMA_SIGNALING'
)

pathways =c(cc_green,isa_blue)

dt = dt[c('ID','NES','pvalue','p.adjust','qvalues')]
dt1 = dt %>% filter(ID %in% pathways)
dt1 = mutate(dt1,
             Category = case_when(ID %in% cc_green ~ 'CelCycle',
                               ID %in% isa_blue ~ 'ImmuneSystem Activation',
                               TRUE ~ NA_character_))
write.table(dt, file = "01_data.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

df = dt1
old_names = df$ID
new_names = NULL
for (i in old_names){
  split_list <- strsplit(i, "_")[[1]]
  modified_string <- paste(tail(split_list, -1), collapse = "_")
  new_names = c(new_names, modified_string)
}
new_names


df$Name = new_names
p = ggplot(df, aes(reorder(Name, NES), NES)) +
  geom_col(aes(fill= Category)) + scale_fill_manual(values = c('#7DB46CFF','#ABD6DFFF','#E7EBE0FF')) +
  coord_flip() + 
  theme(panel.background = element_blank())+#,plot.background = element_blank()) +
  labs(x=NULL, y="Normalized Enrichment Score",title=NULL)

ggsave(plot = p, file = 'pS15.pdf', width = 10, height = 5)