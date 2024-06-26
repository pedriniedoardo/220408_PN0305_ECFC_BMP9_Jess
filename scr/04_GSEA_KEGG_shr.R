# AIM ---------------------------------------------------------------------
# run GSEA unbiased with KEGG annotation

# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(patchwork)

# prepare the dataset with all the annoration needed ---------------------- 
# locate the files 
# file <- dir("out/") %>%
#   str_subset(pattern = "res_") %>%
#   str_subset(pattern = ".txt") %>%
#   str_subset(pattern = "shr",negate = T)
file <- dir("../../out/table/") %>%
  str_subset(pattern = "res_BMP9_vs_Mock_shr") %>%
  str_subset(pattern = ".txt") %>%
  str_subset(pattern = "shr",negate = F)
file 

# load the results 
results <- lapply(paste0("../../out/table/",file),function(x){
  read_tsv(x) 
}) %>%
  setNames(str_remove_all(file,pattern = ".txt"))

# GSEA -------------------------------------------------------------------- 
# use the FC dataset to create the ranked list of genes 
# Symbol or Entrez? 
list_ranks <- lapply(results, function(x){
  
  x <- dplyr::filter(x,!is.na(symbol)) %>%
    # average logFC in case of duplicated genenames
    group_by(symbol) %>%
    summarise(logFC = mean(log2FoldChange))
  
  ranks <- setNames(x$logFC, x$symbol)
  ranks
}) 
glimpse(list_ranks)

# score all the signatures in MsigDB from C2 category ---------------------
# library("msigdbr")
msigdbr_collections() %>%
  print(n=30)
#
# reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
# get all the C2 terms
# cgp_gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
# C2_gene_sets = msigdbr(species = "Mus musculus", category = "C2")
# h_gene_sets <- msigdbr(species = "Mus musculus", category = "H")
head(kegg_gene_sets)

# format in order to be accepted by GSEA
pathways <- split(x = kegg_gene_sets$gene_symbol, f = kegg_gene_sets$gs_name)
# head(pathways)

# RUN GSEA ----------------------------------------------------------------
df_tables_GSEA_all <- lapply(list_ranks, function(x){
  fgsea(pathways, x, minSize=10, maxSize=500)  
}) %>%
  bind_rows(.id = "dataset") %>% 
  # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list) 
  mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
    paste0(x,collapse = "|")
  }))) %>%
  arrange(padj,-abs(NES)) 

dim(df_tables_GSEA_all)

head(df_tables_GSEA_all,n=20) 

# save the whole table
df_tables_GSEA_all %>%
  write_tsv("../../out/table/df_tables_GSEA_res_BMP9_vs_Mock_shr_all_KEGG.tsv")

# explore the correlation between each subset focussing on the ES
# library(GGally)
# test_NES <- df_tables_GSEA_all %>%
#   dplyr::select(dataset,ES,pathway) %>%
#   pivot_wider(names_from = dataset,values_from = ES) %>%
#   column_to_rownames("pathway")
# 
# ggpairs(test_NES)

# COLLAPSE REDUNDANT ------------------------------------------------------
# collapsing the similar pathways 

# split the dataset per type
list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)

names(list_tables_GSEA_all)

list_collapsedPathways <- lapply(names(list_tables_GSEA_all),function(x){
  collapsePathways(list_tables_GSEA_all[[x]], pathways, list_ranks[[x]])
}) %>%
  setNames(names(list_tables_GSEA_all))

# collapsedPathways_MG <- collapsePathways(df_tables_GSEA_all %>%
#                                            filter(dataset %in% "list_df_comparisons_MG_old_young"), pathways, list_ranks$list_df_comparisons_MG_old_young)
# glimpse(collapsedPathways_MG)
str(list_collapsedPathways)

list_mainPathways <- pmap(list(list_tables_GSEA_all,list_collapsedPathways),function(x,y){
  x %>%
    dplyr::filter(pathway %in% y$mainPathways) %>%
    arrange(padj,-abs(NES)) %>%
    pull(pathway) 
})

str(list_mainPathways)

# save list of non redundant terms
# chackt the order of the names is the same
sum(!names(list_tables_GSEA_all) == names(list_mainPathways))

# filter only the non redundant fro each comparison
df_tables_GSEA_all_non_redundant <- 
  pmap(list(list_tables_GSEA_all,list_mainPathways),function(x,y){
    x %>%
      dplyr::filter(pathway %in% y)
  }) %>%
  bind_rows()

# save the table
df_tables_GSEA_all_non_redundant %>%
  write_tsv(file = "../../out/table/df_table_GSEA_res_BMP9_vs_Mock_shr_nonredundant_KEGG.tsv")

test <- df_tables_GSEA_all_non_redundant %>%
  group_by(dataset) %>%
  top_n(wt = padj*(-1),n = 5)

# test plot to show the main terms in each dataset
library(ggrepel)
df_tables_GSEA_all_non_redundant %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() + theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)
ggsave("../../out/plot/GSEA_unbiased_res_BMP9_vs_Mock_shr_nonredundant_KEGG.pdf",width = 10,height = 9)

# library(ggrepel)
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)
  # geom_hline(yintercept = -log(0.05),col="gray",linetype="dashed")
ggsave("../../out/plot/GSEA_unbiased_res_BMP9_vs_Mock_shr_KEGG.pdf",width = 10,height = 9)

# PLOT PROFILE ------------------------------------------------------------ 
# library("patchwork")
plot_pathway <- lapply(list_mainPathways$res_BMP9_vs_Mock_shr[1:9],function(x){ 
  name <- str_sub(x,start = 1,end = -4) 
  # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3) 
  plotEnrichment(pathways[[x]], list_ranks$res_BMP9_vs_Mock_shr) + labs(title = name) 
})

wrap_plots(plot_pathway)
ggsave(filename = "../../out/plot/GSEA_plot_profile_BMP9_vs_Mock_shr_nonredundant_KEGG.pdf",width = 15,height = 10)
# wrap_plots(plot_pathway) + ggsave(filename = "image/GSEA_plot_profile_CGP.pdf",width = 15,height = 10)
# wrap_plots(plot_pathway) + ggsave(filename = "image/GSEA_plot_profile_H.pdf",width = 15,height = 10)
