# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

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
  
  x <- filter(x,!is.na(symbol)) %>%
    # average logFC in case of duplicated genenames
    group_by(symbol) %>%
    summarise(logFC = log2FoldChange)
  
  ranks <- setNames(x$logFC, x$symbol)
  ranks
}) 
glimpse(list_ranks)

# # read in the signatures
# file <- dir("data/signatures/") %>%
#   # str_subset(pattern = "EPC_TEAM|FRIDMAN|ISM|MOSERLE") %>%
#   str_subset(pattern = "review.txt") %>%
#   str_sub(start = 1,end = -5)

# read in the tailored pathways -------------------------------------------
# for reference to the signatures see the file "scr/00_build_siganture.R"
pathways_all <- readRDS("../../data/signatures/list_costume_signatures_pathways.rds")
pathways <- pathways_all[c("HALLMARK_ANGIOGENESIS","GOBP_SPROUTING_ANGIOGENESIS","signature_tip_goveia_2020")]

# RUN GSEA ----------------------------------------------------------------
df_tables_GSEA_all <- lapply(list_ranks, function(x){
  # fgsea(pathways, x, minSize=1, maxSize=500, nperm=1000)
  fgsea(pathways, x, minSize=1, maxSize=500)  
}) %>%
  bind_rows(.id = "dataset") %>% 
  # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list) 
  mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
    paste0(x,collapse = "|")
  }))) %>%
  arrange(padj,-abs(NES)) 

# generate an heatmap for the comparison ofthe nes across datasets and pathways
mat_dataset <- df_tables_GSEA_all %>%
  dplyr::select(dataset,pathway,NES) %>%
  mutate(pathway = str_remove(pathway,pattern = "_review")) %>%
  pivot_wider(names_from = dataset,values_from = NES) %>%
  column_to_rownames(var = "pathway")

ht <- Heatmap(mat_dataset,name = "NES", column_title = "GSEA",
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

# run this to include P26 vs P08
# pdf(file = "images/heatmap_NES_targeted2.pdf", width = 7, height = 5)
pdf(file = "../../out/image/heatmap_NES_targeted_res_BMP9_vs_Mock_shr.pdf", width = 7, height = 7)
draw(ht,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 70), "mm"))
dev.off()

dim(df_tables_GSEA_all)
head(df_tables_GSEA_all,n=20) 

# sasve the table
# run this to include P26 vs P08
# write_csv(df_tables_GSEA_all,file = "out/df_tables_GSEA_signature_targeted2.csv")
write_tsv(df_tables_GSEA_all,file = "../../out/table/df_tables_GSEA_res_BMP9_vs_Mock_shr_all_TARGETED.tsv")

# PLOT PROFILE ------------------------------------------------------------ 
# library("patchwork")

# for the H I did mae non redundant
plot_pathway <- lapply(unique(df_tables_GSEA_all$pathway),function(x){ 
  name <- str_sub(x,start = 1,end = -4) 
  # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3) 
  plotEnrichment(pathways[[x]], list_ranks$res_BMP9_vs_Mock_shr) + labs(title = name) 
})

wrap_plots(plot_pathway)
ggsave(filename = "../../out/image/GSEA_plot_profile_BMP9_vs_Mock_shr_nonredundant_TARGETED.pdf",width = 15,height = 5)

# -------------------------------------------------------------------------
# plot the heatmap of the leading edges for endo cells
df_test <- df_tables_GSEA_all %>%
  filter(dataset == "res_BMP9_vs_Mock_shr")

df_leading_edges <- df_test %>% 
  pull(leadingEdge) %>%
  str_split(pattern = "\\|") %>%
  setNames(df_test$pathway)

# pull the fold change data, I could pull the level of expression for the pseudobulk but I would have just one sample, therefore the row normalization would just make the plot useless
# filter the one in the leading edges that are significant
subset_genes <- results$res_BMP9_vs_Mock_shr %>%
  filter(symbol %in% df_leading_edges$signature_tip_goveia_2020) %>%
  filter(padj<0.05,abs(log2FoldChange)>1) %>%
  pull(symbol)

vds_filter <- readRDS(file = "../../out/object/vds_filter.rds")
# gene_id <- rownames(assay(vds_filter)) %in% df_leading_edges$signature_tip_goveia_2020
gene_id <- rownames(assay(vds_filter)) %in% subset_genes
mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

# build the annotation object  
sample_ordered <- str_extract(colnames(mat2),pattern = c("BMP9|mock"))
column_ha <- HeatmapAnnotation(treat = sample_ordered,  
                               col = list(treat = c("mock" = "green", "BMP9" = "gray"))) 

# row_ha <- rowAnnotation(class = rep(c("Apoptosis","Senescence","Fibrosis","SASP"),c(2,3,4,7)), 
#                         col = list(class = c("Apoptosis" = "violet", "Senescence" = "black","Fibrosis" = "yellow","SASP"="brown"))) 

ht2 <- Heatmap(mat2, 
               name = "exp",
               top_annotation = column_ha, 
               # cluster_rows = F, 
               # col = colorRamp2(c(-2, 0, 1), c("green", "white", "red")),
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7)),
               column_title = "GOVEIA sig leading edges") 

pdf("../../out/image/heatmap_leadingEdges_goveia_shr.pdf",width = 6,height = 6) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# plot heatmap with rank of genes -----------------------------------------
# get the top and the bottom genes of the signatures in the rank
df_rank <- results$res_BMP9_vs_Mock_shr %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(rank = nrow(.):1) %>%
  mutate(rank2 = 1:nrow(.))

# filter the ranks of the gene in the signature

id_UP <- df_rank %>%
  filter(symbol %in% pathways$HALLMARK_ANGIOGENESIS)
head(id_UP)

# library(ComplexHeatmap)
m = matrix(df_rank$rank,ncol = 1)
ha_up = rowAnnotation(foo = anno_mark(at = id_UP$rank2, labels = id_UP$symbol))
hm_up <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha_up)

pdf("../../out/image/heatmap_id_HALLMARK_ANGIOGENESIS_BMP9_shr.pdf",width = 4.5,height = 20) 
draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

#
# get the top and the bottom genes of the signatures in the rank
df_rank <- results$res_BMP9_vs_Mock_shr %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(rank = nrow(.):1) %>%
  mutate(rank2 = 1:nrow(.))

# filter the ranks of the gene in the signature

id_UP <- df_rank %>%
  filter(symbol %in% pathways$GOBP_SPROUTING_ANGIOGENESIS)
head(id_UP)

# library(ComplexHeatmap)
m = matrix(df_rank$rank,ncol = 1)
ha_up = rowAnnotation(foo = anno_mark(at = id_UP$rank2, labels = id_UP$symbol))
hm_up <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha_up)

pdf("../../out/image/heatmap_id_GOBP_SPROUTING_ANGIOGENESIS_BMP9_shr.pdf",width = 4.5,height = 20) 
draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()
