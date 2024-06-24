# AIM ---------------------------------------------------------------------
# run GSEA with a targeted set of genes

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
  
  x <- dplyr::filter(x,!is.na(symbol)) %>%
    # average logFC in case of duplicated genenames
    group_by(symbol) %>%
    summarise(logFC = mean(log2FoldChange))
  
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
pathways <- pathways_all[c("Arterial_vs_mural","Capillary_vs_mural","Tip_vs_mural","Venous_vs_mural","Arterial_vs_all","Capillary_vs_all","Venous_vs_all")]

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
pdf(file = "../../out/plot/heatmap_NES_targeted_BBB_res_BMP9_vs_Mock_shr.pdf", width = 6, height = 7)
draw(ht,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 70), "mm"))
dev.off()

dim(df_tables_GSEA_all)
head(df_tables_GSEA_all,n=20) 

# sasve the table
# run this to include P26 vs P08
# write_csv(df_tables_GSEA_all,file = "out/df_tables_GSEA_signature_targeted2.csv")
write_tsv(df_tables_GSEA_all,file = "../../out/table/df_tables_GSEA_res_BMP9_vs_Mock_shr_all_TARGETED_BBB.tsv")

# PLOT PROFILE ------------------------------------------------------------ 
# library("patchwork")

# for the H I did mae non redundant
plot_pathway <- lapply(unique(df_tables_GSEA_all$pathway),function(x){ 
  name <- str_sub(x,start = 1,end = -1) 
  # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3) 
  plotEnrichment(pathways[[x]], list_ranks$res_BMP9_vs_Mock_shr) + labs(title = name) 
})

wrap_plots(plot_pathway)
ggsave(filename = "../../out/plot/GSEA_plot_profile_BMP9_vs_Mock_shr_nonredundant_TARGETED_BBB.pdf",width = 15,height = 10)

# -------------------------------------------------------------------------
# plot the heatmap of the leading edges for endo cells
df_test <- df_tables_GSEA_all %>%
  dplyr::filter(dataset == "res_BMP9_vs_Mock_shr")

df_leading_edges <- df_test %>% 
  pull(leadingEdge) %>%
  str_split(pattern = "\\|") %>%
  setNames(df_test$pathway)

# pull the fold change data, I could pull the level of expression for the pseudobulk but I would have just one sample, therefore the row normalization would just make the plot useless
# filter the one in the leading edges that are significant
df_GOI <- pmap(list(df_leading_edges,names(df_leading_edges)),function(x,y){
  results$res_BMP9_vs_Mock_shr %>%
    dplyr::filter(symbol %in% x) %>%
    dplyr::filter(padj<0.05,abs(log2FoldChange)>log2(1.5)) %>%
    mutate(dataset = y) %>%
    dplyr::select(symbol,dataset)
}) %>%
  bind_rows() %>%
  mutate(dataset=str_remove_all(dataset,pattern = "_vs_all|_vs_mural")) %>%
  group_by(dataset,symbol) %>%
  summarise()

vds_filter <- readRDS(file = "../../out/object/vds_all_filter.rds")
# gene_id <- rownames(assay(vds_filter)) %in% df_leading_edges$signature_tip_goveia_2020
gene_id <- rownames(assay(vds_filter)) %in% df_GOI$symbol
mat <- assay(vds_filter)[gene_id, ]
# scale the rows and reorder the matrix in otder to match the gene order in teh df_GOI
mat2 <- ((mat - rowMeans(mat))/rowSds(mat))[df_GOI$symbol,]

# build the annotation object  
#
sample_ordered <- data.frame(sample = colnames(mat2)) %>%
  left_join(vds_filter@colData %>%
              data.frame(),by="sample")

# update the column name of the matrix
colnames(mat2) <- sample_ordered$sample_name

column_ha <- HeatmapAnnotation(treat = sample_ordered$treat,
                                   gender = sample_ordered$gender,
                                   col = list(treat = c("mock" = "gray", "BMP9" = "black"),
                                              gender = c("M" = "blue", "F" = "pink"))) 

row_ha <- rowAnnotation(class = df_GOI$dataset, 
                        col = list(class = c("Arterial" = "red","Capillary"="yellow","Tip"="green","Venous"="blue"))) 

ht2 <- Heatmap(mat2, 
               name = "exp",
               top_annotation = column_ha, 
               cluster_rows = F, 
               right_annotation = row_ha, 
               row_split = rep(c(1,2,3,4),c(8,16,3,7)),
               # cluster_rows = F, 
               # col = colorRamp2(c(-2, 0, 1), c("green", "white", "red")),
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7)),
               column_title = "BBB YANG22 sig leading edges (FC1.5 padj<0.05)") 

pdf("../../out/plot/heatmap_leadingEdges_BBB_shr.pdf",width = 6,height = 8) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# plot genes dotplot ------------------------------------------------------
# read in the DESeq2 object 
data <- readRDS("../../out/object/dds_all_filter_DESeq.rds")

# generate the look up table
lut <- data@colData %>%
  data.frame()

# generate the tabel of milion reads
MR <- counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  pivot_longer(names_to = "sample",values_to = "exp",-symbol)%>%
  group_by(sample)%>%
  summarise(MR = sum(exp)/10^6)

counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  dplyr::filter(symbol %in% df_GOI$symbol) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = "sample") %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=treat,y = count_norm_adj,group=clone,col=clone))+geom_point(alpha=0.6)+
  geom_line() +
  facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/plot/scatterplot_GOI_BBB.pdf",width = 15,height = 15) 
