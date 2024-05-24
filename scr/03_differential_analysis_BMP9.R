# libraries ---------------------------------------------------------------
source(file = "scr/00_library.R")

# read in the data --------------------------------------------------------
# ddsHTSeq <- readRDS("../../out/object/ddsHTSeq_BMP9.rds")
design <- readRDS("../../out/object/design_BMP9.rds")
vds_filter <- readRDS("../../out/object/vds_filter.rds")

# Pre-filtering -----------------------------------------------------------
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we increase the speed of the transformation and testing functions within DESeq2. Here we perform a minimal pre-filtering to keep only rows that have at least 11 reads total.

# keep <- rowSums(counts(ddsHTSeq)) >= 11
# table(keep)
# ddsHTSeq_filter <- ddsHTSeq[keep,]

# differential expression analyisis ---------------------------------------
# ddsHTSeq_filter <- DESeq(ddsHTSeq_filter)
# # if needed is possible to check the distributions of the counts before and after the normalizatoin
# # boxplot(log(counts(ddsHTSeq_structure,normalized = T)))
# # boxplot(log(counts(ddsHTSeq_structure,normalized = F)))
# #
# # save the filtered object
# saveRDS(ddsHTSeq_filter,"../../out/object/ddsHTSeq_filter.rds")
#
ddsHTSeq_filter <- readRDS("../../out/object/ddsHTSeq_filter.rds")

# print the contrast
resultsNames(ddsHTSeq_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast <- makeContrasts(BMP9vsMock = treatBMP9,
                          levels = design)

res_BMP9 <- results(ddsHTSeq_filter, contrast=contrast[,"BMP9vsMock"],alpha = 0.05)

summary(res_BMP9)

# add the gene symbols
list_df <- 
  list("BMP9_vs_Mock" = res_BMP9) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue) %>%
      # add the symbol
      mutate(ensembl = mapIds(org.Hs.eg.db,
                             keys = symbol,
                             column = "ENSEMBL",
                             keytype = "SYMBOL",
                             multiVals = "first")) %>%
      mutate(entrez = mapIds(org.Hs.eg.db,
                             keys = symbol,
                             column="ENTREZID",
                             keytype="SYMBOL",
                             multiVals="first"))
    
    df
  })

# save the tables
pmap(list(list_df,names(list_df)),function(x,y){
  name <- paste0("res_",y,".txt")
  # name
  write_tsv(x,file = paste0("../../out/table/",name))
})


# shrink ------------------------------------------------------------------  
res_BMP9_shr <- lfcShrink(ddsHTSeq_filter, res = res_BMP9, type = "ashr")
summary(res_BMP9_shr)

# add the gene symbols
list_df_shr <- 
  list("BMP9_vs_Mock_shr" = res_BMP9_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue) %>%
      # add the symbol
      mutate(ensembl = mapIds(org.Hs.eg.db,
                              keys = symbol,
                              column = "ENSEMBL",
                              keytype = "SYMBOL",
                              multiVals = "first")) %>%
      mutate(entrez = mapIds(org.Hs.eg.db,
                             keys = symbol,
                             column="ENTREZID",
                             keytype="SYMBOL",
                             multiVals="first"))
    
    df
  })

# save the tables
pmap(list(list_df_shr,names(list_df_shr)),function(x,y){
  name <- paste0("res_",y,".txt")
  # name
  write_tsv(x,file = paste0("../../out/table/",name))
})

# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
list_df$BMP9_vs_Mock %>%
  data.frame()%>%
  filter(baseMean>1)%>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  theme_bw()
ggsave("../../out/image/histogram_pvalue_BMP9.pdf",width = 4,height = 3)

# PLOTTING RESULTS --------------------------------------------------------
# add the info of the genename
test_plot <- list_df$BMP9_vs_Mock%>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0))

test_plot %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = test_plot[test_plot$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = test_plot[test_plot$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = test_plot[test_plot$col==1,],
    aes(label = symbol),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_text_BMP9.pdf",width = 12,height = 12)
#
test_plot_shr <- list_df_shr$BMP9_vs_Mock_shr%>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0))

test_plot_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = test_plot[test_plot$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = test_plot[test_plot$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = test_plot[test_plot$col==1,],
    aes(label = symbol),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_text_BMP9_shr.pdf",width = 12,height = 12)

# plotMA(res, ylim = c(-5, 5))
list_df$BMP9_vs_Mock %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) + 
  scale_x_log10() + scale_color_manual(values = c("gray","red")) + theme_bw() + 
  theme(legend.position = "none")
ggsave("../../out/image/MA_plot_BMP9.pdf",width = 4,height = 3)

# using the shrinked values
list_df_shr$BMP9_vs_Mock_shr %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color))+geom_point(alpha=0.2) + 
  scale_x_log10() + 
  scale_color_manual(values = c("gray","red")) + 
  theme_bw() +
  theme(legend.position = "none")
ggsave("../../out/image/MA_plot_shr_BMP9.pdf",width = 4,height = 3)

# the DEGs plot stringent
DEG_1 <- list_df$BMP9_vs_Mock %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  filter(col==1) %>%
  pull(symbol)

mat_filter <- assay(vds_filter) %>%
  data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

mat <- mat_filter[rownames(vds_filter) %in% DEG_1, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)
#

sample_ordered <- str_extract(colnames(mat2),pattern = "BMP9|mock")
column_ha <- HeatmapAnnotation(treat = sample_ordered,  
                               col = list(treat = c("mock" = "green", "BMP9" = "gray"))) 

ht2 <- Heatmap(mat2, 
               name = "exp", 
               column_title = "BMP9",
               row_names_gp = gpar(fontsize = 3),
               top_annotation = column_ha, 
               # cluster_rows = F, 
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7))
               ) 
pdf("../../out/image/heatmap_BMP9_DEG.pdf",width = 4,height = 15) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# the DEGs plot shr
DEG_2 <- list_df_shr$BMP9_vs_Mock_shr %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  filter(col==1) %>%
  pull(symbol)

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

mat_shr <- mat_filter[rownames(vds_filter) %in% DEG_2, ]
mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr)
#

sample_ordered_shr <- str_extract(colnames(mat2_shr),pattern = "BMP9|mock")
column_ha_shr <- HeatmapAnnotation(treat = sample_ordered_shr,  
                               col = list(treat = c("mock" = "green", "BMP9" = "gray"))) 

ht2_shr <- Heatmap(mat2_shr, 
               name = "exp", 
               column_title = "BMP9",
               row_names_gp = gpar(fontsize = 3),
               top_annotation = column_ha_shr
               # cluster_rows = F, 
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7))
) 
pdf("../../out/image/heatmap_BMP9_DEG_shr.pdf",width = 4,height = 15) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# PLOT DISPERSION ---------------------------------------------------------
pdf("../../out/image/ddsHTSeq_filter_dispersion.pdf",width = 5,height = 5) 
plotDispEsts(ddsHTSeq_filter)
dev.off()
