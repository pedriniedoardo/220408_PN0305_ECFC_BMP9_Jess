# AIM ---------------------------------------------------------------------
# plot some genes

# libraries ---------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
# read in the table of DE
df_DE <- read_tsv("../../out/table/res_BMP9_vs_Mock_shr.txt")

# read in the GSEA result
df_GSEA <- read_tsv("../../out/table/df_tables_GSEA_res_BMP9_vs_Mock_shr_all_KEGG.tsv")

# read in the scaled data
vds_filter <- readRDS(file = "../../out/object/vds_all_filter.rds")

# read in the DESeq2 object 
data <- readRDS("../../out/object/dds_all_filter_DESeq.rds")

# wrangling ---------------------------------------------------------------
# read in the output frmo GSEA to identify the target genes of interest
df_GSEA %>%
  dplyr::filter(str_detect(pathway,pattern = "WNT"))
# the pathway do not seems to be specifically enriched from p value

# print the genes involved in the signature
le <- df_GSEA %>%
  dplyr::filter(str_detect(pathway,pattern = "WNT")) %>%
  pull(leadingEdge) %>%
  str_split("\\|") %>%
  unlist()

# print the stat of all the ledding edges
df_DE %>%
  dplyr::filter(symbol %in% le)

# since only two genes would cross the threshold of significance. To have more options I will focus on the genes with a lower FC
df_DE %>%
  dplyr::filter(symbol %in% le) %>%
  dplyr::filter(abs(log2FoldChange)>log2(1.5),padj<0.05)

# subset the gene for plotting
GOI <- df_DE %>%
  dplyr::filter(symbol %in% le) %>%
  dplyr::filter(abs(log2FoldChange)>log2(1.5),padj<0.05) %>%
  pull(symbol)

# plot genes with heatmap -------------------------------------------------
# index of the GOI
gene_id <- rownames(assay(vds_filter)) %in% GOI

# generate the matrix
mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

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

# row_ha <- rowAnnotation(class = rep(c("Apoptosis","Senescence","Fibrosis","SASP"),c(2,3,4,7)), 
#                         col = list(class = c("Apoptosis" = "violet", "Senescence" = "black","Fibrosis" = "yellow","SASP"="brown"))) 

ht2 <- Heatmap(mat2, 
               name = "exp",
               top_annotation = column_ha, 
               # cluster_rows = F, 
               # col = colorRamp2(c(-2, 0, 1), c("green", "white", "red")),
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7)),
               column_title = "KEGG_WNT_SIGNALING_PATHWAY") 

pdf("../../out/plot/heatmap_GOI_WNT.pdf",width = 6,height = 4) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# plot genes dotplot ------------------------------------------------------
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

# plot the data following the methods implemented in the plotCounts funciton from DESeq2
# Normalized counts plus a pseudocount of 0.5 are shown by default.
counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  dplyr::filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = "sample") %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=treat,y = count_norm_adj))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6)+facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/plot/boxplot_GOI_WNT.pdf",width = 4,height = 4) 

counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  dplyr::filter(symbol %in% GOI) %>%
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
ggsave("../../out/plot/scatterplot_GOI_WNT.pdf",width = 4,height = 4) 
