# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(UpSetR)
library(gplots)
library(scales)
library(pals)

# plot genes with heatmap -------------------------------------------------
vds_filter <- readRDS(file = "../../out/object/vds_all.rds")

# # read in the offtargets
# genes_B <- read_csv("data/offtarget_B.csv") %>%
#   pull(Gene)
# 
# genes_D <- read_csv("data/offtarget_D.csv") %>%
#   pull(Gene)
# 
# list_genes <- list("offset_B" = genes_B,"offset_D" = genes_D)
# 
# pdf("out/image/upset_offset_genes.pdf",width = 5,height = 5)
# upset(fromList(list_genes), order.by = "freq") 
# dev.off()
# 
# list_genes2 <- venn(list_genes, show.plot=FALSE)
# list_genes3 <- attr(list_genes2,"intersections")
# 
# # how many genes are in the dataset
# lapply(list_genes3, function(x){
#   sum(rownames(assay(vds_filter)) %in% x)
# })

gene_id <- rownames(assay(vds_filter)) %in% c("XIST","DDX3Y","RPS4Y1","USP9Y")
# gene_id <- rownames(assay(vds_filter)) %in% unlist(list_genes3)
# notice that for both dataset the design didn't affect mucht he normalizationo,therefore the topmost variable genes are the same for both
# I will prduce the plot for just one dataset
# The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene’s average across all samples. Hence, we center each genes’ values across samples, and plot a heatmap (figure below). We provide a data.frame that instructs the pheatmap function how to label the columns.
# mat <- assay(vds_filter)[topVarGenes, ]
mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat,useNames = TRUE)

# change the rownames to match the symbol 
# rownames(mat2) <- rownames(mat2) %>% 
#   data.frame(ensembl=.) %>% 
#   left_join(DEG_1,by = "ensembl") %>% 
#   pull(symbol) 

# build the annotation object  
meta <- colData(vds_filter) %>% 
  data.frame()

LUT_sample_matrix <- data.frame(sample_matrix = colnames(mat2)) %>%
  left_join(meta,by = c("sample_matrix"="sample"))

# sample_ordered <- LUT_sample_matrix$clone

# define a LUT for the clones
color_id <- alphabet(length(unique(LUT_sample_matrix$clone)))
# check the colors
show_col(color_id)

# build the named vector
names(color_id) <- unique(LUT_sample_matrix$clone)

column_ha <- HeatmapAnnotation(clone = LUT_sample_matrix$clone,
                               gender = LUT_sample_matrix$gender,
                               BMP_treat = LUT_sample_matrix$treat,
                               col = list(clone = color_id,
                                          condition = c("M" = "blue", "F" = "pink"),
                                          BMP_treat = c("BMP9" = "yellow", "mock" = "purple"))) 

# change the name of the column in the matrix
colnames(mat2) <- LUT_sample_matrix$sample_name

# row_ha <- rowAnnotation(class = rep(c("common B_D","offtarget B","offset D"),c(22,3,24)),
#                         col = list(class = c("common B_D" = "violet", "offtarget B" = "yellow","offset D"="brown")))

ht2 <- Heatmap(mat2, 
               name = "exp",
               top_annotation = column_ha, 
               # cluster_rows = F, 
               # col = colorRamp2(c(-2, 0, 1), c("green", "white", "red")),
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7)),
               column_title = "gender genes") 

pdf("../../out/plot/heatmap_GOI_gender.pdf",width = 8,height = 5) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()