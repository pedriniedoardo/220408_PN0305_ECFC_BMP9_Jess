# # libraries ---------------------------------------------------------------
# library(tidyverse)
# library(GSVA)
# library(limma)
# library(ComplexHeatmap)
# library(org.Hs.eg.db)
# 
# # read in the data --------------------------------------------------------
# # read in the normalzied expression table geneXsample
# exp <- read_tsv("out/ddsHTSeq_filter_company_tot_counts_norm.tsv") %>%
#   dplyr::select("ensembl",contains("4WK")) %>%
#   column_to_rownames("ensembl") %>%
#   as.matrix()
# 
# colnames(exp) <- str_sub(colnames(exp),start = 18,end = -1)
# 
# # read in the signatures
# file_sig <- dir("signatures/signatures_pietro/") %>%
#   # str_subset(pattern = "EPC_TEAM|FRIDMAN|ISM|MOSERLE") %>%
#   str_subset(pattern = ".gmt") %>%
#   str_sub(start = 1,end = -5)
# 
# pathways_gmt <- lapply(file_sig, function(x){
#   name <- paste0("signatures/signatures_pietro/",x,".gmt")
#   gene <- getGmt(name) %>%
#     geneIds() %>%
#     .[[1]]
#     
#   mapIds(org.Hs.eg.db,
#          keys = gene,
#          column = "ENSEMBL",
#          keytype = "SYMBOL",
#          multiVals = "first")
# }) %>%
#   setNames(file_sig)
# 
# # merge all the pathways in a single list
# pathways <- c(pathways_gmt)
# 
# # -------------------------------------------------------------------------
# # perform the GSEA based on the normalized counts
# es <- gsva(exp,
#            pathways,
#            min.sz=5,
#            max.sz=500,
#            kcdf="Poisson",
#            mx.diff=TRUE,
#            verbose=FALSE,
#            parallel.sz=1)
# 
# es_log <- gsva(log(exp+1),
#                pathways,
#                min.sz=5,
#                max.sz=500,
#                kcdf="Poisson",
#                mx.diff=TRUE,
#                verbose=FALSE,
#                parallel.sz=1)
# 
# # -------------------------------------------------------------------------
# # library(ComplexHeatmap)
# mat <- es
# mat_norm <- es %>%
#   data.frame() %>%
#   rownames_to_column() %>%
#   # filter(rowname %in% rownames(DEgeneSets)) %>%
#   # scale the values rowwise
#   gather(key = sample,value = exp,-rowname) %>%
#   group_by(rowname) %>%
#   mutate(norm = (exp - mean(exp))/sd(exp)) %>%
#   dplyr::select(-exp) %>%
#   spread(key = sample,value = norm) %>%
#   column_to_rownames()
# 
# sample_ordered <- str_extract(colnames(mat_norm),pattern = "DG|CT") 
# sample_ordered
# 
# # build the annotation object 
# column_ha <- HeatmapAnnotation(treat = sample_ordered, 
#                                col = list(treat = c("DG" = "gray", "CT" = "green")))  
# 
# hm <- Heatmap(mat_norm,
#               # add annotation for the columns
#               # hide columns labels
#               # show_column_names = F,
#               # fix width of the lables
#               top_annotation = column_ha,
#               row_names_max_width = max_text_width(
#                 rownames(mat),
#                 gp = gpar(fontsize = 12)
#               ))
# 
# pdf(file = "images/heatmap_GSVA_targeted_es_nonLog_ZScore_DG.pdf", width = 12, height = 4)
# draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm")) 
# dev.off()
