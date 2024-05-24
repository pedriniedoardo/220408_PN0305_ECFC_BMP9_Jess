# # libraries ---------------------------------------------------------------
# library(tidyverse)
# library(GSVA)
# library(limma)
# library(ComplexHeatmap)
# 
# # read in the data --------------------------------------------------------
# # read in the normalzied expression table geneXsample
# exp <- read_tsv("out/table_of_counts_normalzied_pseudobulk.txt") %>%
#   column_to_rownames("gene") %>%
#   as.matrix()
# 
# # read in the list of signatures
# file <- dir("../signatures/") %>%
#   # str_subset(pattern = "EPC_TEAM|FRIDMAN|ISM|MOSERLE") %>%
#   str_subset(pattern = "review.txt") %>%
#   str_sub(start = 1,end = -5)
# 
# library(GSEABase)
# pathways <- lapply(file, function(x){
#   name <- paste0("../signatures/",x,".txt")
#   read_tsv(name) %>% 
#     dplyr::select(human_gene) %>% 
#     drop_na() %>%
#     pull(human_gene)
#   
# }) %>%
#   setNames(file)
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
#            pathways,
#            min.sz=5,
#            max.sz=500,
#            kcdf="Poisson",
#            mx.diff=TRUE,
#            verbose=FALSE,
#            parallel.sz=1)
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
# pdf(file = "images/heatmap_GSVA_targeted_es_nonLog_ZScore.pdf", width = 7, height = 4)
# Heatmap(mat_norm,
#         # add annotation for the columns
#         # hide columns labels
#         # show_column_names = F,
#         # fix width of the lables
#         row_names_max_width = max_text_width(
#         rownames(mat),
#         gp = gpar(fontsize = 12)
#         ))
# dev.off()
# 
# # Heatmap(es,
# #         # add annotation for the columns
# #         # hide columns labels
# #         # show_column_names = F,
# #         # fix width of the lables
# #         row_names_max_width = max_text_width(
# #           rownames(mat),
# #           gp = gpar(fontsize = 12)
# #         ))
# #   