# AIM ---------------------------------------------------------------------
# run some preprocessing and QC of the dataset to understand if there are outliers or global trends.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(vsn)
library(hexbin)
library(viridis)
library(pheatmap)
library(PoiClaClu)
library(AnnotationDbi) 
library(AnnotationHub)
library(GGally)
library(RNAseqQC)
library(edgeR)

# read in the data --------------------------------------------------------
ddsHTSeq <- readRDS("../../out/object/dds_all.rds")

# remove low expressed genes ----------------------------------------------
# To reduce the noise filter out the lowly expressed genes
# plot the total counts per sample
colSums(counts(ddsHTSeq)) %>%
  data.frame(tot_counts=.) %>%
  rownames_to_column("sample") %>% 
  # add the metadata
  left_join(colData(ddsHTSeq) %>% 
              data.frame(),by = "sample") %>% 
  ggplot(aes(x=sample_name,y = tot_counts))+geom_col()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = margin(0.5, 0.5, 2, 2, "cm"))
ggsave("../../out/plot/barplot_tot_count_all.pdf",width = 7,height = 5)

# check the total size of the matrix before filtering fow low expressed genes
nrow(ddsHTSeq)

# remove potential non infirmative genes (lowly expressed genes).
# use the automatic implementation provided by edgeR
ddsHTSeq_filter <- ddsHTSeq[edgeR::filterByExpr(counts(ddsHTSeq), group = colData(ddsHTSeq)$treat),]

# this is the old implementaion 
# ddsHTSeq_filter <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 100, ]

# check the size of the dataset after filtering for lowly expressed genes
nrow(ddsHTSeq_filter)

# save the filtered object
saveRDS(ddsHTSeq_filter,file = "../../out/object/dds_all_filter.rds")

# scaling transformation of the data --------------------------------------
# vsd
vds_filter <- vst(ddsHTSeq_filter, blind = F)
head(assay(vds_filter), 3)

# just of gender identification, generate also the scaled version of the object before filtering
vds_gender <- vst(ddsHTSeq, blind = F)

# rlog
rld_filter <- rlog(ddsHTSeq_filter, blind = F)
head(assay(rld_filter), 3)

# diagnostic plot
meanSdPlot_vsd <- meanSdPlot(assay(vds_filter))
ggsave(plot = meanSdPlot_vsd$gg+theme_bw(),filename = "../../out/plot/meanSdPlot_vsd.pdf",width = 4,height = 4)

meanSdPlot_rlog <- meanSdPlot(assay(rld_filter))
ggsave(plot = meanSdPlot_rlog$gg+theme_bw(),filename = "../../out/plot/meanSdPlot_rlog.pdf",width = 4,height = 4)

ddsHTSeq_filter <- estimateSizeFactors(ddsHTSeq_filter)

# sample distance ---------------------------------------------------------
# using vsd scaling
sampleDists_vsd <- dist(t(assay(vds_filter)))
sampleDists_vsd

head(assay(vds_filter))

sampleDistMatrix_vsd <- as.matrix(sampleDists_vsd)

rownames(sampleDistMatrix_vsd) <- paste(vds_filter$sample_name)
colnames(sampleDistMatrix_vsd) <- NULL

map_colors<-colorRampPalette(viridis(12))(255)

hm_1 <- pheatmap(sampleDistMatrix_vsd,
                 clustering_distance_rows = sampleDists_vsd,
                 clustering_distance_cols = sampleDists_vsd,
                 col = map_colors)

pdf("../../out/plot/heatmap_vsd.pdf",width = 7,height = 5)
hm_1
dev.off()

# using rlog scaling
sampleDists_rld <- dist(t(assay(rld_filter)))
sampleDists_rld

head(assay(rld_filter))

sampleDistMatrix_rld <- as.matrix(sampleDists_rld)

rownames(sampleDistMatrix_rld) <- paste(rld_filter$sample_name)
colnames(sampleDistMatrix_rld) <- NULL

hm_1_2 <- pheatmap(sampleDistMatrix_rld,
                   clustering_distance_rows = sampleDists_rld,
                   clustering_distance_cols = sampleDists_rld,
                   col = map_colors)

pdf("../../out/plot/heatmap_rld.pdf",width = 7,height = 5)
hm_1_2
dev.off()

# using poisson distance
poisd <- PoissonDistance(t(counts(ddsHTSeq_filter,normalized = F)))

samplePoisDistMatrix <- as.matrix(poisd$dd)

rownames(samplePoisDistMatrix) <- paste(ddsHTSeq_filter$sample_name)
colnames(samplePoisDistMatrix) <- NULL

hm_p <- pheatmap(samplePoisDistMatrix,
                 clustering_distance_rows = poisd$dd,
                 clustering_distance_cols = poisd$dd,
                 col = map_colors)

pdf("../../out/plot/heatmap_poisd.pdf",width = 7,height = 5)
hm_p
dev.off()

# plot cluster alternative ------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/plot/heatmap_cluster_all.pdf",width = 7,height = 5)
set.seed(2)
plot_sample_clustering(vds_filter,
                       anno_vars = c("clone","treat","gender"),
                       distance = "euclidean")
dev.off()

# PCA plot ----------------------------------------------------------------
plot_vsd <- plotPCA(vds_filter,
                    intgroup = c("clone","treat","gender")) +
  theme_bw()

plot_vsd$data %>%
  ggplot(aes(x=PC1,y=PC2,label = clone)) +
  # ggplot(aes(x=PC1,y=PC2,col=BMP_treat,shape=condition)) +
  geom_point(aes(col=treat,shape=gender),size =3) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../../out/plot/PCA_vsd_colorBMP9.pdf",width = 5,height = 4)

# pull more PC
rv <- rowVars(assay(vds_filter),useNames = TRUE)
# select the ntop genes by variance
select_var <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
test <- prcomp(t(assay(vds_filter)[select_var,]))$x %>% 
  data.frame() %>% 
  rownames_to_column("sample")

# plot more PC by condition
left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample")) %>%
  # ggpairs(columns = 5:14,ggplot2::aes(colour=condition),upper = "blank")+
  ggpairs(columns = 6:15,ggplot2::aes(colour=treat),upper = "blank") +
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/PCA_panel_vsd_BMP9.pdf",width = 20,height = 20)

# explore pc score by metadata fo the samples
test2 <- left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample"))

test_df1 <- test2 |> 
  dplyr::select(clone,treat,name) |> 
  # dplyr::rename(approx.time = Approx.Time.between.collection.and.processing) |> 
  pivot_longer(names_to = "var_1",values_to = "value_1",-c(name))

test_df2 <- test2 |> 
  dplyr::select(name,PC1:PC9) |>
  pivot_longer(names_to = "var_2",values_to = "value_2",-c(name))

left_join(test_df1,test_df2,by=c("name")) |> 
  mutate(comparison = paste0(var_1,"_vs_",var_2)) |> 
  ggplot(aes(x=value_1,y=value_2)) +
  facet_wrap(~comparison,scales = "free",ncol=9) +
  # facet_grid(var_1~var_2,scales = "free_x",drop = T) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/panel_metadata_PC.pdf",width = 20,height = 4)

# plot PCA also fro rlog scaled data
plot_rld <- plotPCA(rld_filter,
                    intgroup = c("clone","treat","gender")) +
  theme_bw()

plot_rld$data %>%
  ggplot(aes(x=PC1,y=PC2,label = clone)) +
  # ggplot(aes(x=PC1,y=PC2,col=BMP_treat,shape=condition)) +
  geom_point(aes(col=treat,shape=gender),size =3) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../../out/plot/PCA_rld_colorBMP9.pdf",width = 5,height = 4)

# MSD plot ----------------------------------------------------------------
mds_vsd <- as.data.frame(colData(vds_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_vsd))

ggplot(mds_vsd, aes(x = `1`, y = `2`, color = treat)) +
  geom_point(size = 3)+ theme_bw()
ggsave("../../out/plot/MDS_vsd.pdf",width = 5,height = 4)

mds_rld <- as.data.frame(colData(rld_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_rld))

ggplot(mds_rld, aes(x = `1`, y = `2`, color = treat)) +
  geom_point(size = 3) + theme_bw()
ggsave("../../out/plot/MDS_rld.pdf",width = 5,height = 4)

mdsPois <- as.data.frame(colData(ddsHTSeq_filter)) %>%
  cbind(cmdscale(samplePoisDistMatrix))

ggplot(mdsPois, aes(x = `1`, y = `2`, color = treat)) +
  geom_point(size = 3) + theme_bw()
ggsave("../../out/plot/PoissonDistance_scatter.pdf",width = 5,height = 4)

# save the object of interest ---------------------------------------------
saveRDS(vds_filter,file = "../../out/object/vds_all_filter.rds")
saveRDS(vds_gender,file = "../../out/object/vds_all.rds")
saveRDS(rld_filter,file = "../../out/object/rld_all_filter.rds")

# biotype of the dataset --------------------------------------------------
# read in the gtf file
gtf_test <- read.table(gzfile("../../data/gencode.v31.basic.annotation.gtf.gz"), header = FALSE, sep = '\t')
# gtf_test <- read.table("data/GCF_000001405.40_GRCh38.p14_genomic.gtf", header = FALSE, sep = '\t')
head(gtf_test)

# from V9 pull the gene_id and the filter only at the gene level
test <- gtf_test %>%
  dplyr::filter(V3 == "gene")

dim(test)
head(test)
# from the colum 9 pull the biotype and gene name
# test2 <- test %>%
#   separate(V9,into = c("gene_id","transcript_id","db_xref", "db_xref2", "description", "gbkey", "gene", "gene_biotype","pseudo"),sep = ";")

gene_ids <- str_extract(test$V9,pattern = "gene_id .*; gene_type") %>%
  str_remove_all("gene_id |; gene_type")

sum(is.na(gene_ids))
head(gene_ids)

gene_name <- str_extract(test$V9,pattern = "gene_name .*; level") %>%
  str_remove_all("gene_name |; level")

sum(is.na(gene_name))
head(gene_name)

biotype <- str_extract(test$V9,pattern = "gene_type \\w+;") %>%
  str_remove_all("gene_type |;")
sum(is.na(biotype))
head(biotype)


LUT_gtf_file <- data.frame(gene_ids,gene_name,biotype)
sum(is.na(LUT_gtf_file))
write_tsv(LUT_gtf_file,"../../out/table/df_LUT_gtf_genes.tsv")

LUT_gtf_file <- read_tsv("../../out/table/df_LUT_gtf_genes.tsv") %>% 
  group_by(gene_name,biotype) %>% 
  summarise() %>% 
  ungroup()

# confirm there are no NAs
sum(is.na(LUT_gtf_file))
# some gene names have multiple annotation
LUT_gtf_file %>% 
  group_by(gene_name) %>% 
  mutate(n = n()) %>% 
  arrange(desc(n))

# add the annotation to the dataset
df_biotype <- data.frame(gene_name = rownames(ddsHTSeq_filter)) %>% 
  left_join(LUT_gtf_file)

df_biotype %>%
  group_by(biotype) %>% 
  summarise(n=n()) %>% 
  mutate(prop = n/sum(n)) %>% 
  # dplyr::filter(!is.na(GENEBIOTYPE)) %>% 
  mutate(biotype = fct_reorder(biotype,-prop)) %>% 
  ggplot(aes(x=biotype,y=prop))+geom_col(position = "dodge")+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90),strip.background = element_blank())
# facet_wrap(~condition)
ggsave("../../out/plot/barplot_biotype.pdf",width = 6,height = 4)
