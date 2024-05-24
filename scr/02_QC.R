# libraries ---------------------------------------------------------------
source(file = "scr/00_library.R")

# read in the data --------------------------------------------------------
ddsHTSeq <- readRDS("../../out/object/ddsHTSeq_BMP9.rds")

# remove low epressed genes -----------------------------------------------
colSums(counts(ddsHTSeq)) %>%
  data.frame(tot_counts=.) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y = tot_counts))+geom_col()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("../../out/image/barplot_tot_count.pdf",width = 5,height = 5)

nrow(ddsHTSeq)
# remove potential non infirmative genes
ddsHTSeq_filter <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 10, ]

nrow(ddsHTSeq_filter)

# scaling transformation of the data --------------------------------------
#vsd
vds_filter <- vst(ddsHTSeq_filter, blind = F)
vsd_blind <- vst(ddsHTSeq_filter, blind = T)
head(assay(vds_filter), 3)
head(assay(vsd_blind), 3)

# rlog
rld_filter <- rlog(ddsHTSeq_filter, blind = F)
rld_blind <- rlog(ddsHTSeq_filter, blind = T)
head(assay(rld_filter), 3)
head(assay(rld_blind), 3)

meanSdPlot_vsd <- meanSdPlot(assay(vds_filter))
ggsave(plot = meanSdPlot_vsd$gg+theme_bw(),filename = "../../out/image/meanSdPlot_vsd.pdf",width = 4,height = 4)

meanSdPlot_rlog <- meanSdPlot(assay(rld_filter))
ggsave(plot = meanSdPlot_rlog$gg+theme_bw(),filename = "../../out/image/meanSdPlot_rlog.pdf",width = 4,height = 4)

ddsHTSeq_filter <- estimateSizeFactors(ddsHTSeq_filter)

# sample distance ---------------------------------------------------------
sampleDists_vsd <- dist(t(assay(vds_filter)))
sampleDists_vsd

head(assay(vds_filter))

sampleDistMatrix_vsd <- as.matrix(sampleDists_vsd)

rownames(sampleDistMatrix_vsd) <- paste(vds_filter$clone,vds_filter$treat,sep = "_")
colnames(sampleDistMatrix_vsd) <- NULL

map_colors<-colorRampPalette(viridis(12))(255)

hm_1 <- pheatmap(sampleDistMatrix_vsd,
                 clustering_distance_rows = sampleDists_vsd,
                 clustering_distance_cols = sampleDists_vsd,
                 col = map_colors)

pdf("../../out/image/heatmap_vsd.pdf",width = 5,height = 3)
hm_1
dev.off()

sampleDists_rld <- dist(t(assay(rld_filter)))
sampleDists_rld

head(assay(rld_filter))
head(assay(rld_blind))

sampleDistMatrix_rld <- as.matrix(sampleDists_rld)

rownames(sampleDistMatrix_rld) <- paste(rld_filter$clone,rld_filter$treat,sep = "_")
colnames(sampleDistMatrix_rld) <- NULL

hm_1_2 <- pheatmap(sampleDistMatrix_rld,
                   clustering_distance_rows = sampleDists_rld,
                   clustering_distance_cols = sampleDists_rld,
                   col = map_colors)

pdf("../../out/image/heatmap_rld.pdf",width = 5,height = 3)
hm_1_2
dev.off()

poisd <- PoissonDistance(t(counts(ddsHTSeq_filter,normalized = F)))

samplePoisDistMatrix <- as.matrix(poisd$dd)

rownames(samplePoisDistMatrix) <- paste(ddsHTSeq_filter$clone,ddsHTSeq_filter$treat,sep="_")
colnames(samplePoisDistMatrix) <- NULL

hm_p <- pheatmap(samplePoisDistMatrix,
                 clustering_distance_rows = poisd$dd,
                 clustering_distance_cols = poisd$dd,
                 col = map_colors)

pdf("../../out/image/heatmap_poisd.pdf",width = 5,height = 3)
hm_p
dev.off()

# PCA plot ----------------------------------------------------------------
plot_vsd <- plotPCA(vds_filter, intgroup = c("clone","treat")) + theme_bw()

plot_vsd$data %>%
  ggplot(aes(x=PC1,y=PC2,col=treat,label=clone)) +
  geom_point() +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../../out/image/PCA_vsd.pdf",width = 5,height = 4)

plot_rld <- plotPCA(rld_filter, intgroup = c("clone","treat")) + theme_bw()

plot_rld$data %>%
  ggplot(aes(x=PC1,y=PC2,col=treat,label=clone)) +
  geom_point() +
  ggrepel::geom_text_repel(show.legend = F)+
  theme_bw() + ylab(plot_rld$labels[1]) + xlab(plot_rld$labels[2])
ggsave("../../out/image/PCA_rld.pdf",width = 5,height = 4)

# MSD plot ----------------------------------------------------------------
mds_vsd <- as.data.frame(colData(vds_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_vsd))

ggplot(mds_vsd, aes(x = `1`, y = `2`, color = treat)) +
  geom_point(size = 3)+ theme_bw()
ggsave("../../out/image/MDS_vsd.pdf",width = 5,height = 4)

mds_rld <- as.data.frame(colData(rld_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_rld))

ggplot(mds_rld, aes(x = `1`, y = `2`, color = treat)) +
  geom_point(size = 3) + theme_bw()
ggsave("../../out/image/MDS_rld.pdf",width = 5,height = 4)

mdsPois <- as.data.frame(colData(ddsHTSeq_filter)) %>%
  cbind(cmdscale(samplePoisDistMatrix))

ggplot(mdsPois, aes(x = `1`, y = `2`, color = treat)) +
  geom_point(size = 3) + theme_bw()
ggsave("../../out/image/PoissonDistance_scatter.pdf",width = 5,height = 4)

# save the object of interest ---------------------------------------------
saveRDS(vds_filter,file = "../../out/object/vds_filter.rds")
saveRDS(vsd_blind,file = "../../out/object/vsd_blind.rds")
saveRDS(rld_filter,file = "../../out/object/rld_filter.rds")
saveRDS(rld_blind,file = "../../out/object/rld_blind.rds")

# define the gender of the samples ----------------------------------------
gene_id <- rownames(assay(vds_filter)) %in% c("XIST","DDX3Y","RPS4Y1","USP9Y")

mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

anno_structure <- as.data.frame(colData(vds_filter)[, c("clone","treat")])

hm_var <- pheatmap(mat2, annotation_col = anno_structure)

pdf(file = "../../out/image/heatmap_gender_genes.pdf", width = 6, height = 3.5)
hm_var
dev.off()
