# libraries ---------------------------------------------------------------
source(file = "scr/00_library.R")

# read in the data --------------------------------------------------------
ddsHTSeq_gender <- readRDS("../../out/object/ddsHTSeq_BMP9_gender.rds")

# PRE-FILTER DATASET ------------------------------------------------------
nrow(ddsHTSeq_gender)
# remove potential non infirmative genes
ddsHTSeq2_filter <- ddsHTSeq_gender[rowSums(counts(ddsHTSeq_gender)) > 10, ]

nrow(ddsHTSeq2_filter)

# TRANSFORM DATASET -------------------------------------------------------
#vsd2
vsd2_filter <- vst(ddsHTSeq2_filter, blind = F)
vsd2_blind <- vst(ddsHTSeq2_filter, blind = T)
head(assay(vsd2_filter), 3)
head(assay(vsd2_blind), 3)

# rlog
rld2_filter <- rlog(ddsHTSeq2_filter, blind = F)
rld2_blind <- rlog(ddsHTSeq2_filter, blind = T)
head(assay(rld2_filter), 3)
head(assay(rld2_blind), 3)

# save the scatter plot of the variance vs mean expression
meanSdPlot_vsd2 <- meanSdPlot(assay(vsd2_filter))
ggsave(plot = meanSdPlot_vsd2$gg+theme_bw(),filename = "../../out/image/meanSdPlot_vsd_gender.pdf",width = 4,height = 4)

meanSdPlot_rlog <- meanSdPlot(assay(rld2_filter))
ggsave(plot = meanSdPlot_rlog$gg+theme_bw(),filename = "../../out/image/meanSdPlot_rlog_gender.pdf",width = 4,height = 4)

# In the above function calls, we specified blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
# To show the effect of the transformation, in the figure below we plot the first sample against the second, first simply using the log2 function (after adding 1, to avoid taking the log of zero), and then using the VST and rlog-transformed values. For the log2 approach, we need to first estimate size factors to account for sequencing depth, and then specify normalized=TRUE. Sequencing depth correction is done automatically for the vst and rlog.
ddsHTSeq2_filter <- estimateSizeFactors(ddsHTSeq2_filter)

# SAMPLE DISTANCE ---------------------------------------------------------
# A useful first step in an RNA-seq analysis is often to assess overall similarity between samples: Which samples are similar to each other, which are different? Does this fit to the expectation from the experiment’s design?
# We use the R function dist to calculate the Euclidean distance between samples. To ensure we have a roughly equal contribution from all genes, we use it on the VST data. We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.
sampleDists_vsd2 <- dist(t(assay(vsd2_filter)))
# sampleDists_blind_vsd2 <- dist(t(assay(vsd2_blind)))
sampleDists_vsd2
# sampleDists_blind_vsd2

head(assay(vsd2_filter))
# head(assay(vsd2_blind))

# We visualize the distances in a heatmap in a figure below, using the function pheatmap from the pheatmap package.
# library("pheatmap")
# library("RColorBrewer")
# In order to plot the sample distance matrix with the rows/columns arranged by the distances in our distance matrix, we manually provide sampleDists to the clustering_distance argument of the pheatmap function. Otherwise the pheatmap function would assume that the matrix contains the data values themselves, and would calculate distances between the rows/columns of the distance matrix, which is not desired. We also manually specify a blue color palette using the colorRampPalette function from the RColorBrewer package.
sampleDistMatrix_vsd2 <- as.matrix(sampleDists_vsd2)
# sampleDistMatrix_blind_vsd2 <- as.matrix(sampleDists_blind_vsd2)
# colData(vsd2_blind)

rownames(sampleDistMatrix_vsd2) <- paste(vsd2_filter$clone,vsd2_filter$gender,vsd2_filter$treat,sep = "_")
colnames(sampleDistMatrix_vsd2) <- NULL

map_colors<-colorRampPalette(viridis(12))(255)

hm_1 <- pheatmap(sampleDistMatrix_vsd2,
                 clustering_distance_rows = sampleDists_vsd2,
                 clustering_distance_cols = sampleDists_vsd2,
                 col = map_colors)
# pdf("images/heatmap_vsd2_tot.pdf",width = 4,height = 3)
pdf("../../out/image/heatmap_vsd_gender.pdf",width = 5,height = 3)
hm_1
dev.off()

# same report using rlog normalization
sampleDists_rld2 <- dist(t(assay(rld2_filter)))
# sampleDists_blind_rld2 <- dist(t(assay(rld2_blind)))
sampleDists_rld2
# sampleDists_blind_rld2

head(assay(rld2_filter))
head(assay(rld2_blind))

# We visualize the distances in a heatmap in a figure below, using the function pheatmap from the pheatmap package.
sampleDistMatrix_rld2 <- as.matrix(sampleDists_rld2)
# sampleDistMatrix_blind_rld2 <- as.matrix(sampleDists_blind_rld2)

rownames(sampleDistMatrix_rld2) <- paste(rld2_filter$clone,rld2_filter$gender,rld2_filter$treat,sep = "_")
colnames(sampleDistMatrix_rld2) <- NULL

# rownames(sampleDistMatrix_blind_rld2) <- paste(vsd2_blind$var)
# colnames(sampleDistMatrix_blind_rld2) <- NULL

hm_1_2 <- pheatmap(sampleDistMatrix_rld2,
                   clustering_distance_rows = sampleDists_rld2,
                   clustering_distance_cols = sampleDists_rld2,
                   col = map_colors)

#  pdf("images/heatmap_rld2_tot.pdf",width = 4,height = 3)
pdf("../../out/image/heatmap_rld_gender.pdf",width = 5,height = 3)
hm_1_2
dev.off()
# save_pheatmap_pdf(width = 4,height = 3,x = hm_1_2,filename = "images/heatmap_distance_diabetic_rld2.pdf")

# Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples. The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead of columns, so we need to transpose the counts in dds.
# since this thistance rely on the raw counts there is no difference between the structure and diabetic dataset in this case (the sample are the same)
#library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsHTSeq2_filter,normalized = F)))

samplePoisDistMatrix <- as.matrix(poisd$dd)

rownames(samplePoisDistMatrix) <- paste(ddsHTSeq2_filter$clone,ddsHTSeq2_filter$gender,ddsHTSeq2_filter$treat,sep="_")
colnames(samplePoisDistMatrix) <- NULL

hm_p <- pheatmap(samplePoisDistMatrix,
                 clustering_distance_rows = poisd$dd,
                 clustering_distance_cols = poisd$dd,
                 col = map_colors)
#
# pdf("images/heatmap_poisd_tot.pdf",width = 4,height = 3)
pdf("../../out/image/heatmap_poisd_gender.pdf",width = 5,height = 3)
hm_p
dev.off()
# save_pheatmap_pdf(width = 4,height = 3,x = hm_p,filename = "images/heatmap_distance_poisd.pdf")

# PCA plot ----------------------------------------------------------------
# Another way to visualize sample-to-sample distances is a principal components analysis (PCA). In this ordination method, the data points (here, the samples) are projected onto the 2D plane such that they spread out in the two directions that explain most of the differences. The x-axis is the direction that separates the data points the most. The values of the samples in this direction are written PC1. The y-axis is a direction (it must be orthogonal to the first direction) that separates the data the second most. The values of the samples in this direction are written PC2. The percent of the total variance that is contained in the direction is printed in the axis label. Note that these percentages do not add to 100%, because there are more dimensions that contain the remaining variance (although each of these remaining dimensions will explain less than the two that we see).
plot_vsd2 <- plotPCA(vsd2_filter, intgroup = c("clone","treat","gender")) + theme_bw()
# ggsave(plot = plot_vsd2,"image/PCA_vsd2_MG.pdf",width = 4,height = 4)

# plot_vsd2$data %>%
#   # mutate(senescence = factor(senescence,levels = c("EP","LP","ETO"))) %>%
#   ggplot(aes(x=PC1,y=PC2,col=treat,label=clone)) +
#   geom_text() +
#   scale_x_continuous(expand = expansion(mult = 0.1))+
#   # geom_point(alpha=0.8,size=3) +
#   theme_bw() + ylab(plot_vsd2$labels[1]) + xlab(plot_vsd2$labels[2])
#   # scale_color_manual(values = c(viridis(3)))
#   # purple
#   # scale_color_manual(values = c(c("#ffccff","#cc00cc","#4d004d")))

# scale_color_manual(values = hue_pal()(4)[c(2,1,3,4)])
plot_vsd2$data %>%
  # mutate(senescence = factor(senescence,levels = c("EP","LP","ETO"))) %>%
  ggplot(aes(x=PC1,y=PC2,col=gender,shape=treat,label=clone)) +
  geom_point() +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  # geom_point(alpha=0.8,size=3) +
  theme_bw() + ylab(plot_vsd2$labels[1]) + xlab(plot_vsd2$labels[2])
ggsave("../../out/image/PCA_vsd_gender.pdf",width = 5,height = 4)

# plot_vsd2$data %>%
#   left_join(meta,c("name"="sample_name","treat","clone")) %>%
#   # mutate(senescence = factor(senescence,levels = c("EP","LP","ETO"))) %>%
#   ggplot(aes(x=PC1,y=PC2,col=gender,label=clone)) +
#   geom_point() +
#   ggrepel::geom_text_repel(show.legend = F)+
#   scale_x_continuous(expand = expansion(mult = 0.1))+
#   # geom_point(alpha=0.8,size=3) +
#   theme_bw() + ylab(plot_vsd2$labels[1]) + xlab(plot_vsd2$labels[2])

# plot_vsd2$data %>%
#   ggplot(aes(x=PC1,y=PC2,col=class)) + geom_point(alpha=0.8,size=3) + theme_bw() + ylab(plot_vsd2$labels[1]) + xlab(plot_vsd2$labels[2]) + facet_wrap(~stim)
# # scale_color_manual(values = hue_pal()(4)[c(2,1,3,4)])
# ggsave("images/PCA_vsd2_split.pdf",width = 6,height = 4)


plot_rld2 <- plotPCA(rld2_filter, intgroup = c("clone","treat","gender")) + theme_bw()

plot_rld2$data %>%
  ggplot(aes(x=PC1,y=PC2,col=gender,shape=treat,label=clone)) +
  # geom_point(alpha=0.8,size=3) +
  geom_point() +
  ggrepel::geom_text_repel(show.legend = F)+
  theme_bw() + ylab(plot_rld2$labels[1]) + xlab(plot_rld2$labels[2])
# scale_color_manual(values = hue_pal()(4)[c(2,1,3,4)])
ggsave("../../out/image/PCA_rld_gender.pdf",width = 5,height = 4)

# From the PCA plot, we see that the differences between cells (the different plotting shapes) are considerable, though not stronger than the differences due to treatment with dexamethasone (red vs blue color). This shows why it will be important to account for this in differential testing by using a paired design (“paired”, because each dex treated sample is paired with one untreated sample from the same cell line). We are already set up for this design by assigning the formula ~ cell + dex earlier.

# MSD plot ----------------------------------------------------------------
# Another plot, very similar to the PCA plot, can be made using the multidimensional scaling (MDS) function in base R. 
# This is useful when we don’t have a matrix of data, but only a matrix of distances. Here we compute the MDS for the distances calculated from the VST data and plot these in a figure below.
mds_vsd2 <- as.data.frame(colData(vsd2_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_vsd2))

ggplot(mds_vsd2, aes(x = `1`, y = `2`, color = gender)) +
  geom_point(size = 3)+ theme_bw()
# + coord_fixed() 
ggsave("../../out/image/MDS_vsd_gender.pdf",width = 4,height = 4)

mds_rld2 <- as.data.frame(colData(rld2_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_rld2))

ggplot(mds_rld2, aes(x = `1`, y = `2`, color = gender)) +
  geom_point(size = 3) + theme_bw()
# + coord_fixed()
ggsave("../../out/image/MDS_rld_gender.pdf",width = 4,height = 4)

# In a figure below we show the same plot for the PoissonDistance:
# as before the counts are the same for both the datasets (structure, diabetic)
mdsPois <- as.data.frame(colData(ddsHTSeq2_filter)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
# mutate(senescence = factor(senescence,levels = c("EP","LP","ETO")))

ggplot(mdsPois, aes(x = `1`, y = `2`, color = gender)) +
  geom_point(size = 3) + theme_bw()
# + coord_fixed()
# scale_color_manual(values = c(viridis(3)))
# purple
# scale_color_manual(values = c(c("#ffccff","#cc00cc","#4d004d")))

# ggsave("images/PoissonDistance_scatter_senescence_company.pdf",width = 4,height = 4)
# ggsave("images/PoissonDistance_scatter_senescence_company_viridis.pdf",width = 4,height = 4)
ggsave("../../out/image/PoissonDistance_scatter_gender.pdf",width = 4,height = 4)

# save the object of interest ---------------------------------------------
saveRDS(vsd2_filter,file = "../../out/object/vsd_filter_gender.rds")
saveRDS(vsd2_blind,file = "../../out/object/vsd_blind_gender.rds")
saveRDS(rld2_filter,file = "../../out/object/rld_filter_gender.rds")
saveRDS(rld2_blind,file = "../../out/object/rld_blind_gender.rds")
