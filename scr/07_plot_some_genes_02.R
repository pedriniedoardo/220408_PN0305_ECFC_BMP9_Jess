# AIM ---------------------------------------------------------------------
# plot some genes

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
df_DE <- read_tsv("../../out/table/res_BMP9_vs_Mock_shr.txt")

# the leading edges of the goveia signature
vect <- "KIT|CXCR4|N4BP3|INHBB|ARHGAP25|APLN|PGF|PALD1|RGCC|SMAD1|MYO1B|HTRA3|DAPK2|LXN|PCDH12|GNG2|SLC12A2|DLL4|RASSF2|SLC16A6|SH3TC1|ANGPT2|VASH1|SEMA6B|MCAM|PLEKHG1|TNFAIP8L1|DYSF|RAPGEF5|ADAM19|PXDN|USP31|ROBO4|FLT4|LTB|PDGFB|COL4A1|BCL6B|COL4A2|MMP15|NUDT4|LAMA4|PRR5|CTNNBIP1|ARHGAP27|EIF2B2|PLK2|ADAMTS7|GRK5|PRKD2|MARCKS|TCF4|FSCN1|CHST7|MYH10|DGKH|LAMC1|KLF6|UBTD1|TRAF3IP2|DLG1|OSBPL10|COLGALT1|MAP4K4|PPP1R13B|TANK|GRPEL2|ARL4C|ADM|NID1|SNN|FHL3|ST8SIA4|RHBDF1|RALGDS|EGLN1|COL18A1|NOTCH4|OTUD4|WIPF1|GALNT1|GALNT7|PRNP|DOCK6|VIM|GNAI2|SLC4A7|SVIL|RALA|ROCK2|GMFB|GPRIN3|CREB3L2|REV3L|LIMD1|STRN4|GRAP|EPHA2|FAM43A|INKA1|PLXNA2|MYH9|PPP1R18|TNPO1|ITPRIPL2|SIPA1|IGBP1|STX6"

vect2 <- str_split(vect,pattern = "\\|") %>%
  unlist()

GOI <- df_DE %>%
  dplyr::filter(symbol %in% vect2) %>%
  dplyr::filter(abs(log2FoldChange)>1,padj<0.05) %>%
  pull(symbol)

df_DE %>%
  dplyr::filter(symbol %in% c("KDR"))

# plot genes with heatmap -------------------------------------------------
vds_filter <- readRDS(file = "../../out/object/vds_all_filter.rds")

# gene_id <- rownames(assay(vds_filter)) %in% c("XIST","DDX3Y","RPS4Y1","USP9Y")
# gene_id <- rownames(assay(vds_filter)) %in% c("VCAN","CXCL12","ID2","ABI3BP","APLP1","MGP","TIMP3","SPOCK1","THY1","NID2","FBLN5","VEGFA")
# gene_id <- rownames(assay(vds_filter)) %in% c("SELP","C7","CPE")
# gene_id <- rownames(assay(vds_filter)) %in% c("GJA5","SEMA3G","DKK2","FBLN5")
gene_id <- rownames(assay(vds_filter)) %in% GOI

# notice that for both dataset the design didn't affect mucht he normalizationo,therefore the topmost variable genes are the same for both
# I will prduce the plot for just one dataset
# The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene’s average across all samples. Hence, we center each genes’ values across samples, and plot a heatmap (figure below). We provide a data.frame that instructs the pheatmap function how to label the columns.
# mat <- assay(vds_filter)[topVarGenes, ]
mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

# change the rownames to match the symbol 
# rownames(mat2) <- rownames(mat2) %>% 
#   data.frame(ensembl=.) %>% 
#   left_join(DEG_1,by = "ensembl") %>% 
#   pull(symbol) 

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

ht2 <- Heatmap(mat2, 
               name = "exp",
               top_annotation = column_ha, 
               # cluster_rows = F, 
               # col = colorRamp2(c(-2, 0, 1), c("green", "white", "red")),
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7)),
               column_title = "GOI") 

pdf("../../out/plot/heatmap_GOI_02.pdf",width = 6,height = 6) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# plot genes dotplot ------------------------------------------------------
data <- readRDS("../../out/object/dds_all_filter_DESeq.rds")

lut <- data@colData %>%
  data.frame()

# GOI <- c("SMAD1","SMAD9","SMAD6","BMPR1B","SMAD7")
# GOI <- subset_genes
# GOI <- c("VCAN","CXCL12","ID2","ABI3BP","APLP1","MGP","TIMP3","SPOCK1","THY1","NID2","FBLN5","VEGFA")
# GOI <- c("SELP","C7","CPE")
# GOI <- c("GJA5","SEMA3G","DKK2","FBLN5")

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
ggsave("../../out/plot/boxplot_GOI_02.pdf",width = 10,height = 8) 

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
ggsave("../../out/plot/scatterplot_GOI_02.pdf",width = 10,height = 8) 
