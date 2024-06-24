# AIM ---------------------------------------------------------------------
# the aim of the script is to read in the raw table of counts and build the object

# libraries ---------------------------------------------------------------
library("tidyverse")
library("DESeq2")

# read in the expression data ---------------------------------------------
# sample read one file
Cts <-"../../data/all.counts"
test <- read.delim(Cts,header = T)
head(test)

# metadata ----------------------------------------------------------------
# build the annotation besed on the sample metadata
LUT_samples <- read_csv("../../data/LUT_samples.csv")

# extract only the count information
mat_exp <- test %>% 
  dplyr::select(-c("Chr","Start","End","Strand","Length")) %>% 
  column_to_rownames("Geneid") %>% 
  # dplyr::select(contains("PN0495")) %>% 
  # keep all the samples of interest
  dplyr::select(LUT_samples$sample_id) %>% 
  as.matrix()

# match the order of the sample in the matrix with the sample in the sample sheet
# build the metadata for the deseq object
coldata <- data.frame(sample = colnames(mat_exp)) %>% 
  left_join(LUT_samples,by = c("sample"="sample_id")) %>% 
  mutate(rowname = sample) %>% 
  column_to_rownames("rowname")

# save the table
write.csv(coldata,file = "../../data/LUT_samples_final.csv",row.names = T)

# define the model --------------------------------------------------------
clone <- coldata$clone
gender <- coldata$gender
treat <- factor(coldata$treat,levels = c("mock","BMP9"))

# build the design
design <- model.matrix(~ clone + treat)
colnames(design) <- c("intercept","cloneclone_4A","cloneClone_51","cloneClone_70","cloneclone_9C","treatBMP9")

saveRDS(design,file = "../../out/object/design_all.rds")

# build the object --------------------------------------------------------
# is keeping only the objext in the lut_sample
dds <- DESeqDataSetFromMatrix(countData = mat_exp,
                              colData = coldata,
                              design = design)

saveRDS(dds,file = "../../out/object/dds_all.rds")
