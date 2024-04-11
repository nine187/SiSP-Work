#### CMS CLASS FOR SINGLE-SAMPLE PREDICTION: CRC TCGA COHORT == REFERENCE #### 
# AUTHOR: JANTAPPAPA CHANTHERCROB 
# CONTACT: J.CHANTHERCROB@GMAIL.COM

#### CALL LIBRARIES ####
library(magrittr)
library(dplyr)
library(stringr)

library(reticulate)
use_condaenv('r-reticulate')
library(tensorflow)
library(keras)
library(DeepCC)

library(umap)
library(ggplot2)
library(ggrepel)

#### GENE EXPRESSION DATA FOR CMS CLASSIFICATION ####
# LOAD TCGA DATA OF COLORECTAL CANCER THAT DEEPCC USED FOR TRAINING THEIR MODEL
# TRANSPOSE DATA IN FORM OF PATIENT SAMPLES AS ROWS AND GENES AS COLUMNS
# CONVERT MATRIX TYPE INTO DATA FRAME TYPE
tcga.ref <- read.csv("data/GEP-CRC-TCGA_log2TPM_456samples.csv", 
                     row.names = 1, check.names = FALSE) %>% 
  t() %>%
  as.data.frame()

tcga.lab <- read.delim("data/cms_labels_public_all.txt") 
tcga.lab$sample <- paste(tcga.lab$sample, "-01", sep = "")
tcga.lab <- tcga.lab %>% 
  dplyr::filter(dataset == "tcga") %>%
  dplyr::select(sample, dataset, CMS_network) %>% 
  dplyr::filter(sample %in% rownames(tcga.ref)) 

# CHECK ID NAMES FOR BOTH DATA FRAMES
# MUST BE TRUE
all(tcga.lab$sample == row.names(tcga.ref))

#### GET FUNCTIONAL SPECTRA OF THE REFFERENT DATA SET: CRC TCGA DATA SET ####
# RUN FUNCTIONAL SPECTRA FUNCTION FOR ALL SAMPLES
# CHECK NUMBER OF CORES FOR THIS COMPUTER: parallel::detectCores()
#tcga.ref_fs <- getFunctionalSpectra(tcga.ref, geneSets = "MSigDBv7")
# SAVE FUNCTIONAL SPECTRA RESULT AS AN R DATA
#saveRDS(tcga.ref_fs, file = "FS-CRC-TCGA_log2TPM_456samples.rds")

# LOAD FUNCTIONAL SPECTRA RESULTS OF TCGA DATA
tcga.ref_fs <- readRDS(file = "data/FS-CRC-TCGA_log2TPM_456samples.rds")
ifelse(all(class(tcga.ref_fs) == c("matrix", "array")), 
       "Go to get thier deep features", 
       "Do not forget to convert FS data to the maxtrix type")

# LOAD FUNCTIONAL SPECTRA RESULTS OF SIRIRAJ COHORT
si.test_fs <- readRDS(file = "data/FS-CRC-SIRIRAJ_RSEM-TPM_3samples.rds")

#### GET 10 DEEP FEATURES OF ALL SAMPLES IN THE REFFERENT DATA SET: CRC TCGA DATA SET ####
# CALL THE CRC-TCGA CLASSIFIER
model.name <- "data/CRC_TCGA"
cms.model <- load_DeepCC_model(model.name)

# GET 10 FEATURES FROM THE DEEPCC CLASIFIER
tcga.ref_df <- get_DeepCC_features(cms.model, tcga.ref_fs)

si.test_df <- get_DeepCC_features(cms.model, si.test_fs)
si.test_df <- data.frame(si.test_df)

#### CREATE A UMAP TRANFORMATION OF THE REFERENCE DATA: CRC TCGA DATA SET ####
# PERFORM UMAP
set.seed(1234)
custom.config <- umap.defaults
custom.config$n_components <- 3
#custom.config$min_dist <- 0.5
#custom.config$n_neighbors <- 20
tcga.ref_umap <- umap(d = tcga.ref_df, 
                      config = custom.config)

si.test_umap <- data.frame(predict(tcga.ref_umap, si.test_df)) #modify P Arm's code

tcga.ref_umap.dat <- data.frame("umap1" = tcga.ref_umap$layout[, 1], 
                                "umap2" = tcga.ref_umap$layout[, 2],
                                "umap3" = tcga.ref_umap$layout[, 3],
                                "cms.lab" = tcga.lab$CMS_network)

si.test_umap.dat <- data.frame("umap1" = si.test_umap$X1, 
                               "umap2" = si.test_umap$X2,
                               "umap3" = si.test_umap$X3)
rownames(si.test_umap.dat) <- rownames(si.test_umap)

cms.col <- c("#FFA9A9", "#D7BEFF", "#9FE2BF", "#FFE493")
ggplot(data = NULL) + 
  stat_ellipse(data = tcga.ref_umap.dat, 
               aes(x = umap1, y = umap2, color = cms.lab), 
               type = "norm", linetype = 2, level = 0.95) +
  stat_ellipse(data = tcga.ref_umap.dat, 
               aes(x = umap1, y = umap2, color = cms.lab), 
               type = "norm", linetype = 2, level = 0.91) +
  geom_point(data = tcga.ref_umap.dat, 
             aes(x = umap1, y = umap2, color = cms.lab), 
             size = 3, alpha = 0.85, shape = 16) +
  geom_point(data = si.test_umap.dat, 
             aes(x = umap1, y = umap2), 
             color = "black", size = 4, shape = 21, stroke = 1) +
  geom_text_repel(data = si.test_umap.dat, 
                  aes(x = umap1, y = umap2, label = rownames(si.test_umap.dat)), 
                  color = "black", size = 7, box.padding = 1, 
                  family = "Arial") +
  scale_color_manual(values = cms.col) +
  labs(x = "UMAP1", y = "UMAP2", color = "CMS class", shape = "CMS class") + 
  #xlim(-6.5, 8) + ylim(-4, 6) +
  theme(legend.position = c(0.05, 0.05),
        legend.justification = c(0.05, 0.05), 
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 14), 
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 16), 
        axis.title.y = element_text(face = "bold", size = 16), 
        axis.line = element_line(size = 1), 
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14, face = "bold"),
        panel.background = element_rect(color = "white", fill = "white"), 
        strip.background = element_rect(color = "white", fill = "white"),
        panel.grid.major = element_line(color = "gray", linewidth = 0.25), 
        panel.grid.minor = element_line(color = "gray87", linewidth = 0.12))

