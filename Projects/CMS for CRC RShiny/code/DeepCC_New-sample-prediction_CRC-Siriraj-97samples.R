#### CMS CLASS FOR SINGLE-SAMPLE PREDICTION: CRC SIRIRAJ COHORT == REFERENCE #### 
# AUTHOR: JANTAPPAPA CHANTHERCROB 
# CONTACT: J.CHANTHERCROB@GMAIL.COM

#### CALL LIBRARIES ####
library(magrittr)
library(dplyr)
library(stringr)
library(tibble)

library(reticulate)
use_condaenv('r-reticulate')
library(tensorflow)
library(keras)
library(DeepCC)

library(umap)
library(ggplot2)
library(ggrepel)

#### PREPARATION GENE EXPRESSION DATA FOR CMS CLASSIFICATION ####
# LOAD COLORECTAL CANCER GENE EXPRESSION OF SIRIRAJ HOSPITAL
# CUT THE UNWANTED SAMPLME AND COLUMNS OUT
# REMOVE NA OF ENTREZ IDS
si.dat <- read.csv("data/GEP-CRC-SIRIRAJ_RSEM-TPM_100samples.csv ") %>%
  dplyr::select(-c(2, 3)) %>%
  dplyr::filter(entrez_ID != is.na(entrez_ID))

# FIND THE UNIQUE DUPLICATED GENE SYMBOLS
dup.entrez <- unique(si.dat$entrez_ID[duplicated(si.dat$entrez_ID)])

# DELETE DUPLICATED ENTREZ IDS
# REMOVE THE OLD ROW NAMES AND SET THE NEW ROW NAMES INTO THE DATA FRAME
si.dat <- si.dat[-which(si.dat$entrez_ID %in% dup.entrez), ] %>%
  remove_rownames() %>%
  tibble::column_to_rownames(., var = "entrez_ID") 
# DO LOG2-TRANSFORMATION
si.dat <- log2(si.dat + 1)

 #### PREPARE TESTING DATA: 3 SAMPLES ####
# DEFINE THE NEW SAMPLES FOR MAPPING IN THE UMAP PLANE
test.samp <- c("TBi_RNA02", "TBi_RNA46", "TBi_RNA32")
# 3 SAMPLES
# TRANSPOSE DATA IN FORM OF PATIENT SAMPLES AS ROWS AND GENES AS COLUMNS
# CONVERT MATRIX TYPE INTO DATA FRAME TYPE
si.test <- si.dat %>% 
  dplyr::select(all_of(test.samp)) %>% 
  t() %>%
  as.data.frame()

# LOAD COLORECTAL CANCER GENE EXPRESSION OF SIRIRAJ HOSPITAL: 97 SAMPLES
# TRANSPOSE DATA IN FORM OF PATIENT SAMPLES AS ROWS AND GENES AS COLUMNS
# CONVERT MATRIX TYPE INTO DATA FRAME TYPE
si.ref <- si.dat %>% 
  dplyr::select(-all_of(test.samp)) %>% 
  t() %>%
  as.data.frame()

#### GET FUNCTIONAL SPECTRA OF THE REFERENCE DATA SET: CRC SIRIRAJ DATA SET ####
# RUN FUNCTIONAL SPECTRA FUNCTION FOR ALL SAMPLES
# CHECK NUMBER OF CORES FOR THIS COMPUTER: parallel::detectCores()
#si.ref_fs <- getFunctionalSpectra(si.ref, geneSets = "MSigDBv7")
# SAVE FUNCTIONAL SPECTRA RESULT AS AN R DATA
#saveRDS(si.ref_fs, file = "FS-CRC-SIRIRAJ_RSEM-TPM_97samples.rds")

# LOAD FUNCTIONAL SPECTRA RESULTS OF SIRIRAJ COHORT
si.ref_fs <- readRDS(file = "data/FS-CRC-SIRIRAJ_RSEM-TPM_97samples.rds")

#### GET FUNCTIONAL SPECTRA OF THE TESTING DATA SET: 3 SAMPLES ####
#si.test_fs <- getFunctionalSpectra(si.test, , geneSets = "MSigDBv7", cores = 3)
# SAVE FUNCTIONAL SPECTRA RESULT AS AN R DATA
#saveRDS(si.test_fs, file = "FS-CRC-SIRIRAJ_RSEM-TPM_3samples.rds")

# LOAD FUNCTIONAL SPECTRA RESULTS OF SIRIRAJ COHORT
si.test_fs <- readRDS(file = "data/FS-CRC-SIRIRAJ_RSEM-TPM_3samples.rds")

#combine both dataframe for 100 patients
si.all_fs <- rbind(si.ref_fs, si.test_fs)

#save the dataframe 
#write.csv(x = si.all_fs, file = "data/Si.CMSfeat_100.csv", row.names = F)

#### CMS LABEL ####
# CALL THE CRC-TCGA CLASSIFIER
model.name <- "data/CRC_TCGA"
cms.model <- load_DeepCC_model(model.name)
# GET CMS LABELS FOR ALL SAMPLES
si.ref_lab <- get_DeepCC_label(cms.model, si.ref_fs, 
                               cutoff = 0.5, prob_mode = TRUE)
si.ref_na <- which(is.na(si.ref_lab$DeepCC))

#### GET 10 DEEP FEATURES OF ALL SAMPLES IN THE TESTING DATA SET: CRC SIRIRAJ DATA SET ####
# GET 10 FEATURES FROM THE DEEPCC CLASIFIER
si.ref_df <- get_DeepCC_features(cms.model, si.ref_fs)

#get deepfeature of all Si cohort
si.feature_df <- get_DeepCC_features(cms.model, si.all_fs)
si.label <- as.data.frame(get_DeepCC_label(cms.model, si.all_fs, cutoff = 0.5))

#umap visualization code
si_umap <- umap(d = si.feature_df,config = custom.config)
si_umap <- si_umap$layout

#combine the dataframe
si.all_df <- cbind(si_umap, si.label)

si.all_df <- which(is.na(si.all_df$`get_DeepCC_label(cms.model, si.all_fs, cutoff = 0.5)`))
# GET CMS LABELS FOR ALL SAMPLES
si.ref_lab_all <- get_DeepCC_label(cms.model, si.ref_fs, 
                               cutoff = 0.5, prob_mode = TRUE)
si.ref_na <- which(is.na(si.ref_lab$DeepCC))

#### GET 10 DEEP FEATURES OF ALL SAMPLES IN THE REFERENCE DATA SET: 3 SAMPLES ####
# GET 10 FEATURES FROM THE DEEPCC CLASIFIER
si.test_df <- get_DeepCC_features(cms.model, si.test_fs)

#### CREATE A UMAP TRANFORMATION OF THE REFERENCE DATA: CRC TCGA DATA SET ####
# SET SEED FOR REPRODUCIBLE
set.seed(1234)
# DEFINE CONFIGURATION SETTINGS OF UMAP
custom.config <- umap.defaults
custom.config$n_components <- 3
custom.config$min_dist <- 0.4
custom.config$n_neighbors <- 30
# PERFORM UMAP
si.ref_umap <- umap(d = si.ref_df, config = custom.config)
# NEW SAMPLE PREDICTION
si.test_pred <- data.frame(predict(si.ref_umap, si.test_df))

#perform UMAP on the dataframe
si.reference_umap <- umap(d = si.all_fs, config = custom.config)


#### VISUALIZATION ####
# PREPARE DATA FRAMES FOR VIZ: REFERENCE DATA, N = 97
si.ref_umap.dat <- data.frame("umap1" = si.ref_umap$layout[-si.ref_na, 1], 
                              "umap2" = si.ref_umap$layout[-si.ref_na, 2],
                              "umap3" = si.ref_umap$layout[-si.ref_na, 3],
                              "cms.lab" = si.ref_lab$DeepCC[-si.ref_na])
# PREPARE DATA FRAMES FOR VIZ: NEW DATA, N = 3
si.test_pred.dat <- data.frame("umap1" = si.test_pred$X1, 
                               "umap2" = si.test_pred$X2, 
                               "umap3" = si.test_pred$X3)
rownames(si.test_pred.dat) <- rownames(si.test_pred)
# DEFINE CMS SUBYPES COLORS
cms.col <- c("#FFA9A9", "#D7BEFF", "#9FE2BF", "#FFE493")
# UMAP 2 VS UMAP 3 DIMENSION
ggplot(data = NULL) + 
  stat_ellipse(data = si.ref_umap.dat, 
               aes(x = umap2, y = umap3, color = cms.lab), 
               type = "norm", linetype = 2, level = 0.95) +
  stat_ellipse(data = si.ref_umap.dat, 
               aes(x = umap2, y = umap3, color = cms.lab), 
               type = "norm", linetype = 2, level = 0.91) +
  geom_point(data = si.ref_umap.dat, 
             aes(x = umap2, y = umap3, color = cms.lab), 
             size = 3, alpha = 0.85, shape = 16) +
  geom_point(data = si.test_pred.dat, 
             aes(x = umap2, y = umap3), 
             color = "black", size = 4, shape = 21, stroke = 1) +
  geom_text_repel(data = si.test_pred.dat, 
                  aes(x = umap2, y = umap3, label = rownames(si.test_pred.dat)), 
                  color = "black", size = 6, box.padding = 1, 
                  family = "Arial") +
  scale_color_manual(values = cms.col) +
  labs(x = "UMAP2", y = "UMAP3", color = "CMS class", shape = "CMS class") + 
  #xlim(-6.5, 8) + ylim(-4, 6) +
  theme(#legend.position = c(1, 0.05),
        #legend.justification = c(1, 0.05), 
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

