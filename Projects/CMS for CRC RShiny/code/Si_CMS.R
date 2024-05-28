rm(list=ls())
library(DeepCC)
library(keras)

Si_cohort2 <- read.csv("data/Si_cohort.csv")
Si_cohort <- read.csv("data/GEP-CRC-SIRIRAJ_RSEM-TPM_100samples.csv")
#Si_cohort <- Si_cohort[,-c(2,3)]
Si_cohort <- Si_cohort[,-c(1,3)]
Si_cohort <- t(Si_cohort)
#rename col names to 
colnames(Si_cohort) <- Si_cohort[1,]
Si_cohort <- Si_cohort[-c(1),]

load_DeepCC_model <- function(prefix){
  load(file = paste0(prefix, ".RData"))
  classifer <- keras::load_model_hdf5(filepath =paste0(prefix, ".hdf5"))
  list(classifier = classifer, levels = levels)
}
#change file path later
# for CRC model trained with Sample_ID data
CRC_TCGA <- load_DeepCC_model("data/CRC_TCGA")

#transform the data for DeepCC algorithm
#multiply the value by 10^6 and log2-transformed the data 
Si_cohort_TPM <- Si_cohort * (10^6)
#add 1/0.25 to avoid result of log2 of 0
Si_cohort_TPM <- Si_cohort_TPM + 1
Si_cohort_TPM <- log2(Si_cohort_TPM)
#swap column and row for DeepCC data
Si_cohort <- t(Si_cohort)
#get functionalspectra
Si_cohort_TPM_fs <- getFunctionalSpectra(Si_cohort)

# obtain deep features 
deepfeature <- get_DeepCC_features(CRC_TCGA, Si_cohort_TPM_fs)

#get DeepCC label
deeplabel <- get_DeepCC_label(CRC_TCGA, Si_cohort_TPM_fs, cutoff = 0.5, 
                              prob_mode = T, prob_raw = T)
deepprob <- get_DeepCC_prob(CRC_TCGA, Si_cohort_TPM_fs)

deepprob <- round(deepprob,2)

DeepCC_label <- cbind(DeepCC_label, deepprob)

rownames(DeepCC_label) <- rownames(Si_cohort)
# Assuming df1 and df2 are your dataframes
# Extract row names of df1 as a column
patient_ID <- data.frame(rowname = rownames(Si_cohort_fs), stringsAsFactors = FALSE)
print(patient_ID)

# Merge the rownames column with df2
DeepCC_label_TPM <- cbind(patient_ID, DeepCC_label_TPM)

DeepCC_label_TPM$rowname <- gsub("RNA","",DeepCC_label_TPM$rowname)
rownames(DeepCC_label_TPM) <- DeepCC_label_TPM$rowname

#############################################################
Si_cohort_clinical <- read.csv(file="data/110si_crc.csv")

#remove tumor id w/o name
Si_cohort_clinical <- Si_cohort_clinical[8:108,]

tbl_summary(
  Si_cohort_clinical,
  include = c(sex, age, metastasis, stage),
  by = cms_type,
  missing = "no",
) %>%
  add_n() %>%
  add_p() %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()

#rename the Si_cohort_clinical dataframe to match with Si_cohort 
Si_cohort_clinical$tumor_tbi <- gsub("-","_",Si_cohort_clinical$tumor_tbi)
Si_cohort_clinical$tumor_tbi <- gsub("DNA", "", Si_cohort_clinical$tumor_tbi)
Si_cohort_clinical$tumor_tbi <- gsub("I", "i", Si_cohort_clinical$tumor_tbi)

Si_cohort_clinical <- head(Si_cohort_clinical, -5)

rownames(Si_cohort_clinical) <- Si_cohort_clinical$tumor_tbi
Si_cohort_clinical <- Si_cohort_clinical[,-1]

#rename the Si_cohort df
rownames(Si_cohort) <- gsub("RNA", "", rownames(Si_cohort))

#combine the dataframes
combined_Si <- cbind(Si_cohort_clinical, DeepCC_label_TPM)

