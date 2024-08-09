rm(list=ls())

#load the CMS results data
library(readxl)
library(dplyr)

#Si cohort gene expression data
Si_cohort_gene_xpr <- read.csv(file = "data/Si_cohort.csv")

#read si cohort result data, 100 samples
Si_cohort_res <- read.csv(file = "data/Si_cohort_res.csv")

#read si cohort from P'O, 110 samples
Si_cohort_dat <- read_excel(path ="data/110 Si_clinical.xlsx")

new_col <- c("tumor_tbi", "SiSP_CMS", "SiSP_CMS1", "SiSP_CMS2", "SiSP_CMS3", "SiSP_CMS4", "SiSP_CRIS", "SiSP_CRISdist")

colnames(Si_cohort_res) <- new_col

#rename the column to be the same as Si_cohort_dat
Si_cohort_res$tumor_tbi <- gsub("_RNA", "",x = Si_cohort_res$tumor_tbi)
Si_cohort_dat$tumor_tbi <- gsub("-", "", x = Si_cohort_dat$tumor_tbi)
Si_cohort_dat$tumor_tbi <- gsub("I", "i", x= Si_cohort_dat$tumor_tbi)
Si_cohort_dat$tumor_tbi <- gsub("DNA", "", x = Si_cohort_dat$tumor_tbi)

#extract the ID for both dataframe
Si_RNAseq <- Si_cohort_dat$tumor_tbi
Si_clinical <- Si_cohort_res$tumor_tbi

nonoverlap_sample <- Si_RNAseq[!Si_RNAseq %in% Si_clinical]
print(nonoverlap_sample)

Si_cohort_combined <- merge(Si_cohort_dat, Si_cohort_res, by = "tumor_tbi")

#Si_cohort_combined$cris_type <- mapping_

writexl::write_xlsx(Si_cohort_combined, path = "results/Si_cohort.xlsx")

