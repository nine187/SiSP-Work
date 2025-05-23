#### CMS clinical data table #### 
# AUTHOR: PASITH PRAYOONRAT
# CONTACT: PASITH.P@GMAIL.COM


library(gtsummary)

Si_cohort_clinical <- read.csv(file="data/110si_crc.csv")
Si_cohort_RNA <- read.csv(file="data/Si_cohort.csv")

#remove tumor id w/o name
Si_cohort_clinical <- Si_cohort_clinical[8:108,]

#extract first 3 samples in the cohort from both clinical and RNA-seq data
Si_cohort_clinical_3 <- Si_cohort_clinical[c(1:3),]
Si_cohort_RNA_3 <- Si_cohort_RNA[,c(1:4)]


tbl_summary(
  Si_cohort_clinical,
  include = c(sex, age, metastasis, stage, location),
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
rownames(Si_cohort_clinical) <- gsub("RNA", "", rownames(Si_cohort_clinical))
