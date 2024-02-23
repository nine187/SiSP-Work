rm(list=ls())
graphics.off()

library(TCGAbiolinks)

#import the matched csv
matched_tcga <- read.csv(file="data/matched_TCGA.csv")

#load the TCGA COAD and READ dataset
TCGA_COAD <- GDCquery_clinic("TCGA-COAD")
TCGA_READ <- GDCquery_clinic("TCGA-READ")

#look for the matched columns in COAD and READ
which(colnames(TCGA_COAD) %in% c("submitter_id", "vital_status", "days_to_last_follow_up", "days_to_death", "site_of_resection_or_biopsy"))
TCGA_COAD <- TCGA_COAD[,c(2,9,27,39,45)]

which(colnames(TCGA_READ) %in% c("submitter_id", "vital_status", "days_to_last_follow_up", "days_to_death", "site_of_resection_or_biopsy"))
TCGA_READ <- TCGA_READ[,c(2,9,25,107,113)]

#merge the dataframe
TCGA_COADREAD <- rbind(TCGA_COAD,TCGA_READ)

#convert to binary 

TCGA_COADREAD$deceased <- ifelse(TCGA_COADREAD$vital_status =="Alive", FALSE,TRUE)

TCGA_COADREAD$deceased <- as.character(lapply(TCGA_COADREAD$deceased, function(x) ifelse(x,0,1)))

TCGA_COADREAD$overall_survival <- ifelse(TCGA_COADREAD$vital_status == "Alive",
                                     TCGA_COADREAD$days_to_last_follow_up,
                                     TCGA_COADREAD$days_to_death)

TCGA_COADREAD <- TCGA_COADREAD[, !colnames(TCGA_COADREAD) %in% "vital_status"]
TCGA_COADREAD <- TCGA_COADREAD[, !colnames(TCGA_COADREAD) %in% "days_to_death"]
TCGA_COADREAD <- TCGA_COADREAD[, !colnames(TCGA_COADREAD) %in% "days_to_last_follow_up"]

write.csv(TCGA_COADREAD, file="data/test_cohort2.csv", row.names = FALSE)
