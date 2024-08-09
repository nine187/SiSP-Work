rm(list=ls())

#data normalization package 
library(limma)

#read the TPM files
CRC_78samples <- read.csv(file="data/Processed-RNAseq_CRCProject_20240704/Processed-RNAseq_CRCProject_20240704/02.ProcessedData/GEP_TPM/GEP_CRCproject78samp_TPM.csv")
CRC_78samples <- CRC_78samples[,-1]
CRC_10samples <- read.csv(file="data/Processed-RNAseq_CRCProject_20240704/Processed-RNAseq_CRCProject_20240704/02.ProcessedData/GEP_TPM/GEP_CRCproject10samp_TPM.csv")
CRC_10samples <- CRC_10samples[,-1]
CRC_90samples <- read.csv(file="data/Processed-RNAseq_CRCProject_20240704/Processed-RNAseq_CRCProject_20240704/02.ProcessedData/GEP_TPM/GEP_CRCproject90samp_TPM.csv")
CRC_90samples <- CRC_90samples[,-1]
CRC_30samples <- read.csv(file="data/Processed-RNAseq_CRCProject_20240704/Processed-RNAseq_CRCProject_20240704/02.ProcessedData/GEP_TPM/GEP_CRCproject30samp_TPM.csv")

#read the 100 samples gene expression data
#CRC_100samples_2 <- read.csv(file="data/Si_cohort.csv")

#1.1 combine for 100 samples
CRC_100samples <- cbind(CRC_10samples, CRC_90samples)
symbol_col_index <- which(colnames(CRC_100samples) == "symbol")
columns_to_remove <- c(12)
CRC_100samples <- CRC_100samples[,-columns_to_remove]

#check for zero values column
#CRC_100samples <- CRC_100samples[rowSums(CRC_100samples[, -1] == 0) < (ncol(CRC_100samples) - 1), ]

# check for duplicated rows
duplicated_rows <- duplicated(CRC_100samples$symbol) | duplicated(CRC_100samples)

# identify duplicated rows
duplicated_indices <- which(duplicated_rows)

#filter out duplicated rows
CRC_100samples <- CRC_100samples[-duplicated_indices,]

#write.csv(CRC_100samples, file ="data/Si_RNA_100samples.csv", row.names = F)

#78 samples

CRC_78samples <- CRC_78samples[rowSums(CRC_78samples[, -1] == 0) < (ncol(CRC_78samples) - 1), ]

# Check for duplicated rows
duplicated_rows <- duplicated(CRC_78samples$symbol) | duplicated(CRC_78samples)

# Identify duplicated rows
duplicated_indices <- which(duplicated_rows)

CRC_78samples <- CRC_78samples[-duplicated_indices,]
write.csv(CRC_78samples,file = "data/Si_RNA_78samples.csv", row.names = F)

#2. combine for 178 samples
CRC_178samples <- cbind(CRC_78samples,CRC_10samples, CRC_90samples)

#find out which column is the duplicate
symbol_col_index <- which(colnames(CRC_178samples) == "symbol")
columns_to_remove <- c(80, 91)

#remove duplicated column
CRC_178samples <- CRC_178samples[,-columns_to_remove]

#try to remove rows that are all zeroes apart from the first column (symbol)
# Remove rows where all values are zero except the first column
CRC_178samples <- CRC_178samples[rowSums(CRC_178samples[, -1] == 0) < (ncol(CRC_178samples) - 1), ]

# Check for duplicated rows
duplicated_rows <- duplicated(CRC_178samples$symbol) | duplicated(CRC_178samples)

# Identify duplicated rows
duplicated_indices <- which(duplicated_rows)

CRC_178samples <- CRC_178samples[-duplicated_indices,]
#save 178 samples
#write.csv(CRC_178samples,file = "data/Si_RNA_178samples.csv", row.names=F)

#combine for 208 samples
CRC_208samples <- cbind(CRC_78samples,CRC_10samples, CRC_90samples,CRC_30samples)

#find out which column is the duplicate
symbol_col_index <- which(colnames(CRC_208samples) == "symbol")
ensg_col_index <- which(colnames(CRC_208samples) == "ensgene")
columns_to_remove <- c(80,91,182, 183)

#remove duplicated column
CRC_208samples <- CRC_208samples[,-columns_to_remove]

#try to remove rows that are all zeroes apart from the first column (symbol)
# Remove rows where all values are zero except the first column
CRC_208samples <- CRC_208samples[rowSums(CRC_208samples[, -1] == 0) < (ncol(CRC_208samples) - 1), ]

# Check for duplicated rows
duplicated_rows <- duplicated(CRC_208samples$symbol) | duplicated(CRC_208samples)

# Identify duplicated rows
duplicated_indices <- which(duplicated_rows)

CRC_208samples <- CRC_208samples[-duplicated_indices,]
#save 208 samples
write.csv(CRC_208samples,file = "data/Si_RNA_208samples_raw.csv", row.names=F)
