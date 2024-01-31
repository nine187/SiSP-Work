rm(list=ls())
library(Biobase)
library(MmCMS)

#import the dataset
dataset <- crcTCGAsubset
test_dataset <- TestData_gemm

#assay data
dataset_expr <- crcTCGAsubset@assayData$exprs

#prepare the tcga dataset for MmCMS function
#a numeric expression matrix with sample columns, and HGNC symbol rownames.
#Data should be normalized. see the example in TestData_gemm.

MmCMS(crcTCGAsubset)
