##Script name:class.R
##
##Purpose of script:
##
##Author: Pasith Prayoonrat
##
##Date Created: 25-1-2024
##
##Email: pasith.p@gmail.com
rm(list=ls())
graphics.off()

#DeepCC
library(DeepCC)#DeepCC
library(tidyverse)
library(keras)
library(reshape2)
library(janitor)
library(magrittr) #data manipulation
library(caret) #confusion matrix
library(dplyr) #for selecting certain matrix
library(ggplot2) #visualization
library(CRISclassifier) #load features for CRIS

#CRIS
source("code/cris_classifier_mod.R")
#browser()
#source the function
source("code/function.R")

#prepare the dataset
data_prep(input = "data/test_data.txt")

#run CRIS
cris_classifier(input.exp.filename = "results/CRIS_data.txt", output.name = "results/")

#run DeepCC
#1.get functional spectra
DeepCC_FS <- getFunctionalSpectra(output_deepcc)
#2.get deep feature
# 2.for CRC model trained with Sample_ID data
CRC_TCGA <- load_DeepCC_model("data/CRC_TCGA")
# 3.obtain deep features 
deepfeature <- get_DeepCC_features(CRC_TCGA, DeepCC_FS)
# 4.get DeepCC label
deeplabel <- get_DeepCC_label(CRC_TCGA, DeepCC_FS)

DeepCC_label <- as.data.frame(deeplabel)
rownames(DeepCC_label) <- rownames(output_deepcc)

#compile and export the data for visualization
write.csv(DeepCC_label, file="results/DeepCC_result.txt")