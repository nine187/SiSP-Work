rm(list=ls())

graphics.off()

library(dplyr) #manipulate dataframes

#load CRIS packages
source("src/load_libraries.r") 

#load the COADREAD sample text file from data.txt folder set path later
COADREAD <- read.delim(file = "COADREAD.uncv2.mRNAseq_scaled_estimate.txt", header = TRUE)

#clean the dataset for CRIS classifier, first column should represent gene symbol
# Replace everything after the first "|" with an empty string
modified_names <- sub("\\|.*$", "", rownames(COADREAD))

# Make the modified names unique by adding an index
unique_modified_names <- make.names(modified_names, unique = TRUE)

# Assign the unique modified names back to the data frame
rownames(COADREAD) <- unique_modified_names

#*****get this data back later
#*some rows names are not assigned convert to gene id later now remove it
COADREAD <- COADREAD[-c(1:29),]

#try
#COADREAD <- COADREAD * 1000000

#add a new column for symbol
COADREAD$Symbol <- row.names(COADREAD)

#move the column to the first column
COADREAD <- COADREAD[c("Symbol", names(COADREAD)[-ncol(COADREAD)])]
COADREAD <- COADREAD %>% relocate(Symbol)
#change row names to numberf
rownames(COADREAD) <- 1:nrow(COADREAD)

#save the file as .txt file
write.table(COADREAD, file = "COADREAD_CRIS.txt", quote = FALSE, eol = "\n")

#try reading the .txt files
test <- read.table(file = "COADREAD_CRIS.txt")
test_2 <- read.table(file = "demo.txt")

#source the cris_classifier function
source(file = "cris_classifier_mod.R")

#look at the example dataset
demo <- list.files(pattern="txt.gz$", system.file("data",package="CRISclassifier"), full.names=TRUE)
COADREAD_CRIS <- ("C:/Users/Lenovo/Documents/GitHub/CRIS_single-sample/CRIS_single-sample/data/demo.txt/COADREAD_CRIS.txt")

#use cris_classfier function on the dataset
cris_classifier(input.exp.filename = "COADREAD.uncv2.mRNAseq_scaled_estimate.txt", output.name = "COADREAD")
cris_classifier(input.exp.filename = demo, output.name = "cris", nresmpl = 1)