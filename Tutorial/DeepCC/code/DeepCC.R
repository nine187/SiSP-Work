rm(list=ls())
graphics.off()
#setwd
library(tidyverse)
library(DeepCC)
library(keras)
library(reshape2)
library(janitor)
library(magrittr) #data manipulation
library(caret) #confusion matrix
library(dplyr) #for selecting certain matrix
library(ggplot2) #visualization

#import the datasets to the environment
CMS_label <- read.table(file = "cms_labels_public_all.txt", header = TRUE)
COADREAD <- read.delim(file = "COADREAD.uncv2.mRNAseq_scaled_estimate.txt", header = TRUE)
Sample_ID <- read.csv(file = "TCGA_CRC_scale-estimate_SampleID_456samples.csv")

#explore dataset
#head(CMS_label)

################DATA CLEANING##########################################
#1.CMS_label
#remove other dataset except tcga
CMS_label<- subset(CMS_label, (CMS_label$dataset == "tcga"))
#remove dataset-tcga row (maybe keep it?)
CMS_label <- subset(CMS_label, select = -dataset)

#2.COADREAD
#remove gene symbol in COADREAD dataset look for packages to convert
#head(rownames(COADREAD)), use regmatches to remove all the data from |
rownames(COADREAD) <- sub("^[|^]*", "", rownames(COADREAD))

#try to remove the | from the row
rownames(COADREAD) <- str_replace_all(rownames(COADREAD), "[^[:alnum:]]", "")
#https://stackoverflow.com/questions/66910817/how-to-remove-everything-before-first-occurrence-of-comma-in-r?noredirect=1&lq=1

#3.TCGA 
#change the - to . so that colnames are the same

################DATA WRANGLING/DATA PREPERATION#################################

#match COADREAD sample ID with TCGA dataset
#transpose the column names
COADREAD_tcga <- as.data.frame(colnames(COADREAD))
#use double backspace to escape and use the dot
a <- gsub("\\.", "-", COADREAD_tcga$`colnames(COADREAD)`[1:677])
a <- data.frame(ColumnNames = c(a))
COADREAD_tcga$SampleID <- gsub("\\.", "-", COADREAD_tcga$`colnames(COADREAD)`[1:677])
COADREAD_tcga <- subset(COADREAD_tcga, select = 2)
rm(a)
#check if there is the same value in both dataframes
identical(names(COADREAD_tcga$ModifiedColumnNames), names(Sample_ID$SampleID))
#merge the dataframes, check later
filtered_SampleID <- merge(COADREAD_tcga, Sample_ID)
#mutate data to be character

#remove samples that are not in the filtered_SampleID dataset in the COADREAD dataset
#transpose the filtered_SampleID dataset
filtered_SampleID <- t(filtered_SampleID)
#move the values to colnames with janitor package
filtered_SampleID <- filtered_SampleID %>%
                    row_to_names(row_number = 1)

#fix! to do all 3 at the same time later probably with * or something
colnames(filtered_SampleID) <- str_replace_all(colnames(filtered_SampleID), "[^[:alnum:]]", ".")
#colnames(filtered_SampleID) <- sub("-", ".", colnames(filtered_SampleID), fixed = TRUE)
#colnames(filtered_SampleID) <- sub("-", ".", colnames(filtered_SampleID), fixed = TRUE)
#colnames(filtered_SampleID) <- sub("-", ".", colnames(filtered_SampleID), fixed = TRUE)

#change filtered_SampleID from character matrix to dataframe (same type as COADREAD)
filtered_SampleID <- as.data.frame(filtered_SampleID)

#convert COADREAD from double to character
colnames(COADREAD[1:677]) <- as.character(colnames(COADREAD[1:677]))

#create a list of characters for the filtered sample
colnames(COADREAD) %in% colnames(filtered_SampleID)
list_filter<- colnames(COADREAD)[colnames(COADREAD) %in% colnames(filtered_SampleID) ]

# check column name that are in the list, and remove those that are not
COADREAD <- COADREAD %>% select(one_of(list_filter))

#replace . with -
#COADREAD_tcga <- gsub(".", "-", COADREAD_tcga)
#remove column names that doesn't overlap with the intersect function

#transform the data for DeepCC algorithm
#multiply the value by 10^6 and log2-transformed the data 
COADREAD_tf <- COADREAD*(10^6)
#add 1/0.25 to avoid result of log2 of 0
COADREAD_tf <- COADREAD_tf + 1
COADREAD_tf <- log2(COADREAD_tf)

#swap column and row for DeepCC data
COADREAD_tf <- t(COADREAD_tf)

#whitespace
colnames(COADREAD_tf) <- gsub(" ", "", colnames(COADREAD_tf))

###########################MODEL FITTING/MODELLING#############################

#fit DeepCC model to the data
test <- getFunctionalSpectra(COADREAD)
#save the parameter as csv 

#write.csv(test, file = "DeepCC_result.csv")
#labels <- NA

load_DeepCC_model <- function(prefix){
  load(file = paste0(prefix, ".RData"))
  classifer <- keras::load_model_hdf5(filepath =paste0(prefix, ".hdf5"))
  list(classifier = classifer, levels = levels)
}
#change file path later
# for CRC model trained with Sample_ID data
CRC_TCGA <- load_DeepCC_model("CRC_TCGA")

# obtain deep features 
deepfeature <- get_DeepCC_features(CRC_TCGA, test)

#get DeepCC label
deeplabel <- get_DeepCC_label(CRC_TCGA, test)
DeepCC_label <- as.data.frame(deeplabel)
#DeepCC_label <- as.data.frame(table(deeplabel))

#match DeepCC_label with deepfeature CMS label (assuming the CMS result is in the same order)
rownames(DeepCC_label) <- rownames(deepfeature)

#modify the CMS_label to be the same as DeepCC_label dataframe
#change the rownames into tcga
rownames(CMS_label) <- CMS_label[, c(1)]
#left only the final CMS result
CMS_label <- CMS_label[, -c(1:3), drop = FALSE]

#change the format of the rownames to be the same, can be replace with the same function above
rownames(CMS_label) <- sub("-", ".", rownames(CMS_label), fixed = TRUE)
rownames(CMS_label) <- sub("-", ".", rownames(CMS_label), fixed = TRUE)

#remove .xxx in DeepCC label
rownames(DeepCC_label) <- sub(".01", "", rownames(DeepCC_label), fixed = TRUE)

#match the CMS_label with deepCC_lable
rownames(DeepCC_label) %in% rownames(CMS_label)
list_filter <- rownames(DeepCC_label)[rownames(DeepCC_label) %in% rownames(CMS_label) ]

#subset DeepCC label based on the filter
subsetted_DeepCC_label <- DeepCC_label[list_filter, ]

# merge the two data frames based on row names
merged_df <- merge(CMS_label, subsetted_DeepCC_label, by = 0, all.x = TRUE)

merged_df <- merge(CMS_label, DeepCC_label, by.x)
row.names(merged_df) <- rownames(merged_df$Row.names)

#remove rows with NA and NOLBL from the dataframe
df_cleaned <- merged_df[complete.cases(merged_df), ]
df_cleaned <- df_cleaned[df_cleaned$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples != 'NOLBL', ]
df_cleaned <- df_cleaned[df_cleaned$y != 'NOLBL', ]

#save workspace here
#create a confusion matrix of the merged dataset
confusionMatrix(data = merged_df$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples,
                reference = merged_df$y)
cleaned_label <- as.factor(df_cleaned$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples)
cleaned_feature <- as.factor(df_cleaned$y)
cf <- confusionMatrix(cleaned_label, cleaned_feature)
fourfoldplot(as.fawc(cf),color=c("yellow","pink"),main = "Confusion Matrix")

#visualization
cm <- confusionMatrix(factor(y.pred), factor(y.test), dnn = c("DeepCC", "Reference"))

plt <- as.data.frame(cf$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="blue") +
  labs(x = "DeepCC",y = "Prediction") +
  scale_x_discrete(labels=c("CMS4","CMS3","CMS2","CMS1")) +
  scale_y_discrete(labels=c("CMS1","CMS2","CMS3","CMS4"))

#confusion matrix with true positive, false positive,...
#Test with CMS labels
#dat <- as.factor(CMS_label$CMS_RFclassifier)
#ref <- as.factor(CMS_label$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples)
#head(dat)
#head(ref)
#confusionMatrix(dat, ref)