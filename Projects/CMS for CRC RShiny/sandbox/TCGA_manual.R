rm(list=ls())
graphics.off()

#load the library
TCGA_sample <- read.csv(file="data/matched_TCGA.csv")

#load the CRIS and CMS output
merged_df$Row.names <- str_replace_all(merged_df$Row.names, "[^[:alnum:]]", "-")

CRIS_resu <- read.csv(file="CRIS_COADREAD_prediction_result.xls", sep = "")

test_merge <- merge(merged_df, CRIS_resu, by.x = "Row.names", by.y = "sample.names")