#### CRC BATCH CORRECTION CODE #### 
# AUTHOR: PASITH PRAYOONRAT
# CONTACT: PASITH.P@GMAIL.COM

rm(list=ls())
graphics.off()

library(sva) # ComBat

#load Si cohort TPM data
CRC_208samples <- read.csv(file="data/Si_RNA_208samples.csv") #TPM data
CRC_208samples_raw <- read.csv(file="data/Si_RNA_208samples_raw.csv")

row.names(CRC_208samples) <- CRC_208samples[,1]
CRC_208samples <- CRC_208samples[,-1]

#example RNA-seq combat code
count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
batch <- c(rep(1, 4), rep(2, 4))
group <- rep(c(0,1), 4)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group, full_mod=TRUE)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=NULL, full_mod=FALSE)
adjusted_counts <- as.data.frame(adjusted_counts)

#assign batch (known batch - 1. 78 samples, 2. 100 samples, 3. 30 samples)
CRC_batch <- c(rep(1,78), rep(2,100), rep(3,30))

#include condition group variable
CRC_208samples.combat <- as.data.frame(ComBat_seq(CRC_208samples, batch = CRC_batch))

#save output
write.csv(x = CRC_208samples.combat, file = "data/Si_RNA_208_raw_samples_combat.csv", row.names = F)

#load Si cohort raw count data
#CRC_208samples_raw <- read.csv(file="data/")
