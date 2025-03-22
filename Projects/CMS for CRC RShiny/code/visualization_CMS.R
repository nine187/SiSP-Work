#### CMS result visualization #### 
# AUTHOR: PASITH PRAYOONRAT
# CONTACT: PASITH.P@GMAIL.COM

library(ggplot2)
#read the classifieR CMS data output for Siriraj patient
library(ggplot2)
library(dplyr)

#sankey diagram
library(networkD3)

#assuming you've already read the CSV and created `Si_classifier` dataframe
table_classifier <- as.data.frame(table(Si_classifier$RF.predictedCMS))
table_classifier_closest <- as.data.frame(table(Si_classifier$RF.nearestCMS))

#create a piechart
CMScolors <- c("#FFA9A9", "#D7BEFF", "#9FE2BF", "#FFE493", "#bfc0c0")

#plot for the CMS
ggplot(data = table_classifier, aes(x = "", y = Freq, fill = Var1)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), 
            color = "black", fontface = "bold", size = 6) +
  scale_fill_manual(values = CMScolors) +
  guides(fill = guide_legend(title = "CMS class")) +
  theme_void() +
  coord_polar("y", start = 0) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"))

#plot for nearest CMS
ggplot(data = table_classifier_closest, aes(x = "", y = Freq, fill = Var1)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), 
            color = "black", fontface = "bold", size = 6) +
  scale_fill_manual(values = CMScolors) +
  guides(fill = guide_legend(title = "CMS class")) +
  theme_void() +
  coord_polar("y", start = 0) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"))

Si_clinical <- read.csv("data/110si_crc.csv")
table(Si_clinical$cms_type)

#plot for Si_result from RSHiny

# Replace NA with "UNC"
Si_result <- na_if(Si_result, NA)
Si_result[is.na(Si_result)] <- "UNC"
table_Si_res_RShiny <- as.data.frame(table(Si_result$CMS.classification))

#plot for nearest CMS
ggplot(data = table_Si_res_RShiny, aes(x = "", y = Freq, fill = Var1)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), 
            color = "black", fontface = "bold", size = 6) +
  scale_fill_manual(values = CMScolors) +
  guides(fill = guide_legend(title = "CMS class")) +
  theme_void() +
  coord_polar("y", start = 0) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"))

#read cohort result
Si_result_RShiny <- read.csv(file="SiSP Work/Projects/CMS for CRC RShiny/data/Si_cohort_res.csv")
Si_result_classifieR <- read.csv(file="SiSP Work/Projects/CMS for CRC RShiny/data/classifieR_output_prob.csv")

Si_all <- cbind(Si_result_classifieR, Si_result_RShiny)
Si_all <- Si_all[,c(1,6,7,9)]
# Replace NA with "UNC" in each row of the Si_all data frame
Si_all[is.na(Si_all)] <- "UNC"

accuracy_cutoff <- table(Si_all$CMS.classification == Si_all$RF.predictedCMS)
accuracy_nocutoff <- table(Si_all$CMS.classification == Si_all$RF.nearestCMS)

#Si_cohort TPM
Si_cohort <- read.csv(file="data/Si_cohort.csv")
# Multiply all columns except the first one by 10^6
Si_cohort <- Si_cohort %>%
  mutate_at(vars(-1), ~ . * 10^6)

Si_cohort <- Si_cohort %>%
  mutate_at(vars(-1), ~ . + 1)

Si_cohort <- Si_cohort %>%
  mutate_at(vars(-1), ~ log2(.) )
write.csv(Si_cohort, file = "data/Si_cohort_TPM.csv", row.names = F)
test <- read.csv(file="data/Si_cohort_TPM.csv")
