rm(list=ls())
graphics.off()

library(survival)
library(ggplot2)
library(ggsurvfit)
library(lubridate)
library(survminer)
library(dplyr)
library(rstatix)

#load cohort data
cohort_data <- read.csv(file="data/test_cohort_long2.csv")
Data <- data.frame(
  Year = c(rep(c("2006-07", "2007-08", "2008-09", "2009-10"), each = 4)),
  Category = c(rep(c("A", "B", "C", "D"), times = 4)),
  Frequency = c(168, 259, 226, 340, 216, 431, 319, 368, 423, 645, 234, 685, 166, 467, 274, 251)
)

ggplot(Data, aes(Year, Frequency, fill = Category)) +
  geom_col() +
  geom_text(aes(label = Frequency), size = 3, hjust = 0.5, vjust = 3, position = "stack")

########################## location barplot ##################################
table <- as.data.frame(table(cohort_data$site_of_resection_or_biopsy))

#clean the dataset for 
cohort_data$site_of_resection_or_biopsy[cohort_data$site_of_resection_or_biopsy %in% c("Sigmoid colon","Ascending colon", "Descending colon", "Splenic flexure of colon", "Rectosigmoid junction","Rectum, NOS")] <- "left"
cohort_data$site_of_resection_or_biopsy[cohort_data$site_of_resection_or_biopsy %in% c("Cecum","Hepatic flexure of colon")] <- "right"
cohort_data$site_of_resection_or_biopsy[cohort_data$site_of_resection_or_biopsy %in% c("Colon, NOS", "Transverse colon" ,"Unknown primary site", "Connective, subcutaneous and other soft tissues of abdomen", "NA ")] <- NA

# Group by Location and Category, then summarize to calculate the frequency
frequency_df <- cohort_data %>%
  group_by(site_of_resection_or_biopsy, CMS_final_network_plus_RFclassifier_in_nonconsensus_samples) %>%
  summarise(Frequency = n())

ggplot(frequency_df, aes(CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, Frequency, fill = site_of_resection_or_biopsy))+
  geom_col()
  #scale_fill_manual(values=c("#CC0000", "#006600", "#669999", "#00CCCC", 
                             #"#660099", "#CC0066", "#FF9999", "#FF9900", 
                             #"black", "purple", "green", "blue"))

################Overall survival plot################################
# Assuming your data has 'time' and 'status' columns
fit <- survfit(Surv(time, status) ~ 1, data = data)
cohort_vis <- ezfun::set_ccf_palette("contrast")
  
# Plotting survival curve
ggsurvfit::survfit2(Surv(time,status) ~ 1, data = data)%>% 
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Overall survival probability")

# Split the dataframe into a list of dataframes based on the values of the Location column
split_data <- split(cohort_data, cohort_data$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples)

# Assign each element of the list to a separate dataframe
CMS1_df <- split_data[["CMS1"]]
CMS2_df <- split_data[["CMS2"]]
CMS3_df <- split_data[["CMS3"]]
CMS4_df <- split_data[["CMS4"]]
NOLBL_df <- split_data[["NOLBL"]]

# Fit survival curves for each category, check later for al
surv_cms1 <- survfit(Surv(overall_survival, deceased) ~ 1, data = CMS1_df)
surv_cms2 <- survfit(Surv(overall_survival, deceased) ~ 1, data = CMS2_df)
surv_cms3 <- survfit(Surv(overall_survival, deceased) ~ 1, data = CMS3_df)
surv_cms4 <- survfit(Surv(overall_survival, deceased) ~ 1, data = CMS4_df)
surv_nolbl <- survfit(Surv(overall_survival, deceased) ~ 1, data = NOLBL_df)

surv_pvalue(CMS1_df, split_data)
surv_pvalue(log_rank_test, split_data)

log_rank_test <- survdiff(Surv(overall_survival, deceased) ~ CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, data = CMS1_df)

# Perform log-rank test for the difference in survival curves between CMS groups
log_rank_test <- survdiff(Surv(overall_survival, deceased) ~ CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, data = cohort_data)
#look at log_rank_test p-value 

# Combine survival objects into a list
surv_list <- list(CMS1 = surv_cms1, CMS2 = surv_cms2, CMS3 = surv_cms3, CMS4 = surv_cms4)
CMS1_2 <- list(CMS1 = surv_cms1, CMS2 = surv_cms2)
ggsurvplot_combine(surv_list, data = cohort_data, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                   xlab = "Days", ylab = "Overall survival probability",
                   legend.title = "Location", legend.labs = c("CMS1", "CMS2", "CMS3", "CMS4")) +
  labs(title = "Overall Survival Probability by Location")

#combine the df, loop this later
CMS1_2_df <- rbind(CMS1_df, CMS2_df)
CMS1_2_pval <- survdiff(Surv(overall_survival, deceased) ~ CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, data = CMS1_2_df)
CMS1_3_df <- rbind(CMS1_df, CMS3_df)
CMS1_3_pval <- survdiff(Surv(overall_survival, deceased) ~ CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, data = CMS1_3_df)
CMS1_4_df <- rbind(CMS1_df, CMS4_df)
CMS1_4_pval <- survdiff(Surv(overall_survival, deceased) ~ CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, data = CMS1_4_df)
CMS2_3_df <- rbind(CMS2_df, CMS3_df)
CMS2_3_pval <- survdiff(Surv(overall_survival, deceased) ~ CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, data = CMS2_3_df)
CMS2_4_df <- rbind(CMS2_df, CMS4_df)
CMS2_4_pval <- survdiff(Surv(overall_survival, deceased) ~ CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, data = CMS2_4_df)
CMS3_4_df <- rbind(CMS3_df, CMS4_df)
CMS3_4_pval <- survdiff(Surv(overall_survival, deceased) ~ CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, data = CMS3_4_df)

# Extract the p-value from the log-rank test
p_value <- log_rank_test[["pvalue"]]

# Create a dataframe to store results
p_value_df <- data.frame(
  CMS_comparison = c("CMS1_2", "CMS1_3", "CMS1_4", "CMS2_3", "CMS2_4", "CMS3_4"),
  p_value = NA
)

# Assign p-values to the dataframe
p_value_df$p_value[1] <- CMS1_2_pval[["pvalue"]]
p_value_df$p_value[2] <- CMS1_3_pval[["pvalue"]]
p_value_df$p_value[3] <- CMS1_4_pval[["pvalue"]]
p_value_df$p_value[4] <- CMS2_3_pval[["pvalue"]]
p_value_df$p_value[5] <- CMS2_4_pval[["pvalue"]]
p_value_df$p_value[6] <- CMS3_4_pval[["pvalue"]]

# Round p-values to two decimal places
p_value_df$p_value <- round(p_value_df$p_value, 2)

# Print the rounded dataframe as a table
kable(p_value_df, caption = "Comparison of CMS Groups - p-values (rounded to 2 decimal places)")


# Print the dataframe as a table
test <- kable(p_value_df, caption = "Comparison of CMS Groups - p-values")