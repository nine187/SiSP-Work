rm(list=ls())
graphics.off()

library(survival)
library(ggplot2)
library(ggsurvfit)
library(lubridate)
library(survminer)
library(dplyr)

#load cohort data
cohort_data <- read.csv(file="SiSP Work/Projects/CMS for CRC RShiny/data/test_cohort_long2.csv")
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

# Group by Location and Category, then summarize to calculate the frequency
frequency_df <- cohort_data %>%
  group_by(site_of_resection_or_biopsy, CMS_final_network_plus_RFclassifier_in_nonconsensus_samples) %>%
  summarise(Frequency = n())

ggplot(frequency_df, aes(CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, Frequency, fill = site_of_resection_or_biopsy))+
         geom_col() 
         #geom_text(aes(label = Frequency), size = 3, hjust = 0.5, vjust = 3, position = "stack")

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

# Fit survival curves for each category
surv_cms1 <- survfit(Surv(overall_survival, deceased) ~ 1, data = CMS1_df)
surv_cms2 <- survfit(Surv(overall_survival, deceased) ~ 1, data = CMS2_df)
surv_cms3 <- survfit(Surv(overall_survival, deceased) ~ 1, data = CMS3_df)
surv_cms4 <- survfit(Surv(overall_survival, deceased) ~ 1, data = CMS4_df)
surv_nolbl <- survfit(Surv(overall_survival, deceased) ~ 1, data = NOLBL_df)


# Combine survival objects into a list
surv_list <- list(CMS1 = surv_cms1, CMS2 = surv_cms2, CMS3 = surv_cms3, CMS4 = surv_cms4, NOLBL = surv_nolbl)

ggsurvplot_combine(surv_list, data = cohort_data, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                   xlab = "Days", ylab = "Overall survival probability",
                   legend.title = "Location", legend.labs = c("CMS1", "CMS2", "CMS3", "CMS4", "NOLBL"))
ggsurvplot_combine(surv_list, data = cohort_data, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                   xlab = "Days", ylab = "Overall survival probability",
                   legend.title = "Location", legend.labs = c("CMS1", "CMS2", "CMS3", "CMS4", "NOLBL")) +
  labs(title = "Overall Survival Probability by Location")
