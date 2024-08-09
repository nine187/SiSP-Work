rm(list=ls())
graphics.off()

library(ggsurvfit)

#read the Si 100 samples clinical data
Si_100_surv <- read.csv(file="data/Si_100_clinical_surv.csv")
Si_100_surv <- Si_100_surv[1:100,]

fit <- survfit(Surv(time_all, event_death) ~ CMS.classification, data = Si_100_surv)

p <- ggsurvplot(
  fit,
  data = Si_100_surv,
  xlab = "Months",
  ylab = "Overall survival probability",
  palette = c('red', 'blue', 'green', 'yellow'), # Specify colors for each level
  linetype = 1, # Line type for survival curves
  pval = TRUE, # Display p-value
  test.for.trend = TRUE, # Test for trend
  risk.table = FALSE, # Display risk table
  surv.plot.height = 0.25, # Height of the survival plot
  ylim = c(0.85, 1.00) # Set y-axis limits
)

# Customize the plot labels if needed
p$plot <- p$plot + 
  scale_color_manual(values = c('red', 'blue', 'green', 'yellow'), 
                     labels = c('CMS1', 'CMS2', 'CMS3', 'CMS4')) +
  scale_linetype_manual(values = rep(1, 4), 
                        labels = c('CMS1', 'CMS2', 'CMS3', 'CMS4'))

# Display the plot
print(p)

