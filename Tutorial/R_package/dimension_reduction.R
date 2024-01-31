rm(list=ls())
graphics.off()

library(ggplot2)
library(ggrepel)
library(cowplot)

library(plspm)
data("cereals")
head(cereals)