library(magrittr)
samp.id <- data.frame("IDs" = c("A", "B", "C"))
View(samp.id)
samp.dat <- data.frame("sampl" = c("A","B","C","D"))
View(samp.dat)
samp.dat$sampl %in% samp.id$IDs
samp.dat[samp.dat$sampl %in% samp.id$IDs, ]