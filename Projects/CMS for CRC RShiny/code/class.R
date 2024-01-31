##Script name:class.R
##
##Purpose of script:
##
##Author: Pasith Prayoonrat
##
##Date Created: 25-1-2024
##
##Email: pasith.p@gmail.com
#rm(list=ls())
#graphics.off()

library(DeepCC)#DeepCC
library(CRISclassifier)#CRIS 
#source the function
source("code/function.R")

#prepare the dataset
data_prep(input)

#run DeepCC

#run CRIS

#cris example

demo <- list.files(pattern="txt", system.file("data",package="CRISclassifier"), full.names=TRUE)
demo
cris_classifier(input.exp.filename = demo, output.name="cris", nresmpl=1)

data(demo) #example CRIS data
cris_classifier(input.exp.filename = output_CRIS)
#tomorrow: check if both data are the same type, how to make both dataset compatible