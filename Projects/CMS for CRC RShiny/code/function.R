##Script name:function.R
##
##Purpose of script:a script containing different functions to run the CRC web application website
## the data input requires a gene expression dataset (RNA-seq dataset) (change later)
##Author: Pasith Prayoonrat
##
##Date Created: 25-1-2024
##
##Email: pasith.p@gmail.com

library(tidyverse)#data cleaning
library(tools)#check for data type
library(dplyr)#data frame manipulation

#run using the test prepared data: of 

#function for preparing the data
data_prep <- function(input)
                      #,output_deepcc)
  {
  #browser()
  #input_data = column is patient ID, row Entrez gene ID, must indicate file path 
  #data_prep example file
                      #, output_cris)
  #1.check column for genes and sample id data input, file types, etc.
  if (grepl("\\.txt", input) == FALSE){
    stop("### Please use .txt file type! ###")
  }
    #read the file, set more of this later
    input <- read.csv(file = input, header = TRUE, sep="")

  #2.let the user choose which type (genesymbol, entrez, etc.)
  #check the data if it is in the right value type (ex.RNA-seq scaled estimate, 
  # TPM-transformed, raw count(int))
    
  #how to check if the row is gene symbol, placeholder function
  if (is.integer(rownames(input))){stop("### Please use gene symbol as rownames ###")}
  
  #3.check for number of samples & print out the number of samples in the dataset
  print("The number of samples are:")
  print((ncol(input) - 1))
  
  ###DeepCC data preperation

  output_deepcc <- input
  #remove everythin after | in case a dataset have both gene and 
  #rownames(output_deepcc) <- sub("^[|^]*", "", rownames(output_deepcc))
  
  #remove all |
  rownames(output_deepcc) <- str_replace_all(rownames(output_deepcc), "[^[:alnum:]]", "")
  
  #transform the dataset into the format for deepcc
  
  output_deepcc <- input + 1
  output_deepcc <- input * (10^6)
  output_deepcc <- log2(input)
  output_deepcc <- t(output_deepcc)

  #export the data outside the function
  output_deepcc <<- as.data.frame(output_deepcc) #check if this is the correct deepcc output later
  #CRIS-transform the row/column into the format for the cris_classifier column
  
  ###CRIS data preperation
  output_CRIS <- input
  output_CRIS['Symbol'] <- NA
  
  #clean the dataset for CRIS classifier, first column should represent gene symbol
  # Replace everything after the first "|" with an empty string
  modified_names <- sub("\\|.*$", "", rownames(output_CRIS))
  
  # Make the modified names unique by adding an index
  unique_modified_names <- make.names(modified_names, unique = TRUE)
  
  # Assign the unique modified names back to the data frame
  rownames(output_CRIS) <- unique_modified_names
  
  #remove the first few rows (delete this later)
  output_CRIS <- output_CRIS[-c(1:29),]
  
  #move the column to the first column
  output_CRIS <- output_CRIS[c("Symbol", names(COADREAD)[-ncol(COADREAD)])]
  output_CRIS <- output_CRIS %>% relocate(Symbol)

  #put symbol in first column
  output_CRIS[,1] <- rownames(output_CRIS)
  
  #change row names to numberf
  rownames(output_CRIS) <- 1:nrow(output_CRIS)
  
  output_CRIS <<- output_CRIS
  
  #save the file
  write.table(output_CRIS, file = "CRIS_data.txt", quote = FALSE, eol = "\n")
}

##check why the data output is inf tomorrow
#test the function using test_data.txt
data_prep(input="test_data.txt")

#try Deepcc function with the data


#try CRIS classifier with the data

cris_classifier(input.exp.filename = "CRIS_data.txt")