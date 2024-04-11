##Script name:function.R
##
##Purpose of script:a script containing different functions to run the CRC web application website
## the data input requires a gene expression dataset (RNA-seq dataset) (change later)
##Author: Pasith Prayoonrat
##
##Date Created: 25-1-2024
##
##Email: pasith.p@gmail.com

# Call org.Hs.eg.db library from Bioconductor package
library(org.Hs.eg.db)
# Call the libraries for DeepCC
library(reticulate)
use_condaenv("r-reticulate")
library(tensorflow)
library(keras)
library(DeepCC)
library(CRISclassifier)
library(Biobase)
library("AnnotationDbi")


# Load the trained DeepCC model for CMS classification
#prefix <- "~/R_Sisyspharm/CMS/DeepCC/DeepCC_TCGA_scale-estimate_model/DeepCC_CRC_TCGA-scale-estimate-456sam-600epoch-relu_model"
prefix <- "data/CRC_TCGA"
classifier <- load_DeepCC_model(prefix)

# Load the selected feature data set used for CRIS classification
data("features")

## CMS CLASSIFICATION ##
CMSprediction <- function(GEPdata, classifier){
  ## CMS classification utilizing DeepCC ##
  # Call the gene expression data
  #GEPdata <- GEPdataInput()
  # Define the row names of the gene expression data
  rownames(GEPdata) <- GEPdata[, 1]
  # Remove the first column of the 
  GEPdata <- GEPdata[, -1]
  
  #convert to entrez
  # The column names of a gene expression data should be Entrez ID of genes
  # Get column names of gene expression data and map symbol to entrez id
  GeneENT  <- mapIds(org.Hs.eg.db, rownames(GEPdata) , 'ENTREZID', 'SYMBOL')
  
  # After converting symbol to entrez id, there are NA samples in the data
  # This loop for finding NA positions of entrez id data
  GeneENTnapos <- which(is.na(GeneENT))
  # Remove NA positions from gene expression data (some columns will be removed)
  GEPdata <- GEPdata[-as.integer(GeneENTnapos),]
  # Remove some genes containing NA in the Entrez ID genes
  GeneENTnona <- as.data.frame(GeneENT[!is.na(GeneENT)])
  # Change column names from symbol to entrez id
  rownames(GEPdata) <- GeneENTnona$`GeneENT[!is.na(GeneENT)]`
  
  # Log2-transformation of gene expression data before performing functional spectra
  # Transpose the gene expression data set before performing functional spectra
  #TPM option
  GEPdata <- GEPdata*(10^6) #log2(TPM + 1)
  GEPdata <- GEPdata + 1
  GEPdata <- t(log2(GEPdata))
  #no TPM option (already TPM transformed)
  # GEPdata <- t(GEPdata)
  # Classify new data set by utilizing the trained DeepCC model
  print("fs")
  Freqspectra <- getFunctionalSpectra(GEPdata)
  #prefix <- "~/R_Sisyspharm/Coding/Platform/DeepCC_CRC_TCGA-scale-estimate-456sam-600epoch-relu_model"
  #classifier <- load_DeepCC_model(prefix)
  Predlabel <- get_DeepCC_label(classifier, Freqspectra, cutoff = 0.5)
  Predlabel_prob <- get_DeepCC_label(classifier, Freqspectra, cutoff = 0.5, 
                                     prob_mode = T, prob_raw = T)
  Predlabel_prob <- as.data.frame(Predlabel_prob)
  colnames(Predlabel_prob) <- c("CMS1", "CMS2","CMS3", "CMS4")
  print("pd")
  CMSpred <- as.data.frame(as.character(Predlabel), stringsAsFactors = FALSE)
  colnames(CMSpred) <- "CMS classification"
  CMSpred <- cbind(CMSpred, Predlabel_prob)
  # Return the CMS prediction result
  return(CMSpred)
}

## CRIS CLASSIFICATION ##
CRISprediction <- function(GEPdata, GEPsamples){
  # suppressWarnings()
  temp.nn.wt <- "TRUE"
  dist.selection <- "cosine"
  GenePattern.output <- "TRUE"
  rnd.seed <- 7392854
  nresmpl <- 1000
  
  # Advanced setting
  norm.method <- "row.std" # "row.std.ref","ratio.ref"
  within.sig <- "FALSE"
  dchip.output <- "FALSE"
  signature.heatmap <- "TRUE"
  FDR.sample.bar <- 0.2 # NA if not needed
  plot.FDR <- "TRUE"
  col.range <- 3 # SD in heatmap
  heatmap.legend <- signature.heatmap
  histgram.null.dist <- "FALSE" # Histgram of null dist for the distance
  hist.br <- 30
  
  # For GenePattern
  nresmpl <- as.numeric(nresmpl)
  col.range <- as.numeric(col.range)
  hist.br <- as.numeric(hist.br)  # Bin number for resampled dist histgram
  rnd.seed <- as.numeric(rnd.seed)
  if (FDR.sample.bar != "NA"){
    FDR.sample.bar <- as.numeric(FDR.sample.bar)
    if (is.numeric(FDR.sample.bar) == FALSE){
      stop("### Provide numerical value (0~1) for FDR.sample.bar! ###")
    }
  }
  
  # Set random seed
  set.seed(rnd.seed)
  
  # File format check
  if (length(features[1,])!= 3 & length(features[1,])!= 4){
    stop("### Please use features file format! ###")
  }
  if (length(features[1,]) < 4 & temp.nn.wt == "TRUE"){
    temp.nn.wt <- "FALSE"
  }
  third.col <- rownames(table(features[, 3]))
  if (is.na(as.numeric(third.col[1]))){
    stop("### TRUEhe 3rd column of feature file should be numerical! ###")
  }
  feat.col.names <- colnames(features)
  feat.col.names[1:2] <- c("ProbeID","GeneName")
  colnames(features) <- feat.col.names
  num.features <- length(features[, 1])
  num.cls <- length(table(features[, 3]))
  feature.col.num <- length(features[1, ])
  ord <- seq(1:num.features)
  features <- cbind(ord,features)  # Add order column to "features"
  
  # Expression data
  # Call expression data file (.csv or .txt)
  #GEPdata <- GEPdataInput()
  # Change the first column name
  colnames(GEPdata)[1] <- c("GeneName")
  # Other dataset's mean & SD for row normalization (optional)
  ProbeID <- GEPdata[,1]
  gene.names <- GEPdata[,1]
  num.samples <- (length(GEPdata[1,])-1)
  exp.dataset <- GEPdata
  #GEPsamples <- GEPdataInput()
  sample.names <- as.vector(colnames(GEPsamples)[-1])
  print(num.samples)
  
  #rownames
  rownames(exp.dataset) <- exp.dataset[,1]
  exp.dataset <- exp.dataset[,-1]
  
  # Row normalize
  normed.exp.dataset <- exp.dataset
  if (norm.method == "row.std"){
    exp.mean <- apply(exp.dataset, 1, mean, na.rm=TRUE)
    exp.sd <- apply(exp.dataset, 1, sd, na.rm=TRUE)
    normed.exp.dataset <- (exp.dataset - exp.mean)/exp.sd   
  }
  normed.exp.dataset<-cbind(ProbeID, normed.exp.dataset)
  
  # extract features from normed.exp.dataset
  out <- list()
  out[[1]] <- features
  out[[2]] <- normed.exp.dataset
  exp.dataset.extract <- merge(features, normed.exp.dataset, sort = FALSE)
  if (length(exp.dataset.extract[,1]) < 1){
    stop("### No matched probes! ###")
  }
  
  order.extract <- order(exp.dataset.extract[, 2])
  exp.dataset.extract <- exp.dataset.extract[order.extract,]
  order.extract.after <- exp.dataset.extract[, 2]
  exp.dataset.extract <- exp.dataset.extract[-2]
  
  if (temp.nn.wt == "FALSE"){
    features.extract <- exp.dataset.extract[, 1:3]
    if (feature.col.num == 4){
      exp.dataset.extract <- exp.dataset.extract[-4]
    }
    features.extract <- cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)
    num.features.extract <- length(features.extract[, 1])
    
    ProbeID.extract <- as.vector(exp.dataset.extract[, 1])
    exp.dataset.extract <- exp.dataset.extract[-c(1:3)]
    rownames(exp.dataset.extract) <- ProbeID.extract
  }
  
  if (temp.nn.wt == "TRUE" & num.cls == 2){
    features.extract <- exp.dataset.extract[, 1:4]
    features.extract <- cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)
    temp.nn.wt.vector <- as.numeric(as.vector(features.extract[, 5]))
    if (is.numeric(temp.nn.wt.vector) == FALSE){
      stop("# Please use numeric values in 4th column!#")
    }
    num.features.extract <- length(features.extract[, 1])
    ProbeID.extract <- as.vector(exp.dataset.extract[, 1])
    exp.dataset.extract <- exp.dataset.extract[-c(1:4)]
    rownames(exp.dataset.extract) <- ProbeID.extract
  }
  
  # make template
  
  for (i in 1:num.cls){
    temp.temp <- as.numeric(as.vector(features.extract[, 4]))
    temp.temp[temp.temp!=i] <- 0
    temp.temp[temp.temp==i] <- 1
    eval(parse(text = paste("temp.", i, "<-temp.temp", sep="")))
  }
  
  # Weighted template (only for 2cls)
  if (temp.nn.wt == "TRUE" & num.cls == 2){
    temp.1 <- temp.nn.wt.vector
    temp.2 <- -temp.nn.wt.vector
  }
  
  # Compute distance and p-value
  predict.label <- vector(length = num.samples, mode = "numeric")
  dist.to.template <- vector(length = num.samples, mode = "numeric")
  dist.to.cls1 <- vector(length = num.samples, mode = "numeric")
  
  rnd.feature.matrix <- matrix(0, 
                               nrow = num.features.extract, 
                               ncol = nresmpl)
  
  perm.dist.vector <- vector(length = nresmpl*num.cls, mode = "numeric")
  nominal.p <- vector(length = num.samples, mode = "numeric")
  BH.FDR <- vector(length = num.samples, mode = "numeric")
  Bonferroni.p <- vector(length = num.samples, mode = "numeric")
  
  for (i in 1:num.samples){
    print(paste("sample # ", i, sep = ""))
    current.sample <- as.vector(exp.dataset.extract[, i])
    # Compute original distance
    orig.dist.to.all.temp <- vector(length = num.cls, mode = "numeric")
    if (temp.nn.wt == "TRUE"){   # Weight sample data
      current.sample <- current.sample*abs(temp.nn.wt.vector)
    }
    if (dist.selection == "cosine"){
      for (o in 1:num.cls){      # Compute distance to all templates
        current.temp <- vector()
        eval(parse(text = paste("current.temp <- temp.", o, sep = "")))
        orig.dist.to.all.temp[o] <-
          sum(current.temp * current.sample, na.rm = TRUE) /
          (sqrt (sum(current.temp ^ 2, na.rm = TRUE)) * sqrt(sum(current.sample ^ 2, 
                                                                 na.rm = TRUE)))
      }
    }
    if (dist.selection == "correlation"){
      for (o in 1:num.cls){      # Compute distance to all templates
        eval(parse(text = paste("current.temp <- temp.", o, sep = "")))
        orig.dist.to.all.temp[o] <- cor(current.temp, 
                                        current.sample, 
                                        method = "pearson", 
                                        use = "complete.obs")
      }
    }
    
    if (num.cls == 2){           # Find nearest neighbor (2 classes)
      if (orig.dist.to.all.temp[1] >= orig.dist.to.all.temp[2]){
        predict.label[i] <- 1
        dist.to.template[i] <- 1-orig.dist.to.all.temp[1]
        dist.to.cls1[i] <- -(orig.dist.to.all.temp[1] + 1)
      }
      if (orig.dist.to.all.temp[1] < orig.dist.to.all.temp[2]){
        predict.label[i] <- 2
        dist.to.template[i] <- 1-orig.dist.to.all.temp[2]
        dist.to.cls1[i] <- orig.dist.to.all.temp[2] + 1
      }
    }
    
    if (num.cls > 2){
      for (o in 1:num.cls){       # find nearest neighbor (>2 classes)
        if (is.na(orig.dist.to.all.temp[o]) != TRUE){
          if (orig.dist.to.all.temp[o] == max(orig.dist.to.all.temp, 
                                              na.rm = TRUE)){
            predict.label[i] <- o
            dist.to.template[i] <- 1-orig.dist.to.all.temp[o]
            dist.to.cls1[i] <- (1-orig.dist.to.all.temp[o]) + o
          }
        }
      }
    }
    # Permutation test
    if (within.sig == "FALSE"){     # Generate resampled features from all probes
      for (p in 1:nresmpl){
        rnd.feature.matrix[, p] <- sample(normed.exp.dataset[, (i+1)], 
                                          num.features.extract, 
                                          replace = FALSE)
      }
    }
    if (within.sig == "TRUE"){     # Generate resampled features from only signature genes
      for (p in 1:nresmpl){
        rnd.feature.matrix[, p] <- sample(exp.dataset.extract[, i], 
                                          num.features.extract, 
                                          replace = FALSE)
      }
    }
    
    if (temp.nn.wt == "TRUE" & num.cls == 2){
      rnd.feature.matrix <- rnd.feature.matrix * abs(temp.nn.wt.vector)
    }
    
    # Compute distance to all templates
    if (dist.selection == "cosine"){          # Cosine
      for (res in 1:num.cls){
        temp.resmpl <- vector()
        eval(parse(text = paste("temp.resmpl<-temp.", res,sep = "")))
        prod.sum <- sum(t(t(rnd.feature.matrix) * temp.resmpl), na.rm = TRUE)
        data.sq.sum <- sum(rnd.feature.matrix ^ 2, na.rm = TRUE)
        temp.sq.sum <- sum(temp.resmpl ^ 2, na.rm = TRUE)
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)] <-
          (1-(prod.sum/(sqrt(data.sq.sum)*sqrt(temp.sq.sum))))
      }
    }
    if (dist.selection=="correlation"){          # Correlation
      for (res in 1:num.cls){
        eval(parse(text = paste("temp.resmpl<-temp.", res,sep = "")))
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-
          (1-as.vector(cor(rnd.feature.matrix, 
                           temp.resmpl, 
                           method = "pearson", 
                           use="complete.obs")))
      }
    }
    
    # Compute nominal p-value
    combined.stats.rank <- rank(c(dist.to.template[i],perm.dist.vector))
    nominal.p[i] <- combined.stats.rank[1]/length(combined.stats.rank)
    
    # Histgram of combined null distributions
    #  if (histgram.null.dist == "TRUE" & capabilities("png") == TRUE){
    #    png(paste("resampled_", dist.selection, "_dist_histgram_", sample.names[i], ".png",sep=""))
    #    hist(c(dist.to.template[i],perm.dist.vector),br=hist.br,main=paste(sample.names[i],", # resampling: ",nresmpl,sep=""))
    #  }
    
  } # main sample loop END
  print("out of main loop")
  # MCTRUE correction
  BH.FDR <- nominal.p * num.samples/rank(nominal.p)
  Bonferroni.p <- nominal.p * num.samples
  BH.FDR[BH.FDR > 1] <- 1
  Bonferroni.p[Bonferroni.p > 1] < -1
  
  ### output ###
  # prediction results
  #return(dist.to.template)
  dist.to.cls1.rank <- rank(dist.to.cls1)
  predict.label2 <- vector()
  predict.label2[predict.label == 1] <- "CRIS-A"
  predict.label2[predict.label == 2] <- "CRIS-B"
  predict.label2[predict.label == 3] <- "CRIS-C"
  predict.label2[predict.label == 4] <- "CRIS-D"
  predict.label2[predict.label == 5] <- "CRIS-E"
  CRISpred <- as.data.frame(predict.label2, stringsAsFactors = FALSE)
  
  #bind the dist.to.cls1 probabilty to the table
  CRISpred[,2] <- dist.to.cls1
  
  #assign the column name to the dataframe
  colnames(CRISpred) <- c("CRIS classification", "CRIS: distance to cluster")
  
  # Return the CRIS prediction result
  return(CRISpred)
}

############################RShiny function################################
#load file in RShiny
load_file <- function(name, path) {
  ext <- tools::file_ext(name)
  switch(ext,
         csv = vroom::vroom(path, delim = ","),
         tsv = vroom::vroom(path, delim = "\t"),
         validate("Invalid file; Please upload a .csv or .tsv file")
  )
}
