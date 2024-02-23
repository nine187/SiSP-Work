# Call the involving library
library(shiny)
library(shinyjs)
# Call org.Hs.eg.db library from Bioconductor package
library(org.Hs.eg.db)
# Call the libraries for DeepCC
library(reticulate)
use_condaenv("r-reticulate")
library(tensorflow)
library(keras)
library(DeepCC)
library(CRISclassifier)
# Call the libraries for plotting
library(ggplot2)

# Call library for plotting the sankey network
library(networkD3)

# Modify the limitation of file size
options(shiny.maxRequestSize = 30*1024^2)

# Load the trained DeepCC model for CMS classification
#prefix <- "~/R_Sisyspharm/CMS/DeepCC/DeepCC_TCGA_scale-estimate_model/DeepCC_CRC_TCGA-scale-estimate-456sam-600epoch-relu_model"
prefix <- "data/CRC_TCGA"
classifier <- load_DeepCC_model(prefix)

# Load the selected feature data set used for CRIS classification
data("features")

## CMS CLASSIFICATION ##
CMSprediction <- function(GEPdata, classifier){
  browser()
  ## CMS classification utilizing DeepCC ##
  # Call the gene expression data
  #GEPdata <- GEPdataInput()
  # Define the row names of the gene expression data
  rownames(GEPdata) <- GEPdata[, 1]
  # Remove the first column of the 
  GEPdata <- GEPdata[, -1]
  # Log2-transformation of gene expression data before performing functional spectra
  # Transpose the gene expression data set before performing functional spectra
  GEPdata <- t(log2(GEPdata + 1))
  
  # The column names of a gene expression data should be Entrez ID of genes
  # Get column names of gene expression data and map symbol to entrez id
  #GeneENT  <- mapIds(org.Hs.eg.db, colnames(GEPdata) , 
                                   #'ENTREZID', 'SYMBOL')
  
  # After converting symbol to entrez id, there are NA samples in the data
  # This loop for finding NA positions of entrez id data
  #GeneENTnapos <- which(is.na(GeneENT[, 1]))
  # Remove NA positions from gene expression data (some columns will be removed)
  #GEPdata <- GEPdata[, -as.integer(GeneENTnapos)]
  # Remove some genes containing NA in the Entrez ID genes
  #GeneENTnona <- as.data.frame(GeneENT[-as.integer(GeneENTnapos), 1])
  # Change column names from symbol to entrez id
  #colnames(GEPdata) <- GeneENTnona[, 1]
  
  # Classify new data set by utilizing the trained DeepCC model
  Freqspectra <- getFunctionalSpectra(GEPdata)
  print("fs")
  #prefix <- "~/R_Sisyspharm/Coding/Platform/DeepCC_CRC_TCGA-scale-estimate-456sam-600epoch-relu_model"
  #classifier <- load_DeepCC_model(prefix)
  Predlabel <- get_DeepCC_label(classifier, Freqspectra, cutoff = 0.5)
  print("pd")
  CMSpred <- as.data.frame(as.character(Predlabel), stringsAsFactors = FALSE)
  colnames(CMSpred) <- "CMS classification"
  # Return the CMS prediction result
  return(CMSpred)
}

## CRIS CLASSIFICATION ##
CRISprediction <- function(GEPdata, GEPsamples){
  browser()
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
  colnames(CRISpred) <- "CRIS classification"
  
  # Return the CRIS prediction result
  return(CRISpred)
}

server <- function(input, output, session) {
  # Get gene expression data set from user
  GEPdataInput <- reactive({
    # Require an external file from user
    req(input$GEPfile)
    # Check file extension --> csv or txt
    if (length(grep("csv$", input$GEPfile)) > 0){
      GEPdata <- read.csv(input$GEPfile$datapath,
                          header = TRUE)
    }
    if (length(grep("txt$", input$GEPfile)) > 0){
      GEPdata <- read.delim(input$GEPfile$datapath,
                            header = TRUE)
    }
    return(GEPdata)
  })
  
  # Get a value from the check box group as a function
  ClasscheckInput <- reactive({
    length(input$Classcheck)
  })
  # Observe the status of the check box group, and the uploaded input file and be disable/ able for the action button
  observe({
    if (ClasscheckInput() == 0 | is.null(input$GEPfile$datapath)){
      disable("Classprediction")
    }
    if (ClasscheckInput() > 0 & !is.null(input$GEPfile$datapath)) {
      enable("Classprediction")
    }
  })
  
  Subtypeprediction <- eventReactive(input$Classprediction, {
    if (length(input$Classcheck) == 1) {
      if (input$Classcheck == "CMS"){
        Subtype <- CMSprediction(GEPdata = GEPdataInput(), 
                                 classifier = classifier)
      }else{
        Subtype <- CRISprediction(GEPdata = GEPdataInput(), 
                                  GEPsamples = GEPdataInput())
      }
    }
    if (length(input$Classcheck) == 2) {
      CMSsubtype <- CMSprediction(GEPdata = GEPdataInput(), 
                                  classifier = classifier)
      CRISsubtype <- CRISprediction(GEPdata = GEPdataInput(),
                                    GEPsamples = GEPdataInput())
      Subtype <- cbind(CMSsubtype, CRISsubtype, stringsAsFactors = FALSE)
      print(class(Subtype[, 1]))
      print(class(Subtype[, 2]))
      print(class(Subtype[, 1][1]))
    }
    return(Subtype)
  })
  
  # Plot pie chart for showing the proportion of classification labels
  output$CMSpie <- renderPlot({
    # Prepare data for plotting
    Subtype <- Subtypeprediction()
    Subtype <- Subtype$`CMS classification`
    Subtype[is.na(Subtype)] <- "N/A"
    Subtype <- data.frame(table(Subtype))
    # Define color for plotting
    CMScolors <- c("#FFA9A9", "#D7BEFF", "#9FE2BF", "#FFE493", "#bfc0c0")
    if (length(input$Classcheck) == 1) {
      if (input$Classcheck == "CMS"){
        Piechart <- ggplot(data = Subtype, aes(x = 2, y = Freq, fill = Subtype)) +
          geom_col(color = "white") +
          coord_polar(theta = "y") +
          geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), 
                    color = "black", fontface = "bold", size = 6) +
          scale_fill_manual(values = CMScolors) +
          guides(fill = guide_legend(title = "CMS class")) +
          theme_void() +
          xlim(0.5, 2.5) +
          theme(legend.position = "right",
                legend.title = element_text(face = "bold"))
      }
    }
    if (length(input$Classcheck) == 2) {
      Piechart <- ggplot(data = Subtype, aes(x = 2, y = Freq, fill = Subtype)) +
        geom_col(color = "white") +
        coord_polar(theta = "y") +
        geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), 
                  color = "black", fontface = "bold", size = 6) +
        scale_fill_manual(values = CMScolors) +
        guides(fill = guide_legend(title = "CMS class")) +
        theme_void() +
        xlim(0.5, 2.5) +
        theme(legend.position = "right",
              legend.title = element_text(face = "bold"))
    }
    return(Piechart)
  })
  
  # Plot pie chart for showing the proportion of classification labels
  output$CRISpie <- renderPlot({
    # Prepare data for plotting
    Subtype <- Subtypeprediction()
    Subtype <- Subtype$`CRIS classification`
    Subtype[is.na(Subtype)] <- "N/A"
    Subtype <- data.frame(table(Subtype))
    # Define color for plotting
    CRIScolors <- c("#D7263D", "#F46036", "#2E294E", "#1B998B", "#00B7E0", "#bfc0c0")
    if (length(input$Classcheck) == 1) {
      if (input$Classcheck == "CRIS"){
        Piechart <- ggplot(data = Subtype, aes(x = 2, y = Freq, fill = Subtype)) +
          geom_col(color = "white") +
          coord_polar(theta = "y") +
          geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), 
                    color = "white", fontface = "bold", size = 6) +
          scale_fill_manual(values = CRIScolors) +
          guides(fill = guide_legend(title = "CRIS class")) +
          theme_void() +
          xlim(0.5, 2.5) +
          theme(legend.position = "right",
                legend.title = element_text(face = "bold"))
      }
    }
    if (length(input$Classcheck) == 2) {
      Piechart <- ggplot(data = Subtype, aes(x = 2, y = Freq, fill = Subtype)) +
        geom_col(color = "white") +
        coord_polar(theta = "y") +
        geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), 
                  color = "white", fontface = "bold", size = 6) +
        scale_fill_manual(values = CRIScolors) +
        guides(fill = guide_legend(title = "CRIS class")) +
        theme_void() +
        xlim(0.5, 2.5) +
        theme(legend.position = "right",
              legend.title = element_text(face = "bold"))
    }
    return(Piechart)
  })
  
  output$CMSCRISsankey <- renderSankeyNetwork({
    if (length(input$Classcheck) == 2) {
      # Define node
      Subtype <- Subtypeprediction()
      CMSpred <- Subtype$`CMS classification`
      CMSpred[is.na(CMSpred)] <- "N/A"
      CMSnode <- data.frame(table(CMSpred))
      colnames(CMSnode) <- c("Var1", "Freq")
      print(CMSnode)
      
      CRISpred <- Subtype$`CRIS classification`
      CRISpred[is.na(CRISpred)] <- "N/A"
      CRISnode <- data.frame(table(CRISpred))
      colnames(CRISnode) <- c("Var1", "Freq")
      print(CRISnode)
      
      Node <- data.frame(rbind(CMSnode, CRISnode))
      Node <- data.frame(as.character(Node$Var1))
      colnames(Node) <- "name"
      
      # Define link
      Link <- data.frame("source" = CMSpred, 
                         "target" = CRISpred, 
                         "value" = c(1:length(CMSpred)))
      print(Link)
      Sourcepos <- vector(mode = "integer")
      for (i in 1:length(Link$source)) {
        pos <- which(Link$source[i] == CMSnode$Var1) 
        pos <- pos - 1
        Sourcepos <- append(Sourcepos, pos)
      }
      print(Sourcepos)
      Targetpos <- vector(mode = "integer")
      for (i in 1:length(Link$target)) {
        pos <- which(Link$target[i] == CRISnode$Var1) 
        pos <- pos - 1 + length(CMSnode$Var1)
        Targetpos <- append(Targetpos, pos)
      }
      print(Targetpos)
      
      Link$source <- Sourcepos
      Link$target <- Targetpos
      print("ok")
      
      # Combine data in a list
      CRCsubtypes <- list("links" = Link, "nodes" = Node)
      
      # Add a factor column to list in the nodes 
      Group <- as.factor(Node$name)
      CRCsubtypes$nodes$group <-  Group
      print(CRCsubtypes)
      
      node_color <- 'd3.scaleOrdinal() .domain(["CMS1", "CMS2", "CMS3", "CMS4", "CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E", "N/A"]) .range(["#FFA9A9", "#D7BEFF", "#9FE2BF", "#FFE493", "#D7263D", "#F46036", "#2E294E", "#1B998B", "#00B7E0", "#bfc0c0"])'
      # Plot
      sankeyNetwork(Links = CRCsubtypes$links, Nodes = CRCsubtypes$nodes, 
                    Source = 'source', Target = 'target', Value = 'value', 
                    NodeID = 'name', units = 'TWh', fontSize = 12, 
                    nodeWidth = 30, iterations = 0, colourScale = node_color,
                    NodeGroup = "group")
    }
  })      
  
  
  
  # Show the content of prediction as a table
  output$Content <- renderTable({
    GEP <- GEPdataInput()
    GEP <- GEP[, -1]
    Subtype <- Subtypeprediction()
    Predictionresult <- cbind("Sample ID" = colnames(GEP),
                              Subtype)
    return(Predictionresult)
  })
  
}

# Construct UI
ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel("SiSP LAB"),
  # Divide columns of the UI page
  fluidRow(
    column(width = 2, 
           # Input: Select a file ----
           fileInput(inputId = "GEPfile", 
                     label = "Choose Gene Expression File (.CSV/ .TXT file)",
                     multiple = FALSE,
                     accept = c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv")),
           # Horizontal line ----
           tags$hr(),
           # Input: Checkbox group if file has header ----
           checkboxGroupInput(inputId = "Classcheck",
                              label = "Colorectal Cancer Subgruop Classification(s)",
                              selected = "CMS",
                              choices = list("CMS", "CRIS")),
           # Horizontal line ----
           tags$hr(),
           # Construct an action button for prediction
           actionButton(inputId = "Classprediction", 
                        label = "Predict"),
           p("Press the 'Predict' button after the classification(s) is chosen.")
    ),
    column(width = 10,
           fluidRow(
             column(width = 5, 
                    plotOutput(outputId = "CMSpie")), 
             column(width = 5, 
                    plotOutput(outputId = "CRISpie")),
             fluidRow(
               column(width = 10,
                      sankeyNetworkOutput(outputId = "CMSCRISsankey"),
                      tableOutput(outputId = "Content"))
             )
           )
    )
  )
)

shinyApp(ui = ui, server = server)
