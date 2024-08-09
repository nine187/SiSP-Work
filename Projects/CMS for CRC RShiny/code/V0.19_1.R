##Script name:app.R
##
##Purpose of script: A script containing Rshiny app for running a sample CRC CMS Web Application 
##Author: Jantappapa C., Pasith P.
##
##Date Created: 25-1-2024
##
##Email: pasith.p@gmail.com

# Call the involving library
library(shiny)
library(shinyjs)

# Call the libraries for plotting
library(ggplot2)

# Call library for plotting the sankey network and knowledge-based database
library(networkD3)

#survival analysis curve
library(survival)
library(shiny)
library(ggplot2)
library(lubridate)
library(dplyr)

#visNetwork 
library(visNetwork)

#SECTION1 table
library(gtsummary)

#for umap
library(umap)

#json 
library(jsonlite)

# Load required library for Excel writing
library(openxlsx)

source("code/function_1.R")

# Modify the limitation of size
options(shiny.maxRequestSize = 1000*1024^2)

server <- function(input, output, session) {
  ############# UMAP code ###########################
  cms.col <- c("red", "blue", "green", "yellow")
  #load the csv file
  tcga.feature <- reactive({
    read.csv("data/umap_TCGA.csv")
  })
  
  si.feature <- reactive({
    read.csv("data/umap_Si.csv")
  })
  
  #load the config
  custom.config <- umap.defaults
  custom.config$min_dist <- 0.4
  custom.config$n_neighbors <- 30
  
  #umap plot for tcga
  output$umap_TCGA <- renderPlot({
    tcga.ref_umap.dat <- tcga.feature()
    ggplot(tcga.ref_umap.dat, aes(x = umap1, y = umap2)) +
      geom_point() +
      labs(title = "TCGA_UMAP: Dimension 1 vs Dimension 2")+
      stat_ellipse(data = tcga.ref_umap.dat, 
                   aes(x = umap1, y = umap2, color = cms.lab), 
                   type = "norm", linetype = 2, level = 0.95)+
      geom_point(data = tcga.ref_umap.dat, 
                 aes(x = umap1, y = umap2, color = cms.lab), 
                 size = 3, alpha = 0.85, shape = 16)+
      scale_color_manual(values = cms.col)+
      theme(legend.position = c(0.05, 0.05),
            legend.justification = c(0.05, 0.05), 
            legend.text = element_text(size = 14),
            legend.title = element_text(face = "bold", size = 14), 
            text = element_text(family = "Arial"),
            axis.title.x = element_text(face = "bold", size = 16), 
            axis.title.y = element_text(face = "bold", size = 16), 
            axis.line = element_line(size = 1), 
            axis.text.x = element_text(size = 14, vjust = 0.5), 
            axis.text.y = element_text(size = 14),
            strip.text.x = element_text(size = 14, face = "bold"),
            panel.background = element_rect(color = "white", fill = "white"), 
            strip.background = element_rect(color = "white", fill = "white"),
            panel.grid.major = element_line(color = "gray", linewidth = 0.25), 
            panel.grid.minor = element_line(color = "gray87", linewidth = 0.12))})
  
  
  #umap plot for si
  output$umap_si <- renderPlot({
    si.ref_umap.dat <- si.feature()
    ggplot(si.ref_umap.dat, aes(x = umap1, y = umap2)) +
      geom_point() +
      labs(title = "Si_UMAP: Dimension 1 vs Dimension 2")+
      stat_ellipse(data = si.ref_umap.dat, 
                   aes(x = umap1, y = umap2, color = cms.lab), 
                   type = "norm", linetype = 2, level = 0.95)+
      geom_point(data = si.ref_umap.dat, 
                 aes(x = umap1, y = umap2, color = cms.lab), 
                 size = 3, alpha = 0.85, shape = 16)+
      scale_color_manual(values = cms.col)+
      theme(legend.position = c(0.05, 0.05),
            legend.justification = c(0.05, 0.05), 
            legend.text = element_text(size = 14),
            legend.title = element_text(face = "bold", size = 14), 
            text = element_text(family = "Arial"),
            axis.title.x = element_text(face = "bold", size = 16), 
            axis.title.y = element_text(face = "bold", size = 16), 
            axis.line = element_line(size = 1), 
            axis.text.x = element_text(size = 14, vjust = 0.5), 
            axis.text.y = element_text(size = 14),
            strip.text.x = element_text(size = 14, face = "bold"),
            panel.background = element_rect(color = "white", fill = "white"), 
            strip.background = element_rect(color = "white", fill = "white"),
            panel.grid.major = element_line(color = "gray", linewidth = 0.25), 
            panel.grid.minor = element_line(color = "gray87", linewidth = 0.12))})
  
  ##### Data Input ###########
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
  
  #CMS Summary Table
  CMS_Summary <- reactive({
    read.csv("data/CMS Summary Table.csv")
  })
  
  output$CMS_Summary <- renderTable({
    return((CMS_Summary()))
  })
  
  #CRIS Summary Table
  CRIS_Summary <- reactive({
    read.csv("data/CRIS Summary.csv")
  })
  
  output$CRIS_Summary <- renderTable({
    return((CRIS_Summary()))
  })
  
  #Apply the CMS prediction function to the inputted data
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
    CMScolors <- c("red", "blue", "green", "yellow", "#bfc0c0")
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
  
  #sankey diagram
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
  
  # Define a reactiveValues object to store predicted data
  predictedData <- reactiveValues(data = NULL)
  
  # Render the table
  output$Content <- renderTable({
    GEP <- GEPdataInput()
    GEP <- GEP[, -1]
    Subtype <- Subtypeprediction()
    Predictionresult <- cbind("Sample ID" = colnames(GEP), Subtype)
    predictedData$data <- Predictionresult  # Store in reactiveValues
    return(Predictionresult)
  })
  
  # Download handler for downloading the table as an Excel file
  output$downloadPrediction <- downloadHandler(
    filename = function() {
      paste("Prediction_Result", Sys.Date(), ".xlsx", sep = "_")
    },
    content = function(file) {
      # Retrieve data from reactiveValues
      data <- predictedData$data
      write.xlsx(data, file, rowNames = FALSE)
    }
  )
  #Show the content of prediction for CMS - fix this code later
  #output$CMS_prob <- renderTable({
  #prob <- Subtypeprediction
  #GEP <- GEPdataInput()
  #Subtype <- Subtypeprediction()
  #prob_2 <- cbind("Sample ID" = colnames(GEP),
  #Subtype)
  #return(prob_2)
  #})
  
  #function to read the input file
  clinical_data <- reactive({
    req(input$cohort_file)
    cohort_data <- read.csv(input$cohort_file$datapath)
    return(cohort_data)
  })
  
  #bind the dataframe of Subtype prediction to clinical data
  #figure out how to make the input wait for the subtype output
  clinical_data <- reactive({
    req(input$cohort_file, Subtypeprediction())
    df <- read.csv(input$cohort_file$datapath, header = TRUE) # Assuming clinical data is loaded from a file
    df2 <- Subtypeprediction()
    #bind only the elements we wanted from the dataframe
    combined <- cbind(df, df2)
    #use gt summary to compute a table
    combined_table <- tbl_summary(
      combined,
      include = c(sex, age, metastasis, stage, location),
      by = "CMS classification",
      missing = "no",
    ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") %>%
      bold_labels()
    return(combined_table)
  })
  
  #function to render the clinical data 
  output$summary_table <- renderTable({
    clinical_data()
  })
  
  #json table, check this code later
  json_data <- reactive({
    fromJSON("code/pybel.indra.json")})
  
  output$json <-renderTable(
    json_data()
  )
  
  output$networkCMS_output <- renderVisNetwork({
    if(input$networkCMS == "netCMS_pred") { 
  nodes <- data.frame(id = 1:8,
                      label = paste(c("CMS1", 
                                      "pMMR_CRC", 
                                      "dMMR_CRC",
                                      "recurrence rate",
                                      "adjuvant chemotherapy",
                                      "anti PD1 therapy",
                                      "pembrolizumab",
                                      "bevacizumab")),
                      group = c("GrA", "GrB", "GrC", "GrD", "GrE", "GrF", "GrG", "GrH"),
                      value = 1:8,
                      shape = "dot",
                      size = "30",
                      title = paste0("<p><b>", 1:8, "</b><br>Node !</p>"),
                      color = "grey",
                      scaling = list(min=30,max = 30))
  
  edges <- data.frame(from = c(1, 1, 2, 3, 4), 
                      to = c(6, 8, 4, 7, 5),
                      label = c("clinical benefit", "clinical benefit", "association", "clinical benefit", "inhibition"),
                      length = 200,
                      font = list(size = 20))  # Shorter edge length to reduce overlap
  
  visNetwork(nodes, edges, 
             main = "CMS", 
             submain = list(text = "Predictive Biomarker",
                            style = "font-family:Comic Sans MS;color:#ff0000;font-size:15px;text-align:center;"), 
             footer = "database gathered from 10 papers", 
             height = "700px",
             width = "100%") %>%
    visPhysics(stabilization = FALSE) %>%
    visOptions(selectedBy = "group", 
               highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visEdges(smooth = FALSE, color = list(color = "black")) %>%
    visLayout(randomSeed = NA, improvedLayout = TRUE) %>%
    visHierarchicalLayout(direction = "LR", levelSeparation = 500) %>%
    visInteraction(navigationButtons = TRUE) %>%
    visNodes(font = list(size = 20))  # Adjust font size here
    
    # Use improvedLayout for better spreading
    
    }
    else {
      treCMS1nodes <- data.frame(id = 1:12,
                                 
                                 #add labels on node
                                 label = paste(c("CMS1", 
                                                 "CMS1_MSS_BRAF", 
                                                 "five_year_OS",
                                                 "CMS1_TP53_mut",
                                                 "CDK5",
                                                 "overall_survival",
                                                 "BRAF_mCRC_CDX_loss",
                                                 "BRAF_mCRC_CK7_loss",
                                                 "LIMK1",
                                                 "CMS2",
                                                 "CMS3",
                                                 "worse_OS")),
                                 
                                 #add groups on nodes,
                                 group = c("A","B","C","D","E","F","G","H","I","J","K","L"),
                                 
                                 #control shape of nodes
                                 shape = c("circle", "circle", "circle", "circle", "circle" ,"circle"
                                           ,"circle","circle","circle","circle","circle","circle"),
                                 
                                 # tooltip (html or character), when the mouse is above
                                 title = paste0("<p><b>", 1:12,"</b><br>Node !</p>",
                                 color = "black", size = 25)
                                 
      )
      treCMS1edges <- data.frame(from = c(1,1,1,1,2,4,7,8,9,10,11), to = c(2,4,5,12,6,3,6,6,6,9,9),
                                 #add labels to the edge
                                 label = c("Associate", "Associate","Associate","Associate",
                                           "Inhibition", "Inhibition","Inhibition","Inhibition",
                                           "Association", "Upregulation", "Upregulation"),
                                 length = c(500,500,500,500,500,500,500,500,500,500,500)
      )
      visNetwork(treCMS1nodes, treCMS1edges, 
                 main = "CMS", 
                 submain = list(text = "Prognostic Biomarker",
                                style = "font-family:Comic Sans MS;color:#ff0000;font-size:15px;text-align:center;"), 
                 footer = "database gathered from 10 papers", 
                 height = "500px",
                 width = "100%")
    }})
  
  # Define reactive function to read CMS Prognostic data
  CMS_prog <- reactive({
    req(input$CMS_df)
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - RShiny-Prognostic.csv")
  })
  
  # Define reactive function to read CMS Predictive data
  CMS_pred <- reactive({
    req(input$CMS_df)
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - RShiny-Predictive.csv")
  })
  
  # Render table based on selected input
  output$table_CMS <- renderTable({
    if (input$CMS_df == "CMS_prog") {
      return(CMS_prog())
    } else {
      return(CMS_pred())
    }
  })
}

# Construct UI
ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel(h1("SiSP Colorectal Cancer Subtyping Platform V 0.19", align = "center"),
             windowTitle = "SiSP Colorectal Cancer Subtyping Platform V 0.19"),
  # Divide columns of the UI page
  fluidPage(
    column(width = 12, 
           # Input: Select a file ----
           fileInput(inputId = "GEPfile", 
                     label = "Upload gene expression matrix (counts/normalised expression) in .csv format.",
                     multiple = FALSE,
                     accept = c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv")),
           # Horizontal line ----
           tags$hr(),
           # Input: Checkbox group if file has header ----
           checkboxGroupInput(inputId = "Classcheck",
                              label = "Colorectal Cancer Subgroup Classification(s)",
                              selected = "CMS",
                              choices = list("CMS", "CRIS")),
           # Horizontal line ----
           tags$hr(),
           #radioButtons("TPM", "TPM option", TPM_options),
           # Construct an action button for prediction
           actionButton(inputId = "Classprediction", 
                        label = "Button"),
           p("Press the 'Predict' button after the classification(s) is chosen. \n
             The input file must be a .csv file. \n
             The column must be patient ID and row must be Gene Symbol."),
           
           #cohort data file
           fileInput(inputId = "cohort_file", 
                     label = "Upload clinical data file in .csv format (the first column must contain the matching ID with the 
                     first column of the RNA-seq .csv file)",
                     multiple = FALSE,
                     accept = c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv")),
           tableOutput("CMS_Summary"),
           tableOutput("CRIS_Summary"),
           #actionButton(inputId = "Cohort_vis",
           #label = "Visualization"),
           #remove the visualization buttion
           #p("Press the button to visualize the cohort data.
           #The input file must be a .csv file
           #The column must be patient ID, days, status, and location.")
    ),
    ########################################
    p("SECTION 1: Tumor Identity"),
    # Show a table of the parsed data
    mainPanel(
      tableOutput("summary_table"),
      ########################################
      p("SECTION 2: Subtype Classification"),
      column(width = 10,
             fluidRow(
               column(width = 5, 
                      plotOutput(outputId = "CMSpie")), 
               column(width = 5, 
                      plotOutput(outputId = "CRISpie")),
               fluidRow(
                 column(width = 10,
                        sankeyNetworkOutput(outputId = "CMSCRISsankey"),
                        
                        #cohort data visualization
                        tableOutput(outputId = "Content"),
                        #plotOutput(outputId = "CMS_prob")
                        
                        #UMAP
                        plotOutput(outputId = "umap_TCGA"),
                        plotOutput(outputId = "umap_si"),
                        downloadButton("downloadPrediction", "Download Prediction Result"),
                        #################################################################
                        #CMS Summary Table
                        p("SECTION 3: Clinical Link (Prognostic + Predictive Biomarker)"),
                        selectInput(inputId = "networkCMS",
                                    label = "CMS Network:",
                                    choices = c("Predictive Biomarker" = "netCMS_pred", 
                                                "Prognostic Biomarker" = "netCMS_prog"),
                                    selected = "netCMS_pred"),
                        visNetworkOutput("networkCMS_output"),
                        #####################################
                        #button for dataframe choices
                        selectInput(inputId = "CMS_df",
                                    label = "Subtypes:",
                                    choices =  c("CMS predictive biomarker" = "CMS_pred",
                                                 "CMS prognostic biomarker" = "CMS_prog")),
                        tableOutput("table_CMS"),
                 )
               )
             )
      ))))

shinyApp(ui = ui, server = server)

