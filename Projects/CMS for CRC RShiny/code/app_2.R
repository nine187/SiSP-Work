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

# Call library for plotting the sankey network
library(networkD3)

#survival analysis curve
library(survival)
library(shiny)
library(ggplot2)
library(ggsurvfit)
library(lubridate)
library(dplyr)
library(survminer)

#visNetwork 
library(visNetwork)

#SECTION1 table
library(gtsummary)

#for umap
library(umap)

source("code/function_1.R")

# Modify the limitation of size
options(shiny.maxRequestSize = 1000*1024^2)

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
  
  # Show the content of prediction as a table
  output$Content <- renderTable({
    GEP <- GEPdataInput()
    GEP <- GEP[, -1]
    Subtype <- Subtypeprediction()
    Predictionresult <- cbind("Sample ID" = colnames(GEP),
                              Subtype)
    #Add CRIS prediction probability
    return(Predictionresult)
  })
  
  #Show the content of prediction for CMS - fix this code later
  #output$CMS_prob <- renderTable({
  #prob <- Subtypeprediction
  #GEP <- GEPdataInput()
  #Subtype <- Subtypeprediction()
  #prob_2 <- cbind("Sample ID" = colnames(GEP),
  #Subtype)
  #return(prob_2)
  #})
  
  #cohort visualization
  cohort_vis <- ezfun::set_ccf_palette("contrast")
  
  output$surv_plot <- renderPlot({
    req(input$cohort_file)
    
    data <- read.csv(input$cohort_file$datapath)
    
    # Assuming your data has 'time' and 'status' columns
    fit <- survfit(Surv(overall_survival, deceased) ~ 1, data = data)
    
    # Plotting survival curve
    ggsurvfit::survfit2(Surv(overall_survival, deceased) ~ 1, data = data)%>% 
      ggsurvfit() +
      labs(
        x = "Days",
        y = "Overall survival probability")
  })
  
  #Plotting the piechart for CRC location
  output$loc_plot <- renderPlot({
    req(input$cohort_file)
    data <- read.csv(input$cohort_file$datapath)
    loc <- as.data.frame(table(data$site_of_resection_or_biopsy))
    ggplot(loc, aes(x = Var1, y = Freq)) +
      geom_bar(stat = "identity", fill = "blue") +  # Basic bar graph
      labs(title = "Bar Graph", x = "Category", y = "Value") +  # Add title and axis labels
      theme_minimal()  })  # Optional: set a minimal theme
  
  output$loc_plot_CMS <- renderPlot({
    req(input$cohort_file)
    cohort_data <- read.csv(input$cohort_file$datapath)
    
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
  })
  
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
  
  output$surv_plot_CMS <- renderPlot({
    req(input$cohort_file)
    cohort_data <- read.csv(input$cohort_file$datapath)
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
    
    # Combine survival objects into a list
    surv_list <- list(CMS1 = surv_cms1, CMS2 = surv_cms2, CMS3 = surv_cms3, CMS4 = surv_cms4)
    
    ggsurvplot_combine(surv_list, data = cohort_data, risk.table = TRUE, pval = TRUE,
                       xlab = "Days", ylab = "Overall survival probability",
                       legend.title = "Location", legend.labs = c("CMS1", "CMS2", "CMS3", "CMS4")) +
      labs(title = "Overall Survival Probability by Location")
  })
  output$networkCMS_output <- renderVisNetwork({
    if(input$networkCMS == "netCMS1_pred") {
      nodes <- data.frame(id = 1:3,
                          label = paste(c("CMS1 predictive", 
                                          "HSP90 inhibition", 
                                          "Anti-PD1 therapy with significantly longer PFS")),
                          group = c("GrA", "GrB", "GrC"),
                          value = 1:3,
                          shape = c("circle", "square", "square"),
                          title = paste0("<p><b>", 1:3,"</b><br>Node !</p>"),
                          shadow = c(TRUE, TRUE, FALSE))
      
      edges <- data.frame(from = c(1,1), to = c(2,3),
                          label = c("predict response"),
                          length = c(500, 500)
      )
      
      visNetwork(nodes, edges, 
                 main = "CMS1", 
                 submain = list(text = "Predictive Biomarker",
                                style = "font-family:Comic Sans MS;color:#ff0000;font-size:15px;text-align:center;"), 
                 footer = "Fig.1 Example", 
                 height = "500px",
                 width = "100%")
    }
    else {
      treCMS1nodes <- data.frame(id = 1:5,
                                 
                                 #add labels on node
                                 label = paste(c("CMS1 treatment", 
                                                 "Pembrolizumab", 
                                                 "Bevacizumab",
                                                 "oxaliplatin + bevacizumab",
                                                 "fluorouracil")),
                                 
                                 #add groups on nodes,
                                 group = c("PMIDA", "PMIDB", "PMIDC","PMIDD","PMIDE"),
                                 
                                 #control shape of nodes
                                 shape = c("circle", "square", "square", "square", "square"),
                                 
                                 # tooltip (html or character), when the mouse is above
                                 title = paste0("<p><b>", 1:5,"</b><br>Node !</p>")
                                 
      )
      treCMS1edges <- data.frame(from = c(1,1,1,1), to = c(2,3,4,5),
                                 #add labels to the edge
                                 label = c("benefit treatment", "benefit treatment", 
                                           "first line", "doesn't benefit"),
                                 length = c(500,500,500,500)
      )
      visNetwork(treCMS1nodes, treCMS1edges, 
                 main = "CMS1", 
                 submain = list(text = "Treatment",
                                style = "font-family:Comic Sans MS;color:#ff0000;font-size:15px;text-align:center;"), 
                 footer = "Fig.1 Example", 
                 height = "500px",
                 width = "100%")
    }
  })
  # Read CSV file
  CMS1_pred <- reactive({
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - CMS1 - Predictive.csv")
  })
  CMS2_pred <- reactive({
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - CMS2 - Predictive.csv")
  })
  CMS3_pred <- reactive({
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - CMS3- Predictive.csv")
  })
  CMS4_pred <- reactive({
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - CMS4 - Predictive.csv")
  })
  CMS1_prog <- reactive({
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - CMS1 - Prognostic .csv")
  })
  CMS2_prog <- reactive({
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - CMS2 - Prognostic.csv")
  })
  CMS3_prog <- reactive({
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - CMS3 - Prognostic.csv")
  })
  CMS4_prog <- reactive({
    read.csv("data/CMS Knowledge based table summary (credit_ น้องนศพ) - CMS4 - Prognostic.csv")
  })
  output$table_CMS <- renderTable({
    if(input$CMS_df == "CMS1_pred"){
      return(CMS1_pred())
    } 
    if(input$CMS_df == "CMS2_pred"){
      return(CMS2_pred())
    }
    if(input$CMS_df == "CMS3_pred"){
      return(CMS3_pred())
    }
    if(input$CMS_df == "CMS4_pred"){
      return(CMS4_pred())
    }
    if(input$CMS_df == "CMS1_prog"){
      return(CMS1_prog())
    }
    if(input$CMS_df == "CMS2_prog"){
      return(CMS2_prog())
    }
    if(input$CMS_df == "CMS3_prog"){
      return(CMS3_prog())
    }
    else{
      return(CMS4_prog())
    }}
  )
}

TPM_options <- c("TPM", "no TPM")

# Construct UI
ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel(h1("SiSP Colorectal Cancer Subtyping Platform V 0.12", align = "center"),
             windowTitle = "SiSP Colorectal Cancer Subtyping Platform V 0.12"),
  # Divide columns of the UI page
  fluidPage(
    column(width = 12, 
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
                              label = "Colorectal Cancer Subgroup Classification(s)",
                              selected = "CMS",
                              choices = list("CMS", "CRIS")),
           # Horizontal line ----
           tags$hr(),
           radioButtons("TPM", "TPM option", TPM_options),
           # Construct an action button for prediction
           actionButton(inputId = "Classprediction", 
                        label = "Button"),
           p("Press the 'Predict' button after the classification(s) is chosen. \n
             The input file must be a .csv file. \n
             The column must be patient ID and row must be Gene Symbol."),
           
           #cohort data file
           fileInput(inputId = "cohort_file", 
                     label = "Choose cohort file (.CSV/ .TXT file)",
                     multiple = FALSE,
                     accept = c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv")),
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
      tableOutput("CMS_Summary"),
      tableOutput("CRIS_Summary"),
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
                        plotOutput(outputId = "loc_plot"),
                        plotOutput(outputId = "surv_plot"),
                        plotOutput(outputId = "loc_plot_CMS"),
                        plotOutput(outputId = "surv_plot_CMS"),
                        #plotOutput(outputId = "CMS_prob")
                        #CMS Summary Table
                        p("SECTION 3: Clinical Linke (Prognostic + Predictive Biomarker"),
                        selectInput(inputId = "networkCMS",
                                    label = "CMS Network:",
                                    choices = c("CMS1-Predictive Biomarker" = "netCMS1_pred", 
                                                "CMS1-Treatment" = "netCMS2_pred"),
                                    selected = "netCMS1_pred"),
                        visNetworkOutput("networkCMS_output"),
                        #####################################
                        #button for dataframe choices
                        selectInput(inputId = "CMS_df",
                                    label = "Subtypes:",
                                    choices =  c("CMS1 predictive biomarker + Treatment" = "CMS1_pred",
                                                 "CMS2 predictive biomarker + Treatment" = "CMS2_pred",
                                                 "CMS3 predictive biomarker + Treatment" = "CMS3_pred",
                                                 "CMS4 predictive biomarker + Treatment" = "CMS4_pred",
                                                 "CMS1 prognostic biomarker" = "CMS1_prog",
                                                 "CMS2 prognostic biomarker" = "CMS2_prog",
                                                 "CMS3 prognostic biomarker" = "CMS3_prog",
                                                 "CMS4 prognostic biomarker" = "CMS4_prog")),
                        tableOutput("table_CMS"),
                 )
               )
             )
      )
    )))

shinyApp(ui = ui, server = server)

#test on Siriraj 100 RNA and clinical data samples
