library(shiny)
library(tidyverse)
library(shinydashboard)
library(MASS)
library(ggplot2)
library(shinythemes)
library(shiny)
library(tidyverse)
library(shinydashboard)
library(MASS)
library(ggplot2)
library(shinythemes)

ui <- fluidPage(
  theme = shinytheme("lumen"),
  titlePanel("Data Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", multiple = FALSE, accept = ".csv"),
      br(),
      actionButton("analyzeBtn", "Statistical Analysis Results")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data", DT::dataTableOutput('tbl')), # Data dalam tabel
        tabPanel("Summary Data", verbatimTextOutput("summaryText")),
        tabPanel("Box Plot", 
                 fluidPage(
                   sidebarLayout(
                     sidebarPanel(
                       h5("select variables for the boxplot:"),
                       selectInput("var_box", "Variabel:", ""),
                       
                     ),
                     plotOutput("boxPlot")))),
        tabPanel("Bar Plot", 
                 fluidPage(
                   sidebarLayout(
                     sidebarPanel(
                       h5("select variables for the barplot:"),
                       selectInput("var_bar", "Variabel:", "")
                       
                     ),
                     plotOutput("barPlot")))),
        tabPanel("Data and ANOVA Test Results", 
                 fluidPage(
                   titlePanel("Analysis of Variance"),
                   sidebarLayout(
                     sidebarPanel(
                       selectInput("vardipen", label = h3("Dependen"),
                                   list("y" = "y"), selected = 1),
                       
                       selectInput("varindepen", label = h3("Independen"),
                                   list("group" = "group"), selected = 1)
                       
                     ),
                     
                     mainPanel(
                       titlePanel("Data"),
                       DT::dataTableOutput('tbl1'),
                       titlePanel("Result"),
                       verbatimTextOutput("anovaText")
                     )
                   )
                 )
        )
      )
    )
  )
)


server <- function(input, output, session) {
  data1 <- data.frame(
    y = c(2.5,2.7,2.8,2.6,3.0,2.4,2.9,2.5,2.6,2.7,
          3.8,3.5,4.0,3.7,3.9,3.6,4.1,3.4,3.8,3.9,
          3.1,2.9,3.0,3.2,3.3,2.8,3.4,3.1,3.2,3.5),   # variabel dependen
    group = rep(c("Left Sidebar", "Center Page", "Right Sidebar"), each = 10)  # variabel independen (faktor)
  )
  data <- reactive({
    req(input$file)
    read.csv(input$file$datapath, sep = ";")
  })
  
  # Output Data
  output$tbl = DT::renderDataTable({
    DT::datatable(data(), options = list(lengthChange = FALSE))
  })
  
  output$summaryText <- renderPrint({
    if (input$analyzeBtn > 0) {
      summary(data())
    }
  })
  
  output$anovaText <- renderPrint({
    if (input$analyzeBtn > 0) {
      anova_result <- aov(y ~ group, data = data1)
      summary(anova_result)
    }
  })
  output$tbl1 = DT::renderDataTable({
    if (input$analyzeBtn > 0) {
      DT::datatable(data1, options = list(lengthChange = FALSE))
    }
  })
  
  output$boxPlot <- renderPlot({
    if (input$analyzeBtn > 0) {
      req(input$var_box)
      boxplot(t(data()[, input$var_box]), 
              main = paste("Boxplot of", input$var_box,"per Hari"),
              xlab = "Day", ylab = input$var_box)
    }
  })
  # Update pilihan variabel pada sidebar
  observe({
    updateSelectInput(session, "var_box", choices = c("Left.Sidebar", "Center.Page", "Right.Sidebar"))
  })
  
  output$barPlot <- renderPlot({
    if (input$analyzeBtn > 0) {
      req(input$var_bar)
      
      # Gambar bar plot
      barplot(t(data()[, input$var_bar]),
              main = paste("Bar Plot of", input$var_bar, "per Hari"),
              xlab = "Day", ylab = input$var_bar)
    }
  })
  # Update pilihan variabel pada sidebar
  observe({
    updateSelectInput(session, "var_bar", choices = c("Left.Sidebar", "Center.Page", "Right.Sidebar"))
  })
}

shinyApp(ui = ui, server = server)