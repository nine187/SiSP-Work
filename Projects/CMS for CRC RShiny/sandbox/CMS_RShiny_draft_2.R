library(shiny)
#source("data/function.R")

ui <- fluidPage(
  fileInput("upload", NULL, accept = c(".csv", ".tsv", ".txt")),
  #Divide the column of the UI page
  fluidRow(column(width = 2,
                  fileInput(inputID =))))

server <- function(input, output, session) {
  options(shiny.maxRequestSize=1000*1024^2) #set file size limit to 1GB
  
  data <- reactive({
    req(input$upload)
    test <- read.csv(file = input$upload$datapath, header = TRUE)
    return(test)
  })
  
  subtype_pred <- eventReactive(input$Classprediction, {
    Subtype <- deepcc(data, classifier)
    return(Subtype)
  })
  #visualize the dataframe using renderTable
  output$head <- renderTable({
    df <- data()
    df <- df[, -1]
    return(df)
  })
}
shinyApp(ui = ui, server = server)