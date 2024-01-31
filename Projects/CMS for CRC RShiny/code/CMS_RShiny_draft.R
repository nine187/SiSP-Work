library(shiny)
library(shinydashboard)
#library(gridExtra) #arrange RShiny output #ask P'Arm about this
library(networkD3) #data visualization https://stackoverflow.com/questions/48024037/printing-a-sankey-diagram-in-shiny


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ---- CMS for CRC Web Application for SiSP
  selectInput("dataset", label = "Dataset", choices = c("CMS", "CRIS", "DeepCC", "iCMS", "imCMS")),
  verbatimTextOutput("summary"),
  tableOutput("table"),
  dashboardHeader(title = h4("CRC CMS Web Application (V.0.1)", align = "center"),
                  dropdownMenu(messageItem("New User", "Can I get some help?",
                              time = "Today"))),
  fileInput("upload", "Upload a .csv file"),
  dashboardSidebar(   #textOutput("greeting"),
                      sliderInput(inputId = "Freq",
                                     label = "Number of bins:",
                                     min = 1,
                                     max = 50,
                                     value = 30)),
  dashboardBody(
      fluidRow(  
      textOutput("greeting"),
      splitLayout(cellWidths = c("33%", "33%", "33%"), plotOutput("examplePlot"), plotOutput("examplePlot2"), plotOutput("examplePlot3")),
      #plotOutput("examplePlot3"),
      sankeyNetworkOutput("sankey"),
  
  ))
)
# Define server logic required to draw a histogram ----
server <- function(input, output) {
  #greetings
  output$greeting <- renderText(paste0("The total number of CMS is  ", sum(CMS$Freq), "."))
  
  #plot
  #CRIS
  output$examplePlot <- renderPlot(ggplot(CRIS, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 0.5, color = "black") +
    geom_text(aes(label = Freq),
              position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(title = "CRIS Data Distribution"))
  
  #CMS
  output$examplePlot2 <- renderPlot(ggplot(CMS, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 0.5, color = "black") +
    geom_text(aes(label = Freq),
              position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(title = "CMS Data Distribution"))
  
  #DeepCC
  output$examplePlot3 <- renderPlot(ggplot(DeepCC, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 0.5, color = "black") +
    geom_text(aes(label = Freq),
              position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(title = "DeepCC Data Distribution"))
  
  #sankey
  output$sankey <- renderSankeyNetwork({
    sankeyNetwork(Links = CRCsubtypes$links, 
                        Nodes = CRCsubtypes$nodes, 
                        Source = 'source.pos', Target = 'target.pos', Value = 'value', 
                        NodeID = 'name', fontSize = 14, 
                        nodeWidth = 20, iterations = 0,
                        units = 'TWh',
                        NodeGroup = "group", 
                        fontFamily = "arial")
  })}
  
  #use grid.arrange to organize the graphs into
  
  #ptlist <- list(examplePlot,examplePlot2,examplePlot3)
  #wtlist <- c(input$)
  #grid.arrange(grobs=ptlist)
  #https://stackoverflow.com/questions/34384907/how-can-put-multiple-plots-side-by-side-in-shiny-r
  
# Create Shiny app ----
shinyApp(ui = ui, server = server)