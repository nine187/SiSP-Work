library(shiny)
library(survival)
library(ggplot2)
library(ggsurvfit)
library(lubridate)

server <- function(input, output, session) {
  # Cohort visualization
  cohort_vis <- ezfun::set_ccf_palette("contrast")
  
  output$surv_plot <- renderPlot({
    req(input$cohort_file)
    
    data <- read.csv(input$cohort_file$datapath)
    
    # Assuming your data has 'time' and 'status' columns
    fit <- survfit(Surv(time, status) ~ 1, data = data)
    
    # Plotting survival curve
    ggsurvfit::survfit2(Surv(time,status) ~ 1, data = data)%>% 
      ggsurvfit() +
      labs(
        x = "Days",
        y = "Overall survival probability")
  })
}

ui <- fluidPage(
  fileInput(inputId = "cohort_file", 
            label = "Choose cohort file (.CSV/ .TXT file)",
            multiple = FALSE,
            accept = c("text/csv",
                       "text/comma-separated-values,text/plain",
                       ".csv")),
  actionButton(inputId = "Cohort_vis",
               label = "Visualization"),
  p("Press the button to visualize the cohort data"),
  plotOutput(outputId = "surv_plot")
)

shinyApp(ui = ui, server = server)
