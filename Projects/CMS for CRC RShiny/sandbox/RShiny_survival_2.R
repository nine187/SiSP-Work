library(shiny)
library(survival)
library(ggplot2)
library(lubridate)

server <- function(input, output, session) {
  # Cohort visualization
  cohort_vis <- ezfun::set_ccf_palette("contrast")
  
  #use cohort visualization 
  
  Surv(lung$time, lung$status)[1:10]
  
  s1 <- survfit(Surv(time, status) ~ 1, data = lung)
  str(s1)
  
  output$surv_plot <- renderPlot({
    ggsurvfit::survfit2(Surv(time, status) ~ 1, data = lung) %>% 
      ggsurvfit() +
      labs(
        x = "Days",
        y = "Overall survival probability"
      )
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
  #plotOutput(outputId = "surv_plot")
)

shinyApp(ui = ui, server = server)