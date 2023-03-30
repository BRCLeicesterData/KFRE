#####
# KRFE app
# App developed by Naomi Bradbury
# nvm4 at leicester.ac.uk
#####

# Libraries
library(ggplot2)

source("plots.R") # File containing plot code
source("util.R") # File containing utility functions


# Load data
data_sets <- readRDS("data_sets.rds") 
# List of patient data subset by ethnicity
# List order: [1]k (all), [2]sa (south asian), [3]w(white)

results <- readRDS("results.rds")
# Data frame of analysis results

ui <- fluidPage(
  
  fluidRow(
    style = "padding-left:10px; padding-right:10px;",
    HTML("This R Shiny app is designed to accompany the article: Maher", 
         "<em>", "et al", "</em>", "Using the Kidney Failure Risk Equation to predict 
         end-stage kidney disease in CKD patients of South Asian ethnicity: an external 
         validation study (in submission) [link to be added]. For more information please 
         consult the User Guide [link to be added].")
  ),
  
  hr(style = 'border-left: 1px solid grey'), 
  
  fluidRow(

    # Left column
    column(6, align="center",
           column(6,
                  # Cohort selector 
                  selectInput("cohort_left", label = h5("Choose cohort:"), 
                              choices = list("All" = "All", "South Asian" = "South_Asian", 
                                             "White" = "White"), 
                              selected = "South_Asian")
           ),
           
           column(6,
                  # Model selector
                  selectInput("model_left", label = h5("Choose model:"),
                              choices = list("Model 2" = 2, "Model 3" = 3,
                                             "Model 4" = 4, "Model 5" = 5),
                              selected = 5)
           ),
           
           fluidRow(
             column(8,
                    # Axes truncation
                    sliderInput("trunc_left", label = h5("Axes truncation"), min = 10, 
                                max = 100, value = 60, step = 10, post = "%"),
             ),
             
             column(4, 
                    textOutput("warning_left"), # Truncation warning text
             ),
           ),
           
           plotOutput("calibration_left"), # Calibration plot
           plotOutput("distribution_left") # Distribution plot
    ),
    
    # Right column      
    column(6, align="center",
           column(6,
                  # Cohort selector 
                  selectInput("cohort_right", label = h5("Choose cohort:"), 
                              choices = list("All" = "All", "South Asian" = "South_Asian", 
                                             "White" = "White"), 
                              selected = "White")
           ),
           
           column(6,
                  # Model selector
                  selectInput("model_right", label = h5("Choose model:"),
                              choices = list("Model 2" = 2, "Model 3" = 3,
                                             "Model 4" = 4, "Model 5" = 5),
                              selected = 5)
           ),
           
           fluidRow(
             column(8,
                    # Axes truncation
                    sliderInput("trunc_right", label = h5("Axes truncation"), min = 10, 
                                max = 100, value = 60, step = 10, post = "%"),
             ),
             
             column(4,
                    textOutput("warning_right"), # Truncation warning text
             ),
           ),
           
           plotOutput("calibration_right"), # Calibration plot
           plotOutput("distribution_right"), # Distribution plot
           style = 'border-left: 1px solid grey' # Separating line between columns
    ),
    
    hr(style = 'border-left: 1px solid grey'),  
    
    fluidRow(
      style = "padding-left:20px; padding-right:10px;",
      HTML("This study was funded by the National Institute for Health and Care Research (NIHR) 
      Applied Research Collaboration East Midlands (ARC EM) and Leicester NIHR Biomedical 
      Research Centre (BRC). The views expressed are those of the author(s) and not necessarily 
      those of the NIHR or the Department of Health and Social Care.")
    )
  )
)

server <- function(input, output) {
  
  # Truncation warning text left - function in util.R
  output$warning_left <- renderText({
    warning_text(input$trunc_left)
  })
  
  # Truncation warning text right - function in util.R
  output$warning_right <- renderText({
    warning_text(input$trunc_right)
  })
  
  # Calibration plot left - function in plot.R
  output$calibration_left <- renderPlot({
    
    # Select correct subset of data
    res <- results[(results$cohort == input$cohort_left) & (results$model == input$model_left),] 
    
    # Call plot function (in plots.R)
    calibration_plot(res, input$cohort_left, input$trunc_left) 

  })
  
  
  # Calibration plot right - function in plot.R
  output$calibration_right <- renderPlot({
    
    # Select correct subset data
    res <- results[(results$cohort == input$cohort_right) & (results$model == input$model_right),] 
    
    # Call plot function (in plots.R)
    calibration_plot(res, input$cohort_right, input$trunc_right) 
    
  })
  
  # Distribution plot left - function in plot.R
  output$distribution_left <- renderPlot({
    
    # Select correct cohort 
    if (input$cohort_left == "All") {data <- data_sets[[1]] 
    } else if (input$cohort_left == "South_Asian") {data <- data_sets[[2]]
    } else if (input$cohort_left == "White") {data <- data_sets[[3]]
    }
    
    # Select correct kfre column based on model, then add to df
    if (input$model_left == 2) {select_kfre <- data$kfre2
    } else if (input$model_left == 3) {select_kfre <- data$kfre3
    } else if (input$model_left == 4) {select_kfre <- data$kfre4
    } else if (input$model_left == 5) {select_kfre <- data$kfre5
    }

    data <- cbind(data, select_kfre)
    
    # Call plot function (in plots.R)
    distribution_plot(data, input$trunc_left) 
    
  })
  
  # Distribution plot right - function in plot.R
  output$distribution_right <- renderPlot({
    
    # Select correct cohort 
    if (input$cohort_right == "All") {data <- data_sets[[1]] 
    } else if (input$cohort_right == "South_Asian") {data <- data_sets[[2]]
    } else if (input$cohort_right == "White") {data <- data_sets[[3]]
    }
    
    # Select correct kfre column based on model, then add to df
    if (input$model_right == 2) {select_kfre <- data$kfre2
    } else if (input$model_right == 3) {select_kfre <- data$kfre3
    } else if (input$model_right == 4) {select_kfre <- data$kfre4
    } else if (input$model_right == 5) {select_kfre <- data$kfre5
    }
    
    data <- cbind(data, select_kfre)
    
    # Call plot function (in plots.R)
    distribution_plot(data, input$trunc_right) 
    
  })
  
}

shinyApp(ui = ui, server = server)