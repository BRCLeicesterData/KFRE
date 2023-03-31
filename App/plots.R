# Calibration plot

# INPUTS:
# results: df containing subset of results.rds to be plotted
# cohort: cohort input selection
# model: model input selection
# max: max axis limit from axis truncation slider as a percentage
#
# OUTPUT: calibration plot as ggplot

calibration_plot <- function(results, cohort, model, maximum) {
  
  # Make correct title
  if (cohort == "All") {t1 <- "Overall - "
  } else if (cohort == "South_Asian") {t1 <- "South Asian cohort - "
  } else if (cohort == "White") {t1 <- "White cohort - "}
  
  if (model == 2) {t2 <- "Model 2"
  } else if (model == 3) {t2 <- "Model 3" 
  } else if (model == 4) {t2 <- "Model 4" 
  } else if (model == 5) {t2 <- "Model 5"
  }
  
  title <- paste(t1, t2, sep = "")
  
  # Convert max to a proportion
  max <- maximum/100
  
  plot <- ggplot(results, aes(pred.mean, risk1)) +
    geom_errorbar(aes(ymin = risk.l, ymax = risk.u, width=0.02), colour="blue") +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") + # reference line
    labs(title = title, y = "Observed risk", x = "Predicted risk") +
    scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, max)) +  # truncate scale
    scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, max)) +
    theme(axis.text = element_text(face="bold", size=12),
          axis.title = element_text(face="bold", size=14))
  
  return(plot)
}

# Distribution plot (using histogram) 

# INPUTS:
# data: df containing correct cohort from data_sets.rds to be plotted with an additional column
# named 'select_kfre' that is the kfre value associated with the model to be plotted
# max: max axis limit from axis truncation slider as a percentage
#
# OUTPUT: distribution histogram as ggplot

distribution_plot <- function(data, maximum) {
  
  # Convert max to a proportion
  max <- maximum/100
  
  plot <- ggplot(data, aes(x = select_kfre)) +
    geom_histogram(fill = "coral1", bins = 200) +
    scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, max)) + # truncate scale
    xlab("Predicted risk") +
    ylab("Count") +
    theme_minimal() +
    scale_y_continuous() +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size=12),
          axis.title = element_text(size=14))

  return(plot)
}
