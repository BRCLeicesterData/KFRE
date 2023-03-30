warning_text <- function(trunc) {
  ifelse(trunc < 60,
         "Warning: some data points may not be displayed when the axes are truncated below 60% ",
         '')
}