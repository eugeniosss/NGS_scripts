library(plotly)
library(tidyverse)

# Define a function to create the plot
create_plot <- function(df) {
  fig <- plot_ly(data = df, x = ~Number_of_reads, y = ~Mean_Coverage, color = ~RUN, type = 'scatter', mode = 'lines') %>%
    layout(
      xaxis = list(
        title = "Number of raw sequenced reads"),
      yaxis = list(
        title = "Expected mean mean coverage")
    )
  return(fig)
}

# Main function to execute the script
main <- function(input_file, output_file) {
  # Read the input data
  df <- read.csv(input_file, sep = "\t", header = TRUE, fill = TRUE)
  
  # Create the plot
  fig <- create_plot(df)
  
  # Save the plot as HTML widget
  htmlwidgets::saveWidget(as_widget(fig), output_file)
}

# Call the main function with provided arguments
main(commandArgs(trailingOnly = TRUE)[1], commandArgs(trailingOnly = TRUE)[2])
