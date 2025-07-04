library(plotly)
library(tidyverse)

# Function to read data from a file
read_data <- function(file_path) {
  read.csv(file_path, sep = "\t", header = TRUE, fill = TRUE)
}

# Function to create the plot
create_plot <- function(input_file, output_file, log_file) {
  cat("Input file:", input_file, "\n")
  cat("Output file:", output_file, "\n")
  cat("Log file:", log_file, "\n")
  
  # Read data
  df <- read_data(input_file)
  
  # Pivot data
  data2 <- df %>%
    pivot_longer(
      cols = c("EXPECTED_DISTINCT", "distinct_reads"),
      names_to = "condition",
      values_to = "Value"
    )
  
  # Create plot
  fig <- plot_ly(data = data2, x = ~TOTAL_READS, y = ~Value, linetype = ~condition,
                 marker = list(size = 2)) %>% 
    add_trace(data = subset(data2, condition == "EXPECTED_DISTINCT"), y = ~Value, mode = 'lines', 
              opacity = 0.25, name = "Reads projected to be unique", line = list(color = 'blue')) %>%
    add_trace(data = subset(data2, condition == "distinct_reads"), y = ~Value, mode = 'lines', 
              opacity = 1, name = "Unique", line = list(color = 'blue')) %>%
    layout(
      xaxis = list(title = "Number of raw reads", range = c(0, 100000000)),
      yaxis = list(title = "Number of reads")
    )
  
  # Save plot to HTML file
  htmlwidgets::saveWidget(as_widget(fig), output_file)
}

# Call the create_plot function with command-line arguments
args <- commandArgs(trailingOnly = TRUE)
create_plot(args[1], args[2], args[3])

