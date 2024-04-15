# Load required libraries
library(plotly)
library(tidyverse)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_html <- args[2]

# Read input data
df <- read.csv(input_file, sep = "\t", header = TRUE, fill = TRUE)

# Perform data transformation
data2 <- df %>%
  pivot_longer(
    cols = c("EXPECTED_DISTINCT", "distinct_reads"),
    names_to = "condition",
    values_to = "Value"
  )

# Create the Plotly figure
fig <- plot_ly(data = data2, x = ~TOTAL_READS, y = ~Value, color = ~RUN, type = 'scatter', mode = 'lines',
               linetype = ~condition) %>%
  layout(
    xaxis = list(
      title = "Number of raw reads",
      range = c(0, 10000000)
    ),
    yaxis = list(
      title = "Number of reads"
    ),
    legend = list(
      tracegroupgap = 100
    )
  )

# Save the plot as HTML widget
htmlwidgets::saveWidget(as_widget(fig), output_html)

# Modify the HTML file
system2(command = "sed", args = c("-i", "'s/EXPECTED_DISTINCT/: reads projected to be unique/g'", output_html))
system2(command = "sed", args = c("-i", "'s/distinct_reads/: actually unique reads/g'", output_html))

