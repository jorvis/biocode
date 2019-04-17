#!/usr/bin/env Rscript

library(ggplot2)
library(readr)

args <- commandArgs(TRUE)
file_path <- args[1]

df <- read_tsv(file_path,
  # First row is not the name of columns
  col_names=FALSE)

bar_plot <- ggplot(df) +
  geom_col(
    # Aesthetic mappings describe how variables
    # in the data are mapped to visual properties
    mapping=aes(
        x=X3, # GO Labels
        y=X4  # counts
      )
    ) +
  # Subplot into the 3 categories of X1
  facet_wrap(~X1,
    # X1's categories are mutually exclusive
    scales = "free") +
  # Flip axes so X3 labels are horizontal
  coord_flip() +
  # Remove axis titles
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank())
    
ggsave("plot.png", plot=bar_plot, width=15, height=5)

# ggsave will sometimes auto-generate Rplots.pdf when being run on the command line
# https://stackoverflow.com/questions/17348359/how-to-stop-r-from-creating-empty-rplots-pdf-file-when-using-ggsave-and-rscript
if(file.exists("Rplots.pdf")) {
  # Silently delete file if it exists
  invisible(file.remove("Rplots.pdf"))
}
