# Code to create a Chord Diagram related to the frequency of gene conversion events.
# This code was adapted for Horse IGHV functional and pseudogenes. You should change the gene names as necessary for your dataset.

# -------------------------------------------------------------
# 1. LIBRARY LOADING
# -------------------------------------------------------------

# Load necessary libraries. If you already have these libraries installed, you can just comment with #.
install.packages("circlize")
install.packages("dplyr")

library(circlize)
library(dplyr)

# -------------------------------------------------------------
# 2. READ DATA DIRECTORY PATH FROM COMMAND-LINE ARGUMENT
# -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- "." # Default to current directory

if (length(args) > 0) {
  DATA_DIR <- args[1]
  cat(paste("Using data directory:", DATA_DIR, "\n"))
} else {
  cat("No data path provided. Using current directory (./). Ensure the input file is here.\n")
}

# -------------------------------------------------------------
# 3. FILE PATH DEFINITION
# -------------------------------------------------------------
INPUT_FILE <- file.path(DATA_DIR, "final_result_filtered.csv")

# -------------------------------------------------------------
# 4. MAIN LOGIC
# -------------------------------------------------------------

data <- read.csv(INPUT_FILE, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Select only needed columns (Note: 'df_combined_circus' is used in the original logic 
# and is assumed to be defined elsewhere or is a placeholder for 'data' preprocessing.)
data <- df_combined_circus %>% select(gene, V.GENE.and.allele) 

# Create a frequency table with the combination of the columns gene and V.GENE.and.allele
freq_table <- table(data$gene, data$V.GENE.and.allele)

# Defining gene names and colors, note that you have to change the IGHV gene names depending on the genes you have on you dataset.
genes_cores <- c("IGHV(II)-43" = "blue", "IGHV1-31" = "red", "IGHV1-79" = "green", "IGHV3-24" = "orange",
                 "IGHV3-26" = "yellow", "IGHV3-32" = "purple", "IGHV3-36" = "cyan", "IGHV3-39" = "magenta",
                 "IGHV3-40" = "brown", "IGHV4-20N" = "darkgreen", "IGHV4-23" = "darkblue", "IGHV4-25" = "darkred",
                 "IGHV4-28" = "darkcyan", "IGHV4-30" = "darkmagenta", "IGHV4-33" = "darkgray", "IGHV4-35" = "darkorange",
                 "IGHV4-38" = "#330033", "IGHV4-81" = "#99CC00", "IGHV4-83" = "#FF3300", "IGHV7-34" = "#006600",
                 "IGHV1-41" = "#00CCFF", "IGHV1-5" = "#FFCC00", "IGHV1-70" = "#FF6600", "IGHV3-78" = "#CC00FF",
                 "IGHV4-21" = "#FF9933", "IGHV4-22" = "#009999", "IGHV4-29" = "#9900CC",
                 "IGHV4-29" = "#336699", "IGHV4-37" = "#FF9966", "IGHV4-65" = "#66CCCC",
                 "IGHV4-65" = "#6699FF", "IGHV4-82" = "#993366", "IGHV4S1" = "#99FF99",
                 "IGHV4-53" = "blue", "IGHV3-46" = "red")

# Create the Chord Diagram
chordDiagram(freq_table, grid.col = genes_cores,  annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), link.sort = TRUE, link.decreasing = TRUE)

# Adding gene names and adjusting font size
circos.track(track.index = 1, panel.fun = function(x, y) {
  gene_name <- CELL_META$sector.index
  # Defining the size
  font_size <- 0.8
  
  # Increasing the font size for a few genes of interest
  if (gene_name %in% c("IGHV4-21", "IGHV4-22", "IGHV4-29", "IGHV4-35", "IGHV4-53")) {
    font_size <- 2.0 
  }
  
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], gene_name, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = font_size)
}, bg.border = NA)
