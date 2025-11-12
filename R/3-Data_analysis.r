# Code for data analysis from gene order filter

# -------------------------------------------------------------
# 1. LIBRARY LOADING
# -------------------------------------------------------------

# Loading needed libraries. If you already have these libraries installed, you can just comment these lines with #.
# install.packages("Biostrings")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("parallel")
# install.packages("psych")
# install.packages("ggplot2")
# install.packages("tidyverse")

#library(Biostrings)
library(dplyr)
library("plyr")
library(parallel)
library(psych)
library(ggplot2)
library(tidyverse)

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
# Importing data filtered from Gene_Order_Filter.r
INPUT_FILE <- file.path(DATA_DIR, "final_result_filtered.csv")

# -------------------------------------------------------------
# 4. MAIN ANALYSIS LOGIC
# -------------------------------------------------------------

# Importing data filtered from Gene_Order_Filter.r
data1 <- read.csv(INPUT_FILE, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Counting the repetitions from each column
cat("\n--- Repetitions per Column ---\n")
cat("\n--- distance_to_AID_motif ---\n")
table(data1$distance_to_AID_motif)
cat("\n--- Functional gene ---\n")
table(data1$V.GENE.and.allele)
cat("\n--- edit_distance ---\n")
table(data1$edit_distance)
cat("\n--- nearest_AID_motif ---\n")
table(data1$nearest_AID_motif)

# Seeing the most repeated pseudogenes
cat("\n--- Most Repeated Pseudogenes ---\n")
names_P <- data1$gene
repetitions_P <- table(unlist(names_P))
repetitions_P <- sort(repetitions_P, decreasing = TRUE)

# Printing the names and repetitions
for (name in names(repetitions_P)) {
  print(paste(name, repetitions_P[name]))
}

# Seeing the most repeated functional genes
cat("\n--- Most Repeated Functional Genes ---\n")
names_F <- data1$V.GENE.and.allele
repetitions_F <- table(unlist(names_F))
repetitions_F <- sort(repetitions_F, decreasing = TRUE)

# Printing the names and repetitions
for (name in names(repetitions_F)) {
  print(paste(name, repetitions_F[name]))
}

# Counting seq_event repetitions
# Calculating the size (length) of seq_event
seq_event_size = nchar(data1$seq_event)

# Counting Seq_event size repetitions
cat("\n--- Seq_event Size Repetitions ---\n")
repetitions_seq_event <- table(seq_event_size)
print(repetitions_seq_event)
