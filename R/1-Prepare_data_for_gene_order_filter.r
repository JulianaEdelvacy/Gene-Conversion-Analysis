# Script to prepare data for gene order filter (Updated to use command-line arguments)

# Loading needed libraries. If you already have these libraries installed, you can just comment these lines with #.
install.packages("Biostrings")
install.packages("dplyr")
install.packages("plyr")
install.packages("parallel")
install.packages("psych")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("openxlsx") # Required for saving results

library(Biostrings)
library(dplyr)
library("plyr")
library(parallel)
library(psych)
library(ggplot2)
library(tidyverse)
library(openxlsx)

# -------------------------------------------------------------
# 1. READ DATA DIRECTORY PATH FROM COMMAND-LINE ARGUMENT
# -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- "." # Default to current directory

if (length(args) > 0) {
  DATA_DIR <- args[1]
  cat(paste("Using data directory:", DATA_DIR, "\n"))
} else {
  cat("No data path provided. Using current directory (./). Ensure input files are here.\n")
}

# -------------------------------------------------------------
# 2. FILE PATH DEFINITIONS
# -------------------------------------------------------------

# Use file.path() to correctly construct full file paths for cross-platform compatibility
IMGT_ANNOTATION_FILE <- file.path(DATA_DIR, "2_IMGT-gapped-nt-sequences.txt")
BREPCONVERT_RESULT_FILE <- file.path(DATA_DIR, "Input_test_Horse_Gene_Conversion.csv")
OUTPUT_FILE <- "final_result.csv" # The output file will be written to the script's execution directory

# -------------------------------------------------------------
# 3. MAIN SCRIPT LOGIC
# -------------------------------------------------------------

### Importing IMGT annotation table 
dataset <- read.csv(
  IMGT_ANNOTATION_FILE, 
  sep = "\t", 
  stringsAsFactors = FALSE, 
  header = TRUE
)

#### Using only columns of interest to work with. 
dataset <- dataset %>% select(Sequence.ID, V.D.J.REGION, V.REGION, V.GENE.and.allele)

# Leaving only unique sequences on the dataset
dataset <- dataset[!duplicated(dataset$V.REGION),]

### Counting total rows on the table
TamanhoTotal <- nrow(dataset)

# Counting the number of lines with empty spaces (for report/debug purposes)
datasetvazios <- nrow(dataset[dataset == '',])

# Deleting lines with empty spaces or NA values
dataset <- dataset[!(dataset$V.REGION == "" | is.na(dataset$V.REGION)), ]
dataset <- dataset[!(dataset$V.D.J.REGION == "" | is.na(dataset$V.D.J.REGION)), ]
dataset <- dataset[!(dataset$Sequence.ID == "" | is.na(dataset$Sequence.ID)), ]

# Counting again the number of lines with empty spaces, the expected result is zero.
datasetcont <- nrow(dataset[dataset == '',])

# Reading result file from Brepconvert
resultado1 <- read.csv(
  BREPCONVERT_RESULT_FILE, 
  sep = ",", 
  stringsAsFactors = FALSE, 
  header = TRUE
)

# Making sure that the size of the event has the minimum size of three.
filtro <- nchar(resultado1$seq_event) >= 3
new_resultado1 <- resultado1[filtro, ]

# Adding SeqID values to the table
colnames(dataset) <- c("SeqID", colnames(dataset) [2:ncol(dataset)])
new_resultado1 <- left_join(new_resultado1, dataset, by = "SeqID")

# Joining tables
tabela_juntada <- new_resultado1

# Modifying column V Gene and Allele
tabela_juntada$V.GENE.and.allele <- gsub("Equcab", "", tabela_juntada$V.GENE.and.allele)
tabela_juntada$V.GENE.and.allele <- gsub(" ", "-", tabela_juntada$V.GENE.and.allele, fixed = TRUE, perl = TRUE)
tabela_juntada$V.GENE.and.allele <- gsub("-IG", "IG", tabela_juntada$V.GENE.and.allele, fixed = TRUE, perl = TRUE)

### Keeping in the column only the first annotation and deleting everything after comma 
tabela_juntada$V.GENE.and.allele <- c(sub(",(.*)", "",tabela_juntada$V.GENE.and.allele))

# Removing NA values
tabela_juntada1 <- na.omit(tabela_juntada)

# Saving the results
write.table(tabela_juntada1, OUTPUT_FILE, sep = "\t", row.names = FALSE)
