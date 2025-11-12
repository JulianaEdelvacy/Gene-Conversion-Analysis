# Code to select sequences based on the locus gene order in IMGT.
# The rationale is that if the gene has a direct or opposite orientation, it doesn't select anything that comes after it.

# First, get the functional gene, check if it's present in the direct or opposite dictionary.
# And get the index of i.
# If it's direct, exclude everything after i
# If it's opposite, exclude everything before i
# Get the pseudogene and check if it's present in the dictionary if it is, store its value.

# -------------------------------------------------------------
# 1. LIBRARY LOADING
# -------------------------------------------------------------
library(openxlsx) 

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
# 3. FILE PATH DEFINITIONS
# -------------------------------------------------------------
INPUT_FILE <- file.path(DATA_DIR, "final_result.csv")
OUTPUT_FILE <- "final_result_filtered.csv" # Output written to the script's execution directory

# -------------------------------------------------------------
# 4. HORSE IMGT GENE ORDER DICTIONARY (direto)
# -------------------------------------------------------------
# Dictionary with horse IMGT gene order, you could also load the file.
direto <- c(
  "IGHV4S1" = 0, "IGHV1-41N" = 1, "IGHV3-40N" = 2, "IGHV3-39N" = 3,
  "IGHV4-38N" = 4, "IGHV4-37N" = 5, "IGHV1-84" = 6, "IGHV4-83" = 7,
  "IGHV4-82" = 8, "IGHV4-81" = 9, "IGHV3-80" = 10, "IGHV1-79" = 11,
  "IGHV3-78" = 12, "IGHV2-77" = 13, "IGHV4-33D" = 14, "IGHV7-34D" = 15,
  "IGHV4-35D" = 16, "IGHV3-76" = 17, "IGHV4-75" = 18, "IGHV4-74" = 19,
  "IGHV3-73" = 20, "IGHV4-72" = 21, "IGHV5-71" = 22, "IGHV1-70" = 23,
  "IGHV3-69" = 24, "IGHV4-68" = 25, "IGHV5-67" = 26, "IGHV9-66" = 27,
  "IGHV4-65" = 28, "IGHV3-64" = 29, "IGHV2-63" = 30, "IGHV3-62" = 31,
  "IGHV1-61" = 32, "IGHV3-60" = 33, "IGHV4-59" = 34, "IGHV3-58" = 35,
  "IGHV4-57" = 36, "IGHV3-56" = 37, "IGHV4-55" = 38, "IGHV(II)-54" = 39,
  "IGHV4-53" = 40, "IGHV4-52" = 41, "IGHV4-51" = 42, "IGHV3-50" = 43,
  "IGHV4-49" = 44, "IGHV4-48" = 45, "IGHV(II)-43N" = 46, "IGHV4-42N" = 47,
  "IGHV3-47" = 48, "IGHV3-46" = 49, "IGHV1-45" = 50, "IGHV4-42D" = 51,
  "IGHV(II)-43D" = 52, "IGHV3-44" = 53, "IGHV(II)-43" = 54, "IGHV4-42" = 55,
  "IGHV1-41D" = 56, "IGHV3-40D" = 57, "IGHV3-39D" = 58, "IGHV4-38D" = 59,
  "IGHV4-37D" = 60, "IGHV1-41" = 61, "IGHV3-40" = 62, "IGHV3-39" = 63,
  "IGHV4-38" = 64, "IGHV4-37" = 65, "IGHV3-36" = 66, "IGHV4-35" = 67,
  "IGHV7-34" = 68, "IGHV4-33" = 69, "IGHV3-32" = 70, "IGHV1-31" = 71,
  "IGHV4-30" = 72, "IGHV4-29" = 73, "IGHV4-28" = 74, "IGHV4-20N" = 75,
  "IGHV3-27" = 76, "IGHV3-26" = 77, "IGHV4-25" = 78, "IGHV3-24" = 79,
  "IGHV4-23" = 80, "IGHV4-22" = 81, "IGHV4-20D" =	82, "IGHV4-21" =	83,
  "IGHV4-20" =	84, "IGHV3-19" =	85, "IGHV3-18" =	86, "IGHV4-17" =	87,
  "IGHV3-16" =	88, "IGHV(II)-15" =	89, "IGHV4-14" =	90, "IGHV3-13" =	91,
  "IGHV(II)-12" =	92, "IGHV4-11" =	93, "IGHV4-10" =	94, "IGHV3-9" =	95,
  "IGHV4-8" =	96, "IGHV3-7" =	97, "IGHV4-6" =	98, "IGHV1-5" =	99,
  "IGHV4-4" =	100, "IGHV4-3" =	101, "IGHV4-2" =	102, "IGHV3-1" =	103)

### Import Brepconvert output file
data <- read.csv(INPUT_FILE, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Prepare data
data$gene <- c(gsub("\\*(.*)", "",data$gene)) ### Keep in the column only the first annotation and delete everything after *
data$V.GENE.and.allele <- gsub("\\*(.*)", "", data$V.GENE.and.allele) ### Keep in the column only the first annotation and delete everything after *
pseudogenes_unique <- unique(data$gene) # Create a list with unique copy of each pseudogene
funcional_unique <- unique(data$V.GENE.and.allele) # Create a list with unique copy of each functional gene
dic_direto <- list() # Create empty directory

# Gene Order Filter
result_data <- data.frame() # initialize an empty data frame to store the result

for (i in funcional_unique) {
  # First, get the functional gene, check if it's present in the direct or opposite dictionary.
  if (i %in% names(direto)) {
    # And get the index of i.
    j = direto[[i]]
    # If it's direct, exclude everything after i (by taking everything less than or equal to j)
    dic_direto <- direto[direto <= j]
    # Get the pseudogene and check if it's present in the dictionary; if it is, store its value.
    pseudogenes_to_keep <- intersect(names(dic_direto), pseudogenes_unique)
    
    # Get the functional and pseudogene that are on the same line
    data_to_keep <- data[data$gene %in% pseudogenes_to_keep & data$V.GENE.and.allele == i, ]
    
    # Join to the final result
    result_data <- rbind(result_data, data_to_keep)
  }
}

# Saving the results
write.table(result_data, OUTPUT_FILE, sep = "\t", row.names = FALSE)
cat(paste("Filtering complete. Results saved to:", OUTPUT_FILE, "\n"))
