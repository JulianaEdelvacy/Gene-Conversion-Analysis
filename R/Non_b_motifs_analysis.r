# Result analysis from Non Beta motifs search tool
# The tool is available at this link: https://nonb-abcc.ncifcrf.gov/apps/nBMST/default/

# -------------------------------------------------------------
# 1. LIBRARY LOADING
# -------------------------------------------------------------

# Loading needed libraries. If you already have these libraries installed, you can just comment these lines with #.
install.packages("data.table")
# install.packages("ggplot2")
# install.packages("dplyr")
install.packages("ggnewscale")

library(data.table)
library(ggplot2)
library(dplyr)
library(ggnewscale)

# -------------------------------------------------------------
# 2. READ DATA DIRECTORY PATH FROM COMMAND-LINE ARGUMENT
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
# 3. FILE PATH DEFINITION
# -------------------------------------------------------------
DR_FILE <- file.path(DATA_DIR, "Direct_Repeats_motifs.tsv")
DATA1_FILE <- file.path(DATA_DIR, "final_result_filtered.csv")
# Caminhos para arquivos intermediários e finais (dentro da pasta de dados)
RESULTS_CSV <- file.path(DATA_DIR, "results_dt_DR_unique.csv")
PLOT_OUTPUT_FILE <- file.path(DATA_DIR, "Motifs_Scatter_Plot.png")

# -------------------------------------------------------------
# 4. MAIN ANALYSIS LOGIC
# -------------------------------------------------------------

# Counting and plotting the repetitions of Direct Repeats motifs.
# Loading Direct repeat result file - USANDO O DR_FILE CORRETO
DR <- read.csv(DR_FILE, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Load your gene conversion results from Brepconvert - USANDO O DATA1_FILE CORRETO
data1 <- read.csv(DATA1_FILE, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
df_filtered_selecionado <- as.data.table(data1)

# Use data.table to filter and clean the Sequence column
DR_filtered <- as.data.table(DR)
DR_filtered[, Sequence := trimws(sub(".*seq\\s*", "", Sequence))]

# Initialize results list with a maximum size
results <- vector("list", nrow(DR_filtered) * nrow(df_filtered_selecionado))
result_index <- 1

# Iterate over the DR sequences
for (dr_seq in DR_filtered$Sequence) {
  # Find all rows in df_filtered_selecionado where DR sequence is substring in V.REGION
  matches_dt <- df_filtered_selecionado[grepl(dr_seq, V.REGION)]
  
  if (nrow(matches_dt) > 0) {
    # For each matched row, find start and end positions of DR sequence inside the V.REGION sequence
    start_positions <- regexpr(dr_seq, matches_dt$V.REGION)
    
    # Filter out non-matching sequences
    valid_matches <- start_positions != -1
    
    if (any(valid_matches)) {
      # Calculate start and end positions for valid matches
      start_pos <- start_positions[valid_matches]
      end_pos <- start_pos + nchar(dr_seq) - 1
      
      # Append to results list
      for (row_index in which(valid_matches)) {
        results[[result_index]] <- list(
          SeqID = matches_dt$Sequence.ID[row_index],
          Sequence = dr_seq,
          Start = start_pos[row_index],
          End = end_pos[row_index]
        )
        result_index <- result_index + 1
      }
    }
  }
}

# Trim the results list to the actual size
results <- results[1:(result_index - 1)]
# Combine results into a data.table
results_dt_DR <- rbindlist(results)
results_dt_DR_unique <- unique(results_dt_DR)
results_dt_DR_unique <- results_dt_DR_unique[!is.na(Sequence) & Sequence != "", ]

# Save end_counts as a CSV file - USANDO O CAMINHO CORRETO
write.csv(results_dt_DR_unique, RESULTS_CSV, row.names = FALSE)


# Ploting the results
# Carrega o CSV que acabou de ser salvo (do caminho correto)
DR_unique <- read.csv(RESULTS_CSV, sep = ",", stringsAsFactors = FALSE, header = TRUE)

# Convert results_dt_DR_unique to data.frame
DR_unique <- as.data.frame(DR_unique)

# Prepare both datasets
ir_data <- DR_unique %>%
  mutate(source = "Direct Repeats and Slipped Motifs")

# A variável 'data1' já foi carregada no início
gc_data <- data1 %>%
  mutate(source = "Gene Conversion")

# Calculate frequency for Gene Conversion events
gc_data_freq <- gc_data %>%
  group_by(start, end) %>%
  summarise(frequency = n(), .groups = 'drop') %>% 
  mutate(source = "Gene Conversion")

# Calculate frequency for Direct Repeats and Slipped Motifs
ir_data_freq <- ir_data %>%
  group_by(Start, End) %>%
  summarise(frequency = n(), .groups = 'drop') %>%
  mutate(source = "Direct Repeats and Slipped Motifs") %>%
  # rename columns to match gc_data_freq
  rename(start = Start, end = End)

# Combine for axis limits calculation
all_starts <- c(gc_data_freq$start, ir_data_freq$start)
all_ends <- c(gc_data_freq$end, ir_data_freq$end)

# Combine the data for plotting
combined_data <- rbind(
  mutate(gc_data_freq, source = "Gene Conversion"),
  mutate(ir_data_freq, source = "Direct Repeats and Slipped Motifs")
)

plot <- ggplot(combined_data, aes(x = start, y = end, size = frequency, color = source)) +
  # Points with size representing frequency and color representing source
  geom_point(alpha = 0.8) +
  
  # Manual color scale for the sources
  scale_color_manual(values = c("Gene Conversion" = "lightblue", 
                                "Direct Repeats and Slipped Motifs" = "red"),
                     name = "Event Type") +
  
  scale_size_continuous(range = c(2, 10),
                        breaks = seq(1, 15001, by = 2500),
                        name = "Frequency") + # Adjust size range as needed
  
  labs(x = "Initial Position",
       y = "Final Position",
       title = "Analysis of Gene Conversion vs. Direct Repeats Motifs") + # Adicionei um título
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.80, 0.3),
    legend.box.background = element_rect(color = "black", linewidth = 0.5),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.text = element_text(size = 9)
  ) +
  
  scale_x_continuous(
    limits = c(0, max(all_starts, na.rm = TRUE) + 25),
    breaks = seq(0, max(all_starts, na.rm = TRUE), by = 25)
  ) +
  scale_y_continuous(
    limits = c(0, max(all_ends, na.rm = TRUE) + 25),
    breaks = seq(0, max(all_ends, na.rm = TRUE), by = 25)
  )

# -------------------------------------------------------------
# 5. SAVE PLOT TO PNG FILE
# -------------------------------------------------------------

# Salva o gráfico gerado ('plot') como PNG no caminho correto
ggsave(PLOT_OUTPUT_FILE, plot = plot, width = 8, height = 6, units = "in", dpi = 300)

cat(paste("Plot saved successfully to:", PLOT_OUTPUT_FILE, "\n"))
