# Result analysis from Non Beta motifs search tool

# Loading needed libraries, if you already have these libraries you can just comment these lines with #.
install.packages("data.table")
install.packages("ggplot2")
install.packages("data.table")
install.packages("dplyr")
install.packages("ggnewscale")

library(data.table)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggnewscale)

#Couting and plotting the repetitions of Direct Repeats motifs.
# Loading Direct repeat result file
DR <- read.csv("Direct_Repeats_motifs.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Load your gene conversion results from Brepconvert
data1 <- read.csv("final_result_filtered.cs.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
df_filtered_selecionado <- as.data.table(data1)

# Use data.table to filter and clean the Sequence column
DR_filtered <- as.data.table(DR)
DR_filtered[, Sequence := trimws(sub(".*seq\\s*", "", Sequence))]

# Initialize results list with a maximum size
results <- vector("list", nrow(DR_filtered) * nrow(df_filtered_selecionado))
result_index <- 1

# Iterate over the IR sequences
for (ir_seq in DR_filtered$Sequence) {
  # Find all rows in df_combined_case where IR sequence is substring in V.REGION
  matches_dt <- df_filtered_selecionado[grepl(ir_seq, V.REGION)]
  
  if (nrow(matches_dt) > 0) {
    # For each matched row, find start and end positions of IR sequence inside the V.REGION sequence
    start_positions <- regexpr(ir_seq, matches_dt$V.REGION)
    
    # Filter out non-matching sequences
    valid_matches <- start_positions != -1
    
    if (any(valid_matches)) {
      # Calculate start and end positions for valid matches
      start_pos <- start_positions[valid_matches]
      end_pos <- start_pos + nchar(ir_seq) - 1
      
      # Append to results list
      for (row_index in which(valid_matches)) {
        results[[result_index]] <- list(
          SeqID = matches_dt$Sequence.ID[row_index],
          Sequence = ir_seq,
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
# Remove single column
results_dt_DR <- rbindlist(results)
results_dt_DR_unique <- unique(results_dt_DR)
results_dt_DR_unique <- results_dt_DR_unique[!is.na(Sequence) & Sequence != "", ]
# Save end_counts as a CSV file
write.csv(results_dt_DR_unique, "results_dt_DR_unique.csv", row.names = FALSE)


#Ploting the results
DR_unique <- read.csv("results_dt_DR_unique_SeqID.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
# Convert results_dt_IR to data.frame
DR_unique <- as.data.frame(results_dt_DR_unique)
# Prepare both datasets
ir_data <- DR_unique %>%
  mutate(source = "Direct Repeats and Slipped Motifs")
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

ggplot(combined_data, aes(x = start, y = end, size = frequency, color = source)) +
  # Points with size representing frequency and color representing source
  geom_point(alpha = 0.8) +
  
  # Manual color scale for the sources
  scale_color_manual(values = c("Gene Conversion" = "lightblue", 
                                "Direct Repeats and Slipped Motifs" = "red"),
                     name = "Event Type") +
  
  scale_size_continuous(range = c(2, 10),
                        breaks = seq(1, 15001, by = 2500),
                        name = "Frequency") +  # Adjust size range as needed
  
  labs(x = "Initial Position",
       y = "Final Position") +
  
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

