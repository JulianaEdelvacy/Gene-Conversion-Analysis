#This code aims to see if gene conversion positions overlap with AID hotspots.

library(ggplot2)
library(dplyr)
library(reshape2)

# Import Gene Conversion Data
# You can obtain this result from previous analysis (R files 1-3).
data1 <- read.csv("final_result_filtered_834.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
data2 <- read.csv("final_result_filtered_839.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
data3 <- read.csv("final_result_filtered_843.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
data4 <- read.csv("final_result_filtered_844.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
data5 <- read.csv("final_result_filtered_1N.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
data6 <- read.csv("final_result_filtered_2N.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
data7 <- read.csv("final_result_filtered_3N.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
data8 <- read.csv("final_result_filtered_4N.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

df_combined_case <- rbind(data1,data2,data3, data4,
                          data5, data6, data7, data8)

# see if gene conversion positions overlap with AID hotspots

library(dplyr)
library(data.table)
library(stringr)
library(tidyr)

# 1. Reading, Combining, and Filtering AID Hotspot Files
# ==============================================================================

# Creates a list of dataframes (hotspots)
# You can obtain these files from the IMGT Hight V Quest Annotation Files
aid_hotspots_list <- list(
  read.csv("1N_10_V-REGION-mutation-hotspots.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE),
  read.csv("2N_10_V-REGION-mutation-hotspots.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE),
  read.csv("3N_10_V-REGION-mutation-hotspots.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE),
  read.csv("4N_10_V-REGION-mutation-hotspots.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE),
  read.csv("834_10_V-REGION-mutation-hotspots.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE),
  read.csv("839_10_V-REGION-mutation-hotspots.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE),
  read.csv("843_10_V-REGION-mutation-hotspots.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE),
  read.csv("844_10_V-REGION-mutation-hotspots.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
)

# Combine all dataframes
df_hotspots_all <- rbindlist(aid_hotspots_list, fill = TRUE) %>%
  as_tibble()

# Filter by functionality and ensure the SeqID name
df_hotspots_filtered <- df_hotspots_all %>%
  filter(V.DOMAIN.Functionality == "productive") %>%
  # Rename Sequence.ID 
  dplyr::rename(SeqID = `Sequence.ID`) 


# 2. Definition of the Position Extraction Function 
# ==============================================================================

# Adaptation of the logic to extract position ranges (e.g., “103-106”)
# from the format: “NNNN,POS_START-POS_END(...)|...”
extract_hotspot_positions <- function(text_vector) {
  sapply(text_vector, function(x) {
    if (is.na(x) || x == "") return(NA)
    
    # Pattern: Captures the group of numbers separated by hyphens after the comma.
    # Match example: “agca,103-106(CDR1)” -> Captures “103-106”.
    matches <- str_extract_all(x, pattern = ",(\\d+-\\d+)\\(")
    
    # Extracts the captured text and cleans up punctuation characters (commas and parentheses)
    positions <- unlist(matches) %>%
      str_remove_all(pattern = ",|\\(")
    
    if (length(positions) == 0) return(NA)
    
    # Join all positions found into a single string separated by semicolons
    paste(positions, collapse = ";")
    
  }, USE.NAMES = FALSE)
}

#3. Function Application and Final Preparation
# ==============================================================================

df_hotspots_processed <- df_hotspots_filtered %>%
  mutate(
    # Apply the function to extract the position intervals for Colunm 1
    AID_Hotspot1_Pos = extract_hotspot_positions(`X.a.g.g.c.t..a.t.`),
    
    # Apply the function to extract the position intervals for Colunm 2
    AID_Hotspot2_Pos = extract_hotspot_positions(`X.a.t..a.g.c.c.t.`)
  ) %>%
  # Keep only the columns that are essential for the next step.
  select(
    SeqID,
    V.DOMAIN.Functionality,
    AID_Hotspot1_Raw = `X.a.g.g.c.t..a.t.`, # Keep the original column for reference
    AID_Hotspot2_Raw = `X.a.t..a.g.c.c.t.`, # Keep the original column for reference
    AID_Hotspot1_Pos,
    AID_Hotspot2_Pos
  )

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# ==============================================================================
#1. Preparation of AID Hotspot Data
# ==============================================================================
df_hotspots_aid_long <- df_hotspots_processed %>%
  select(SeqID, AID_Hotspot1_Pos, AID_Hotspot2_Pos) %>%
  
  pivot_longer(
    cols = c(AID_Hotspot1_Pos, AID_Hotspot2_Pos),
    names_to = "hotspot_motif",
    values_to = "hotspot_intervals",
    values_drop_na = TRUE 
  ) %>%
  
  mutate(hotspot_motif = str_remove(hotspot_motif, "_Pos")) %>%
  
  separate_rows(hotspot_intervals, sep = ";") %>%
  
  filter(hotspot_intervals != "") %>%
  
  separate(
    col = hotspot_intervals,
    into = c("start_aid", "end_aid"),
    sep = "-",
    convert = TRUE,
    remove = TRUE
  )

# ==============================================================================
#2. Joining gene conversion and AID hotspot tables
# ==============================================================================
df_overlap_aid <- df_combined_case %>%
  mutate(tam_CG = end - start + 1) %>%
  select(SeqID, 
         start, 
         end, 
         tam_CG) %>%
  left_join(df_hotspots_aid_long, by = "SeqID", relationship = "many-to-many")

# ==============================================================================
#3. Overlap Calculation and Classification
# ==============================================================================

df_aid_analysis_final <- df_overlap_aid %>%
  mutate(
    overlap_start = pmax(start, start_aid),
    #overlap can only go as far as the “shortest” interval allows.
    overlap_end = pmin(end, end_aid),
    overlap_length = pmax(0, overlap_end - overlap_start + 1),
    
    # 2. Verification 
    is_contained_in_aid = (start >= start_aid) & (end <= end_aid)
  ) %>% 
  mutate(
    overlap_type_aid = case_when(
      #1: Remove those without a hotspot (NAs)
      is.na(start_aid) ~ "No Hotspot in Sequence",
      
      # 2: If the Boolean variable is TRUE (The cluster is contained)
      is_contained_in_aid ~ "Full Overlap",
      
      # 3rd: If it is not contained, but close (Partial)
      overlap_length > 0 ~ "Partial Overlap",
      
      # 4: If the Boolean variable is FALSE and there is no overlap
      !is_contained_in_aid & overlap_length == 0 ~ "No Overlap",
      
      TRUE ~ "Others"
    )
  )

write.csv(df_aid_analysis_final, "CG_insideAIDMOTIF.csv", row.names = FALSE)

df_aid_analysis_final <- read.csv("CG_insideAIDMOTIF.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)

# ==============================================================================
# 5. Making Grafic
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Generating the summary table based on all rows in the dataframe
table_porcentagem_simple <- df_aid_analysis_final %>%
  group_by(overlap_type_aid) %>%
  summarise(
    n_ocorrencias = n(), # Count all lines/hotspots
    .groups = 'drop'
  ) %>%
  mutate(
    total_geral = sum(n_ocorrencias),
    porcentagem = (n_ocorrencias / total_geral) * 100,
    porcentagem_arredondada = round(porcentagem, 2)
  )

# Display the result in the console
print("--- SIMPLE PERCENTAGE SUMMARY (ALL OCCURRENCES) ---")
print(table_porcentagem_simple)


library(ggplot2)

# Preparing data for the publication standard
plot_data_pub <- table_porcentagem_simple %>%
  mutate(
    overlap_category = case_when(
      overlap_type_aid == "Full Overlap",
      overlap_type_aid == "Partial Overlap",
      overlap_type_aid == "No Overlap",
      overlap_type_aid == "No Hotspot",
      TRUE ~ overlap_type_aid
    ),
    # Sets the exact order of the bars in the chart
    overlap_category = factor(overlap_category, 
                              levels = c("Full Overlap", "Partial Overlap", "No Overlap", "No Hotspot"))
  )

# Colors
colors_pub <- c(
  "Full Overlap"    = "#1B4F72", 
  "Partial Overlap" = "#5DADE2", 
  "No Overlap"      = "#D5D8DC",
  "No Hotspot"      = "#EBEDEF"
)

# Making Graphic
final_plot <- ggplot(plot_data_pub, aes(x = overlap_category, y = porcentagem, fill = overlap_category)) +
  geom_bar(stat = "identity", color = "black", width = 0.65, linewidth = 0.4) +
  geom_text(
    aes(label = paste0(round(porcentagem, 1), "%\n(n=", n_ocorrencias, ")")),
    vjust = -0.4, 
    size = 3.8, 
    fontface = "bold", 
    lineheight = 0.85
  ) +
  # Application of colors and scales
  scale_fill_manual(values = colors_pub) +
  scale_y_continuous(
    limits = c(0, max(plot_data_pub$porcentagem) * 1.3), # Espaço extra para o texto
    expand = c(0, 0) # Remove o espaço entre a barra e o eixo X
  ) +
  labs(x = "AID Motif Overlap Type", y = "Frequency (%)") +
  # Theme (Classic)
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.line = element_line(linewidth = 0.6)
  )

# Visualization
print(final_plot)

# Save in high resolution (TIFF 600 DPI)
ggsave(
  filename = "Figure_AID_Overlap_Publication.tiff", 
  plot = final_plot, 
  device = "tiff", 
  width = 9, 
  height = 5, 
  units = "in", 
  dpi = 1000, 
  compression = "lzw"
)

# Save in JPEG format 
ggsave(
  filename = "Figure_AID_Overlap_Publication.jpg", 
  plot = final_plot, 
  device = "jpeg", 
  width = 9, 
  height = 5, 
  units = "in", 
  dpi = 1000,   # Alta densidade de pixels
  quality = 95 # Qualidade máxima da imagem (0-100)
)