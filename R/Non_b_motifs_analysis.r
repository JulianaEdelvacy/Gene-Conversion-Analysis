# Result analysis from Non Beta motifs search tool

# Loading needed libraries, if you don't have any of these libraries just use the command install.packages(your_library).
library(data.table)

#Couting and plotting the repetitions of direct repeats motifs.
# Loading Direct repeat result file
DR <- read.csv("687514b02f75d_DR.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Load your gene conversion results from Brepconvert
data1 <- read.csv("final_filtered_result.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
df_filtered_selecionado <- as.data.table(df_filtered_selecionado)
