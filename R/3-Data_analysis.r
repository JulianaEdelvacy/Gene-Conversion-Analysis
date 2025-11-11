# Code for data analysis from gene order filter

# Loading needed libraries, if already have these libraries you can just comment these lines with #.
install.packages("Biostrings")
install.packages("dplyr")
install.packages("plyr")
install.packages("parallel")
install.packages("psych")
install.packages("ggplot2")
install.packages("tidyverse")

library(Biostrings)
library(dplyr)
library("plyr")
library(parallel)
library(psych)
library(ggplot2)
library(tidyverse)

# Importing data filtered from Gene_Order_Filter.r
data1 <- read.csv("final_result_filtered.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Couting the repetitions from each column
table(data1$distance_to_AID_motif)
table(data$V.GENE.and.allele)
table(data1$edit_distance)
table(data1$nearest_AID_motif)

# Seeing the most repetated pseudogenes
names_P <- data1$gene
repetitions_P <- table(unlist(names_P))
repetitions_P <- sort(repetitions_P, decreasing = TRUE)

# Printing the names and repetitions
for (name in names(repetitions_P)) {
  print(paste(name, repetitions_P[name]))
}

# Seeing the most repetated functional genes
names_F <- data1$V.GENE.and.allele
repetitions_F <- table(unlist(names_F))
repetitions_F <- sort(repetitions_F, decreasing = TRUE)

# Printing the names and repetitions
for (name in names(repetitions_F)) {
  print(paste(name, repetitions_F[name]))
}

#Counting seq_event repetitions
#Tirando a média de seq event
seq_event_size = nchar(seq_event)
#Couting Seq_event repetitions
repetitions_seq_event <- table(seq_event_size)
print(repetitions_seq_event)
