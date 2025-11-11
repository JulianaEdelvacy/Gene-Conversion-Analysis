#Script to prepare data for gene order filter

#Loading needed libraries, if you already have these libraries intalled you can just coment these lines with #.
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

### Importing IMGT annotation table 
dataset <- read.csv("2_IMGT-gapped-nt-sequences.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

#### Using only columns of interest to work with. 
dataset <- dataset %>% select(Sequence.ID, V.D.J.REGION, V.REGION, V.GENE.and.allele)

#Leaving only unique sequences on the dataset
dataset <- dataset[!duplicated(dataset$V.REGION),]

#length(unique(dataset$V.REGION))

### Counting nrwos on the table
TamanhoTotal <- nrow(dataset)

# Counting the number of lines with empty spaces
datasetvazios <- nrow(dataset[dataset == '',])

# Deleting lines with empty spaces or NA values
dataset <- dataset[!(dataset$V.REGION == "" | is.na(dataset$V.REGION)), ]
dataset <- dataset[!(dataset$V.D.J.REGION == "" | is.na(dataset$V.D.J.REGION)), ]
dataset <- dataset[!(dataset$Sequence.ID == "" | is.na(dataset$Sequence.ID)), ]

# Counting again the number of line with empty spaces, the expected result is to be zero.
datasetcont <- nrow(dataset[dataset == '',])

# Reading result file from Brepconvert
resultado1 <- read.csv("Input_test_Horse_Gene_Conversion.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)

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

#Saving the results
library(openxlsx)
write.table(tabela_juntada1, "final_result.csv", sep = "\t", row.names = FALSE)
