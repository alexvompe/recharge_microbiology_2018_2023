## Author: Alex Vompe
## Date: 3/8/2024
## Title: Preparing a fasta file for tree building to add to ps object

# Load the libraries====
library(Biostrings)
library(here)
library(tidyverse)

# Read in all sequences and make fasta file====
df_seqs = read.csv(here::here("./sequencing data and dada2/allrunsseqs.csv"), header = TRUE)
df_seqs = column_to_rownames(df_seqs, var = "X")

seqs = rownames(df_seqs)
length(seqs) #44606

seqs_strings = DNAStringSet(seqs)
names(seqs_strings) = paste0("ASV", seq(length(seqs_strings)))
writeXStringSet(x = seqs_strings, "seqs.fa")

