## Author: Alex Vompe
## Date: 2/11/26
## Title:  Script to generate needed phyloseq objects using the filtered phyloseq object

# Load the libraries----
library(tidyverse)
library(phyloseq)
library(here)

# Load the filtered object----
ps_recharge_filtered = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filtered.rds"))

# Generate rarefied object at 1000ASVs/sample based on rarefaction curves----
ps_recharge_rarefied = rarefy_even_depth(ps_recharge_filtered, 
                                         rngseed=1, 
                                         sample.size=1000, 
                                         replace=F)

saveRDS(ps_recharge_rarefied, "ps_recharge_rarefied.RDS")

# Generate rarefied object agglomerated to family----
ps_recharge_rarefied_families = tax_glom(ps_recharge_rarefied,
                                         taxrank = "Family",
                                         NArm = FALSE)
saveRDS(ps_recharge_rarefied_families, "ps_recharge_rarefied_families.rds")

# Generate filtered phyloseq object with phylogeny from FastTree----
tree_filt = read_tree(here::here("./sequencing data and dada2/seqs.tre"))
phy_tree(ps_recharge_filtered) = tree_filt

#tree linux code:
#1) used mothur align seqs with silva 138_1 reference, 30 processors
#2) used mothur filter seqs, 30 processors
#3) used FastTreeMP -nt -gtr to construct tree

saveRDS(ps_recharge_filtered, "ps_recharge_filt_tree.rds")

# Filtered physeq w/ tree agglomerated to family----
ps_filtered_tree_fams = tax_glom(ps_filtered, taxrank = "Family", 
                                 NArm = FALSE)

saveRDS(ps_filtered_tree_fams, "ps_recharge_filt_tree_fams.rds")

# Adding experiment stages to the above phyloseq object----
ps_beta = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree_fams.rds"))
df_sample = as.data.frame(sample_data(ps_beta))
df_sample$Stage = NA

df_sample$Stage[df_sample$Date=="Jul18" |
                  df_sample$Date=="Nov18"] = "Pre-MHWs"
df_sample$Stage[df_sample$Date=="Mar19" |
                  df_sample$Date=="Aug19" |
                  df_sample$Date=="Nov19" |
                  df_sample$Date=="Mar20"] = "MHWs"
df_sample$Stage[df_sample$Date=="Aug20" |
                  df_sample$Date=="May21" |
                  df_sample$Date=="Aug21" |
                  df_sample$Date=="Nov21" |
                  df_sample$Date=="Apr22" |
                  df_sample$Date=="Jul22"] = "MHW Recovery"
df_sample$Stage[df_sample$Date=="Nov22" |
                  df_sample$Date=="Apr23" |
                  df_sample$Date=="Jul23"] = "Enrichment Recovery"

df_sample = sample_data(df_sample)

sample_data(ps_beta) = df_sample
saveRDS(ps_beta, here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree_fams_stages.rds"))