## Author: Alex Vompe
## Date: 2/11/26
## Title: generation of data frames and csvs used in the analysis
## NOTE: some data frames used in the analysis may be generated in other scripts in this repository

# Load the libraries----
library(tidyverse)
library(phyloseq)
library(vegan)
library(here)
library(Rmisc)

# Microbiome beta dispersion with WUniFrac----
ps_beta = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree_fams.rds"))

ord_families = ordinate(ps_beta, "PCoA", "unifrac", weighted=TRUE)
saveRDS(ord_families, "beta_div_ordination_object.rds")

#Distances for each group
dist_aret_nut_low = phyloseq::distance(subset_samples(ps_beta,
                                                      Coral=="Aret" &
                                                        Nutrients=="Enriched" &
                                                        cp_1=="low"), method="unifrac", 
                                       weighted=TRUE)

dist_aret_nut_high = phyloseq::distance(subset_samples(ps_beta,
                                                       Coral=="Aret" &
                                                         Nutrients=="Enriched" &
                                                         cp_1=="high"), method="unifrac", 
                                        weighted=TRUE)

dist_aret_amb_low = phyloseq::distance(subset_samples(ps_beta,
                                                      Coral=="Aret" &
                                                        Nutrients=="Ambient" &
                                                        cp_1=="low"), method="unifrac", 
                                       weighted=TRUE)

dist_aret_amb_high = phyloseq::distance(subset_samples(ps_beta,
                                                       Coral=="Aret" &
                                                         Nutrients=="Ambient" &
                                                         cp_1=="high"), method="unifrac", 
                                        weighted=TRUE)

dist_plob_nut_low = phyloseq::distance(subset_samples(ps_beta,
                                                      Coral=="Plob" &
                                                        Nutrients=="Enriched" &
                                                        cp_1=="low"), method="unifrac", 
                                       weighted=TRUE)

dist_plob_nut_high = phyloseq::distance(subset_samples(ps_beta,
                                                       Coral=="Plob" &
                                                         Nutrients=="Enriched" &
                                                         cp_1=="high"), method="unifrac", 
                                        weighted=TRUE)

dist_plob_amb_low = phyloseq::distance(subset_samples(ps_beta,
                                                      Coral=="Plob" &
                                                        Nutrients=="Ambient" &
                                                        cp_1=="low"), method="unifrac", 
                                       weighted=TRUE)

dist_plob_amb_high = phyloseq::distance(subset_samples(ps_beta,
                                                       Coral=="Plob" &
                                                         Nutrients=="Ambient" &
                                                         cp_1=="high"), method="unifrac", 
                                        weighted=TRUE)

dist_poc_nut_low = phyloseq::distance(subset_samples(ps_beta,
                                                     Coral=="Poc" &
                                                       Nutrients=="Enriched" &
                                                       cp_1=="low"), method="unifrac", 
                                      weighted=TRUE)

dist_poc_nut_high = phyloseq::distance(subset_samples(ps_beta,
                                                      Coral=="Poc" &
                                                        Nutrients=="Enriched" &
                                                        cp_1=="high"), method="unifrac", 
                                       weighted=TRUE)

dist_poc_amb_low = phyloseq::distance(subset_samples(ps_beta,
                                                     Coral=="Poc" &
                                                       Nutrients=="Ambient" &
                                                       cp_1=="low"), method="unifrac", 
                                      weighted=TRUE)

dist_poc_amb_high = phyloseq::distance(subset_samples(ps_beta,
                                                      Coral=="Poc" &
                                                        Nutrients=="Ambient" &
                                                        cp_1=="high"), method="unifrac", 
                                       weighted=TRUE)

#Sample data for each group
sample_aret_nut_low = as.data.frame(sample_data(subset_samples(ps_beta,
                                                               Coral=="Aret" &
                                                                 Nutrients=="Enriched" &
                                                                 cp_1=="low")))

sample_aret_nut_high = as.data.frame(sample_data(subset_samples(ps_beta,
                                                                Coral=="Aret" &
                                                                  Nutrients=="Enriched" &
                                                                  cp_1=="high")))

sample_aret_amb_low = as.data.frame(sample_data(subset_samples(ps_beta,
                                                               Coral=="Aret" &
                                                                 Nutrients=="Ambient" &
                                                                 cp_1=="low")))

sample_aret_amb_high = as.data.frame(sample_data(subset_samples(ps_beta,
                                                                Coral=="Aret" &
                                                                  Nutrients=="Ambient" &
                                                                  cp_1=="high")))

sample_plob_nut_low = as.data.frame(sample_data(subset_samples(ps_beta,
                                                               Coral=="Plob" &
                                                                 Nutrients=="Enriched" &
                                                                 cp_1=="low")))

sample_plob_nut_high = as.data.frame(sample_data(subset_samples(ps_beta,
                                                                Coral=="Plob" &
                                                                  Nutrients=="Enriched" &
                                                                  cp_1=="high")))

sample_plob_amb_low = as.data.frame(sample_data(subset_samples(ps_beta,
                                                               Coral=="Plob" &
                                                                 Nutrients=="Ambient" &
                                                                 cp_1=="low")))

sample_plob_amb_high = as.data.frame(sample_data(subset_samples(ps_beta,
                                                                Coral=="Plob" &
                                                                  Nutrients=="Ambient" &
                                                                  cp_1=="high")))

sample_poc_nut_low = as.data.frame(sample_data(subset_samples(ps_beta,
                                                              Coral=="Poc" &
                                                                Nutrients=="Enriched" &
                                                                cp_1=="low")))

sample_poc_nut_high = as.data.frame(sample_data(subset_samples(ps_beta,
                                                               Coral=="Poc" &
                                                                 Nutrients=="Enriched" &
                                                                 cp_1=="high")))

sample_poc_amb_low = as.data.frame(sample_data(subset_samples(ps_beta,
                                                              Coral=="Poc" &
                                                                Nutrients=="Ambient" &
                                                                cp_1=="low")))

sample_poc_amb_high = as.data.frame(sample_data(subset_samples(ps_beta,
                                                               Coral=="Poc" &
                                                                 Nutrients=="Ambient" &
                                                                 cp_1=="high")))

#Beta dispersion calculations
disp_date_aret_nut_low = betadisper(dist_aret_nut_low,
                                    sample_aret_nut_low$Date,
                                    type = "median")
disp_date_aret_nut_high = betadisper(dist_aret_nut_high,
                                     sample_aret_nut_high$Date,
                                     type = "median")
disp_date_aret_amb_low = betadisper(dist_aret_amb_low,
                                    sample_aret_amb_low$Date,
                                    type = "median")
disp_date_aret_amb_high = betadisper(dist_aret_amb_high,
                                     sample_aret_amb_high$Date,
                                     type = "median")
disp_date_plob_nut_low = betadisper(dist_plob_nut_low,
                                    sample_plob_nut_low$Date,
                                    type = "median")
disp_date_plob_nut_high = betadisper(dist_plob_nut_high,
                                     sample_plob_nut_high$Date,
                                     type = "median")
disp_date_plob_amb_low = betadisper(dist_plob_amb_low,
                                    sample_plob_amb_low$Date,
                                    type = "median")
disp_date_plob_amb_high = betadisper(dist_plob_amb_high,
                                     sample_plob_amb_high$Date,
                                     type = "median")
disp_date_poc_nut_low = betadisper(dist_poc_nut_low,
                                   sample_poc_nut_low$Date,
                                   type = "median")
disp_date_poc_nut_high = betadisper(dist_poc_nut_high,
                                    sample_poc_nut_high$Date,
                                    type = "median")
disp_date_poc_amb_low = betadisper(dist_poc_amb_low,
                                   sample_poc_amb_low$Date,
                                   type = "median")
disp_date_poc_amb_high = betadisper(dist_poc_amb_high,
                                    sample_poc_amb_high$Date,
                                    type = "median")

#Make the data frames
df_aret_nut_low = data.frame(dispersion=disp_date_aret_nut_low$distances,
                             sample_aret_nut_low)
df_aret_nut_high = data.frame(dispersion=disp_date_aret_nut_high$distances,
                              sample_aret_nut_high)
df_aret_amb_low = data.frame(dispersion=disp_date_aret_amb_low$distances,
                             sample_aret_amb_low)
df_aret_amb_high = data.frame(dispersion=disp_date_aret_amb_high$distances,
                              sample_aret_amb_high)
df_plob_nut_low = data.frame(dispersion=disp_date_plob_nut_low$distances,
                             sample_plob_nut_low)
df_plob_nut_high = data.frame(dispersion=disp_date_plob_nut_high$distances,
                              sample_plob_nut_high)
df_plob_amb_low = data.frame(dispersion=disp_date_plob_amb_low$distances,
                             sample_plob_amb_low)
df_plob_amb_high = data.frame(dispersion=disp_date_plob_amb_high$distances,
                              sample_plob_amb_high)
df_poc_nut_low = data.frame(dispersion=disp_date_poc_nut_low$distances,
                            sample_poc_nut_low)
df_poc_nut_high = data.frame(dispersion=disp_date_poc_nut_high$distances,
                             sample_poc_nut_high)
df_poc_amb_low = data.frame(dispersion=disp_date_poc_amb_low$distances,
                            sample_poc_amb_low)
df_poc_amb_high = data.frame(dispersion=disp_date_poc_amb_high$distances,
                             sample_poc_amb_high)

#Make the analysis data frame
df_dispersion_1 = rbind(df_aret_nut_low,
                        df_aret_nut_high,
                        df_aret_amb_low,
                        df_aret_amb_high,
                        df_plob_nut_low,
                        df_plob_nut_high,
                        df_plob_amb_low,
                        df_plob_amb_high,
                        df_poc_nut_low,
                        df_poc_nut_high,
                        df_poc_amb_low,
                        df_poc_amb_high)
write.csv(df_dispersion_1, "main dispersion analysis object.csv")

# Shannon diversity analyses----
df_alpha = data.frame(sample_data(ps_recharge_rarefied_families),
                      estimate_richness(ps_recharge_rarefied_families, 
                                        measures = "Shannon"))
write.csv(df_alpha, "micro alpha analyses.csv")
#Add time bins for improved sample sizes in excel

# Ordination objects----
#By time point
ps_beta = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree_fams.rds"))

ord_families = ordinate(ps_beta, "PCoA", "unifrac", weighted=TRUE)
saveRDS(ord_families, "beta_div_ordination_object.rds")

#By stage
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

ps_beta_relative = transform_sample_counts(ps_beta, function (x) {x/sum(x)})
ordination = ordinate(ps_beta_relative, method = "PCoA",
                      distance = "unifrac", weighted = TRUE)
saveRDS(ordination, here::here("./analysis data/data frames and csvs/beta_div_ordination_object_stages.rds"))