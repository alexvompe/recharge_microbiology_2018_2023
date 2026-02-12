## Author: Alex Vompe
## Date: 03/05/24
## Title: Setting up unfiltered and filtered phyloseq objects
## for Recharge 16S analysis

# Load the libraries====
library(phyloseq)
library(readxl)
library(dplyr)
library(here)
library(decontam)

# Save data objects to large excel files for each run====
taxa = readRDS(here("./sequencing data and dada2/dada2 output/taxa_species.rds"))
asvs = readRDS(here("./sequencing data and dada2/dada2 output/seq_table.rds"))
asvs = t(asvs)
write.csv(taxa, "tax_table.csv")
write.csv(asvs, "asv_table.csv")

# Load phyloseq excel file and make Run 1 ps object====
asv_mat1 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_1.xlsx"), sheet = "asv table")
tax_table1 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_1.xlsx"), sheet = "tax table")
samples_df1 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_1.xlsx"), sheet = "sample data")

asv_mat1 = asv_mat1 %>%
  tibble::column_to_rownames("asv")
tax_table1 = tax_table1 %>% 
  tibble::column_to_rownames("asv")
samples_df1 = samples_df1 %>% 
  tibble::column_to_rownames("sample")
asv_mat1 = as.matrix(asv_mat1)
tax_table1 = as.matrix(tax_table1)
OTU1 = otu_table(asv_mat1, taxa_are_rows = TRUE)
TAX1 = tax_table(tax_table1)
samples1 = sample_data(samples_df1)
ps_1 = phyloseq(OTU1, TAX1, samples1)
ps_1

# Load phyloseq excel file and make Run 2 ps object====
asv_mat2 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_2.xlsx"), sheet = "asv table")
tax_table2 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_2.xlsx"), sheet = "tax table")
samples_df2 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_2.xlsx"), sheet = "sample data")

asv_mat2 = asv_mat2 %>%
  tibble::column_to_rownames("asv")
tax_table2 = tax_table2 %>% 
  tibble::column_to_rownames("asv")
samples_df2 = samples_df2 %>% 
  tibble::column_to_rownames("sample")
asv_mat2 = as.matrix(asv_mat2)
tax_table2 = as.matrix(tax_table2)
OTU2 = otu_table(asv_mat2, taxa_are_rows = TRUE)
TAX2 = tax_table(tax_table2)
samples2 = sample_data(samples_df2)
ps_2 = phyloseq(OTU2, TAX2, samples2)
ps_2

# Load phyloseq excel file and make Run 3 ps object====
asv_mat3 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_3.xlsx"), sheet = "asv table")
tax_table3 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_3.xlsx"), sheet = "tax table")
samples_df3 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_3.xlsx"), sheet = "sample data")

asv_mat3 = asv_mat3 %>%
  tibble::column_to_rownames("asv")
tax_table3 = tax_table3 %>% 
  tibble::column_to_rownames("asv")
samples_df3 = samples_df3 %>% 
  tibble::column_to_rownames("sample")
asv_mat3 = as.matrix(asv_mat3)
tax_table3 = as.matrix(tax_table3)
OTU3 = otu_table(asv_mat3, taxa_are_rows = TRUE)
TAX3 = tax_table(tax_table3)
samples3 = sample_data(samples_df3)
ps_3 = phyloseq(OTU3, TAX3, samples3)
ps_3

# Load phyloseq excel file and make Run 4 ps object====
asv_mat4 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_4.xlsx"), sheet = "asv table")
tax_table4 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_4.xlsx"), sheet = "tax table")
samples_df4 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_4.xlsx"), sheet = "sample data")

asv_mat4 = asv_mat4 %>%
  tibble::column_to_rownames("asv")
tax_table4 = tax_table4 %>% 
  tibble::column_to_rownames("asv")
samples_df4 = samples_df4 %>% 
  tibble::column_to_rownames("sample")
asv_mat4 = as.matrix(asv_mat4)
tax_table4 = as.matrix(tax_table4)
OTU4 = otu_table(asv_mat4, taxa_are_rows = TRUE)
TAX4 = tax_table(tax_table4)
samples4 = sample_data(samples_df4)
ps_4 = phyloseq(OTU4, TAX4, samples4)
ps_4

# Load phyloseq excel file and make Run 5 ps object====
asv_mat5 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_5.xlsx"), sheet = "asv table")
tax_table5 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_5.xlsx"), sheet = "tax table")
samples_df5 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_5.xlsx"), sheet = "sample data")

asv_mat5 = asv_mat5 %>%
  tibble::column_to_rownames("asv")
tax_table5 = tax_table5 %>% 
  tibble::column_to_rownames("asv")
samples_df5 = samples_df5 %>% 
  tibble::column_to_rownames("sample")
asv_mat5 = as.matrix(asv_mat5)
tax_table5 = as.matrix(tax_table5)
OTU5 = otu_table(asv_mat5, taxa_are_rows = TRUE)
TAX5 = tax_table(tax_table5)
samples5 = sample_data(samples_df5)
ps_5 = phyloseq(OTU5, TAX5, samples5)
ps_5

# Load phyloseq excel file and make Run 6 ps object====
asv_mat6 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_6.xlsx"), sheet = "asv table")
tax_table6 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_6.xlsx"), sheet = "tax table")
samples_df6 = read_excel(here::here("./analysis data/initial phyloseq assembly/seq_run_6.xlsx"), sheet = "sample data")

asv_mat6 = asv_mat6 %>%
  tibble::column_to_rownames("asv")
tax_table6 = tax_table6 %>% 
  tibble::column_to_rownames("asv")
samples_df6 = samples_df6 %>% 
  tibble::column_to_rownames("sample")
asv_mat6 = as.matrix(asv_mat6)
tax_table6 = as.matrix(tax_table6)
OTU6 = otu_table(asv_mat6, taxa_are_rows = TRUE)
TAX6 = tax_table(tax_table6)
samples6 = sample_data(samples_df6)
ps_6 = phyloseq(OTU6, TAX6, samples6)
ps_6

# Assemble composite unfiltered ps object ====
ps_recharge_unfiltered = merge_phyloseq(ps_1, ps_2, ps_3,
                                        ps_4, ps_5, ps_6)
#Set date order as a factor
sample_data(ps_recharge_unfiltered)$Date = factor(sample_data(ps_recharge_unfiltered)$Date, 
                                                          levels = c("Jul18", "Nov18", "Mar19", 
                                                                     "Aug19", "Nov19", "Mar20", 
                                                                     "Aug20", "May21", "Aug21", 
                                                                     "Nov21", "Apr22", "Jul22", 
                                                                     "Nov22", "Apr23", "Jul23",
                                                                     "control"))

#Save for easy loading
saveRDS(ps_recharge_unfiltered, "ps_recharge_unfiltered.rds")

# QC and filter the ps object by sample data and taxonomy====
ps_recharge_unfiltered = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_unfiltered.rds"))

#Remove mitochondria, chloroplasts, any other non-bacterial/archaeal seqs
ps_recharge_unfiltered = ps_recharge_unfiltered %>% subset_taxa(Family!= "Mitochondria" | is.na(Family))
ps_recharge_unfiltered = ps_recharge_unfiltered %>% subset_taxa(Order!= "Chloroplast" | is.na(Order))
ps_recharge_unfiltered = ps_recharge_unfiltered %>% subset_taxa(Kingdom!= "Eukaryota" | is.na(Kingdom))
ps_recharge_unfiltered = ps_recharge_unfiltered %>% subset_taxa(Kingdom!= "NA" | is.na(Kingdom))
ps_recharge_unfiltered

#Remove NA problem samples
ps1 = subset_samples(ps_recharge_unfiltered, Nutrients!="NA")
NoSampleNA = subset_samples(ps1, Herbivory!="NA")
NoSampleNA = subset_samples(NoSampleNA, Run!="NA")
NoProblem = subset_samples(NoSampleNA, Nutrients!="multi")
OnlyCorals = subset_samples(NoProblem, Coral!="Hal")

#Remove the pos ctrl
OnlyCorals = subset_samples(OnlyCorals, Coral!="pos")

#Next check for contaminants using prevalence and threshold 0.5 (more conservative)
sample_data(OnlyCorals)$is.neg = sample_data(OnlyCorals)$Coral == "control"
contamdf.prev = isContaminant(OnlyCorals, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant) #135 contam ASVs
head(which(contamdf.prev$contaminant))

#Remove contaminants
physeq.noncont = prune_taxa(!contamdf.prev$contaminant, OnlyCorals)

#Remove controls as we have used them for decontam
physeq.noncont = subset_samples(physeq.noncont, Coral!="control")

#Remove singletons
NoSingle = prune_taxa(taxa_sums(physeq.noncont) > 1, physeq.noncont)

#Let's get the final numbers and remove low-read samples
NoSingleover1000 = prune_samples(sample_sums(NoSingle)
                                 >=1000, NoSingle)
NoSingleover1000 #1381 samples, 44606 ASVs

#Extract the sequences of the filtered object to make the tree
write.csv(tax_table(NoSingleover1000), "allrunsseqs.csv")

#Change ASV names for more simple work
taxa_names(NoSingleover1000) = paste0("ASV", seq(ntaxa(NoSingleover1000)))
tax_table(NoSingleover1000) = cbind(tax_table(NoSingleover1000),
                            rownames(tax_table(NoSingleover1000)))
head(taxa_names(NoSingleover1000))

#Export filtered phyloseq object
saveRDS(NoSingleover1000, "ps_recharge_filtered.rds")