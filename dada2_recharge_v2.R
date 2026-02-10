## Author: Alex Vompe
## Date: 2/29/24
## Title: DADA2 pipeline for complete Recharge microbes dataset

# Load the libraries====
library(ggplot2)
library(reshape2)
library(dada2)
library(here)

# Let's start by inspecting the quality profiles for each run====
#Set path to input files for each MiSeq run
path = here("./path/to/seqs")
list.files(path)
forwardReads = sort(list.files(path, pattern="_R1_001.fastq.gz",
                               full.names = TRUE))
reverseReads = sort(list.files(path, pattern="_R2_001.fastq.gz", 
                               full.names = TRUE))

sample.names = sapply(strsplit(basename(forwardReads), "_"), `[`, 1)
sample.names

#Assess read quality distribution for the run
plotQualityProfile(forwardReads,aggregate = TRUE)
plotQualityProfile(reverseReads,aggregate = TRUE)

# Quality filter and trim reads for each run====
forwardReads_filt = file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
reverseReads_filt = file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

filteredFastq = filterAndTrim(forwardReads, forwardReads_filt, 
                              reverseReads, reverseReads_filt,
                              truncLen = c(225,185),
                              maxN = 0, 
                              rm.phix=TRUE, maxEE=c(2,2),
                              compress=TRUE, 
                              multithread=FALSE)
filteredFastq
filter_fun = data.frame(filteredFastq)
filter_fun
ratio = sum(filter_fun$reads.out)/sum(filter_fun$reads.in)
ratio #~78%

#Filtered Fastq Quality. Make sure filter works for all runs!!
path = here("./sequencing data/cutadapted sequence data/MiSeq 5/filtered")
list.files(path)
forwardReads_filt = sort(list.files(path, pattern="_F_filt.fastq.gz",
                               full.names = TRUE))
reverseReads_filt = sort(list.files(path, pattern="_R_filt.fastq.gz", 
                               full.names = TRUE))

sample.names = sapply(strsplit(basename(forwardreads_filt), "_"), `[`, 1)
sample.names

plotQualityProfile(forwardReads_filt, aggregate = TRUE)
plotQualityProfile(reverseReads_filt, aggregate = TRUE)

# Proceed with dada2 for each run with same parameters====
path = here("./path/to/seqs/filtered")
list.files(path)
forwardReads_filt = sort(list.files(path, 
                                    pattern="_F_filt.fastq.gz",
                                    full.names = TRUE))
reverseReads_filt = sort(list.files(path, 
                                    pattern="_R_filt.fastq.gz", 
                                    full.names = TRUE))

sample.names = sapply(strsplit(basename(forwardReads_filt), "_"), `[`, 1)
sample.names

#Learn error rates
error_forward = learnErrors(forwardReads_filt,multithread=TRUE)
error_reverse = learnErrors(reverseReads_filt,multithread=TRUE)

#De-replicate reads
forwardReads_filt_derep = derepFastq(forwardReads_filt, 
                                     verbose=TRUE)
head(forwardReads_filt_derep)
reverseReads_filt_derep = derepFastq(reverseReads_filt, 
                                     verbose=TRUE)
head(reverseReads_filt_derep)

names(forwardReads_filt_derep) = sample.names
names(reverseReads_filt_derep) = sample.names

#Sample Inference
dadaForward = dada(forwardReads_filt_derep, err=error_forward, 
                   multithread=TRUE)
dadaReverse = dada(reverseReads_filt_derep, err=error_reverse, 
                   multithread=TRUE)

#Generate contigs
contigs = mergePairs(dadaForward, forwardReads_filt_derep, 
                     dadaReverse, reverseReads_filt_derep)

#Create an ASV table
seq_table = makeSequenceTable(contigs)
table(nchar(getSequences(seq_table)))

#Identify and remove chimeras
seq_table_nochimeri = removeBimeraDenovo(seq_table, 
                                        method="consensus", 
                                        multithread=TRUE, 
                                        verbose=TRUE)
dim(seq_table_nochimeri)

#Sanity check and summary table
getN = function(x) sum(getUniques(x))
track = cbind(sapply(dadaForward, getN), 
              sapply(dadaReverse, getN), 
              sapply(contigs, getN), 
              rowSums(seq_table_nochimeri))
colnames(track) = c("denoisedF", "denoisedR", "merged", 
                    "nonchim")
rownames(track) = sample.names
track

X = melt(track,varnames = c("Sample","Stage"))

write.csv(track, "tracking_samples.csv")

#Save seq table, clear workspace, and read it back in to save
#memory for taxonomic assignment
saveRDS(seq_table_nochimeri, "seq_table.rds")
seq_table_nochimeri = readRDS("seq_table.rds")

#Classify sequences
taxa = assignTaxonomy(seq_table_nochimeri, 
                      "silva_nr99_v138.1_train_set.fa.gz",
                      multithread=TRUE)

saveRDS(taxa, "taxa.rds")

#save and re-load here: RAM for species assignment
taxa = readRDS("taxa.rds")

taxa = addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa, "taxa_species.rds")