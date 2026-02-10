## Author: Alex Vompe
## Date: 3/11/2024
## Title: Microbiome composition throughout Recharge

# Load the libraries====
library(phyloseq)
library(tidyverse)
library(here)
library(readxl)
library(ggtree)
library(microbiome)
library(ggpubr)
library(patchwork)
library(ggh4x)

# Read in rarefied ps families object====
ps_recharge_rarefied_families = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_rare_tree_fams.rds"))

# Stacked bar with top 20 microbial families data processing====
Top20Fams = names(sort(taxa_sums(ps_recharge_rarefied_families), TRUE)[1:20])
Top20Fams

#Replace the taxonomy of non-top families with "Other" manually to get 20 taxa w Other
# write.csv(tax_table(ps_recharge_rarefied_families), "taxother_rarefied.csv")

#Assemble new phyloseq object with "other" in taxtable
tax_table_other = read_excel(here::here("./analysis data/data frames and csvs/taxother_rarefied.xlsx"))
tax_table_other = tax_table_other %>% 
  tibble::column_to_rownames("asv")
tax_table_other = as.matrix(tax_table_other)
TAX = tax_table(tax_table_other)
tax_table(ps_recharge_rarefied_families) = TAX
Rarefied_w_Other = ps_recharge_rarefied_families

#Save as RDS for downstream access and check file
saveRDS(Rarefied_w_Other, here::here("./analysis data/phyloseq objects/ps_for_stacked_bar.rds"))
Rarefied_w_Other = readRDS(here::here("./analysis data/phyloseq objects/ps_for_stacked_bar.rds"))

#merge samples by covariates of interest
variable1 = as.character(get_variable(Rarefied_w_Other, "Coral"))
variable2 = as.character(get_variable(Rarefied_w_Other, "Date"))
sample_data(Rarefied_w_Other)$Coralbydate <- mapply(paste0, variable1, variable2, 
                                                    collapse = "_")
merge = merge_samples(Rarefied_w_Other, "Coralbydate")

#Transform sample counts with the merged object
relative = transform_sample_counts(merge, function(x) {x/sum(x)})

#For plots that need faceting, make manual modifications
x = data.frame(sample_data(relative))
x$Coralbydate = NULL
x = tibble::rownames_to_column(x, "Coralbydate")
x$data = x$Coralbydate
x = tibble::column_to_rownames(x, "Coralbydate")
x = x %>% separate(data, c("Coral", "Date"), 
                   sep = "(?<=[a-z])(?=[A-Z])")
x$Date = factor(x$Date, levels = c("Jul18", "Nov18", "Mar19", "Aug19",
                                   "Nov19", "Mar20", "Aug20", "May21",
                                   "Aug21","Nov21", "Apr22", "Jul22", 
                                   "Nov22","Apr23", "Jul23"))
x = sample_data(x)
sample_data(relative) = x

# Make the plot====
cbPalette = c("lightgray", "#56B4E9", "#E69F00", "#0055CC",
              "#922B21", "#196F3D", "#7A604B", "#C5B5D4", 
              "#009E73", "#0072B2", "purple", 
              "#CC79A7", "pink", "#FF468F", "#89472F", 
              "#F0E442", "#FF4040", "#66CCCC", "darkorange", 
              "#B4CEFF", "darkblue")

p1=plot_bar(relative3, fill = "Family", 
           x="Date")+
  theme_classic()+
  facet_grid2(~Coral, strip=strip)+
  geom_bar(stat="identity") + 
  scale_fill_manual(values=cbPalette) + 
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  xlab("Date") + 
  ylab("Relative Abundance")
p1$data$Family = factor(p$data$Family)
p1$data$Family = relevel(p$data$Family, "Other")

# Add a tree to the rarefied object====
ps = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree.rds"))

ps_rare = rarefy_even_depth(ps, 
                            rngseed=1, 
                            sample.size=1000, 
                            replace=F)
saveRDS(ps_rare, "ps_recharge_rare_tree.rds")

ps_rare_fams = tax_glom(ps_rare, taxrank = "Family", NArm = FALSE)

saveRDS(ps_rare_fams, "ps_recharge_rare_tree_fams.rds")

# Read in the rarefied phylogenetic family object====
ps = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_rare_tree_fams.rds"))

# merge by coral and date====
variable1 = as.character(get_variable(ps, "Coral"))
variable2 = as.character(get_variable(ps, "Date"))
sample_data(ps)$Coralbydate = mapply(paste0, variable1, variable2,
                                     collapse = "_")
merge = merge_samples(ps, "Coralbydate")

# transform to relative abundance and fix sample data====
relative = transform_sample_counts(merge, function(x) {x/sum(x)})

x = data.frame(sample_data(relative))
x$Coralbydate = NULL
x = tibble::rownames_to_column(x, "Coralbydate")
x$data = x$Coralbydate
x = tibble::column_to_rownames(x, "Coralbydate")
x = x %>% separate(data, c("Coral", "Date"), 
                   sep = "(?<=[a-z])(?=[A-Z])")
x$Date = factor(x$Date, levels = c("Jul18", "Nov18", "Mar19", "Aug19",
                                   "Nov19", "Mar20", "Aug20", "May21",
                                   "Aug21","Nov21", "Apr22", "Jul22", 
                                   "Nov22","Apr23", "Jul23"))
x = sample_data(x)
sample_data(relative) = x

# filter the ps to only contain the families in the relative abundance plot====
ps_phylo = subset_taxa(relative, Family=="Alteromonadaceae" |
                         Family=="Amoebophilaceae" |
                         Family=="Cyanobiaceae" |
                         Family=="Cyclobacteriaceae" |
                         Family=="Endozoicomonadaceae" |
                         Family=="Entomoplasmatales inc. sed." |
                         Family=="Flavobacteriaceae" |
                         Family=="Francisellaceae" |
                         Family=="NA_Alphaproteobacteria" |
                         Family=="NA_Bacilli" |
                         Family=="NA_Bacteria" |
                         Family=="NA_Bacteroidia" |
                         Family=="NA_Campylobacterales" |
                         Family=="Oxalobacteraceae" |
                         Family=="Pirellulaceae" |
                         Family=="Rhodobacteraceae" |
                         Family=="Simkaniaceae" |
                         Family=="Sphingomonadaceae" |
                         Family=="Vibrionaceae" |
                         Family=="Xenococcaceae")

# Make the plot====
melt_simple <- psmelt(ps_phylo) %>% 
  select(Family, OTU, Abundance, Coral, Date)

g = ggtree(ps_phylo, branch.length = "none") + scale_x_reverse()+
  theme(plot.margin = unit(c(0,0,0,0),"pt"))

strip = strip_themed(background_x = elem_list_rect(fill = c("black", 
                                                            "#E69F00", 
                                                            "#56B4E9")),
                     text_x = elem_list_text(color = c("white", "black", 
                                                       "black")))

p = ggplot(melt_simple, aes(x=Date, y=Family, fill=Abundance))+
  theme_classic()+
  facet_grid2(~Coral, strip=strip)+
  geom_tile()+
  scale_fill_gradient2(low = "darkblue", mid = "red", high = "gold",
                       midpoint = 0.5, "Relative Abundance")+
  ylab(NULL)+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust = 1),
        legend.position = "top",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,0,0,0),"pt"))+
  theme(axis.text.x = element_text(color = c(rep("black",2), rep("red", 1),
                                             rep("black",2), rep("red", 1),
                                             rep("black",9))))

p$data$Family = factor(p$data$Family, levels = rev(c("Xenococcaceae",
                                                     "Cyanobiaceae",
                                                     "Entomoplasmatales inc. sed.",
                                                     "NA_Bacilli",
                                                     "Cyclobacteriaceae",
                                                     "Flavobacteriaceae",
                                                     "Amoebophilaceae",
                                                     "NA_Bacteroidia",
                                                     "Pirellulaceae",
                                                     "Rhodobacteraceae",
                                                     "Sphingomonadaceae",
                                                     "NA_Alphaproteobacteria",
                                                     "Simkaniaceae",
                                                     "Vibrionaceae",
                                                     "Alteromonadaceae",
                                                     "NA_Campylobacterales",
                                                     "Oxalobacteraceae",
                                                     "Endozoicomonadaceae",
                                                     "Francisellaceae",
                                                     "NA_Bacteria")))
p2 = p + g + plot_layout(widths = c(6,1))

ggsave(plot=p2, "phylogenetic community composition figure.tiff", 
       units = "mm", 
       height = 185, width = 300, scale = 0.8, dpi = 1000)