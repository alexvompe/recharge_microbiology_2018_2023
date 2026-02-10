## Author: Alex Vompe
## Date: 3/6/24
## Title: Rarefaction curves for microbial alpha diversity

# Load the libraries====
library(phyloseq)
library(tidyverse)
library(vegan)
library(here)
library(ggpubr)

# Load the phyloseq object====
ps_recharge_filtered = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filtered.rds"))

# Samples have highly uneven libraries. Inspect rarefaction curves====
set.seed(1)

calculate_rarefaction_curves <- function(ps_recharge_filtered, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(ps_recharge_filtered, measures, depth) {
    if(max(sample_sums(ps_recharge_filtered)) < depth) return()
    ps_recharge_filtered <- prune_samples(sample_sums(ps_recharge_filtered) >= depth, ps_recharge_filtered)
    
    rarified_microbes.filt <- rarefy_even_depth(ps_recharge_filtered, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_microbes.filt, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, ps_recharge_filtered = ps_recharge_filtered, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(ps_recharge_filtered, c('Observed', 'Shannon'), rep(c(1, 10, 100, 1:100 * 1000), each = 10))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

sampledf = data.frame(sample_data(ps_recharge_filtered))
sampledf = tibble::rownames_to_column(sampledf, "Sample")
sampledf$Sample = gsub('-','.',sampledf$Sample)
head(sampledf)

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, sampledf, by = "Sample")
rarefaction_curve_data_summary_verbose_shan = subset(rarefaction_curve_data_summary_verbose, Measure=="Shannon")
rarefaction_curve_data_summary_verbose_obs = subset(rarefaction_curve_data_summary_verbose, Measure=="Observed")

p1=ggplot(
  data = rarefaction_curve_data_summary_verbose_shan,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    group=Sample)) + 
  geom_point()+
  geom_line()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 1000, linetype="dashed", color="darkred")+
  xlim(0,10000)+
  ylim(0,6)+
  ylab("Shannon Diversity")+
  ggtitle("Shannon Diversity by Random Read \n Depth Subsampling by Sample")

p2=ggplot(
  data = rarefaction_curve_data_summary_verbose_shan,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean)) + 
  stat_summary(fun=median, colour="red", geom="line") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylim(0,6)+
  xlim(0, 10000)+
  geom_vline(xintercept = 1000, linetype="dashed", color="darkred")+
  ylab("Shannon Diversity")+
  ggtitle("Average Shannon Diversity by \n Random Read Depth Subsampling")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p=ggarrange(p1,p2, ncol=2, nrow=1, labels = c("(a)","(b)"),
            align = "hv")
ggsave(plot = p, "shannon rarefaction curves.tiff",
       units = "mm", height = 150, width = 200,
       scale = 1, dpi = 1000)