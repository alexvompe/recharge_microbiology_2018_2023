## Author: Alex Vompe
## Date: 6/12/24
## Title: Integration of coral success data with microbiome data

# load the libraries====
library(here)
library(phyloseq)
library(vegan)
library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)
library(Rmisc)
library(ggh4x)

# load the filtered phyloseq object and export sample data====
ps = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree.rds"))
write.csv(data.frame(sample_data(ps)),
          "ps_recharge_filt_tree_sample_data_integration.csv")
# create an 'integration' column in both coral success.xlsx and the sample
# data with date_ID (unique identifier)

# load coral success.xlsx and sample data and merge by 'integration'====
sample = read.csv(here::here("./analysis data/data frames and csvs/ps_recharge_filt_tree_sample_data_integration.csv"),
                  header = TRUE)

full_df = read_excel(here::here("./analysis data/data frames and csvs/coral success.xlsx"))
full_df = data.frame(full_df)

sample_integrated = merge(sample, full_df)

sample_integrated = column_to_rownames(sample_integrated, 'X')

# QC the integrated sample data====
sample_integrated$Date = factor(sample_integrated$Date, levels = c("Jul18", 
                                                                   "Nov18",
                                                                   "Mar19",
                                                                   "Aug19", 
                                                                   "Nov19",
                                     "Aug20", "May21", "Aug21",
                                     "Nov21", "Apr22", "Jul22",
                                     "Nov22", "Apr23", "Jul23"))

sample_integrated$percent_bleached = as.numeric(sample_integrated$percent_bleached)
sample_integrated$percent_dead = as.numeric(sample_integrated$percent_dead)

# coral success over time paired with micro samples====
df = sample_integrated

df_bleaching = df
df_bleaching$percent_dead = NULL
df_bleaching = na.omit(df_bleaching)
summary_bleaching = summarySE(df_bleaching, measurevar="percent_bleached",
                              groupvars=c("Coral","date_bin",
                                          "Nutrients","cp_1"))
summary_bleaching = dplyr::rename(summary_bleaching, N_bleaching=N,
                                  sd_bleaching=sd, se_bleaching=se,
                                  ci_bleaching=ci)

df_mortality = df
df_mortality$percent_bleached = NULL
df_mortality = na.omit(df_mortality)

df_mortality$date_bin = factor(df_mortality$date_bin,
                               levels = c("pre-MHWs", "MHWs", "MHW recovery",
                                          "MHW + nutrient recovery"))

summary_mortality = summarySE(df_mortality, measurevar="percent_dead",
                              groupvars=c("Coral","date_bin",
                                          "Nutrients","cp_1"))
summary_mortality = dplyr::rename(summary_mortality, N_mortality=N,
                                  sd_mortality=sd, se_mortality=se,
                                  ci_mortality=ci)

cp_1.labs = c("high cpl", "low cpl")
names(cp_1.labs) = c("high","low")

strip = strip_themed(background_x = elem_list_rect(fill = c("black",
                                                            "#E69F00",
                                                            "#56B4E9")),
                     text_x = elem_list_text(color = c("white", "black",
                                                       "black")),
                     background_y = elem_list_rect(fill = c("darkred",
                                                            "darkblue")),
                     text_y = elem_list_text(color = c("white", "white")))

p_trt = ggplot(subset(summary_mortality, date_bin!="pre-MHWs"),
               aes(x=date_bin,y=percent_dead,color=Nutrients,
                   linetype=Nutrients))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_dead-ci_mortality,
                               ymax=percent_dead+ci_mortality))+
  geom_line(aes(group=Nutrients))+
  facet_grid2(cp_1 ~ Coral, strip=strip,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  labs(x="Stage", y="% of Colony Dead")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
        legend.position = "right",panel.spacing.x = unit(0, "lines"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

ggsave(plot=p_trt, "coral success draft_sampled corals.tiff", units="mm", height=185, width = 300,
       dpi=1000, scale=0.8)

# make the new filtered ps object with host data====
sample_integrated = sample_data(sample_integrated)

sample_data(ps) = sample_integrated

saveRDS(ps, "ps_filt_tree_host.rds")

# generate remaining ps objects====

ps_fams = tax_glom(ps, taxrank = "Family", NArm = FALSE)

saveRDS(ps_fams, "ps_filt_tree_fams_host.rds")

ps_rare = rarefy_even_depth(ps, 
                            rngseed=1, 
                            sample.size=1000, 
                            replace=F)

saveRDS(ps_rare, "ps_rare_tree_host.rds")

ps_rare_fams = tax_glom(ps_rare, taxrank = "Family", NArm = FALSE)

saveRDS(ps_rare_fams, "ps_rare_tree_fams_host.rds")