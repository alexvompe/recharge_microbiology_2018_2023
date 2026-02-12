## Author: Alex Vompe
## Date: 5/24/24
## Title: Changes in different microbiome metrics by date

# Load the libraries====
library(phyloseq)
library(tidyverse)
library(here)
library(vegan)
library(Rmisc)
library(rstatix)
library(ggpubr)
library(ggh4x)
library(patchwork)
library(pairwiseAdonis)

# community changes by date====
ps_beta = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree_fams.rds"))
ord_families = readRDS(here::here("./analysis data/data frames and csvs/beta_div_ordination_object.rds"))
ordination_df = plot_ordination(ps_beta, ord_families, type="samples", 
                                color="Date", shape="Coral", justDF = TRUE)

#Axis.1 [63.5%] Axis.2 [5.8%]
strip = strip_themed(background_x = elem_list_rect(fill = c("white", 
                                                            "white", 
                                                            "darkred",
                                                            "white",
                                                            "white",
                                                            "#FF9A98",
                                                            rep("white", 7),
                                                            "#fbd9d3",
                                                            "white")),
                     text_x = elem_list_text(color = c("black", "black", 
                                                       "white", 
                                                       rep("black", 12))))
p1=ggplot(ordination_df, aes(Axis.1, Axis.2, color=Coral)) +
  theme_classic()+
  facet_grid2(~Date, strip = strip)+
  geom_point(alpha=0.4)+
  stat_ellipse(type="norm")+
  stat_ellipse(geom = "point", type="euclid", 
               aes(color=Coral), size=4, level=0.001)+
  scale_color_manual(values = c("black", "#E69F00", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        legend.position = "right",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing = unit(0, "lines"))+
  xlab("PCoA Axis 1 [63.5%]")+
  ylab("PCoA Axis 2 [5.8%]")

#Pairwise adonis by date:
ps_beta_relative = transform_sample_counts(ps_beta, function (x) {x/sum(x)})

#Jul18
Jul18 = subset_samples(ps_beta_relative, Date=="Jul18")
distances = distance(Jul18, method = "unifrac", weighted = TRUE)
df_Jul18 = pairwise.adonis(distances, sample_data(Jul18)$Coral)
df_Jul18$Date = "Jul18"

#Nov18
Nov18 = subset_samples(ps_beta_relative, Date=="Nov18")
distances = distance(Nov18, method = "unifrac", weighted = TRUE)
df_Nov18 = pairwise.adonis(distances, sample_data(Nov18)$Coral)
df_Nov18$Date = "Nov18"

#Mar19
Mar19 = subset_samples(ps_beta_relative, Date=="Mar19")
distances = distance(Mar19, method = "unifrac", weighted = TRUE)
df_Mar19 = pairwise.adonis(distances, sample_data(Mar19)$Coral)
df_Mar19$Date = "Mar19"

#Aug19
Aug19 = subset_samples(ps_beta_relative, Date=="Aug19")
distances = distance(Aug19, method = "unifrac", weighted = TRUE)
df_Aug19 = pairwise.adonis(distances, sample_data(Aug19)$Coral)
df_Aug19$Date = "Aug19"

#Nov19
Nov19 = subset_samples(ps_beta_relative, Date=="Nov19")
distances = distance(Nov19, method = "unifrac", weighted = TRUE)
df_Nov19 = pairwise.adonis(distances, sample_data(Nov19)$Coral)
df_Nov19$Date = "Nov19"

#Mar20
Mar20 = subset_samples(ps_beta_relative, Date=="Mar20")
distances = distance(Mar20, method = "unifrac", weighted = TRUE)
df_Mar20 = pairwise.adonis(distances, sample_data(Mar20)$Coral)
df_Mar20$Date = "Mar20"

#Aug20
Aug20 = subset_samples(ps_beta_relative, Date=="Aug20")
distances = distance(Aug20, method = "unifrac", weighted = TRUE)
df_Aug20 = pairwise.adonis(distances, sample_data(Aug20)$Coral)
df_Aug20$Date = "Aug20"

#May21
May21 = subset_samples(ps_beta_relative, Date=="May21")
distances = distance(May21, method = "unifrac", weighted = TRUE)
df_May21 = pairwise.adonis(distances, sample_data(May21)$Coral)
df_May21$Date = "May21"

#Aug21
Aug21 = subset_samples(ps_beta_relative, Date=="Aug21")
distances = distance(Aug21, method = "unifrac", weighted = TRUE)
df_Aug21 = pairwise.adonis(distances, sample_data(Aug21)$Coral)
df_Aug21$Date = "Aug21"

#Nov21
Nov21 = subset_samples(ps_beta_relative, Date=="Nov21")
distances = distance(Nov21, method = "unifrac", weighted = TRUE)
df_Nov21 = pairwise.adonis(distances, sample_data(Nov21)$Coral)
df_Nov21$Date = "Nov21"

#Apr22
Apr22 = subset_samples(ps_beta_relative, Date=="Apr22")
distances = distance(Apr22, method = "unifrac", weighted = TRUE)
df_Apr22 = pairwise.adonis(distances, sample_data(Apr22)$Coral)
df_Apr22$Date = "Apr22"

#Jul22
Jul22 = subset_samples(ps_beta_relative, Date=="Jul22")
distances = distance(Jul22, method = "unifrac", weighted = TRUE)
df_Jul22 = pairwise.adonis(distances, sample_data(Jul22)$Coral)
df_Jul22$Date = "Jul22"

#Nov22
Nov22 = subset_samples(ps_beta_relative, Date=="Nov22")
distances = distance(Nov22, method = "unifrac", weighted = TRUE)
df_Nov22 = pairwise.adonis(distances, sample_data(Nov22)$Coral)
df_Nov22$Date = "Nov22"

#Apr23
Apr23 = subset_samples(ps_beta_relative, Date=="Apr23")
distances = distance(Apr23, method = "unifrac", weighted = TRUE)
df_Apr23 = pairwise.adonis(distances, sample_data(Apr23)$Coral)
df_Apr23$Date = "Apr23"

#Jul23
Jul23 = subset_samples(ps_beta_relative, Date=="Jul23")
distances = distance(Jul23, method = "unifrac", weighted = TRUE)
df_Jul23 = pairwise.adonis(distances, sample_data(Jul23)$Coral)
df_Jul23$Date = "Jul23"

#Combine df with all permanova comparisons and adjust the p-values
df_permanova = rbind(df_Jul18,df_Nov18,df_Mar19,df_Aug19,df_Nov19,df_Mar20,
                     df_Aug20,df_May21, df_Aug21,df_Nov21,df_Apr22,df_Jul22,
                     df_Nov22,df_Apr23,df_Jul23)
df_permanova$p.adjusted = NULL
df_permanova$p.adjusted = p.adjust(df_permanova$p.value, method = "holm")

#Export
write_csv(df_permanova, "PERMANOVA.csv")

# PERMANOVA by stages----
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
saveRDS(ps_beta, here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree_fams_stages.rds"))
saveRDS(ordination, here::here("./analysis data/data frames and csvs/beta_div_ordination_object_stages.rds"))
ps_beta = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree_fams_stages.rds"))
ordination = readRDS(here::here("./analysis data/data frames and csvs/beta_div_ordination_object_stages.rds"))
ps_beta_relative = transform_sample_counts(ps_beta, function (x) {x/sum(x)})

ordination_df = plot_ordination(ps_beta_relative, ordination, type="samples", 
                                color="Stage", shape="Coral", justDF = TRUE)
ordination_df$Stage = factor(ordination_df$Stage, levels = c("Pre-MHWs",
                                                             "MHWs",
                                                             "MHW Recovery",
                                                             "Enrichment Recovery"))

#Axis.1 [63.5%] Axis.2 [5.8%]
strip = strip_themed(background_x = elem_list_rect(fill = c("white", "darkred",
                                                            "white", "white")),
                     text_x = elem_list_text(color = c("black", "white",
                                                       "black", "black")))

p_stages=ggplot(ordination_df, aes(Axis.1, Axis.2, color=Coral)) +
  theme_classic()+
  facet_grid2(~Stage, strip = strip)+
  geom_point(alpha=0.4)+
  stat_ellipse(type="norm")+
  stat_ellipse(geom = "point", type="euclid", 
               aes(color=Coral), size=4, level=0.001)+
  scale_color_manual(values = c("black", "#E69F00", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        legend.position = "right",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing = unit(0, "lines"))+
  xlab("PCoA Axis 1 [63.5%]")+
  ylab("PCoA Axis 2 [5.8%]")

ggsave(plot = p_stages, "pcoa by stage.tiff", units = "mm", height = 185,
       width = 300, scale = 0.8, dpi = 600)

#Pairwise adonis by stage:
pre_mhw = subset_samples(ps_beta_relative, Stage=="Pre-MHWs")
mhw = subset_samples(ps_beta_relative, Stage=="MHWs")
mhw_recovery = subset_samples(ps_beta_relative, Stage=="MHW Recovery")
enr_recovery = subset_samples(ps_beta_relative, Stage=="Enrichment Recovery")

#pre-mhws
distances = distance(pre_mhw, method = "unifrac", weighted = TRUE)
df_premhw = pairwise.adonis(distances, sample_data(pre_mhw)$Coral)
df_premhw$Stage = "Pre-MHWs"

#mhws
distances = distance(mhw, method = "unifrac", weighted = TRUE)
df_mhw = pairwise.adonis(distances, sample_data(mhw)$Coral)
df_mhw$Stage = "MHWs"

#mhw recovery
distances = distance(mhw_recovery, method = "unifrac", weighted = TRUE)
df_mhwrec = pairwise.adonis(distances, sample_data(mhw_recovery)$Coral)
df_mhwrec$Stage = "MHW Recovery"

#enrichment recovery
distances = distance(enr_recovery, method = "unifrac", weighted = TRUE)
df_enrrec = pairwise.adonis(distances, sample_data(enr_recovery)$Coral)
df_enrrec$Stage = "Enrichment Recovery"

#Combine df with all permanova comparisons and adjust the p-values
df_permanova = rbind(df_premhw, df_mhw, df_mhwrec, df_enrrec)
df_permanova$p.adjusted = NULL
df_permanova$p.adjusted = p.adjust(df_permanova$p.value, method = "holm")

write_csv(df_permanova, "PERMANOVA_stages.csv")

# shannon and beta dispersion changes by date====
df_dispersion = read.csv(here::here("./analysis data/data frames and csvs/main dispersion analysis object.csv"),
                         header = TRUE)
df_dispersion$Date = factor(df_dispersion$Date, levels = c("Jul18", "Nov18", "Mar19",
                                                           "Aug19", "Nov19", "Mar20",
                                                           "Aug20", "May21", "Aug21",
                                                           "Nov21", "Apr22", "Jul22",
                                                           "Nov22", "Apr23", "Jul23"))
df_alpha = read.csv(here::here("./analysis data/data frames and csvs/micro alpha analyses.csv"),
                    header = TRUE)
df_alpha$Date = factor(df_alpha$Date, levels = c("Jul18", "Nov18", "Mar19",
                                                 "Aug19", "Nov19", "Mar20",
                                                 "Aug20", "May21", "Aug21",
                                                 "Nov21", "Apr22", "Jul22",
                                                 "Nov22", "Apr23", "Jul23"))

summary_disp = summarySE(df_dispersion, measurevar="dispersion", 
                         groupvars=c("Coral","Date"))
summary_disp = dplyr::rename(summary_disp, value = dispersion)
summary_disp$Metric = "Beta Dispersion"

#multiply dispersion values by 4, to match the second y-axis scale
#(dispersion will thus be plotted untransformed)
summary_disp$value = summary_disp$value*4
summary_disp$sd = summary_disp$sd*4
summary_disp$se = summary_disp$se*4
summary_disp$ci = summary_disp$ci*4

#Experiment sample sizes
summary_n_ungrouped = summarySE(df_alpha, measurevar="Shannon", 
                                groupvars=c("Coral","Date","Herbivory",
                                            "Nutrients"))
summary_n_grouped = summarySE(df_alpha, measurevar="Shannon", 
                                groupvars=c("Coral","date_bin","cp_1",
                                            "Nutrients"))
write.csv(summary_n_ungrouped, "unbinned sample sizes.csv")
write.csv(summary_n_grouped, "binned sample sizes.csv")

#alpha summary
summary_alpha = summarySE(df_alpha, measurevar="Shannon", 
                          groupvars=c("Coral","Date"))
summary_alpha = dplyr::rename(summary_alpha, value = Shannon)
summary_alpha$Metric = "Shannon Diversity"

df_microbes = rbind(summary_alpha, summary_disp)
df_microbes$Metric = factor(df_microbes$Metric,
                             levels = c("Shannon Diversity",
                                        "Beta Dispersion"))

lvls = levels(df_microbes$Date)
vline.level.1 = 'Mar19'
vline.level.2 = 'Mar20'

p2=ggplot(df_microbes, aes(x=Date, y=value, color=Coral, linetype=Metric))+ 
  theme_classic()+
  geom_line(aes(group=interaction(Coral, Metric)))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=0.3)+
  geom_vline(xintercept = which(lvls==vline.level.1)+9.5,
             linetype = "dotdash", color="#117c13")+
  scale_y_continuous(sec.axis = sec_axis(transform = ~./4,
                                         name = "Beta Dispersion"),
                     "Shannon Diversity Index")+
  annotate("rect", xmin = which(lvls==vline.level.1)-0.5, 
           xmax = which(lvls==vline.level.1)+0.5, ymin = -Inf, ymax = 3.2,
           alpha = .4,fill = "darkred")+
  annotate("rect", xmin = which(lvls==vline.level.2)-0.5, 
           xmax = which(lvls==vline.level.2)+0.5, ymin = -Inf, ymax = 3.2,
           alpha = .2,fill = "darkred")+
  annotate("rect", xmin = which(lvls==vline.level.2)+7.5, 
           xmax = which(lvls==vline.level.2)+8.5, ymin = -Inf, ymax = 3.2,
           alpha = .08,fill = "darkred")+
  annotate(geom = "label", label = "Nutrient Enrichment Period",
           x="Mar20", y=4.5)+
  annotate(geom = "label", label = "No Enrichment",
           x="Apr23", y=4.5)+
  annotate(geom = "segment", x = "Jul18", xend = "Nov18", y = 4,
           yend = 4, arrow = arrow(ends = "both", angle = 90,
                                     length = unit(.2,"cm")))+
  annotate(geom = "segment", x = "Mar19", xend = "Mar20", y = 4,
           yend = 4, arrow = arrow(ends = "both", angle = 90,
                                     length = unit(.2,"cm")))+
  annotate(geom = "segment", x = "Aug20", xend = "Jul22", y = 4,
           yend = 4, arrow = arrow(ends = "both", angle = 90,
                                     length = unit(.2,"cm")))+
  annotate(geom = "segment", x = "Nov22", xend = "Jul23", y = 4,
           yend = 4, arrow = arrow(ends = "both", angle = 90,
                                     length = unit(.2,"cm")))+
  annotate("text", x = which(lvls==vline.level.2)-4.5,
           y = 3.6, label = "pre-MHWs", fontface=2)+
  annotate("text", x = which(lvls==vline.level.2)-1.5,
           y = 3.6, label = "MHWs", fontface=2)+
  annotate("text", x = which(lvls==vline.level.2)+3.5,
           y = 3.6, label = "MHW recovery", fontface=2)+
  annotate("text", x = "Apr23",
           y = 3.45, label = "enrichment\nrecovery", fontface=2)+
  scale_color_manual(values = c("black", "#E69F00", "#56B4E9"),
                     guide = "none")+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

# bleaching and mortality changes by date====
ps = readRDS(here::here("./analysis data/phyloseq objects/ps_rare_tree_fams_host.rds"))
df = data.frame(sample_data(ps))
df = na.omit(df)

df$Date = factor(df$Date, levels = c("Jul18", "Nov18", "Mar19", "May19",
                                     "Aug19", "Nov19", "Mar20", "Aug20",
                                     "May21", "Aug21", "Nov21", "Apr22",
                                     "Jul22", "Nov22", "Apr23", "Jul23"))

summary_mortality = summarySE(df, measurevar = "percent_dead",
                              groupvars = c("Date", "Coral"))
summary_mortality = dplyr::rename(summary_mortality, value = percent_dead)
summary_mortality$Metric = "Mortality"
summary_mortality = rbind(summary_mortality,
                          list('Mar20', 'Aret', 0, 'NA','NA','NA','NA',
                               'Mortality'))
summary_mortality = rbind(summary_mortality,
                          list('Mar20', 'Plob', 0, 'NA','NA','NA','NA',
                               'Mortality'))
summary_mortality = rbind(summary_mortality,
                          list('Mar20', 'Poc', 0, 'NA','NA','NA','NA',
                               'Mortality'))

summary_bleaching = summarySE(df, measurevar = "percent_bleached",
                              groupvars = c("Date", "Coral"))
summary_bleaching = dplyr::rename(summary_bleaching,
                                  value = percent_bleached)
summary_bleaching$Metric = "Bleaching"

summary_bleaching = rbind(summary_bleaching,
                          list('Mar20', 'Aret', 0, 'NA','NA','NA','NA',
                               'Bleaching'))
summary_bleaching = rbind(summary_bleaching,
                          list('Mar20', 'Plob', 0, 'NA','NA','NA','NA',
                               'Bleaching'))
summary_bleaching = rbind(summary_bleaching,
                          list('Mar20', 'Poc', 0, 'NA','NA','NA','NA',
                               'Bleaching'))

df_host = rbind(summary_mortality, summary_bleaching)

df_host$Date = factor(df_host$Date, levels = c("Jul18", "Nov18", "Mar19",
                                     "Aug19", "Nov19", "Mar20", "Aug20",
                                     "May21", "Aug21", "Nov21", "Apr22",
                                     "Jul22", "Nov22", "Apr23", "Jul23"))
df_host$value = as.numeric(df_host$value)
df_host$ci = as.numeric(df_host$ci)

df_host$Metric = factor(df_host$Metric, levels = c("Mortality",
                                                   "Bleaching"))

p3=ggplot(df_host, aes(x=Date, y=value, color=Coral, linetype=Metric))+ 
  theme_classic()+
  geom_line(aes(group=interaction(Coral, Metric)))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0.3)+
  geom_vline(xintercept = which(lvls==vline.level.1)+9.5,
             linetype = "dotdash", color="#117c13")+
  ylab("% of Colony Affected")+
  annotate("rect", xmin = which(lvls==vline.level.1)-0.5, 
           xmax = which(lvls==vline.level.1)+0.5, ymin = -Inf, ymax = Inf,
           alpha = .4,fill = "darkred")+
  annotate("rect", xmin = which(lvls==vline.level.2)-0.5, 
           xmax = which(lvls==vline.level.2)+0.5, ymin = -Inf, ymax = Inf,
           alpha = .2,fill = "darkred")+
  annotate("rect", xmin = which(lvls==vline.level.2)+7.5, 
           xmax = which(lvls==vline.level.2)+8.5, ymin = -Inf, ymax = Inf,
           alpha = .08,fill = "darkred")+
  scale_color_manual(values = c("black", "#E69F00", "#56B4E9"),
                     guide = "none")+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

# make the figure====
p_full = p1 + p2 + p3 + plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold'))


ggsave(plot=p_full, "changes by date.tiff", units="mm",
       height = 370, width = 300, scale=0.7, dpi=1000)