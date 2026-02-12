## Author: Alex Vompe
## Date: 4/15/24
## Title: Trajectories of Endozoicomonas and Parendozoicomonas in P. lobata

# Load the libraries====
library(phyloseq)
library(here)
library(tidyverse)
library(Rmisc)
library(patchwork)
library(ggh4x)
library(lme4)
library(lmerTest)

# Load the data====
ps_rare = readRDS(here::here("./analysis data/phyloseq objects/ps_rare_tree_host.RDS"))

# Transform sample counts to relative abundance====
relative = transform_sample_counts(ps_rare, function(x) {x/sum(x)})

# Find all ASVs in Family Endozoicomonadaceae====
endos = subset_taxa(relative, Genus=="Endozoicomonas")
top20endos = names(sort(taxa_sums(endos), TRUE)[1:20])
#blast all of these ASVs and see if they are endo or parendo.
#ASVs 17, 20, and 26 are parendozoicomonas!

# Make ps object with only top 20 endo/parendo ASVs====
top_endos = prune_taxa(relative, taxa = top20endos)

# Merge count data with sample data====
df_otu = as.matrix(otu_table(top_endos))
df_otu = t(df_otu)
df_otu = as.data.frame(df_otu)

df_sample = data.frame(sample_data(top_endos))

df_sample = tibble::rownames_to_column(df_sample, "sample")
df_otu = tibble::rownames_to_column(df_otu, "sample")

df_endo_data = merge(df_otu, df_sample, by="sample")

# Plot trajectories over time in Porites lobata====
#let's get the most abundant Endozoicomonas and Parendozoicomonas ASVs

df_endo_plob = subset(df_endo_data, Coral=="Plob")
#Most abundant Endo: ASV1
#Most abundant Parendo: ASV17
df_endo_plob = na.omit(df_endo_plob)

#Endo data frame:
endo_df = summarySE(df_endo_plob, measurevar = "ASV1",
                    groupvars = c("Coral","Date"))
endo_df = dplyr::rename(endo_df, relabund = ASV1)
endo_df$Symbiont = "E. ascidiicola"

parendo_df = summarySE(df_endo_plob, measurevar = "ASV17",
                       groupvars = c("Coral","Date"))
parendo_df = dplyr::rename(parendo_df, relabund = ASV17)
parendo_df$Symbiont = "P. haliclonae"

df_plob = rbind(endo_df, parendo_df)

mortality_df = summarySE(df_endo_plob, measurevar = "percent_dead",
                         groupvars = c("Coral","Date"))

lvls = levels(df_plob$Date)
vline.level.1 = 'Mar19'
vline.level.2 = 'Aug20'

strip = strip_themed(background_x = elem_list_rect(fill = c("#E69F00")))

#Dual-axis plot====
p = ggplot(df_plob, aes(x=Date, y=relabund, color=Symbiont))+
  theme_classic()+
  facet_grid2(~Coral, strip = strip)+
  geom_line(aes(group=Symbiont))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=relabund-ci, ymax=relabund+ci),
                width=0.2)+
  scale_color_manual(values = c("darkolivegreen", "black"),
                     labels = c(expression(paste(italic("E. ascidiicola"),
                                                 " (ASV 1)")),
                                expression(paste(italic("P. haliclonae"),
                                                 " (ASV 17)"))))+
  geom_line(data=mortality_df, color="darkred",
            aes(y=percent_dead/100, group=1))+
  geom_point(data=mortality_df, color="darkred",
             aes(y=percent_dead/100), size=2)+
  geom_errorbar(data=mortality_df, color="darkred", 
                aes(y=percent_dead/100,
                    ymin=(percent_dead-ci)/100,
                    ymax=(percent_dead+ci)/100),
                width=0.2)+
  scale_y_continuous(sec.axis = sec_axis(transform = ~.*100,
                                         name = "% of Colony Dead"),
                     "Symbiont Relative Abundance")+
  annotate("rect", xmin = which(lvls==vline.level.1)-0.5, 
           xmax = which(lvls==vline.level.1)+0.5, ymin = -Inf,
           ymax = Inf,
           alpha = .4,fill = "darkred")+
  annotate("rect", xmin = which(lvls==vline.level.2)-1, 
           xmax = which(lvls==vline.level.2), ymin = -Inf,
           ymax = Inf,
           alpha = .2,fill = "darkred")+
  annotate("rect", xmin = which(lvls==vline.level.2)+6.5, 
           xmax = which(lvls==vline.level.2)+7.5, ymin = -Inf,
           ymax = Inf,
           alpha = .08,fill = "darkred")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
        axis.title.y.right = element_text(color = "darkred"),
        axis.text.y.right = element_text(color = "darkred"),
        axis.ticks.y.right = element_line(color = "darkred"),
        axis.line.y.right = element_line(color="darkred"),
        legend.position = "right",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

#Plot with ordinations
ps_beta = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree_fams_stages.rds"))

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

p_full = p_stages + p +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') +
  plot_layout(ncol = 1) &
  theme(legend.position = "bottom",
        legend.box.background = element_rect(colour = "black"),
        plot.tag = element_text(face = 'bold'))

ggsave(plot=p_full, "Plob_main text.tiff", units="mm",
       scale=0.7, height=370, width=300, dpi=600)

# Analyses----
model = lmer(percent_dead ~ ASV17:Date + (1|ID), data = df_endo_plob)
summary(model)
anova(model) #p=9.797e-05

cor.test(formula=~percent_dead + ASV17, alternative = "less",
         data = df_endo_plob)

model = lm(percent_dead ~ ASV17, data = df_endo_plob)
summary(model)
anova(model)

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

#Pairwise adonis by stage:
ps_beta = readRDS(here::here("./analysis data/phyloseq objects/ps_recharge_filt_tree_fams_stages.rds"))
ordination = readRDS(here::here("./analysis data/data frames and csvs/beta_div_ordination_object_stages.rds"))
ps_beta_relative = transform_sample_counts(ps_beta, function (x) {x/sum(x)})

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