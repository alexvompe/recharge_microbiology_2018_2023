## Author: Alex Vompe
## Date: 6/11/24
## Title: Temporal focal coral bleaching and mortality across the experiment

# load the libraries====
library(tidyverse)
library(here)
library(readxl)
library(Rmisc)
library(ggpubr)
library(ggh4x)
library(rstatix)
library(phyloseq)
library(patchwork)

# load and QC the data====
df = read_excel(here::here("./analysis data/data frames and csvs/coral success.xlsx"))

df$Date = factor(df$Date, levels = c("Jul18", "Nov18", "Mar19", "May19",
                                                 "Aug19", "Nov19",
                                                 "Aug20", "May21", "Aug21",
                                                 "Nov21", "Apr22", "Jul22",
                                                 "Nov22", "Apr23", "Jul23"))
df$date_bin = replace(df$date_bin,
                      df$date_bin=="MHW + nutrient recovery",
                      "enrichment recovery")
df$date_bin = factor(df$date_bin,
                               levels = c("pre-MHWs", "MHWs", "MHW recovery",
                                          "enrichment recovery"))

df$percent_bleached = as.numeric(df$percent_bleached)
df$percent_dead = as.numeric(df$percent_dead)
df = na.omit(df)

df_trt = subset(df, Date != "Aug20" & Date != "Jul18")

#Quick summary mortality stats on all corals
pairwise_wilcox_test(df_trt, percent_dead~Coral)
pairwise_wilcox_test(df_trt, percent_bleached~Coral)

summary_bleaching = summarySE(df_trt, measurevar="percent_bleached",
                              groupvars=c("Coral","date_bin",
                                          "Nutrients","cp_1"))

summary_mortality = summarySE(df_trt, measurevar="percent_dead",
                              groupvars=c("Coral","date_bin",
                                          "Nutrients","cp_1"))

# plot bleaching and mortality experiment stages====
cp_1.labs = c("High CPL", "Low CPL")
names(cp_1.labs) = c("high","low")

#bleaching panel
strip_aret = strip_themed(background_x = elem_list_rect(fill = "black"),
                          text_x = elem_list_text(color = "white"))
strip_plob = strip_themed(background_x = elem_list_rect(fill = "#E69F00"),
                          text_x = elem_list_text(color = "black"))
strip_poc = strip_themed(background_x = elem_list_rect(fill = "#56B4E9"),
                         text_x = elem_list_text(color = "black"),
                         background_y = elem_list_rect(fill = c("darkred",
                                                                "#4a2b7a")),
                         text_y = elem_list_text(color = c("white", "white")))
strip_poc2 = strip_themed(background_y = elem_list_rect(fill = c("darkred",
                                                                 "#4a2b7a")),
                          text_y = elem_list_text(color = c("white", "white")))

p1 = ggplot(subset(summary_bleaching, Coral=="Aret"), aes(x=date_bin,
                                   y=percent_bleached,
                                   color=Nutrients,))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_bleached-ci,
                               ymax=percent_bleached+ci))+
  geom_line(aes(group=Nutrients))+
  facet_grid2(cp_1 ~ Coral, strip=strip_aret,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  labs(x="Stage", y="% of Colony Bleached")+
  ylim(-2,30)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none",panel.spacing.x = unit(0, "lines"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

p2 = ggplot(subset(summary_bleaching, Coral=="Plob"), aes(x=date_bin,
                                                          y=percent_bleached,
                                                          color=Nutrients,))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_bleached-ci,
                               ymax=percent_bleached+ci))+
  geom_line(aes(group=Nutrients))+
  facet_grid2(cp_1 ~ Coral, strip=strip_plob,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  labs(x="Stage", y="% of Colony Bleached")+
  ylim(-2,30)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none",panel.spacing.x = unit(0, "lines"),
        legend.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.box.background = element_rect(colour = "black"))

p3 = ggplot(subset(summary_bleaching, Coral=="Poc"), aes(x=date_bin,
                                                          y=percent_bleached,
                                                          color=Nutrients,))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_bleached-ci,
                               ymax=percent_bleached+ci))+
  geom_line(aes(group=Nutrients))+
  facet_grid2(cp_1 ~ Coral, strip=strip_poc,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  labs(x="Stage", y="% of Colony Bleached")+
  ylim(-2,30)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "right",panel.spacing.x = unit(0, "lines"),
        legend.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.box.background = element_rect(colour = "black"))

#mortality panel
stat.test.aret = subset(df_trt, Coral=="Aret" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(percent_dead ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")
stat.test.aret$y.position = stat.test.aret$y.position - 5

stat.test.plob = subset(df_trt, Coral=="Plob" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(percent_dead ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")
stat.test.plob$y.position = stat.test.plob$y.position - 5

stat.test.poc = subset(df_trt, Coral=="Poc" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(percent_dead ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")

p4 = ggplot(subset(summary_mortality, Coral=="Aret"), aes(x=date_bin,
                                             y=percent_dead,
                                             color=Nutrients))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_dead-ci,
                               ymax=percent_dead+ci))+
  geom_line(aes(group=Nutrients))+
  facet_grid(cp_1 ~ Coral,labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  stat_pvalue_manual(stat.test.aret,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype = 0)+
  ylim(-10, 110)+
  labs(x="Stage", y="% of Colony Dead")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_blank())

p5 = ggplot(subset(summary_mortality, Coral=="Plob"), aes(x=date_bin,
                                                          y=percent_dead,
                                                          color=Nutrients))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_dead-ci,
                               ymax=percent_dead+ci))+
  geom_line(aes(group=Nutrients))+
  facet_grid(cp_1 ~ Coral,labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  stat_pvalue_manual(stat.test.plob,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype = 0)+
  ylim(-10, 110)+
  labs(x="Stage", y="% of Colony Dead")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
        legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p6 = ggplot(subset(summary_mortality, Coral=="Poc"), aes(x=date_bin,
                                                          y=percent_dead,
                                                          color=Nutrients))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_dead-ci,
                               ymax=percent_dead+ci))+
  geom_line(aes(group=Nutrients))+
  facet_grid2(cp_1 ~ Coral, strip = strip_poc2,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  stat_pvalue_manual(stat.test.poc,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype = 0)+
  ylim(-10, 110)+
  labs(x="Stage", y="% of Colony Dead")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
        legend.position = "none",
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p = p1 + p2 + p3 + p4 + p5 + p6 + 
  plot_layout(ncol = 3, nrow = 2)+
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(plot=p, "coral success by trt_ALL CORALS.tiff",
       units="mm", height = 370, width = 300, scale=0.6, dpi=1000)

#Version with no bleaching (was not significant)
p4 = ggplot(subset(summary_mortality, Coral=="Aret"), aes(x=date_bin,
                                                          y=percent_dead,
                                                          color=Nutrients))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_dead-ci,
                               ymax=percent_dead+ci))+
  geom_line(aes(group=Nutrients))+
  facet_grid2(cp_1 ~ Coral, strip=strip_aret,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  stat_pvalue_manual(stat.test.aret,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype = 0)+
  ylim(-10, 110)+
  labs(x="Stage", y="% of Colony Dead")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text.y = element_blank())

p5 = ggplot(subset(summary_mortality, Coral=="Plob"), aes(x=date_bin,
                                                          y=percent_dead,
                                                          color=Nutrients))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_dead-ci,
                               ymax=percent_dead+ci))+
  geom_line(aes(group=Nutrients))+
  facet_grid2(cp_1 ~ Coral, strip = strip_plob,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  stat_pvalue_manual(stat.test.plob,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype = 0)+
  ylim(-10, 110)+
  labs(x="Stage", y="% of Colony Dead")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
        legend.position = "none",
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p6 = ggplot(subset(summary_mortality, Coral=="Poc"), aes(x=date_bin,
                                                         y=percent_dead,
                                                         color=Nutrients))+
  geom_point(size=2)+
  geom_errorbar(width=0.2, aes(ymin=percent_dead-ci,
                               ymax=percent_dead+ci))+
  geom_line(aes(group=Nutrients))+
  facet_grid2(cp_1 ~ Coral, strip = strip_poc,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  scale_color_manual(values = c("blue", "#117c13"))+
  stat_pvalue_manual(stat.test.poc,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype = 0)+
  ylim(-10, 110)+
  labs(x="Stage", y="% of Colony Dead")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p = p4 + p5 + p6 + 
  plot_layout(ncol = 3)+
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(plot=p, "coral mortality by trt_ALL CORALS.tiff",
       units="mm", height = 185, width = 300, scale=0.7, dpi=600)

#save stats
stat.save = rbind(stat.test.aret, stat.test.plob, stat.test.poc)
stat.save$groups=NULL
write.csv(stat.save, "mortality wrst results.csv")