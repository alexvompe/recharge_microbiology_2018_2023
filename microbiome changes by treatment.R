## Author: Alex Vompe
## Date: 5/24/24
## Title: Changes in different microbiome metrics by Recharge treatments

# load the libraries====
library(phyloseq)
library(tidyverse)
library(here)
library(vegan)
library(Rmisc)
library(rstatix)
library(ggpubr)
library(ggh4x)
library(lme4)
library(lmerTest)
library(patchwork)
library(readxl)

# changes by experimental treatment====
#shannon panel
df_alpha = read.csv(here::here("./analysis data/data frames and csvs/micro alpha analyses.csv"),
                    header = TRUE)
df_alpha$Date = factor(df_alpha$Date, levels = c("Jul18", "Nov18", "Mar19",
                                                 "Aug19", "Nov19", "Mar20",
                                                 "Aug20", "May21", "Aug21",
                                                 "Nov21", "Apr22", "Jul22",
                                                 "Nov22", "Apr23", "Jul23"))
df_alpha$date_bin = replace(df_alpha$date_bin,
                            df_alpha$date_bin=="MHW + nutrient recovery",
                            "enrichment recovery")
df_alpha$date_bin = factor(df_alpha$date_bin,
                           levels = c("pre-MHWs", "MHWs", "MHW recovery",
                                      "enrichment recovery"))

df_alpha_trt = subset(df_alpha, Date != "Aug20" & Date != "Jul18")

#stats
shapiro.test(df_alpha_trt$Shannon)#non-normal
pairwise_wilcox_test(df_alpha_trt, Shannon ~ Coral)

stat.test.aret = subset(df_alpha_trt,
                        date_bin != "pre-MHWs" & Coral=="Aret") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(Shannon ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")
stat.test.aret$x = stat.test.aret$x - 1
stat.test.aret$xmin = stat.test.aret$xmin - 1
stat.test.aret$xmax = stat.test.aret$xmax - 1


stat.test.plob = subset(df_alpha_trt,
                        date_bin != "pre-MHWs" & Coral=="Plob") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(Shannon ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")
stat.test.plob$x = stat.test.plob$x - 1
stat.test.plob$xmin = stat.test.plob$xmin - 1
stat.test.plob$xmax = stat.test.plob$xmax - 1

stat.test.poc = subset(df_alpha_trt,
                       date_bin != "pre-MHWs" & Coral=="Poc") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(Shannon ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")
stat.test.poc$x = stat.test.poc$x - 1
stat.test.poc$xmin = stat.test.poc$xmin - 1
stat.test.poc$xmax = stat.test.poc$xmax - 1

stat.save = rbind(stat.test.aret, stat.test.plob, stat.test.poc)
stat.save$groups=NULL
write.csv(stat.save, "shannon wrst results.csv")

cp_1.labs = c("high cpl", "low cpl")
names(cp_1.labs) = c("high","low")

alpha_summary = summarySE(df_alpha_trt, measurevar = "Shannon",
                          groupvars = c("Coral","date_bin",
                                        "Nutrients",
                                        "cp_1"))

#plot
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

p1=ggplot(subset(alpha_summary, date_bin != "pre-MHWs" & Coral=="Aret"),
          aes(x=date_bin, y=Shannon, color=Nutrients))+ 
  facet_grid2(cp_1~Coral, strip = strip_aret,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=Shannon-ci,
                               ymax=Shannon+ci))+
  scale_color_manual(values = c("blue", "#117c13"))+
  stat_pvalue_manual(stat.test.aret,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype=0)+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))+
  ylab("Shannon Diversity Index")+
  xlab("Stage")+
  ylim(0,3.7)+
  theme(strip.text.y = element_blank(),
        legend.position = "none",panel.spacing.x = unit(0, "lines"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p2=ggplot(subset(alpha_summary, date_bin != "pre-MHWs" & Coral=="Plob"),
          aes(x=date_bin, y=Shannon, color=Nutrients))+ 
  facet_grid2(cp_1~Coral, strip = strip_plob,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=Shannon-ci,
                               ymax=Shannon+ci))+
  scale_color_manual(values = c("blue", "#117c13"))+
  stat_pvalue_manual(stat.test.plob,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype=0)+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))+
  ylab("Shannon Diversity Index")+
  xlab("Stage")+
  ylim(0,3.7)+
  theme(strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none",panel.spacing.x = unit(0, "lines"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p3=ggplot(subset(alpha_summary, date_bin != "pre-MHWs" & Coral=="Poc"),
          aes(x=date_bin, y=Shannon, color=Nutrients))+ 
  facet_grid2(cp_1~Coral, strip = strip_poc,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=Shannon-ci,
                               ymax=Shannon+ci))+
  scale_color_manual(values = c("blue", "#117c13"))+
  stat_pvalue_manual(stat.test.poc,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype=0)+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))+
  ylab("Shannon Diversity Index")+
  xlab("Stage")+
  ylim(0,3.7)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "right",panel.spacing.x = unit(0, "lines"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

#dispersion panel
df_dispersion = read.csv(here::here("./analysis data/data frames and csvs/main dispersion analysis object.csv"),
                         header = TRUE)
df_dispersion$Date = factor(df_dispersion$Date, levels = c("Jul18", "Nov18", "Mar19",
                                                           "Aug19", "Nov19", "Mar20",
                                                           "Aug20", "May21", "Aug21",
                                                           "Nov21", "Apr22", "Jul22",
                                                           "Nov22", "Apr23", "Jul23"))
# df_dispersion$date_bin = replace(df_dispersion$date_bin,
#                                  df_dispersion$date_bin=="MHW + nutrient recovery",
#                                  "enrichment recovery")
df_dispersion$date_bin = factor(df_dispersion$date_bin,
                                levels = c("pre-MHWs", "MHWs", "MHW recovery",
                                           "enrichment recovery"))

df_dispersion_trt = subset(df_dispersion, Date != "Aug20" & Date != "Jul18")

strip = strip_themed(background_y = elem_list_rect(fill = c("darkred",
                                                            "#4a2b7a")),
                     text_y = elem_list_text(color = c("white", "white")))
dispersion_summary = summarySE(df_dispersion_trt, measurevar = "dispersion",
                               groupvars = c("Coral","date_bin",
                                             "Nutrients",
                                             "cp_1"))
#QC problem Plob sample
dispersion_summary = dispersion_summary[-c(17), ]

#stats
shapiro.test(df_dispersion_trt$dispersion)#non-normal
stat.test.aret = subset(df_dispersion_trt,
                        date_bin != "pre-MHWs" & Coral=="Aret") %>%
  group_by(Coral, cp_1, date_bin) %>%
  wilcox_test(dispersion ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun = "mean_ci")
stat.test.aret$x = stat.test.aret$x - 1
stat.test.aret$xmin = stat.test.aret$xmin - 1
stat.test.aret$xmax = stat.test.aret$xmax - 1

stat.test.plob = subset(df_dispersion_trt,
                        date_bin != "pre-MHWs" & Coral=="Plob") %>%
  group_by(Coral, cp_1, date_bin) %>%
  wilcox_test(dispersion ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun = "mean_ci")
stat.test.plob$x = stat.test.plob$x - 1
stat.test.plob$xmin = stat.test.plob$xmin - 1
stat.test.plob$xmax = stat.test.plob$xmax - 1

stat.test.poc = subset(df_dispersion_trt,
                        date_bin != "pre-MHWs" & Coral=="Poc") %>%
  group_by(Coral, cp_1, date_bin) %>%
  wilcox_test(dispersion ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun = "mean_ci")
stat.test.poc$x = stat.test.poc$x - 1
stat.test.poc$xmin = stat.test.poc$xmin - 1
stat.test.poc$xmax = stat.test.poc$xmax - 1

stat.save = rbind(stat.test.aret, stat.test.plob, stat.test.poc)
stat.save$groups=NULL
write.csv(stat.save, "dispersion wrst results.csv")

cp_1.labs = c("High CPL", "Low CPL")
names(cp_1.labs) = c("high","low")

#plot
p4=ggplot(subset(dispersion_summary, date_bin!="pre-MHWs" & Coral=="Aret"),
          aes(x=date_bin,y=dispersion,color=Nutrients))+ 
  facet_grid(cp_1~Coral, labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=dispersion-ci,
                               ymax=dispersion+ci))+
  stat_pvalue_manual(stat.test.aret,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype=0)+
  scale_color_manual(values = c("blue", "#117c13"))+
  ylab("Beta Dispersion")+
  xlab("Stage")+
  ylim(0,0.5)+
  theme(legend.position = "none",
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p5=ggplot(subset(dispersion_summary, date_bin!="pre-MHWs" & Coral=="Plob"),
          aes(x=date_bin,y=dispersion,color=Nutrients))+ 
  facet_grid(cp_1~Coral, labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=dispersion-ci,
                               ymax=dispersion+ci))+
  stat_pvalue_manual(stat.test.plob,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype=0)+
  scale_color_manual(values = c("blue", "#117c13"))+
  ylab("Beta Dispersion")+
  xlab("Stage")+
  ylim(0,0.5)+
  theme(legend.position = "none",
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p6=ggplot(subset(dispersion_summary, date_bin!="pre-MHWs" & Coral=="Poc"),
          aes(x=date_bin,y=dispersion,color=Nutrients))+ 
  facet_grid2(cp_1~Coral, strip = strip_poc2,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=dispersion-ci,
                               ymax=dispersion+ci))+
  stat_pvalue_manual(stat.test.poc,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype=0)+
  scale_color_manual(values = c("blue", "#117c13"))+
  ylab("Beta Dispersion")+
  xlab("Stage")+
  ylim(0,0.5)+
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

#plot with just microbiome data
p = p1 + p2 + p3 + p4 + p5 + p6 +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(plot=p, "microbiome by experimental treaments.tiff",
       units="mm", height = 370, width = 300, scale=0.6, dpi=1000)

#just dispersion for all experiment stages (shannon not significant)
#stats
shapiro.test(df_dispersion_trt$dispersion)#non-normal
stat.test.aret = subset(df_dispersion_trt, Coral=="Aret") %>%
  group_by(Coral, cp_1, date_bin) %>%
  wilcox_test(dispersion ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun = "mean_ci")

stat.test.plob = subset(df_dispersion_trt,
                        date_bin != "pre-MHWs" & Coral=="Plob") %>%
  group_by(Coral, cp_1, date_bin) %>%
  wilcox_test(dispersion ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun = "mean_ci")
stat.test.plob$x = stat.test.plob$x - 1
stat.test.plob$xmin = stat.test.plob$xmin - 1
stat.test.plob$xmax = stat.test.plob$xmax - 1

stat.test.poc = subset(df_dispersion_trt, Coral=="Poc") %>%
  group_by(Coral, cp_1, date_bin) %>%
  wilcox_test(dispersion ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun = "mean_ci")

cp_1.labs = c("High CPL", "Low CPL")
names(cp_1.labs) = c("high","low")

#plot
p4=ggplot(subset(dispersion_summary, Coral=="Aret"),
          aes(x=date_bin,y=dispersion,color=Nutrients))+ 
  facet_grid2(cp_1~Coral, strip = strip_aret, labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=dispersion-ci,
                               ymax=dispersion+ci))+
  stat_pvalue_manual(stat.test.aret,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype=0)+
  scale_color_manual(values = c("blue", "#117c13"))+
  ylab("Beta Dispersion")+
  xlab("Stage")+
  ylim(-0.05,0.65)+
  theme(legend.position = "none",
        panel.spacing.x = unit(0, "lines"),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p5=ggplot(subset(dispersion_summary, Coral=="Plob"),
          aes(x=date_bin,y=dispersion,color=Nutrients))+ 
  facet_grid2(cp_1~Coral, strip = strip_plob,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=dispersion-ci,
                               ymax=dispersion+ci))+
  stat_pvalue_manual(stat.test.plob,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype=0)+
  scale_color_manual(values = c("blue", "#117c13"))+
  ylab("Beta Dispersion")+
  xlab("Stage")+
  ylim(-0.05,0.65)+
  theme(legend.position = "none",
        panel.spacing.x = unit(0, "lines"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p6=ggplot(subset(dispersion_summary, Coral=="Poc"),
          aes(x=date_bin,y=dispersion,color=Nutrients))+ 
  facet_grid2(cp_1~Coral, strip = strip_poc,
              labeller = labeller(cp_1 = cp_1.labs))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=dispersion-ci,
                               ymax=dispersion+ci))+
  stat_pvalue_manual(stat.test.poc,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype=0)+
  scale_color_manual(values = c("blue", "#117c13"))+
  ylab("Beta Dispersion")+
  xlab("Stage")+
  ylim(-0.05,0.65)+
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p = p4 + p5 + p6 +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(plot=p, "dispersion by experimental treaments.tiff",
       units="mm", height = 185, width = 300, scale=0.7, dpi=600)