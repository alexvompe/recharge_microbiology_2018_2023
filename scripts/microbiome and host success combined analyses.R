## Author: Alex Vompe
## Date: 6/13/24
## Title: Correlating host success and microbiome data

# load the libraries====
library(tidyverse)
library(phyloseq)
library(here)
library(ggpubr)
library(ggh4x)
library(vegan)
library(Rmisc)
library(rstatix)
library(readxl)
library(patchwork)

# trends by date====
#import new alpha ps object and make new alpha data frame
ps_alpha = readRDS(here::here("./analysis data/phyloseq objects/ps_rare_tree_fams_host.rds"))

df_alpha = data.frame(sample_data(ps_alpha),
                      estimate_richness(ps_alpha,measures = "Shannon"))

#import new filtered ps object and make new dispersion data frame
ps_beta = readRDS(here::here("./analysis data/phyloseq objects/ps_filt_tree_fams_host.rds"))

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

df_dispersion = rbind(df_aret_nut_low,
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
saveRDS(df_dispersion, "dispersion data frame with host info.rds")

df_dispersion = readRDS(here::here("./analysis data/data frames and csvs/dispersion data frame with host info.rds"))

#Plot of general correlations of host success and microbiomes over time
df_full = merge(df_alpha, df_dispersion)

strip = strip_themed(background_y = elem_list_rect(fill = c("black",
                                                             "#E69F00",
                                                             "#56B4E9")),
                      text_y = elem_list_text(color = c("white", "black",
                                                        "black")))
lvls = levels(df_full$Date)
vline.level.1 = 'Mar19'
vline.level.2 = 'Aug20'

p = ggplot(df_full, aes(x=Date))+
  theme_classic()+
  facet_grid2(Coral~., strip=strip)+
  stat_summary(geom = "point", size=2, fun="mean", aes(y=Shannon),
               color="black")+
  stat_summary(geom = "line", fun="mean", aes(y=Shannon, group=1),
               color="black")+
  stat_summary(geom = "point", size=2, fun="mean", aes(y=dispersion*2),
               color="black")+
  stat_summary(geom = "line", fun="mean", aes(y=dispersion*2, group=1),
               color="black", linetype=2)+
  stat_summary(geom = "point", size=2, fun="mean", aes(y=percent_dead/30),
               color="darkred", shape=17)+
  stat_summary(geom = "line", fun="mean", aes(y=percent_dead/30, group=1),
               color="darkred")+
  stat_summary(geom = "point", size=2, fun="mean",
               aes(y=percent_bleached/30),
               color="darkred", shape=17)+
  stat_summary(geom = "line", fun="mean", aes(y=percent_bleached/30,
                                              group=1),
               color="darkred", linetype=2)+
  scale_y_continuous(sec.axis = sec_axis(transform = ~.*30,
                                         "% of Colony Bleached (dashed) and Dead (solid)"),
                     "Microbiome Shannon Diversity Index")+
  annotate("rect", xmin = which(lvls==vline.level.1)-0.5, 
           xmax = which(lvls==vline.level.1)+0.5, ymin = 0,
           ymax = Inf,
           alpha = .4,fill = "darkred")+
  annotate("rect", xmin = which(lvls==vline.level.2)-1, 
           xmax = which(lvls==vline.level.2), ymin = 0,
           ymax = Inf,
           alpha = .2,fill = "darkred")+
  theme(axis.line.y.right = element_line(color="darkred"),
        axis.ticks.y.right = element_line(color = "darkred"),
        axis.text.y.right = element_text(color = "darkred"),
        axis.title.y.right = element_text(color = "darkred"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p_axis = ggplot(df_full, aes(x=Date))+
  theme_classic()+
  facet_grid(Coral~.)+
  stat_summary(geom = "point", size=2, fun="mean", aes(y=Shannon),
               color="white")+
  stat_summary(geom = "point", size=2, fun="mean", aes(y=dispersion*2),
               color="white")+
  stat_summary(geom = "point", size=2, fun="mean", aes(y=percent_dead/30),
               color="white")+
  stat_summary(geom = "point", size=2, fun="mean",
               aes(y=percent_bleached/30),
               color="white")+
  scale_y_continuous(sec.axis = sec_axis(transform = ~.*30,
                                         "Host Success (% of Colony)"),
                     "Microbiome Beta Dispersion x 2")+
  theme(axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_blank(),
        axis.line.y = element_line(linetype = 2))

p_full = p_axis + p + plot_layout(widths = c(0.1, 5))
ggsave(plot=p_full, "general microbiome and success trends.tiff",
       units="mm", height = 185, width=300, scale=0.8, dpi=1000)

#Survivors
# subset of corals that were alive throughout the study
df_success = read_excel(here::here("./analysis data/data frames and csvs/coral success.xlsx"))
df_success = na.omit(df_success) #remove corals with no data
#keep only IDs present in every year
df_success = df_success %>% group_by(ID) %>%
  filter(all(2018:2023 %in% year)) %>% ungroup()
#filter to only IDs that survived until the end
df_survivors = subset(df_success, Date=="Jul23" & percent_dead < 100)

survivors = df_survivors$ID #90 corals

df_alpha_survivors = dplyr::filter(df_alpha, ID %in% survivors)
df_alpha_survivors = na.omit(df_alpha_survivors)
df_alpha_survivors$date_bin = factor(df_alpha_survivors$date_bin,
                                     levels = c("pre-MHWs","MHWs",
                                                "MHW recovery",
                                                "MHW + nutrient recovery"))

df_dispersion_survivors = dplyr::filter(df_dispersion, ID %in% survivors)
df_dispersion_survivors = na.omit(df_dispersion_survivors)
df_dispersion_survivors$date_bin = factor(df_dispersion_survivors$date_bin,
                                     levels = c("pre-MHWs","MHWs",
                                                "MHW recovery",
                                                "MHW + nutrient recovery"))

summary_alpha_survivors = summarySE(df_alpha_survivors,
                                    measurevar = "Shannon",
                                    groupvars = c("date_bin",
                                                  "Coral",
                                                  "Nutrients",
                                                  "cp_1"))
summary_alpha_survivors$date_bin = factor(summary_alpha_survivors$date_bin,
                                          levels=c("pre-MHWs","MHWs",
                                          "MHW recovery",
                                          "MHW + nutrient recovery"))

summary_dispersion_survivors = summarySE(df_dispersion_survivors,
                                    measurevar = "dispersion",
                                    groupvars = c("date_bin",
                                                  "Coral",
                                                  "Nutrients",
                                                  "cp_1"))
summary_dispersion_survivors$date_bin = factor(summary_dispersion_survivors$date_bin,
                                          levels=c("pre-MHWs","MHWs",
                                                   "MHW recovery",
                                                   "MHW + nutrient recovery"))

summary_mortality_survivors = summarySE(df_alpha_survivors,
                                    measurevar = "percent_dead",
                                    groupvars = c("date_bin","Coral",
                                                  "Nutrients",
                                                  "cp_1"))

summary_mortality_survivors$date_bin = factor(summary_mortality_survivors$date_bin,
                                          levels=c("pre-MHWs","MHWs",
                                                   "MHW recovery",
                                                   "MHW + nutrient recovery"))

#stats layers
#alpha
stat.test.aret = subset(df_alpha_survivors,
                        Coral=="Aret" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(Shannon ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")

stat.test.plob = subset(df_alpha_survivors,
                        Coral=="Plob" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(Shannon ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")

stat.test.poc = subset(df_alpha_survivors,
                        Coral=="Poc" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(Shannon ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")
stat.test = rbind(stat.test.aret, stat.test.plob, stat.test.poc)
stat.test$x = stat.test$x - 1
stat.test$xmin = stat.test$xmin - 1
stat.test$xmax = stat.test$xmax - 1

#dispersion
stat.test.aret = subset(df_dispersion_survivors,
                        Coral=="Aret" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(dispersion ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")

stat.test.plob = subset(df_dispersion_survivors,
                        Coral=="Plob" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(dispersion ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")

stat.test.poc = subset(df_dispersion_survivors,
                       Coral=="Poc" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(dispersion ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")
stat.test = rbind(stat.test.aret, stat.test.plob, stat.test.poc)
stat.test$x = stat.test$x - 1
stat.test$xmin = stat.test$xmin - 1
stat.test$xmax = stat.test$xmax - 1

#mortality
stat.test.aret = subset(df_alpha_survivors,
                        Coral=="Aret" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(percent_dead ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")

stat.test.plob = subset(df_alpha_survivors,
                        Coral=="Plob" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(percent_dead ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")

stat.test.poc = subset(df_alpha_survivors,
                       Coral=="Poc" & date_bin != "pre-MHWs") %>%
  group_by(Coral, cp_1, date_bin) %>%
  pairwise_wilcox_test(percent_dead ~ Nutrients) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "date_bin", fun="mean_ci")
stat.test = rbind(stat.test.aret, stat.test.plob, stat.test.poc)
stat.test$x = stat.test$x - 1
stat.test$xmin = stat.test$xmin - 1
stat.test$xmax = stat.test$xmax - 1

#plots
strip1 = strip_themed(background_x = elem_list_rect(fill = c("black",
                                                             "#E69F00",
                                                             "#56B4E9")),
                      text_x = elem_list_text(color = c("white", "black",
                                                        "black")),
                      background_y = elem_list_rect(fill = c("darkred",
                                                             "darkblue")),
                      text_y = elem_list_text(color = c("white", "white")))

cp_1.labs = c("high cpl", "low cpl")
names(cp_1.labs) = c("high","low")

p_alpha = ggplot(subset(summary_alpha_survivors,
                        date_bin!="pre-MHWs"),
                 aes(x=date_bin, y=Shannon, color=Nutrients))+
  theme_classic()+
  facet_grid2(cp_1~Coral, strip=strip1,
              labeller = labeller(cp_1 = cp_1.labs))+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=Shannon-ci,
                               ymax=Shannon+ci))+
  scale_color_manual(values = c("blue", "darkgreen"))+
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype = 0)+
  labs(y="Shannon Diversity Index",x="Stage")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

strip2 = strip_themed(background_y = elem_list_rect(fill = c("darkred",
                                                             "darkblue")),
                      text_y = elem_list_text(color = c("white", "white")))

p_dispersion = ggplot(subset(summary_dispersion_survivors,
                        date_bin!="pre-MHWs"),
                 aes(x=date_bin, y=dispersion, color=Nutrients))+
  theme_classic()+
  facet_grid2(cp_1~Coral, strip=strip2,
              labeller = labeller(cp_1 = cp_1.labs))+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=dispersion-ci,
                               ymax=dispersion+ci))+
  scale_color_manual(values = c("blue", "darkgreen"))+
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype = 0)+
  labs(y="Dispersion",x="Stage")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p_mortality = ggplot(subset(summary_mortality_survivors,
                             date_bin!="pre-MHWs"),
                      aes(x=date_bin, y=percent_dead, color=Nutrients))+
  theme_classic()+
  facet_grid2(cp_1~Coral, strip=strip2,
              labeller = labeller(cp_1 = cp_1.labs))+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=percent_dead-ci,
                               ymax=percent_dead+ci))+
  scale_color_manual(values = c("blue", "darkgreen"))+
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", hide.ns = "p.adj",
                     size = 6, fontface = 2, linetype = 0)+
  labs(y="% of Colony Dead",x="Stage")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p1 = ggarrange(p_alpha,NULL,p_dispersion,NULL,p_mortality,ncol=1,
              labels=c("(a)","","(b)","","(c)"),
              heights=c(1,-0.35,1,-0.35,1),align = "hv")

ggsave(plot=p, "survivor micro and success correlation_stats.tiff",
       units="mm", height = 370, width = 300, scale=0.6, dpi=1000)

# Panel with date bins of complete mortality
#Those who died during the MHW stage
df_success = read_excel(here::here("./analysis data/data frames and csvs/coral success.xlsx"))
df_success = na.omit(df_success) #remove corals with no data

mhw_dead = subset(df_success, Date == "Aug20" & percent_dead == 100)

mhw_dead = mhw_dead$ID #29 corals

df_alpha_mhw_dead = dplyr::filter(df_alpha, ID %in% mhw_dead)
df_alpha_mhw_dead = na.omit(df_alpha_mhw_dead)
df_dispersion_mhw_dead = dplyr::filter(df_dispersion, ID %in% mhw_dead)
df_dispersion_mhw_dead = na.omit(df_dispersion_mhw_dead)

df_alpha_mhw_dead$death_bin = "MHWs"
df_dispersion_mhw_dead$death_bin = "MHWs"

#Those who died during the MHW recovery stage
mhw_survived = subset(df_success, Date == "Aug20" & percent_dead < 100)
mhw_survived_ids = mhw_survived$ID
dead_by_rec_end = subset(df_success, Date=="Jul22" & percent_dead==100)

died_by_rec_end_after_MHWs = dplyr::filter(dead_by_rec_end,
                                           ID %in% mhw_survived_ids)
recovery_bin_dead = died_by_rec_end_after_MHWs$ID #29 corals

df_alpha_mhw_rec_dead = dplyr::filter(df_alpha, ID %in% recovery_bin_dead)
df_alpha_mhw_rec_dead = na.omit(df_alpha_mhw_rec_dead)
df_dispersion_mhw_rec_dead = dplyr::filter(df_dispersion, ID %in% recovery_bin_dead)
df_dispersion_mhw_rec_dead = na.omit(df_dispersion_mhw_rec_dead)

df_alpha_mhw_rec_dead$death_bin = "MHW recovery"
df_dispersion_mhw_rec_dead$death_bin = "MHW recovery"

#Those who died in the MHW + nutrients recovery stage
mhw_rec_survived = subset(df_success, Date == "Jul22" & percent_dead < 100)
mhw_rec_survived_ids = mhw_rec_survived$ID
dead_by_end = subset(df_success, Date=="Jul23" & percent_dead==100)

died_by_end_after_MHW_rec = dplyr::filter(dead_by_end,
                                           ID %in% mhw_rec_survived_ids)
recovery_nut_bin_dead = died_by_end_after_MHW_rec$ID #36 corals

df_alpha_mhw_nut_rec_dead = dplyr::filter(df_alpha,
                                          ID %in% recovery_nut_bin_dead)
df_alpha_mhw_nut_rec_dead = na.omit(df_alpha_mhw_nut_rec_dead)
df_dispersion_mhw_nut_rec_dead = dplyr::filter(df_dispersion,
                                          ID %in% recovery_nut_bin_dead)
df_dispersion_mhw_nut_rec_dead = na.omit(df_dispersion_mhw_nut_rec_dead)

df_alpha_mhw_nut_rec_dead$death_bin = "MHW + nutrient recovery"
df_dispersion_mhw_nut_rec_dead$death_bin = "MHW + nutrient recovery"

#combine into new data frame
df_death = rbind(df_alpha_mhw_dead, df_alpha_mhw_rec_dead,
                 df_alpha_mhw_nut_rec_dead)
df_death_dispersion = rbind(df_dispersion_mhw_dead,
                            df_dispersion_mhw_rec_dead,
                            df_dispersion_mhw_nut_rec_dead)

df_death$death_bin = factor(df_death$death_bin,
                            levels = c("MHWs",
                                       "MHW recovery",
                                       "MHW + nutrient recovery"))
df_death_dispersion$death_bin = factor(df_death_dispersion$death_bin,
                            levels = c("MHWs",
                                       "MHW recovery",
                                       "MHW + nutrient recovery"))
death_alpha_summary = summarySE(df_death, measurevar = "Shannon",
                                groupvars = c("death_bin",
                                              "Coral",
                                              "Nutrients",
                                              "cp_1"))
death_dispersion_summary = summarySE(df_death_dispersion,
                                     measurevar = "dispersion",
                                groupvars = c("death_bin",
                                              "Coral",
                                              "Nutrients",
                                              "cp_1"))

#stats
#nothing is significant

#plots
p_alpha_m = ggplot(death_alpha_summary,
                 aes(x=death_bin, y=Shannon, color=Nutrients))+
  theme_classic()+
  facet_grid2(cp_1~Coral, strip=strip1,
              labeller = labeller(cp_1 = cp_1.labs),
              scales = "free_y",
              independent = "y")+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=Shannon-ci,
                               ymax=Shannon+ci))+
  scale_color_manual(values = c("blue", "darkgreen"))+
  labs(y="Shannon Diversity Index",x="Death Stage")+
  theme(legend.position = "right",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p_dispersion_m = ggplot(death_dispersion_summary,
                      aes(x=death_bin, y=dispersion, color=Nutrients))+
  theme_classic()+
  facet_grid2(cp_1~Coral, strip=strip2,
              labeller = labeller(cp_1 = cp_1.labs),
              scales = "free_y",
              independent = "y")+
  geom_point(size=2)+
  geom_line(aes(group=Nutrients))+
  geom_errorbar(width=0.2, aes(ymin=dispersion-ci,
                               ymax=dispersion+ci))+
  scale_color_manual(values = c("blue", "darkgreen"))+
  labs(y="Dispersion",x="Death Stage")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust=1, hjust = 1),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank())

p2 = ggarrange(p_alpha_m,NULL,p_dispersion_m,ncol=1,labels=c("(d)","","(e)"),
              heights = c(1,-0.25,1), align = "hv")

p = p1+p2+plot_layout(widths = c(0.9, 1))

ggsave(plot=p, "micro and host correlation figure.tiff",
       units="mm", height = 185, width = 300, scale=1.1, dpi=1000)