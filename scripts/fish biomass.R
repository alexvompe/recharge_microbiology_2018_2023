## Author: Alex Vompe
## Date: 5/23/25
## Title: Fish functional group biomass by treatment

# Load the packages----
library(tidyverse)
library(rstatix)
library(here)
library(ggh4x)

# Load and QC the data----
df = read_csv(here::here("./analysis data/data frames and csvs/fish biomass.csv"))
df$Stage = factor(df$Stage, levels = c("pre-MHWs","MHWs","MHW recovery",
                                       "enrichment recovery"))
df = subset(df, functionalgroup!="Browser")

#Figure
strip = strip_themed(background_x = elem_list_rect(fill = c("blue","darkgreen")),
                     text_x = elem_list_text(color = c(rep("white", 2))),
                     background_y = elem_list_rect(fill = c("darkred",
                                                            "#4a2b7a")),
                     text_y = elem_list_text(color = c("white", "white")))

nutrient.labs = c("Ambient", "Enriched")
names(nutrient.labs) = c("ambient","enriched")

p = ggplot(df, aes(x=Stage, y=log(biomass),
                       color=functionalgroup))+
  theme_classic()+
  facet_grid2(CPL~nutrients, strip = strip,
              labeller = labeller(nutrients = nutrient.labs))+
  geom_boxplot(position = position_dodge(width=1))+
  stat_summary(geom="point", fun = "mean",
               position = position_dodge(width=1), size = 3)+
  labs(x="Experiment Stage", y="ln(Biomass (g))")+
  scale_color_manual(values = c("turquoise","red","navyblue","#eb6841"),
                     "Fish Functional Group")+
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1))

ggsave(plot = p, "fish biomass figure.tiff", dpi=600,
       units = "mm", height = 185, width = 300, scale = 0.8)

#Stats
#High CPL
df_highcpl = subset(df, CPL=="High CPL")

df_highcpl_corallivore = subset(df_highcpl,
                                functionalgroup=="Corallivore")
df_highcpl_detritivore = subset(df_highcpl,
                               functionalgroup=="Detritivore")
df_highcpl_grazer = subset(df_highcpl,
                           functionalgroup=="Grazer")
df_highcpl_scex = subset(df_highcpl,
                         functionalgroup=="Scraper/Excavator")

test1 = data.frame(tukey_hsd(aov(biomass~nutrients*Stage,
                                 data = df_highcpl_corallivore)))
write_csv(test1, "highcpl_corallivore.csv")

test2 = data.frame(tukey_hsd(aov(biomass~nutrients*Stage,
                                 data = df_highcpl_detritivore)))
write_csv(test2, "highcpl_detritivore.csv")

test3 = data.frame(tukey_hsd(aov(biomass~nutrients*Stage,
                                 data = df_highcpl_grazer)))
write_csv(test3, "highcpl_grazer.csv")

test4 = data.frame(tukey_hsd(aov(biomass~nutrients*Stage,
                                 data = df_highcpl_scex)))
write_csv(test4, "highcpl_scex.csv")

#Low CPL
df_lowcpl = subset(df, CPL=="Low CPL")

df_lowcpl_corallivore = subset(df_lowcpl,
                                functionalgroup=="Corallivore")
df_lowcpl_detritivore = subset(df_lowcpl,
                                functionalgroup=="Detritivore")
df_lowcpl_grazer = subset(df_lowcpl,
                           functionalgroup=="Grazer")
df_lowcpl_scex = subset(df_lowcpl,
                         functionalgroup=="Scraper/Excavator")

test5 = data.frame(tukey_hsd(aov(biomass~nutrients*Stage,
                                 data = df_lowcpl_corallivore)))
write_csv(test5, "lowcpl_corallivore.csv")

test6 = data.frame(tukey_hsd(aov(biomass~nutrients*Stage,
                                 data = df_lowcpl_detritivore)))
write_csv(test6, "lowcpl_detritivore.csv")

test7 = data.frame(tukey_hsd(aov(biomass~nutrients*Stage,
                                 data = df_lowcpl_grazer)))
write_csv(test7, "lowcpl_grazer.csv")

test8 = data.frame(tukey_hsd(aov(biomass~nutrients*Stage,
                                 data = df_lowcpl_scex)))
write_csv(test8, "lowcpl_scex.csv")

#Evaluate cage success for pre-MHWs and MHWs
df_preMHWs = subset(df, Stage=="pre-MHWs")
tukey_hsd(aov(biomass~CPL, data = df_preMHWs))#p = 0.00578

df_MHWs = subset(df, Stage=="MHWs")
tukey_hsd(aov(biomass~CPL, data = df_MHWs))#p = 0.0268

df_MHWrec = subset(df, Stage=="MHW recovery")
tukey_hsd(aov(biomass~CPL, data = df_MHWrec))

df_enrrec = subset(df, Stage=="enrichment recovery")
tukey_hsd(aov(biomass~CPL, data = df_enrrec))

df_preMHWsCorallivore = subset(df, Stage=="pre-MHWs" &
                                 functionalgroup=="Corallivore")
tukey_hsd(aov(biomass~CPL, data = df_preMHWsCorallivore))#ns

df_preMHWsDetritivore = subset(df, Stage=="pre-MHWs" &
                                 functionalgroup=="Detritivore")
tukey_hsd(aov(biomass~CPL, data = df_preMHWsDetritivore))#ns

df_preMHWsGrazer = subset(df, Stage=="pre-MHWs" &
                                 functionalgroup=="Grazer")
tukey_hsd(aov(biomass~CPL, data = df_preMHWsGrazer))#ns

df_preMHWsScex = subset(df, Stage=="pre-MHWs" &
                                 functionalgroup=="Scraper/Excavator")
tukey_hsd(aov(biomass~CPL, data = df_preMHWsScex))#p = 0.0401


df_MHWsCorallivore = subset(df, Stage=="MHWs" &
                                 functionalgroup=="Corallivore")
tukey_hsd(aov(biomass~CPL, data = df_MHWsCorallivore))#ns

df_MHWsDetritivore = subset(df, Stage=="MHWs" &
                                 functionalgroup=="Detritivore")
tukey_hsd(aov(biomass~CPL, data = df_MHWsDetritivore))#p = 0.0301

df_MHWsGrazer = subset(df, Stage=="MHWs" &
                            functionalgroup=="Grazer")
tukey_hsd(aov(biomass~CPL, data = df_MHWsGrazer))#ns

df_MHWsScex = subset(df, Stage=="MHWs" &
                          functionalgroup=="Scraper/Excavator")
tukey_hsd(aov(biomass~CPL, data = df_MHWsScex))#ns
