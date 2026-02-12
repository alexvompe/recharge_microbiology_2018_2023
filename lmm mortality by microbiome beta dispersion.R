##Author: Alex Vompe
##Date: 5/14/25
##Title: Adding microbiome beta dispersion as a model term for host mortality

# Load the packages----
library(lme4)
library(tidyverse)
library(lmerTest)
library(phyloseq)
library(vegan)

# Load the data and subset to recovery----
df = readRDS(here::here("./analysis data/data frames and csvs/dispersion data frame with host info.rds"))
df = subset(df, date_bin=="MHW recovery" |
              date_bin=="MHW + nutrient recovery")

# Models----
#Full models by Experiment Stage
mod_full_aret = lmer(percent_dead ~ date_bin*Nutrients*cp_1*dispersion + (1|Plot/ID), 
                      data = subset(df, Coral=="Aret"),
                      na.action=na.omit)
summary(mod_full_aret)
anova(mod_full_aret)

mod_full_plob = lmer(percent_dead ~ date_bin*Nutrients*cp_1*dispersion + (1|Plot/ID), 
                     data = subset(df, Coral=="Plob"),
                     na.action=na.omit)
summary(mod_full_plob)
anova(mod_full_plob)

mod_full_poc = lmer(percent_dead ~ date_bin*Nutrients*cp_1*dispersion + (1|Plot/ID), 
                     data = subset(df, Coral=="Poc"),
                     na.action=na.omit)
summary(mod_full_poc)
anova(mod_full_poc)

#Direct effects of dispersion and alpha div on host mortality
mod_aret = lmer(percent_dead ~ dispersion + (1|Plot/ID), 
                     data = subset(df, Coral=="Aret"),
                     na.action=na.omit)
summary(mod_aret)
anova(mod_aret)#ns

mod_plob = lmer(percent_dead ~ dispersion + (1|Plot/ID), 
                     data = subset(df, Coral=="Plob"),
                     na.action=na.omit)
summary(mod_plob)
anova(mod_plob)#ns

mod_poc = lmer(percent_dead ~ dispersion + (1|Plot/ID), 
                    data = subset(df, Coral=="Poc"),
                    na.action=na.omit)
summary(mod_poc)
anova(mod_poc)#significant

# Exploratory plots----
test = cor.test(subset(df, Coral=="Poc")$dispersion,
         subset(df, Coral=="Poc")$percent_dead,
         method = c("pearson", "kendall", "spearman"))