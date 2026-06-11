##Author: Alex Vompe
##Date: 5/14/25
##Title: Regressing microbiome beta dispersion with host tissue loss

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
anova(mod_full_aret) #dispersion ns

mod_full_plob = lmer(percent_dead ~ date_bin*Nutrients*cp_1*dispersion + (1|Plot/ID), 
                     data = subset(df, Coral=="Plob"),
                     na.action=na.omit)
summary(mod_full_plob)
anova(mod_full_plob) #dispersion ns

mod_full_poc = lmer(percent_dead ~ date_bin*Nutrients*cp_1*dispersion + (1|Plot/ID), 
                     data = subset(df, Coral=="Poc"),
                     na.action=na.omit)
summary(mod_full_poc)
anova(mod_full_poc) #dispersion p = 0.000990, F = 11.02