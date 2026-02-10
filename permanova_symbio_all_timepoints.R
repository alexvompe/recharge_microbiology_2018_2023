# Author: Connor Draney
# Date: 12/5/2024
# Title: PERMANOVA on all timepoints

library(tidyverse)
library(vegan)
library(pairwiseAdonis)

early_profs_abs <- read_rds(here::here("./analysis data/early_profs_abs.rds"))
early_profs <- early_profs_abs[-(c(3,4,7,8))]
early_profs = mutate(early_profs, .after = "cpl",
                   cp_1 = case_when(cpl == "1x1" ~ "low",
                                    cpl == "2x2" ~ "low",
                                    cpl == "3x3" ~ "high",
                                    cpl == "open" ~ "high"))

prof_abs <- read_rds(here::here("./prof_abs.rds"))
late_prof <- subset(prof_abs, prof_abs$intact == "intact")
cpl <- c()
for (i in late_prof$id){
  cpl <- c(cpl, str_sub(i, 4, 7))
}
cpl <- gsub("_", "", cpl)
late_prof <- add_column(late_prof, cpl, .after = "nutrient", .name_repair = "minimal")
late_prof <- late_prof[-c(3,6)]
late_prof = mutate(late_prof, .after = "cpl",
                   cp_1 = case_when(cpl == "1x1" ~ "low",
                                    cpl == "2x2" ~ "low",
                                    cpl == "3x3" ~ "high",
                                    cpl == "open" ~ "high"))

comb_profs_abs <- merge(early_profs, late_prof, all = T)
comb_profs_abs$cpl = NULL

# subset to summer timepoints for analysis
comb_profs_abs = subset(comb_profs_abs, date=="Jul18" | 
                          date=="Aug19" | date=="Aug20" | 
                          date=="Jul22")

# permanova on profiles ====
prof_comb <- comb_profs_abs[5:67]
prof_comb <- mutate_all(prof_comb, function(x) as.numeric(as.character(x)))
prof_comb <- prof_comb %>% replace(is.na(.), 0)
tot_comb <- rowSums(prof_comb)
prof_comb <- prof_comb[tot_comb > 0, ]
profile_comb_2 <- comb_profs_abs[tot_comb > 0, ]

set.seed(123) #reproducible p-values with permutations
prof_perm <- adonis2(prof_comb ~ date+nutrient+cp_1, data = profile_comb_2, method = "bray", perms = 999)
prof_perm #no effect by nutrients of cpl, try pairwise adonis by date

set.seed(123)
pair_prof <- pairwise.adonis2(prof_comb ~ date, data = profile_comb_2, method = "bray", perms = 999)
pair_prof # all singificantly different