## Authors: Alex Vompe and Connor Draney
## Date: 2/10/26
## Title: Symbiodiniaceae figure

# Load the libraries----
library(tidyverse)
library(here)
library(readxl)

# Import data----
df = read_excel(here::here("./analysis data/symbiodiniaceae/relabund_summary.xlsx"))
df$Date = factor(df$Date, levels = c("Jul18", "Aug19", "Aug20", "Jul22"))
df_assigned = subset(df, profile != "unassigned") #only plot assigned type profiles

# profiles plot----

#20 most abundant dominant types
top_20_profiles = df_assigned %>%
  group_by(profile) %>%
  summarize(total_abundance = sum(Relative_Abundance, na.rm = TRUE)) %>%
  slice_max(order_by = total_abundance, n = 20) %>%
  pull(profile)

df_top20 = df_assigned %>%
  mutate(profile = if_else(profile %in% top_20_profiles, 
                                  as.character(profile), 
                                  "Other"))

p_profiles = ggplot(df_top20, aes(x=Date, y=Relative_Abundance,
                                     fill=profile))+
  theme_classic()+
  geom_bar(position = "fill", stat = "identity")+
  scale_fill_manual(values = c("green4", "#03081e", "#0a1738", "#112552", "#18376c", 
                               "#1f4885", "#285b9e", "#356fb3", "#4582c3", 
                               "#5796d1", "#6aabde", "#80bfe6", "#98d2ed", 
                               "#b3e4f3", "#cef2f8", "#e3f8fc", "#f2fcfe", "#D95F02",
                               "#F1A340", "#FEE08B", "darkgrey"),
                    "Top 20 Dominant Types")+
  ylab("Average Relative Abundance")+
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

ggsave(plot=p_profiles, "corrected_symbio_figure.tiff", units = "mm",
       height = 300, width = 300, scale = 0.6, dpi = 400)