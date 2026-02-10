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
p_profiles = ggplot(df_assigned, aes(x=Date, y=Relative_Abundance,
                                     fill=profile))+
  theme_classic()+
  geom_bar(position = "fill", stat = "identity")+
  scale_fill_manual(values = c("darkgreen", "green4","green3","#03045E","#023E8A",
                               "royalblue4", "blue3",
                               "royalblue3","royalblue","royalblue2",
                               "#0077B6","royalblue1", 
                               "#0096C7", "#00B4D8", "#48CAE4",
                               "#90E0EF", "lightblue","#ADE8F4",
                               "#CAF0F8", "darkorange3","orange3","darkorange",
                               "orange2", "orange", "darkgrey"),
                    "Dominant Type Profile")+
  ylab("Relative Abundance")+
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

ggsave(plot=p_profiles, "corrected_symbio_figure.tiff", units = "mm",
       height = 185, width = 300, scale = 0.8, dpi = 600)