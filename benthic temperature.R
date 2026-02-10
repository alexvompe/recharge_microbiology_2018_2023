## Code adapted with permission from Kelly Speare

library(lubridate)
library(scales)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(here)

# Data --------------------------------------------------------------------
# data from Moorea Coral Reef LTER core time series, bottom mounted termistors:
# http://mcrlter.msi.ucsb.edu/cgi-bin/showDataset.cgi?docid=knb-lter-mcr.1035
# Leichter, J, K. Seydel and C. Gotschalk of Moorea Coral Reef LTER. 2018. MCR LTER: Coral Reef: Benthic Water Temperature, ongoing since 2005. knb-lter-mcr.1035.11

LTER0=read.csv(here::here("./analysis data/benthic temperature data/MCR_LTER00_BottomMountThermistors_20230323.csv"), header=TRUE)
LTER0_2023 = read.csv(here::here("./analysis data/benthic temperature data/LTER00_20220730_20240126.csv"), header=TRUE)
LTER2=read.csv(here::here("./analysis data/benthic temperature data/MCR_LTER02_BottomMountThermistors_20230323.csv"), header=TRUE)
LTER2_2023 = read.csv(here::here("./analysis data/benthic temperature data/LTER02_20230108_20240125.csv"), header=TRUE)

LTER0 = rbind(LTER0, LTER0_2023)
LTER2 = rbind(LTER2, LTER2_2023)

temperature = rbind(LTER0, LTER2)
saveRDS(temperature, "temperature data.rds")

# QC====
temperature = readRDS(here::here("./analysis data/benthic temperature data/temperature data.rds"))

temperature$sensor_depth_m=factor(as.factor(temperature$sensor_depth_m))

temperature = temperature %>%
  subset(reef_type_code=='FOR') %>%
  subset(sensor_depth_m=="10")

temperature$time_use = ymd_hms(temperature$time_local)
temperature$day = format(temperature$time_use, '%Y-%m-%d')

temperature$day=as.Date(temperature$day, '%Y-%m-%d')

temperature$day=ymd(temperature$day)

# mean daily temperature====
lter.day = temperature %>% 
  group_by(day) %>% 
  summarise(temp_c = mean(temperature_c)) %>%
  ungroup()

lter.day$day=as.Date(lter.day$day, '%Y-%m-%d')

lter.day$day = ymd(lter.day$day)

lter.day$day = as.Date(lter.day$day, '%Y-%m-%d')

lter.time.seq = data.frame(day=unique(lter.day$day),ind=seq(1:length(unique(lter.day$day))))
lter.time.seq$year = year(lter.time.seq$day)
lter.time.seq$week.num = week(lter.time.seq$day)
lter.time.seq.week = lter.time.seq %>% group_by(week.num,year) %>% summarise(day=min(day)) %>% ungroup()

lter.day.temp=left_join(lter.day, lter.time.seq, by='day')
lter.week.temp= lter.day.temp %>% group_by(year,week.num) %>% summarise(temp_c = mean(temp_c)) %>% ungroup()
lter.week.temp = lter.week.temp[with(lter.week.temp, order(year,week.num)),]

mma_ref = 29

# calculate accumulated heat stress
lter.week.temp$hotspot = lter.week.temp$temp_c - mma_ref
lter.week.temp$hotspot[lter.week.temp$hotspot < 0] = 0

lter.week.temp$cumstress = NA

# 12 week running sum 
lter.out = lter.week.temp 
for(i in 13:nrow(lter.out)){
  lter.out$cumstress[i] = sum(lter.out$hotspot[(i-12):i],na.rm=T)
}

lter.out$cumstress[is.na(lter.out$cumstress)] = 0
lter.out = left_join(lter.out, lter.time.seq.week, by=c('year','week.num'))

lter.out_hist_subset = subset(lter.out, day >= '2020-01-01')
max(lter.out_hist_subset$cumstress)#4.26
lter.out_hist_subset = subset(lter.out, day >= '2010-01-01')
max(lter.out_hist_subset$cumstress)#6.10

# Summary daily average temperature panel====
p1=ggplot(lter.out_hist_subset, aes(x=day, y=temp_c))+
  geom_line()+
  geom_hline(yintercept = 29, color="darkred", linetype = "longdash")+
  theme_classic()+
  geom_vline(xintercept = as.Date("2018-07-02"), color="darkblue")+
  geom_vline(xintercept = as.Date("2023-08-02"), color="darkblue")+
  scale_x_date(breaks = date_breaks("years"),labels = date_format("%Y"))+
  annotate(geom = "segment", x = as.Date("2018-07-02"),
           xend = as.Date("2023-08-02"), y = 30.5,
           yend = 30.5, arrow = arrow(ends = "both", angle = 90,
                                     length = unit(.2,"cm")),
           color="darkblue")+
  annotate(geom = "text", x = as.Date("2021-02-01"),y = 30.2,
           label = "Experiment", color="darkblue")+
  labs(x="Year", y=expression("Weekly Avg. Temperature"~(degree*C)))+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey",
                                          linetype = "dashed"))

# Summary heat stress panel====
p2=ggplot(lter.out_hist_subset, aes(x=day, y=cumstress))+
  geom_line()+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = "longdash", color = "darkred")+
  geom_vline(xintercept = as.Date("2018-07-02"), color="darkblue")+
  geom_vline(xintercept = as.Date("2023-08-02"), color="darkblue")+
  scale_x_date(breaks = date_breaks("years"), labels = date_format("%Y"))+
  labs(x="Year", y="Acc. Heat Stress (°C-weeks)")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        panel.grid.major.x = element_line(color = "grey",
                                          linetype = "dashed"))

p = p1 + p2 + plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(plot = p, "recharge temperature.tiff", units = "mm",
       scale = 0.8, height = 185, width = 300, dpi = 1000)

#2023 MHW
lter.out_hist_subset = subset(lter.out, day >= '2023-01-01')
max(lter.out_hist_subset$cumstress)#1.13
p2=ggplot(lter.out_hist_subset, aes(x=day, y=cumstress))+
  geom_line()+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = "longdash", color = "darkred")+
  scale_x_date(breaks = date_breaks("months"),
               labels = date_format("%b%Y"))+
  labs(x="Year", y="Acc. Heat Stress")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        panel.grid.major.x = element_line(color = "grey",
                                          linetype = "dashed"))