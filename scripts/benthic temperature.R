## Code adapted with permission from Kelly Speare

library(lubridate)
library(scales)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(here)
library(heatwaveR)

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
  geom_hline(yintercept = 4, linetype = "longdash", color = "darkred")+
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
       scale = 0.8, height = 185, width = 300, dpi = 600)

#2023 MHW severity
lter.out_hist_subset = subset(lter.out, day >= '2023-01-01')
max(lter.out_hist_subset$cumstress)#1.13

# Marine Heatwave Criteria exploration
#Testing
lter.out_test = subset(temperature, '2023-03-17' <= day &
                         day <= '2023-04-22')
p1=ggplot(lter.out_test, aes(x=time_use, y=temperature_c))+
  geom_line()+
  geom_hline(yintercept = 29, color="darkred", linetype = "longdash")+
  theme_bw()+
  scale_x_datetime(breaks = breaks_width("days"),labels = date_format("%D"))+
  labs(x="Day", y=expression("Daily Temperature"~(degree*C)))+
  theme(axis.text.x = element_text(angle = 90))

#But is it a MWH? Needs to have >5 days above 90th percentile
temperature.0=subset(temperature, site=="LTER00")
temperature.0$day=as.Date(temperature.0$day, '%Y-%m-%d')
temperature.0$day=ymd(temperature.0$day)

thermTemp.0=temperature.0

# format time and date
#left the time_local and time_utc columns alone, unaltered. created new date column to reformat
thermTemp.0$date_use = ymd_hms(thermTemp.0$time_local)
thermTemp.0$day_mo_yr = format(thermTemp.0$date_use, '%Y-%m-%d')
thermTemp.0$day = format(thermTemp.0$date_use, '%d') #making new factor for each day (number but as a factor)
thermTemp.0$month = format(thermTemp.0$date_use, '%m') #making new factor for each month (number but as a factor)
thermTemp.0$year = format(thermTemp.0$date_use, '%y') 

thermTemp.0$day=as.factor(thermTemp.0$day)
thermTemp.0$month=as.factor(thermTemp.0$month)
thermTemp.0$year=as.factor(thermTemp.0$year)
thermTemp.0$date_use=as.Date(thermTemp.0$date_use, '%Y-%m-%d')

# subsetting data from August 2018 to July 2020. this will become the line for the high thermal stress year
thermTemp_2018_2020=subset(thermTemp.0, date_use>="2018-07-01" & date_use<="2020-08-31")
thermTemp_2018_2020_mean = thermTemp_2018_2020 %>%
  group_by(day, month, year) %>%
  summarize(
    mean_daily_temp = mean(temperature_c),
    quant = NA,
    .groups = "drop"
  )
thermTemp_2018_2020_mean$timeframe=as.factor("2018to2020")

thermTemp_2018_2020_mean$date= as.Date(with(thermTemp_2018_2020_mean, paste(month, day, year,sep="-")), format="%m-%d-%Y")
years_2018_2020_mean=thermTemp_2018_2020_mean[c(1:3,7,6,4,5)]
#making the upper and lower (mean +- sd) columns
years_2018_2020_mean$temp_upper=years_2018_2020_mean$mean_daily_temp + years_2018_2020_mean$quant
years_2018_2020_mean$temp_lower=years_2018_2020_mean$mean_daily_temp - years_2018_2020_mean$quant
colnames(years_2018_2020_mean)[6]="temp"

# creating a new df with only data up to Dec 31 2017. this will become the mean line and SD
# is >700 rows long because each temp value is listed twice so that it plots over 2 years
thermTemp_toDec2017=subset(thermTemp.0, date_use<"2017-12-31")
thermTemp_toDec2017_mean = thermTemp_toDec2017 %>%
  group_by(day, month) %>%
  summarize(
    mean_daily_temp = mean(temperature_c),
    quant = quantile(temperature_c, probs = 0.9),
    .groups = "drop"
  )

#creating a factor column for this year
thermTemp_toDec2017_mean$timeframe=as.factor("mean")

thermTemp_toDec2017_mean=merge(thermTemp_2018_2020_mean[c(1:3,7)], thermTemp_toDec2017_mean, by=c("day", "month"))
years_toDec2017_mean=thermTemp_toDec2017_mean[c(1:4,7,5,6)] #reordering columns
#making the upper and lower (mean +- sd) columns
years_toDec2017_mean$temp_upper=years_toDec2017_mean$mean_daily_temp + years_toDec2017_mean$quant
years_toDec2017_mean$temp_lower=years_toDec2017_mean$mean_daily_temp - years_toDec2017_mean$quant
colnames(years_toDec2017_mean)[6]="temp"

data=rbind(years_toDec2017_mean, years_2018_2020_mean)

data=merge(data, thermTemp.0, by=c("day", "month", "year"))

##Figure 2 Panel 3
p3=ggplot(data)+
  geom_ribbon(aes(x=date_use, ymin=temp_lower, ymax=temp_upper, fill=timeframe), alpha=0.2)+
  geom_hline(yintercept=29, linetype=2, color="black", linewidth = 1.2)+
  geom_line(aes(x=date_use, y=temp, color=timeframe))+
  scale_color_manual(values = c("purple", "black")) +
  scale_fill_manual(values = c("purple", "pink"))+
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%b%y"), 
               limits=as.Date(c('2018-07-01', '2020-08-31')))+
  labs(x="Date", y=expression("Temperature " ( degree*C)))+
  theme_bw()+
  theme(axis.text.x = element_text(colour="black", angle=45, vjust=1, hjust=1), 
        axis.text.y = element_text(colour="black"))+
  theme(text = element_text(size = 30))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")


# Test Hobday thresholds for MHWs during the study
lter.day = lter.day %>% 
  rename(t = day, temp = temp_c)
lter.day = na.omit(lter.day)

clim = ts2clm(lter.day, climatologyPeriod = c("2004-12-31", "2024-01-01"))

mhw = detect_event(clim)
write.csv(mhw$event, "Hobday MHWs Recharge.csv")
#notes: since there isn't 30 years of history, the 90th percentile
#is overinflated, resulting in 16 MHWs, which is too loose.
#I will push back against using the Hobday definition.