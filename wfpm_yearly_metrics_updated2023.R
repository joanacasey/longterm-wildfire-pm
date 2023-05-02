####***********************
#### Code Description ####
# Author: Alex
# Date: 12/2/2022
# Last updated: 3/8/2023
# Goal: PM2.5 Yearly Metrics
####**********************

# Read libraries 

library(tidyverse)
library(here)
library(sf)
library(data.table)
library(lubridate)
library(zoo)

# Read in data 
wfpm06_10 <- fread(here("data","wfpm_updated2023","raw","wfpm25_CT_2006to2010_updated2023.csv")) %>% 
  mutate(geoid = as.character(geoid))

wfpm11_16 <- fread(here("data","wfpm_updated2023","raw","wfpm25_CT_2011to2016_updated2023.csv")) %>% 
  mutate(geoid = as.character(geoid))

wfpm16_20 <- fread(here("data","wfpm_updated2023","raw","wfpm25_CT_2017to2020_updated2023.csv")) %>% 
  mutate(geoid = as.character(geoid))

# Bind together data frames 

wfpm <- rbind(wfpm06_10, wfpm11_16, wfpm16_20) %>%
  ungroup()

# Create epiweek and year 

wfpm <- wfpm %>% 
  mutate(week = epiweek(date)) %>% 
  mutate(year = year(date)) %>% 
  mutate(month = month(date))

########## Frequency ###########

# How often exposed?
# 1.
# Number of days with wildfires PM2.5 > 5 ug/m3

# Weekly

weekly_exposure <- wfpm %>% 
  group_by(geoid,week,year, month) %>% 
  summarize(weekly_pm = mean(wf_pm25_imp)) %>%
  ungroup() %>% 
  mutate(exposed = ifelse(weekly_pm > 5, 1,0)) 

# Yearly 

wfpm_freq_y <- weekly_exposure %>% 
  group_by(geoid, year) %>% 
  summarize(pm_freq = sum(exposed)) %>% 
  ungroup() %>% 
  select(geoid,year,pm_freq)

########## Duration ##########

# How long exposed?
# 2. 
# Number of days with non-zero wildfire pm2.5

wfpm_days <- wfpm %>% 
  mutate(non_zero = ifelse(wf_pm25_imp > 0, 1, 0))

wfpm_days_y <- wfpm_days %>% 
  group_by(geoid, year) %>% 
  summarize(non_zero_days = sum(non_zero)) %>% 
  ungroup() %>% 
  select(geoid,year,non_zero_days)

# Number of smoke waves 
# 3. 
# Defined as >=2 days over the study area with wf pm2.5 >25 ug/m3

wfpm_sw <- wfpm %>% 
  mutate(flag = ifelse(wf_pm25_imp > 25, 1,0)) # Change threshold here, as needed  

smoke_wave_y <- wfpm_sw %>% 
  group_by(geoid, group = rleid(flag > 0)) %>% 
  mutate(smoke_wave = ifelse(sum(flag) > 1, 1, 0)) %>%
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  group_by(year, geoid) %>% 
  summarize(smoke_waves = sum(smoke_wave)) %>% 
  ungroup() 

########## Concentration ###########

# To what level exposed?
# Peak 
# 4.
# Average daily wf pm2.5 during peak week 

# Monthly

pm_weekly <- wfpm %>% 
  group_by(geoid, month,week,year) %>% 
  summarize(week_pm = mean(wf_pm25_imp)) %>% 
  ungroup()

peak_pm_y <- pm_weekly %>% 
  group_by(geoid, year) %>% 
  summarize(peak_pm = max(week_pm)) %>% 
  ungroup() %>% 
  select(geoid,year, peak_pm)

# Cumulative
# 5.
# Sum of wf pm2.5

# Yearly sum

wfpm_cmltv_y <- wfpm %>% 
  group_by(geoid, year) %>% 
  summarize(cmltv_pm = sum(wf_pm25_imp)) %>% 
  select(geoid,year,cmltv_pm)

# Average
# 6.
# Annual average wildfire wf pm2.5

# Monthly
# Note, this will have the same values as the yearly values, but in a monthly time series 

wfpm_mean_y <- wfpm %>% 
  mutate(nonwfpm = (ml_pm25-wf_pm25_imp)) %>%
  group_by(geoid, year) %>% 
  summarize(ann_wfpm_avg = mean(wf_pm25_imp),
            ann_nonwfpm = mean(nonwfpm)) %>% 
  ungroup() %>%
  select(geoid,year,ann_wfpm_avg,ann_nonwfpm)

### Aggregate into monthly and yearly datasets 

# Yearly 

pm_metrics_yearly <- wfpm_freq_y %>% 
  left_join(wfpm_days_y, by = c("geoid","year")) %>% 
  left_join(smoke_wave_y, by = c("geoid","year")) %>% 
  left_join(peak_pm_y, by = c("geoid","year")) %>% 
  left_join(wfpm_cmltv_y, by = c("geoid","year")) %>% 
  left_join(wfpm_mean_y, by = c("geoid","year")) %>% 
  mutate(smoke_waves = ifelse(is.na(smoke_waves), 0, smoke_waves)) %>% 
  arrange(year)

saveRDS(pm_metrics_yearly, here("data","wfpm_updated2023","aggregate","wfpm_metrics_yearly.RData"))