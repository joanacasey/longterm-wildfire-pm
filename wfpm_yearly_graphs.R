####***********************
#### Code Description ####
# Author: Alex
# Date: 12/2/2022
# Updated: 3/08/2023
# Goal: PM2.5 Yearly Metrics
####**********************

# Read libraries 

options(scipen = 999)
library(tidyverse)
library(here)
library(sf)
library(data.table)
library(lubridate)
library(zoo)
library(patchwork)
library(MetBrewer)

# From TIGRIS file 

cact <- read_sf(here("data","shapes","ca_cts_google_drive")) %>% 
  janitor::clean_names() %>%
  mutate(geoid = str_remove(geoid, "^0+")) %>% 
  rename(area = shape_area) %>% 
  st_drop_geometry() %>% 
  select(geoid, area) 

cact <- cact %>% 
  mutate(area_prop = area/sum(cact$area)) 

# Pre-aggregated data

pm_metrics_yearly <- readRDS(here("data","wfpm_updated2023","aggregate","wfpm_metrics_yearly.RData"))

# Population data 

ct_pop <- read_csv(here("data","ca_census_tracts_2010_population",
                        "nhgis0028_ds172_2010_tract.csv")) %>% 
  janitor::clean_names()

ca_ct_pop <- ct_pop %>% 
  filter(statea == "06") %>%
  mutate(geoid = paste0(statea, countya, tracta)) %>%
  mutate(geoid = str_sub(geoid, 2)) %>%
  select(geoid, h7v001) %>%
  rename(population = h7v001)

# Set up population-weighted CTs

pm_metrics_yearly <- pm_metrics_yearly %>% 
  left_join(ca_ct_pop, by = "geoid") 

sum(ca_ct_pop$population) # 37253956

pm_metrics_yearly <- pm_metrics_yearly %>%
  mutate(population = ifelse(is.na(population), 0, population)) %>%
  mutate(pop_prop = population/37253956) 

# Set up area-weighted CTs

pm_metrics_yearly <- pm_metrics_yearly %>% 
  left_join(cact, by = "geoid") 

pm_metrics_yearly <- pm_metrics_yearly %>% 
  mutate(area = ifelse(is.na(area), 0, area),
         area_prop = ifelse(is.na(area_prop), 0, area_prop))

########## Frequency ###########

# How often exposed?
# 1.
# Number of weeks with wildfires PM2.5 > 5 ug/m3

weekly_wf <- pm_metrics_yearly %>% 
  group_by(year) %>% 
  summarize(pm_freq = sum(pm_freq))

plot_freq <- ggplot(weekly_wf, aes(x=year, y=pm_freq, fill = pm_freq)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Year") +
  ylab(expression("N census tract weeks")) +
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(
    text=element_text(family="Gulliver"))
  
########## Duration ###########

# How long exposed?
# 2. 
# Number of days with non-zero wildfire pm2.5

non_zero_days <- pm_metrics_yearly %>% 
  group_by(year) %>% 
  summarize(duration = sum(non_zero_days)) 

plot_non_zero <- ggplot(non_zero_days, aes(x=year, y=duration, fill = duration)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Year") +
  ylab(expression("N census tract days")) +
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver"))
  
# Number of smoke waves 
# 3. 
# Defined as >=2 days over the study area >25 ug/m3 of wildfire PM2.5 

smoke_waves <- pm_metrics_yearly %>% 
  group_by(year) %>% 
  summarize(smoke_wave = sum(smoke_waves)) 

plot_sw <- ggplot(smoke_waves, aes(x=year, y=smoke_wave, fill = smoke_wave)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Year") +
  ylab(expression("N census tract smoke waves")) +
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver"))

########## Concentration ###########

# To what level exposed?
# Peak 
# 4.
# Average daily wf pm2.5 during peak week 

# Non-weighted 

peak <- pm_metrics_yearly %>% 
  group_by(year) %>% 
  summarize(peak_pm = mean(peak_pm)) 

plot_peak <- ggplot(peak, aes(x=year, y=peak_pm, fill = peak_pm)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Year") +
  ylab(expression("Average peak week mean wildfire PM"[2.5])) + 
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver"))

# Population-weighted 

peak_pw <- pm_metrics_yearly %>% 
  group_by(year) %>% 
  summarize(pw_peak = sum((pop_prop * peak_pm)))

plot_peak_pw <- ggplot(peak_pw, aes(x=year, y=pw_peak, fill = pw_peak)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Year") +
  ylab(expression("Population-weighted peak week mean wildfire PM"[2.5])) + 
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver"))

# Areal-weighted 

peak_aw <- pm_metrics_yearly %>% 
  group_by(year) %>% 
  summarize(aw_peak = sum((area_prop * peak_pm)))

plot_peak_aw <- ggplot(peak_aw, aes(x=year, y=aw_peak, fill = aw_peak)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Year") +
  ylab(expression("Area-weighted peak week mean wildfire PM"[2.5])) + 
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver"))

# Test the three different weighting schemes 

plot_peak + plot_peak_pw + plot_peak_aw

# Average
# 6.
# Annual average wildfire wf pm2.5

# Non-weighted 

annual <- pm_metrics_yearly %>% 
  group_by(year) %>% 
  summarize(ann_wfpm_avg = mean(ann_wfpm_avg)) 

plot_average <- ggplot(annual, aes(x=year, y=ann_wfpm_avg, fill = ann_wfpm_avg)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Year") +
  ylab(expression("Average annual mean wildfire PM"[2.5])) + 
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver"))

# Population-weighted 

average_pw <- pm_metrics_yearly %>% 
  group_by(year) %>% 
  summarize(pw_avg = sum((pop_prop * ann_wfpm_avg)))

plot_average_pw <- ggplot(average_pw, aes(x=year, y=pw_avg, fill = pw_avg)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Year") +
  ylab(expression("Population-weighted average annual mean wildfire PM"[2.5])) + 
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver"))

# Areal-weighted 

average_aw <- pm_metrics_yearly %>% 
  group_by(year) %>% 
  summarize(aw_avg = sum((area_prop * ann_wfpm_avg)))

plot_average_aw <- ggplot(average_aw, aes(x=year, y=aw_avg, fill = aw_avg)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Year") +
  ylab(expression("Area-weighted average annual mean wildfire PM"[2.5])) + 
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver"))

# Test the three different weighting schemes 

plot_average + plot_average_pw + plot_average_aw

# Create layout and graphs to print 

layout <- 'ABCDE'

# Plot final output

wrap_plots(A = plot_freq, B = plot_non_zero,C=plot_peak_pw,D=plot_sw, E=plot_average_pw, design = layout)
  
ggsave("fig1_panelB.png", dpi=300, height=5, width=20, units="in")




