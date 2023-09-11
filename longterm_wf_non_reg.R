#Code for running all analyses in CA long-term wildfire metric paper
#Non-regression code
#Updated 07 Sept 2023

#libraries
library(sf)
library(tidyverse)
library(GGally)
library(lme4)
library(maps)
library(tidycensus)
library(patchwork)
library(MetBrewer)
library(naniar)
library(biscale)
library(patchwork)
library(cowplot)
library(dotwhisker)
library(calendR)
library(openair)
library(mgcv)
library(splines)
library(marginaleffects)
library(glue)
library(viridis)
library(spdep)
library(skimr)
options(tigris_use_cache = TRUE)

#data set up
#CES 3--uses 2010 census tract boundaries
ces3 <- st_read("longterm-pm/calenviroscreen-3.0-shapefile/CES3Results.shp")

#CES 4--uses 2010 census tract boundaries
ces4 <- st_read("longterm-pm/calenviron/calenviroscreen40shpf2021shp/CES4 Final Shapefile.shp")

glimpse(ces3)
glimpse(ces4)

#restrict to just sub-scales and tracts and counties
ces4 <- ces4 %>% dplyr::select(Tract, County, CIscore, CIscoreP, PollBurd, PolBurdSc, PolBurdP, PopChar,PopCharSc,PopCharP, PM2_5, PM2_5_P,Shape_Area)
ces3 <- ces3 %>% dplyr::select(Tract_1, Percentile, CIscore, CIscoreP, PollutionS,Poll_pctl,PopCharSco,Pop_Pctl, Shape_Area, PM2_5, PM_2_5_Pct)                               

#add years that we will use in the study: CalEnvironScreen 3 (2006-2012) and 4.0 (2013-2019)
ces3 <- ces3 %>% mutate(ces_year="ces3")
ces4 <- ces4 %>% mutate(ces_year="ces4")

#so they have the same tract id name
ces3 <- ces3 %>% dplyr::rename( "Tract"="Tract_1")

#What are these scores? CIscore = total score; PollutionS = pollution score; PopCharSco = population charact. score

#Bring in race/ethnicity data
race <- read_csv("longterm-pm/race_eth_2010_ct_ca_counts.csv")
#numeric tract variable to match ces
race <- race %>% mutate(Tract=as.numeric(tractid))

#Add race to both ces score dataframes, we will conduct separate analysis from 2006-2012 and 2013-2019
ces3 <- left_join(ces3, race, by=c("Tract"="Tract"))
ces4 <- left_join(ces4, race, by=c("Tract"="Tract"))

#what is missing from CES?
summary(ces3$CIscoreP)
summary(ces4$CIscoreP)

#summartize counts of -999
ces3 <- ces3 %>% mutate(missing_ces=ifelse(CIscoreP==-999,1,0))
ces4 <- ces4 %>% mutate(missing_ces=ifelse(CIscoreP==-999,1,0))

#put ces3 and ces4 together
glimpse(ces3)
glimpse(ces4)
#rename columns in ces3 to match ces4
ces3 <- ces3 %>% dplyr::select(-Percentile)
ces4 <- ces4 %>% dplyr::select(-County)
ces4 <- ces4 %>% dplyr::select(-PollBurd)
ces4 <- ces4 %>% dplyr::select(-PopChar)

ces3 <- ces3 %>% dplyr::rename( "PolBurdSc"="PollutionS")
ces3 <- ces3 %>% dplyr::rename( "PolBurdP"="Poll_pctl")
ces3 <- ces3 %>% dplyr::rename( "PopCharSc"="PopCharSco")
ces3 <- ces3 %>% dplyr::rename( "PopCharP"="Pop_Pctl")
ces3 <- ces3 %>% dplyr::rename( "PM2_5_P"="PM_2_5_Pct")
summary(ces3$CIscore)
summary(ces4$CIscore)
n_distinct(ces4$Tract)

ces <- rbind(ces3, ces4)

#How many of missing ces overlap with missing race/ethnicity
ces <- ces %>% mutate(missing_race=ifelse(is.na(gisjoin)==T,1,0))
table(ces$missing_race, ces$missing_ces) #all those missing population are also missing CES
table(ces3$missing_ces, ces4$missing_ces) #93 of the same tracts, 10 from ces4, 13 from ces3

#Remove tracts that missing ces information 16070 units pre filter 
ces_missing <- ces %>% filter(missing_ces==1 | missing_race==1)
ces_missing <- ces_missing %>% dplyr::select(gisjoin)
ces_missing <- ces_missing %>% group_by(gisjoin) %>% slice(1)
ces_missing <- ces_missing %>% mutate(missing=1)
st_geometry(ces_missing) <- NULL

#Join missing tract identifier to main dataset
ces <- left_join(ces, ces_missing, by=c("gisjoin"="gisjoin"))
ces <- ces %>% mutate(missing=ifelse(is.na(missing),0,missing))

#Remove tracts that missing ces information 16070 units pre filter, 15838 post (drop N=116 tracts) 
table(ces$missing)

#save ces data for future use
st_write(ces, "longterm-pm/ces.shp")

#Bring in annual wildfire data
load("longterm-pm/wfpm_metrics_yearly_updated_09062023.RData")
wf <- pm_metrics_yearly
n_distinct(wf$geoid)
#To link annual wildfire exposures
wf <- wf %>% mutate(ces_year = case_when(year==2006 ~ "ces3",
                                         year==2007 ~ "ces3",
                                         year==2008 ~ "ces3",
                                         year==2009 ~ "ces3",
                                         year==2010 ~ "ces3",
                                         year==2011 ~ "ces3",
                                         year==2012 ~ "ces3",
                                         year==2013 ~ "ces4",
                                         year==2014 ~ "ces4",
                                         year==2015 ~ "ces4",
                                         year==2016 ~ "ces4",
                                         year==2017 ~ "ces4",
                                         year==2018 ~ "ces4",
                                         year==2019 ~ "ces4",
                                         year==2020 ~ "ces4")
)

wf <- wf %>% mutate(geoid=as.numeric(geoid))
wf <- left_join(wf, ces, by=c("geoid"="Tract", "ces_year"="ces_year"))

#Remove missing tracts
wf <- wf %>% filter(missing!=1)
n_distinct(wf$geoid) #7919 is correct

#Add rural/indigenous
tribal <- read_csv("longterm-pm/tl_2019_us_aiannh copy/tribal_tract.csv")
wf <- left_join(wf, tribal, by = c("geoid"="Tract"))

####Supplementary Figure 1
#Map of urban, rural, and tribal tracts
ces4 <- left_join(ces4, tribal,  by = c("Tract"="Tract"))
ces4 <- ces4 %>% mutate(tract_cat = case_when((rural==0 & indigenous==0) ~ "Urban alone",
                                              (rural==0 & indigenous==1) ~ "Urban Tribal",
                                              (rural==1 & indigenous==0) ~ "Rural alone",
                                              (rural==1 & indigenous==1) ~ "Rural Tribal"))

ces4 %>% 
  ggplot() + 
  geom_sf( aes(fill=tract_cat), color=NA) + 
  scale_fill_manual("Tract category", values=rev(met.brewer("Isfahan1", 4)), na.value="grey50")+  
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) 
ggsave("longterm-pm/urban-rural/tract-map.png", dpi=300, height=6, width=8, units="in" )


####Supplementary Figure 2
# ###Summary across the study period###
wf_reg <- readRDS(file = "longterm-pm/analysis/regression/regression_data.rds")

wf_2 <- wf_reg %>%
  group_by(geoid) %>%
  mutate(
    total_5ug = sum(pm_freq),
    total_non_zero = sum(non_zero_days),
    mean_peak_week = mean(peak_pm),
    mean_wf_pm = mean(ann_wfpm_avg),
    total_smoke_waves = sum(smoke_waves))
total_exposure <- wf_2 %>% group_by(geoid) %>% slice(1)

total_exposure <-
  total_exposure %>% dplyr::select(total_5ug,
                                   total_non_zero,
                                   mean_peak_week,
                                   mean_wf_pm,
                                   total_smoke_waves,
                                   geoid)

total_exposure <- left_join(total_exposure, tribal, by = c("geoid"="Tract") )
st_geometry(total_exposure) <- NULL

total_exposure_box <- total_exposure %>% dplyr::select(total_5ug,
                                                       total_non_zero,
                                                       mean_peak_week,
                                                       mean_wf_pm,
                                                       total_smoke_waves,
                                                       tract_cat) %>%
  pivot_longer(!tract_cat, names_to = "category", values_to = "level")
total_exposure_box <- total_exposure_box %>% dplyr::filter(!str_detect(category, "geoid"))

total_exposure_box$tract_cat <- factor(total_exposure_box$tract_cat, levels=c("Urban alone", "Urban Tribal", "Rural alone", "Rural Tribal"))
box1<-total_exposure_box %>% dplyr::filter(category=="total_5ug") %>% ggplot(aes(y = level, fill = tract_cat)) +
  geom_boxplot() +  scale_fill_manual("Tract category",
                                      values = rev(met.brewer("Isfahan1", 4)),
                                      na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        text=element_text(family="Gulliver")) + scale_y_continuous("Exposure value") +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3))

box1<-total_exposure_box %>% dplyr::filter(category=="total_5ug") %>% ggplot(aes(y = level, fill = tract_cat)) +
  geom_boxplot() +  scale_fill_manual("Tract category",
                                      values = rev(met.brewer("Isfahan1", 4)),
                                      na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        text=element_text(family="Gulliver")) + scale_y_continuous("Exposure value") +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3))

box1<-total_exposure_box %>% dplyr::filter(category=="total_5ug") %>% ggplot(aes(y = level, fill = tract_cat)) +
  geom_boxplot() +  scale_fill_manual("Tract category",
                                      values = rev(met.brewer("Isfahan1", 4)),
                                      na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        text=element_text(family="Gulliver"),
        legend.text=element_text(size=rel(1.5)),
        legend.title = element_text(size=rel(1.5)))+ scale_y_continuous("Exposure value") +
  ggtitle(bquote(Days~with~wildfire~PM[2.5]~">5"~"µg/m"^3))

box2<-total_exposure_box %>% dplyr::filter(category=="total_non_zero") %>% ggplot(aes(y = level, fill = tract_cat)) +
  geom_boxplot() +  scale_fill_manual("Tract category",
                                      values = rev(met.brewer("Isfahan1", 4)),
                                      na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        text=element_text(family="Gulliver"),
        legend.text=element_text(size=rel(1.5)),
        legend.title = element_text(size=rel(1.5)))+ scale_y_continuous("Exposure value") +
  ggtitle(bquote(Weeks~with~wildfire~PM[2.5]~">0"))

box3 <- total_exposure_box %>% dplyr::filter(category=="mean_peak_week") %>% ggplot(aes(y = level, fill = tract_cat)) +
  geom_boxplot() +  scale_fill_manual("Tract category",
                                      values = rev(met.brewer("Isfahan1", 4)),
                                      na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        text=element_text(family="Gulliver"),
        legend.text=element_text(size=rel(1.5)),
        legend.title = element_text(size=rel(1.5)))+ scale_y_continuous("Exposure value") +
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5]))

box4<-total_exposure_box %>% dplyr::filter(category=="total_smoke_waves") %>% ggplot(aes(y = level, fill = tract_cat)) +
  geom_boxplot() +  scale_fill_manual("Tract category",
                                      values = rev(met.brewer("Isfahan1", 4)),
                                      na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        text=element_text(family="Gulliver"),
        legend.text=element_text(size=rel(1.5)),
        legend.title = element_text(size=rel(1.5)))+ scale_y_continuous("Exposure value") +
  ggtitle("Smoke waves")

box5<-total_exposure_box %>% dplyr::filter(category=="mean_wf_pm") %>% ggplot(aes(y = level, fill = tract_cat)) +
  geom_boxplot() +  scale_fill_manual("Tract category",
                                      values = rev(met.brewer("Isfahan1", 4)),
                                      na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        text=element_text(family="Gulliver"),
        legend.text=element_text(size=rel(1.5)),
        legend.title = element_text(size=rel(1.5)))+ scale_y_continuous("Exposure value") +
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5]))

combo_box <- box1 + box2 + box3 + box4 + box5 + guide_area() + plot_layout(ncol = 3) + plot_layout(guides = "collect")
combo_box + plot_annotation(tag_levels = 'A')
ggsave("longterm-pm/urban-rural/exposure_boxplots.png", dpi=300, height=8, width=11, units="in" )

####Figure 2 -- map of wildfire PM exposure over the entire period, Panel A
##Need to use all census tracts, including those with missing data
ces <- st_read("longterm-pm/ces.shp")
wf <- readRDS("longterm-pm/wfpm_metrics_yearly.RData")
n_distinct(wf$geoid)
#To link annual wildfire exposures
wf <- wf %>% mutate(ces_year = case_when(year==2006 ~ "ces3",
                                         year==2007 ~ "ces3",
                                         year==2008 ~ "ces3",
                                         year==2009 ~ "ces3",
                                         year==2010 ~ "ces3",
                                         year==2011 ~ "ces3",
                                         year==2012 ~ "ces3",
                                         year==2013 ~ "ces4",
                                         year==2014 ~ "ces4",
                                         year==2015 ~ "ces4",
                                         year==2016 ~ "ces4",
                                         year==2017 ~ "ces4",
                                         year==2018 ~ "ces4",
                                         year==2019 ~ "ces4",
                                         year==2020 ~ "ces4")
)

wf <- wf %>% mutate(geoid=as.numeric(geoid))
wf <- left_join(wf, ces, by=c("geoid"="Tract", "ces_year"="ces_yer"))


wf_2 <- wf %>% 
  group_by(geoid) %>%
  mutate(
    total_5ug = sum(pm_freq),
    total_non_zero = sum(non_zero_days),
    mean_peak_week = mean(peak_pm),
    mean_wf_pm = mean(ann_wfpm_avg),
    total_smoke_waves = sum(smoke_waves))
total_exposure <- wf_2 %>% group_by(geoid) %>% slice(1)

total_exposure <-
  total_exposure %>% dplyr::select(total_5ug,
                                   total_non_zero,
                                   mean_peak_week,
                                   mean_wf_pm,
                                   total_smoke_waves,
                                   geoid, missing)

#Adding spatial data back to main dataframe
wf_spatial <- wf %>% dplyr::select(geoid, geometry)
wf_spatial <- wf_spatial %>% group_by(geoid) %>% slice(1)
total_exposure <- left_join(wf_spatial, total_exposure)
total_exposure <- total_exposure %>% group_by(geoid) %>% slice(1)

#replace wildfire values with 0 if missing = 1
total_exposure <- total_exposure  %>%
  mutate_at(vars(total_5ug,
                 total_non_zero,
                 mean_peak_week,
                 mean_wf_pm,
                 total_smoke_waves), ~ifelse(missing==1, NA, .))

####MAPS OF EXPOSURE WITH MISSING VALUES SO 8089 TRACTS INCLUDED
study_period_map1 <- total_exposure %>% 
  ggplot() + 
  geom_sf( aes(fill=total_5ug, geometry = geometry), color=NA) + 
  scale_fill_gradientn("N weeks", 
                       colors=met.brewer("OKeeffe2"), na.value="grey50",
                       breaks=c(1,50,100, 150), limits=c(0,150))+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = c(0.75, 0.75)) +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10)) +
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3))+
  theme(       plot.title = element_text(size = 12))


study_period_map2 <- total_exposure %>% 
  ggplot() + 
  geom_sf( aes(fill=total_non_zero, geometry = geometry), color=NA) + 
  scale_fill_gradientn("N days", colors=met.brewer("OKeeffe2"), na.value="grey50")+
  # trans="log1p", breaks=c(0, 100, 300, 1000))+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = c(0.75, 0.75)) +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10)) +
  ggtitle(bquote(Wildfire~PM[2.5]~">0"))+ 
  theme(       plot.title = element_text(size = 12))


study_period_map3 <- total_exposure %>% 
  ggplot() + 
  geom_sf( aes(fill=mean_peak_week, geometry = geometry), color=NA) + 
  scale_fill_gradientn(bquote(µg/m^3),   
                       colors=met.brewer("OKeeffe2"), na.value="grey50",
                       breaks=c(1,75,150, 225)
  )+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = c(0.75, 0.75)) +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10))+
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5]))+ 
  theme(       plot.title = element_text(size = 12))

study_period_map4 <- total_exposure %>% 
  ggplot() + 
  geom_sf( aes(fill=total_smoke_waves, geometry = geometry), color=NA) + 
  scale_fill_gradientn("N smoke \nwaves", colors=met.brewer("OKeeffe2"), na.value="grey50"
  )+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = c(0.75, 0.75)) +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10)) +
  ggtitle("Smoke waves")+ 
  theme(       plot.title = element_text(size = 12))

study_period_map5 <- total_exposure %>% 
  ggplot() + 
  geom_sf(aes(fill=mean_wf_pm, geometry = geometry), color=NA) + 
  scale_fill_gradientn(bquote(µg/m^3), 
                       colors=met.brewer("OKeeffe2"), na.value="grey50"
  )+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = c(0.75, 0.75)) +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10)) +
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5]))+   
  theme(       plot.title = element_text(size = 12))

#Figure 1, Panel A
study_period_map1 + study_period_map2 + study_period_map3 + 
  study_period_map4 + study_period_map5 +  plot_layout(ncol = 5) +
  plot_annotation(tag_levels = 'A')
ggsave("longterm-pm/analysis/maps/fig1_panelA_v2.png", dpi=300, height=5, width=20, units="in" )

####Figure 2, Panel B
########## Frequency ###########
# How often exposed?
# 1.
# Number of weeks with wildfires PM2.5 > 5 ug/m3
weekly_wf <- wf %>% 
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
non_zero_days <- wf %>% 
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
smoke_waves <- wf %>% 
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
peak <- wf %>% 
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


# Average
# 6.
# Annual average wildfire wf pm2.5
# Non-weighted 
annual <- wf %>% 
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

# Plot final output
plot_freq + plot_non_zero  + plot_peak + plot_sw+plot_average +   plot_layout(ncol = 5) +
  plot_annotation(tag_levels = 'A')
ggsave("longterm-pm/analysis/maps/fig2_panelB_v2.png", dpi=300, height=5, width=20, units="in")

####Supplementary Figure 3
study_period_map1_sup <- wf %>% 
  ggplot() + 
  geom_sf( aes(fill=pm_freq, geometry = geometry), color=NA) + 
  scale_fill_gradientn("N weeks", 
                       colors=met.brewer("OKeeffe2"), na.value="grey50")+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = "right") +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10)) +
  facet_wrap(~year)
study_period_map1_sup
ggsave("longterm-pm/analysis/maps/supfig1_panelA.png", dpi=300, height=10, width=7, units="in" )

study_period_map2_sup <- wf %>% 
  ggplot() + 
  geom_sf( aes(fill=non_zero_days, geometry = geometry), color=NA) + 
  scale_fill_gradientn("N days", colors=met.brewer("OKeeffe2"), na.value="grey50",
                       trans="sqrt", breaks=c(0,25,100,200,300))+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = "right") +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10)) +
  facet_wrap(~year)
study_period_map2_sup
ggsave("longterm-pm/analysis/maps/supfig1_panelB.png", dpi=300, height=10, width=7, units="in" )

study_period_map3_sup <- wf %>% 
  ggplot() + 
  geom_sf( aes(fill=peak_pm, geometry = geometry), color=NA) + 
  scale_fill_gradientn(bquote(µg/m^3),   
                       colors=met.brewer("OKeeffe2"), na.value="grey50"
  )+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = "right") +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10))+
  facet_wrap(~year)
study_period_map3_sup
ggsave("longterm-pm/analysis/maps/supfig1_panelC.png", dpi=300, height=10, width=7, units="in" )


study_period_map4_sup <- wf %>% 
  ggplot() + 
  geom_sf( aes(fill=smoke_waves, geometry = geometry), color=NA) + 
  scale_fill_gradientn("N smoke \nwaves", colors=met.brewer("OKeeffe2"), na.value="grey50",
                       breaks=c(0,3,6,9,12)
  )+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = "right") +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10)) +
  facet_wrap(~year)
study_period_map4_sup
ggsave("longterm-pm/analysis/maps/supfig1_panelD.png", dpi=300, height=10, width=7, units="in" )

study_period_map5_sup <- wf %>% 
  ggplot() + 
  geom_sf(aes(fill=ann_wfpm_avg, geometry = geometry), color=NA) + 
  scale_fill_gradientn(bquote(µg/m^3), 
                       colors=met.brewer("OKeeffe2"), na.value="grey50",
                       trans="log1p", breaks=c(1,3,10,50)
  )+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = "right") +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        legend.text=element_text(size=10)) +
  facet_wrap(~year)
study_period_map5_sup
ggsave("longterm-pm/analysis/maps/supfig1_panelE2.png", dpi=300, height=10, width=7, units="in")

####Supplementary Figure 4
# BRING IN MONTHLY DATA
pm_metrics_monthly <- readRDS(here("data","wfpm_updated2023","aggregate","wfpm_metrics_monthly.RData"))
########## Frequency ###########
# How often exposed?
# 1.
# Number of weeks with wildfires PM2.5 > 5 ug/m3
weekly_wf <- pm_metrics_monthly %>% 
  group_by(month) %>% 
  summarize(pm_freq = sum(pm_freq))

plot_freq <- ggplot(weekly_wf, aes(x=month, y=pm_freq, fill = pm_freq)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Month") +
  ylab(expression("N census tract weeks")) +
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_discrete(limits=month.abb, expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver", size = 12)) +
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3)) 

########## Duration ###########
# How long exposed?
# 2. 
# Number of days with non-zero wildfire pm2.5
non_zero_days <- pm_metrics_monthly %>% 
  group_by(month) %>% 
  summarize(duration = sum(non_zero_days)) 

plot_non_zero <- ggplot(non_zero_days, aes(x=month, y=duration, fill = duration)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Month") +
  ylab(expression("N census tract days")) +
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_discrete(limits=month.abb, expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver", size = 12)) +
  ggtitle(bquote(Wildfire~PM[2.5]~">0")) 

# Number of smoke waves 
# 3. 
# Defined as >=2 days over the study area >15 ug/m3 of wildfire PM2.5 
smoke_waves <- pm_metrics_monthly %>% 
  group_by(month) %>% 
  summarize(smoke_wave = sum(smoke_waves)) 
plot_sw <- ggplot(smoke_waves, aes(x=month, y=smoke_wave, fill = smoke_wave)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Month") +
  ylab(expression("N census tract smoke waves")) +
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_discrete(limits=month.abb, expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver", size = 12)) +
  ggtitle(expression("Smoke waves"))

########## Concentration ###########
# To what level exposed?
# Peak 
# 4.
# Average daily wf pm2.5 during peak week 
peak <- pm_metrics_monthly %>% 
  group_by(month) %>% 
  summarize(peak_pm = mean(peak_pm)) 
plot_peak <- ggplot(peak, aes(x=month, y=peak_pm, fill = peak_pm)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Month") +
  ylab(expression("Average peak week mean wildfire PM"[2.5])) + 
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_discrete(limits=month.abb, expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver", size = 12)) +
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5]))

# Average
# 6.
# Annual average wildfire wf pm2.5
# Monthly
# Note, this will have the same values as the yearly values, but in a monthly time series 
monthly <- pm_metrics_monthly %>% 
  group_by(month) %>% 
  summarize(mon_wfpm_avg = mean(mon_wfpm_avg)) 
plot_average <- ggplot(monthly, aes(x=month, y=mon_wfpm_avg, fill = mon_wfpm_avg)) +
  geom_bar(stat = "identity") + 
  scale_fill_met_c(name = "OKeeffe2") + 
  xlab("Month") +
  ylab(expression("Average monthly mean wildfire PM"[2.5])) + 
  ggtitle("") + 
  theme_classic() + 
  guides(fill="none") +
  scale_x_discrete(limits=month.abb, expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(text=element_text(family="Gulliver", size = 12)) +
  ggtitle(expression("Average Monthly Mean Wildfire PM"[2.5]))

plot_freq + plot_non_zero + plot_peak_pw + plot_sw + plot_average_pw  +   plot_annotation(tag_levels = "A") +   plot_layout(ncol = 5)
ggsave("Supplement2.png", dpi=300, height=20, width=15, units="in")

####Supplementary Figure 5
##Correlations between the 5 at total study area
#for correlation plot
total_corr <- total_exposure %>% dplyr::select(total_5ug, total_non_zero, mean_peak_week, total_smoke_waves, mean_wf_pm)
total_corr <- total_corr[2:6]
total_corr<-ggcorr(total_corr, method=c("pairwise", "spearman"),
                   nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                   label_size = 3, limits=c(0.25,1), mid="#3B9AB2",
                   color = "grey50")

##Correlations between the 5 measures by year
wf_corr <- wf %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg)
wf_corr_2006 <- wf %>% filter(year==2006) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2007 <- wf %>% filter(year==2007) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2008 <- wf %>% filter(year==2008) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2009 <- wf %>% filter(year==2009) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2010 <- wf %>% filter(year==2010) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2011 <- wf %>% filter(year==2011) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2012 <- wf %>% filter(year==2012) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2013 <- wf %>% filter(year==2013) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2014 <- wf %>% filter(year==2014) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2015 <- wf %>% filter(year==2015) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2016 <- wf %>% filter(year==2016) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2017 <- wf %>% filter(year==2017) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2018 <- wf %>% filter(year==2018) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2019 <- wf %>% filter(year==2019) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
wf_corr_2020 <- wf %>% filter(year==2020) %>% dplyr::select(pm_freq, non_zero_days, peak_pm,smoke_waves, ann_wfpm_avg) 
annual_corr_2006<-ggcorr(wf_corr_2006, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2006")
annual_corr_2007<-ggcorr(wf_corr_2007, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2007")
annual_corr_2008<-ggcorr(wf_corr_2008, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2008")
annual_corr_2009<-ggcorr(wf_corr_2009, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2009")
annual_corr_2010<-ggcorr(wf_corr_2010, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2010")
annual_corr_2011<-ggcorr(wf_corr_2011, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2011")
annual_corr_2012<-ggcorr(wf_corr_2012, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2012")
annual_corr_2013<-ggcorr(wf_corr_2013, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2013")
annual_corr_2014<-ggcorr(wf_corr_2014, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2014")
annual_corr_2015<-ggcorr(wf_corr_2015, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2015")
annual_corr_2016<-ggcorr(wf_corr_2016, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2016")
annual_corr_2017<-ggcorr(wf_corr_2017, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2017")
annual_corr_2018<-ggcorr(wf_corr_2018, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2018")
annual_corr_2019<-ggcorr(wf_corr_2019, method=c("pairwise", "spearman"),
                         nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0.25,1),mid="#3B9AB2",
                         color = "grey50")+ labs(title="2019")
annual_corr_2020<-ggcorr(wf_corr_2020, method=c("pairwise", "spearman"),
                         nbreaks = 10, digits = 3, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
                         label_size = 3, limits=c(0,1), mid="#3B9AB2",
                         color = "grey50") + labs(title="2020")
annual_corr_2006 + annual_corr_2007 + annual_corr_2008 + annual_corr_2009 +
  annual_corr_2010 + annual_corr_2011 + annual_corr_2012 + annual_corr_2013 +
  annual_corr_2014 + annual_corr_2015 + annual_corr_2016 + annual_corr_2017+
  annual_corr_2018 + annual_corr_2019 + annual_corr_2020 + plot_layout(guides = "collect")
total_corr + annual_corr + plot_annotation(tag_levels = 'A')  

####Supplementary Figure 6
#Make race/ethnicity maps for supplement
ces_race <- ces %>%  
  dplyr::select(hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, Tract, ces_year) 
ces_race <- ces_race %>% pivot_longer(-c(Tract, ces_year, geometry), names_to="race", values_to="percent")
summary(ces_race$percent)
ces_race$race <- factor(ces_race$race, levels=c("hisp_p", "nha_p", "nhain_p", "nhb_p", "nh2ormr", "nhw_p"),
                        labels=c("Hispanic", "NH Asian", "NH American Indian", "NH Black", "NH 2+", "NH white"))

#CES by year
colours_manual <- c("#D3D3D3" , "#A89EB9", "#7E6A9F", "#553687") 
map_ces <- ces_map %>% 
  ggplot() + 
  geom_sf( aes(fill=CIscoreP), color=NA) + 
  scale_fill_gradientn("CES score", colors=colours_manual, na.value="grey50")+
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  facet_wrap(~ces_year) +
  theme(legend.position = "right") +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        plot.title = element_text(size = 12),
        legend.text=element_text(size=10))
map_ces
ggsave("longterm-pm/analysis/maps/ces_maps.png", dpi=300, height=5, width=8, units="in" )

#Label race/ethnicity
map_race <- ces_race %>% dplyr::filter(ces_year=="ces3") %>%
  ggplot() + 
  geom_sf( aes(fill=percent), color=NA) + 
  scale_fill_gradientn("Percent racial/ethnic\ncomposition", colors=rev(met.brewer("Hokusai1")), na.value="grey50")+
  facet_wrap(~race, nrow=3) +
  theme_void(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  theme(legend.position = "right") +
  theme(text=element_text(family="Gulliver"),  legend.title=element_text(size=11), 
        plot.title = element_text(size = 12),
        legend.text=element_text(size=10))
map_race
ggsave("longterm-pm/analysis/maps/race_maps.png", dpi=300, height=14, width=8, units="in" )

####Supplementary Figure 7
#Quartiles of CES score by year
#Remove missing tracts
wf <- wf %>% filter(missing!=1)
n_distinct(wf$geoid) #7919 is correct

wf <- wf %>% 
  group_by(year) %>% 
  mutate(ces_q = ntile(CIscorP, 4))

wf <- wf %>% 
  group_by(year) %>% 
  mutate(ces_bi = ifelse(ces_q==4,1,0))

wf_vars <- wf %>% dplyr::select(ces_bi, year, pm_freq, non_zero_days, smoke_waves, peak_pm, ann_wfpm_avg)
wf_vars <- pivot_longer(wf_vars, -c(ces_bi, year), names_to="exposure_type", values_to="exposure_level")
wf_vars_summary <- wf_vars %>% 
  group_by(ces_bi, year, exposure_type) %>% 
  summarize(p10 = quantile(exposure_level, 0.1),
            p25 = quantile(exposure_level, 0.25),
            p50 = quantile(exposure_level, 0.5),
            p75 = quantile(exposure_level, 0.75),
            p90 = quantile(exposure_level, 0.9))
wf_vars_summary <- pivot_longer(wf_vars_summary, -c(ces_bi, year, exposure_type), names_to="quantile", values_to="quantile_value")

wf_vars_summary %>% filter(exposure_type==5) 
wf_vars_summary$exposure_type <- factor(wf_vars_summary$exposure_type, 
                                        levels=c("pm_freq","non_zero_days","peak_pm", "smoke_waves","ann_wfpm_avg"),
                                        labels=c("1","2","3","4","5"))
ggplot(data=wf_vars_summary, aes(x=year, y=quantile_value, group=factor(ces_bi), fill=factor(ces_bi))) + 
  geom_bar(stat="identity", position = 'dodge') + 
  facet_grid(exposure_type~quantile, scales="free") +
  theme_classic(base_size = 12)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(
    text=element_text(family="Gulliver")) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) +
  ylab("Exposure value") + scale_x_continuous("", breaks=c(2010,2015,2020)) + 
  scale_fill_manual("CES quartile", values = c("#A89EB9","#553687"), labels=c( "Q1-Q3","Q4"))
ggsave("longterm-pm/analysis/supp_fig7.png", dpi=300, height=7, width=10, units="in")

####Supplementary Figure 8
#Remove missing tracts
wf <- wf %>% filter(missing!=1)
n_distinct(wf$geoid) #7919 is correct

#Generate mean CES score across two periods for mapping etc.
ces_score_mean <- wf %>%  
  group_by(geoid, ces_year) %>%
  mutate(ces_mean_p=mean(CIscorP, na.rm=TRUE)) %>%
  dplyr::select(geoid, ces_year, ces_mean_p)
ces_score_mean <- ces_score_mean %>% group_by(geoid) %>% slice(1)

wf <- wf %>% 
  group_by(geoid, year) %>%
  mutate(
    total_5ug = sum(pm_freq),
    total_non_zero = sum(non_zero_days),
    mean_peak_week = mean(peak_pm),
    mean_wf_pm = mean(ann_wfpm_avg),
    total_smoke_waves = sum(smoke_waves))
year_exposure <- wf %>% dplyr::select(geoid, pm_freq, non_zero_days, peak_pm, ann_wfpm_avg, smoke_waves, year, geometry)
year_exposure <- left_join(year_exposure, ces_score_mean)

# create classes for weeks over 5 and ces
year_exposure_2006 <- year_exposure %>% filter(year==2006)
year_exposure_2007 <- year_exposure %>% filter(year==2007)
year_exposure_2008 <- year_exposure %>% filter(year==2008)
year_exposure_2009 <- year_exposure %>% filter(year==2009)
year_exposure_2010 <- year_exposure %>% filter(year==2010)
year_exposure_2011 <- year_exposure %>% filter(year==2011)
year_exposure_2012 <- year_exposure %>% filter(year==2012)
year_exposure_2013 <- year_exposure %>% filter(year==2013)
year_exposure_2014 <- year_exposure %>% filter(year==2014)
year_exposure_2015 <- year_exposure %>% filter(year==2015)
year_exposure_2016 <- year_exposure %>% filter(year==2016)
year_exposure_2017 <- year_exposure %>% filter(year==2017)
year_exposure_2018 <- year_exposure %>% filter(year==2018)
year_exposure_2019 <- year_exposure %>% filter(year==2019)
year_exposure_2020 <- year_exposure %>% filter(year==2020)
#2006
fig_map_bi_2006 <- bi_class(year_exposure_2006, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)
#Legend for bivariate map
bi_legend_2006 <- bi_legend(pal = "PurpleOr",
                            dim = 4,
                            xlab = "Higher average wf pm",
                            ylab = "Higher CES score",
                            size = 8)

#Bivariate map
bi_plot_2006 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2006,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme() + theme(text=element_text(family="Gulliver")) + ggtitle("2006")



#2007
fig_map_bi_2007 <- bi_class(year_exposure_2007, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2007 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2007,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2007")


#2008
fig_map_bi_2008 <- bi_class(year_exposure_2008, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2008 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2008,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme() + theme(text=element_text(family="Gulliver")) + ggtitle("2008")


#2009
fig_map_bi_2009 <- bi_class(year_exposure_2009, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2009 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2009,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2009")



#2010
fig_map_bi_2010 <- bi_class(year_exposure_2010, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2010 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2010,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2010")

#2011
fig_map_bi_2011 <- bi_class(year_exposure_2011, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2011 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2011,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme() + theme(text=element_text(family="Gulliver")) + ggtitle("2011")


#2012
fig_map_bi_2012 <- bi_class(year_exposure_2012, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2012 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2012,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme() + theme(text=element_text(family="Gulliver")) + ggtitle("2012")


#2013
fig_map_bi_2013 <- bi_class(year_exposure_2013, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2013 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2013,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme() + theme(text=element_text(family="Gulliver")) + ggtitle("2013")


#2014
fig_map_bi_2014 <- bi_class(year_exposure_2014, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2014 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2014,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme() + theme(text=element_text(family="Gulliver")) + ggtitle("2014")

#2015
fig_map_bi_2015 <- bi_class(year_exposure_2015, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2015 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2015,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2015")

#2016
fig_map_bi_2016 <- bi_class(year_exposure_2016, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2016 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2016,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2016")



#2017
fig_map_bi_2017 <- bi_class(year_exposure_2013, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2017 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2017,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2017")

#2018
fig_map_bi_2018 <- bi_class(year_exposure_2018, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2018 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2018,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2018")

#2019
fig_map_bi_2019 <- bi_class(year_exposure_2019, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2019 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2019,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2019")


#2020
fig_map_bi_2020 <- bi_class(year_exposure_2020, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)
#Legend for bivariate map
bi_legend_2020 <- bi_legend(pal = "PurpleOr",
                            dim = 4,
                            xlab = "Higher concentration",
                            ylab = "Higher CES score",
                            size = 20)+ theme(text=element_text(family="Gulliver")) 

#Bivariate map
bi_plot_2020 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2020,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2020")

bi_plot_2006+ bi_plot_2007+ bi_plot_2008+ bi_plot_2009+bi_plot_2010+bi_plot_2011+bi_plot_2012+
  bi_plot_2013+bi_plot_2014+bi_plot_2015+bi_plot_2016+bi_plot_2017 + bi_plot_2018 + bi_plot_2019 +bi_plot_2020 + bi_legend_2020
ggsave("longterm-pm/analysis/maps/supp_fig8.png", dpi=300, height=36, width=20, units="in" )

####Figure 4
##Panel A
peak_violin <- wf %>% drop_na(peak_pm, ces_bi) %>%
  ggplot(aes(x=peak_pm, y=factor(ces_bi), fill=factor(ces_bi))) + 
  geom_violin() + 
  scale_x_continuous(bquote(Concentration~(µg/m^3)), trans="log1p", breaks=c(0, 1, 4, 12, 36, 108, 324))+
  scale_fill_manual("CES quartile", values = c("#A89EB9", "#553687"), labels=c("Q1-Q3", "Q4"),
                    guide = guide_legend(reverse = TRUE) ) +
  scale_y_discrete("") + xlab("Concentration") +
  facet_wrap(~year, ncol=1, strip.position="left") +
  theme_minimal(base_size=12)  +
  theme(
    axis.text = element_text(size = 12),
    axis.ticks.y=element_blank(),
    axis.text.y = element_blank(),
    plot.title=element_text(size=14))  + 
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5])) +   
  theme(text=element_text(family="Gulliver"))
peak_violin
ggsave("longterm-pm/analysis/reardon/fig4_panela_peak_week_v2.png", dpi=300, height=7.5, width=4, units="in" )

##Figure 4, Panel B
#2009
fig_map_bi_2009 <- bi_class(year_exposure_2009, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)

#Bivariate map
bi_plot_2009 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2009,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2009")


#2020
fig_map_bi_2020 <- bi_class(year_exposure_2020, x = peak_pm, y = ces_mean_p, style = "quantile", dim = 4)
#Legend for bivariate map
bi_legend_2020 <- bi_legend(pal = "PurpleOr",
                            dim = 4,
                            xlab = "Higher concentration",
                            ylab = "Higher CES score",
                            size = 20)+ theme(text=element_text(family="Gulliver")) 

#Bivariate map
bi_plot_2020 <- ggplot() +
  geom_sf(
    data = fig_map_bi_2020,
    aes(fill = bi_class, geometry = geometry),
    lwd=0,
    color=NA,
    show.legend = FALSE
  ) +
  bi_scale_fill(pal = "PurpleOr", dim = 4) +
  bi_theme()  + theme(text=element_text(family="Gulliver")) + ggtitle("2020")

bi_plot_2009 + theme(text=element_text(family="Gulliver"))
ggsave("longterm-pm/analysis/maps/fig4_panelb_bivariate_2009_peakweak_v2.png", dpi=300, height=6, width=6, units="in" )

bi_plot_2020
ggsave("longterm-pm/analysis/maps/fig4b_panelb_bivariate_2020_peakweak_v2.png", dpi=300, height=6, width=6, units="in" )

####Figure 4
quantile(wf$smoke_waves, probs = .999)
wave_seq = seq(0, 8, 1)

reardon_plot1 <- wf %>%
  dplyr::select(geoid, year, smoke_waves, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(year) %>% 
  mutate(smoke_cat = cut(smoke_waves, wave_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(smoke_cat2 = fct_explicit_na(smoke_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(smoke_cat2, year, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(year,smoke_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(smoke_cat2, year, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat <- factor(reardon_plot1_group$smoke_cat2)
glimpse(reardon_plot1_group)

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(smoke_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat2 <- unclass(reardon_plot1_group$smoke_cat2)
reardon_plot1_group_overall$smoke_cat2 <- unclass(reardon_plot1_group_overall$smoke_cat2)
#Replace high category as 0 because it's really 0
reardon_plot1_group <- reardon_plot1_group %>% mutate(smoke_cat2 = ifelse(smoke_cat2==9,0,smoke_cat2))
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(smoke_cat2 = ifelse(smoke_cat2==9,0,smoke_cat2))

#Distribution of California population plot
overall_percents <- overall_percents %>% mutate(prop_race = as.numeric(sprintf("%0.2f", round(mean_race, digits = 3))),
                                                year="Overall composition")
plot1<-overall_percents %>% 
  ggplot(aes(y=mean_race, x="", fill = factor(race) )) +
  geom_col() +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  xlab(NULL) +
  scale_y_continuous("Census tract racial/ethnic composition (%)", expand=c(0,0), limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  scale_x_discrete(expand=c(0,0)) 

#Overall average for the state across study period
#Across all years
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(year="2006-2020 average")
plot2<-reardon_plot1_group_overall %>% 
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  facet_wrap(~year,  nrow=1) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Smoke waves (n)",expand = c(0,0), limits=c(0,8), breaks=c(2,4,6,8)
                     )

#2016
plot3<- reardon_plot1_group %>% filter(year==2016) %>%
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Smoke waves (n)",expand = c(0,0), limits=c(0,8), breaks=c(2,4,6,8)
                     )
#2018
plot4<- reardon_plot1_group %>% filter(year==2018) %>%
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 12)  +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Smoke waves (n)",expand = c(0,0), limits=c(0,8), breaks=c(2,4,6,8)
                     )

#2020
plot5<- reardon_plot1_group %>% filter(year==2020) %>%
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Smoke waves (n)",expand = c(0,0), limits=c(0,8), breaks=c(2,4,6,8)
                     )

plot1 + plot2 + plot3  + plot4 + plot5 + guide_area()   + plot_layout(ncol=6, guides = "collect")
ggsave("longterm-pm/analysis/reardon/figure_3.png", dpi=300, height=5, width=20, units="in" )

####Supplementary Figure Figure 9--smoke waves
quantile(wf$smoke_waves, probs = .999)
wave_seq = seq(0, 8, 1)

reardon_plot1 <- wf %>%
  dplyr::select(geoid, year, smoke_waves, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(year) %>% 
  mutate(smoke_cat = cut(smoke_waves, wave_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(smoke_cat2 = fct_explicit_na(smoke_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(smoke_cat2, year, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(year,smoke_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(smoke_cat2, year, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat <- factor(reardon_plot1_group$smoke_cat2)
glimpse(reardon_plot1_group)

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(smoke_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat2 <- unclass(reardon_plot1_group$smoke_cat2)
reardon_plot1_group_overall$smoke_cat2 <- unclass(reardon_plot1_group_overall$smoke_cat2)
#Replace high category as 0 because it's really 0
reardon_plot1_group <- reardon_plot1_group %>% mutate(smoke_cat2 = ifelse(smoke_cat2==9,0,smoke_cat2))
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(smoke_cat2 = ifelse(smoke_cat2==9,0,smoke_cat2))


#Plot
reardon_allyears<- reardon_plot1_group %>%
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 12)  +
  facet_wrap(~year,  nrow=3) +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),  breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + 
  scale_x_continuous("Smoke waves (n)",expand = c(0,0))

reardon_allyears
ggsave("longterm-pm/analysis/reardon/supp_fig9.png", dpi=300, height=10, width=10, units="in" )

####Supplementary Figure 10
###Reardon annual average for race/ethnicity
quantile(wf$ann_wfpm_avg, probs = .999)
pm_seq = seq(0, 7, 1)

reardon_plot1 <- wf %>%
  dplyr::select(geoid, year, ann_wfpm_avg, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(year) %>% 
  mutate(pm_cat = cut(ann_wfpm_avg, pm_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(pm_cat, year, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(year,pm_cat), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(pm_cat = fct_explicit_na(pm_cat, na_level="0"))

#Overall percents
overall_percents <- reardon_plot1 %>% dplyr::group_by(race) %>%
  dplyr::summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(pm_cat, year, race) %>%
  dplyr::summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$pm_cat <- factor(reardon_plot1_group$pm_cat)
glimpse(reardon_plot1_group)
# reardon_plot1_group <- reardon_plot1_group %>% group_by(pm_cat, year) %>%
#   mutate(mean_prop=mean_p/sum(mean_p))

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(pm_cat, race) %>%
  dplyr::summarize(mean_p=mean(percent, na.rm=T))

#Fix highest category, should be 0
reardon_plot1_group$pm_cat <- unclass(reardon_plot1_group$pm_cat)
reardon_plot1_group <- reardon_plot1_group %>% mutate(pm_cat = ifelse(pm_cat==8,0,pm_cat))
reardon_plot1_group_overall$pm_cat <- unclass(reardon_plot1_group_overall$pm_cat)
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(pm_cat = ifelse(pm_cat==8,0,pm_cat))

#Distribution of California population plot
overall_percents <- overall_percents %>% mutate(prop_race = as.numeric(sprintf("%0.2f", round(mean_race, digits = 3))),
                                                year="Overall composition")
plot1<-overall_percents %>% 
  ggplot(aes(y=mean_race, x="", fill = factor(race) )) +
  geom_col() +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  xlab(NULL) +
  scale_y_continuous("Census tract racial/ethnic composition (%)", expand=c(0,0), limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  scale_x_discrete(expand=c(0,0)) 

#Overall average for the state across study period
#Across all years
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(year="2006-2020 average")
plot2<-reardon_plot1_group_overall %>% 
  ggplot(aes(x=pm_cat, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  facet_wrap(~year,  nrow=1) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'), expand=c(0,0),breaks=c(1,2,3,4,5,6,7)) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))



#2018
plot4<- reardon_plot1_group %>% filter(year==2018) %>%
  ggplot(aes(x=pm_cat, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),  breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'),expand = c(0,0), 
                                                                       breaks=c(1,2,3,4,5,6,7)
                     )

#2016
plot3<- reardon_plot1_group %>% filter(year==2016) %>%
  ggplot(aes(x=pm_cat-1, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'),expand = c(0,0), breaks=c(1,2,3,4,5,6,7),
                                                                       limits=c(0,7)) 

#2020
plot5<- reardon_plot1_group %>% filter(year==2020) %>%
  ggplot(aes(x=pm_cat, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(0,0,0,0, "pt"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'),expand = c(0,0),
                                                                       breaks=c(1,2,3,4,5,6,7)
                     )

plot1 + plot2 + plot3  + plot4 + plot5 + guide_area()   + plot_layout(ncol=6, guides = "collect") 
ggsave("longterm-pm/analysis/reardon/annual_avg_v2.png", dpi=300, height=5, width=20, units="in" )

###Reardon plot for race/ethnicity by smoke waves
quantile(wf$smoke_waves, probs = .999)
wave_seq = seq(0, 8, 1)

reardon_plot1 <- wf %>%
  dplyr::select(geoid, year, smoke_waves, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(year) %>% 
  mutate(smoke_cat = cut(smoke_waves, wave_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(smoke_cat2 = fct_explicit_na(smoke_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(smoke_cat2, year, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(year,smoke_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(smoke_cat2, year, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat <- factor(reardon_plot1_group$smoke_cat2)
glimpse(reardon_plot1_group)
# reardon_plot1_group <- reardon_plot1_group %>% group_by(pm_cat, year) %>%
#   mutate(mean_prop=mean_p/sum(mean_p))

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(smoke_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))
# reardon_plot1_group_overall <- reardon_plot1_group_overall %>% group_by(pm_cat) %>%
#   mutate(mean_prop=mean_p/sum(mean_p))

reardon_plot1_group$smoke_cat2 <- unclass(reardon_plot1_group$smoke_cat2)
reardon_plot1_group_overall$smoke_cat2 <- unclass(reardon_plot1_group_overall$smoke_cat2)
#Replace high category as 0 because it's really 0
reardon_plot1_group <- reardon_plot1_group %>% mutate(smoke_cat2 = ifelse(smoke_cat2==9,0,smoke_cat2))
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(smoke_cat2 = ifelse(smoke_cat2==9,0,smoke_cat2))

#Distribution of California population plot
overall_percents <- overall_percents %>% mutate(prop_race = as.numeric(sprintf("%0.2f", round(mean_race, digits = 3))),
                                                year="Overall composition")
plot1<-overall_percents %>% 
  na.omit() %>%  ggplot(aes(y=mean_race, x="", fill = factor(race) )) +
  geom_col() +
  facet_wrap(~year)+
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  xlab(NULL) +
  scale_y_continuous("Census tract racial/ethnic composition (%)", expand=c(0,0), limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  scale_x_discrete(expand=c(0,0)) 

#Overall average for the state across study period
#Across all years
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(year="2006-2020 average")
plot2<-reardon_plot1_group_overall %>% 
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  facet_wrap(~year,  nrow=1) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous("Smoke waves (n)", expand=c(0,0), breaks=c(2,4,6,8)) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))

#2018
plot4<- reardon_plot1_group %>% filter(year==2018) %>%
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +  facet_wrap(~year,  nrow=1) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),  breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Smoke waves (n)",expand = c(0,0), breaks=c(2,4,6,8)
                     )

#2016
plot3<- reardon_plot1_group %>% filter(year==2016) %>%
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Smoke waves (n)",expand = c(0,0), breaks=c(2,4,6,8), limits=c(0,8)
                     )

#2020
plot5<- reardon_plot1_group %>% filter(year==2020) %>%
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Smoke waves (n)",expand = c(0,0), breaks=c(2,4,6,8)
                     )

plot1 + plot2 + plot3  + plot4 + plot5 + guide_area()   + plot_layout(ncol=6, guides = "collect")
ggsave("longterm-pm/analysis/reardon/smoke_waves_v2.png", dpi=300, height=5, width=20, units="in" )

#####Avg during peak week
quantile(wf$peak_pm, probs = .999)
peak_seq = seq(0, 115, 5)

reardon_plot1 <- wf %>%
  dplyr::select(geoid, year, peak_pm, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(year) %>% 
  mutate(peak_cat = cut(peak_pm, peak_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(peak_cat2 = fct_explicit_na(peak_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(peak_cat2, year, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(year,peak_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(peak_cat2, year, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$peak_cat2 <- factor(reardon_plot1_group$peak_cat2)
glimpse(reardon_plot1_group)
# reardon_plot1_group <- reardon_plot1_group %>% group_by(pm_cat, year) %>%
#   mutate(mean_prop=mean_p/sum(mean_p))

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(peak_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))
# reardon_plot1_group_overall <- reardon_plot1_group_overall %>% group_by(pm_cat) %>%
#   mutate(mean_prop=mean_p/sum(mean_p))

#Multiplying by 10 because we grouped by three above when generating categories
reardon_plot1_group$peak_cat2 <- unclass(reardon_plot1_group$peak_cat2)*5
reardon_plot1_group_overall$peak_cat2 <- unclass(reardon_plot1_group_overall$peak_cat2)*5
#Replace highest group as 0 because it's really 0
reardon_plot1_group <- reardon_plot1_group %>% mutate(peak_cat2 = ifelse(peak_cat2==120,0,peak_cat2))
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(peak_cat2 = ifelse(peak_cat2==120,0,peak_cat2))

#Distribution of California population plot
overall_percents <- overall_percents %>% mutate(prop_race = as.numeric(sprintf("%0.2f", round(mean_race, digits = 3))),
                                                year="Overall composition")
plot1<-overall_percents %>% 
  ggplot(aes(y=prop_race, x="", fill = factor(race) )) +
  geom_col() +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  xlab(NULL) +
  scale_y_continuous("Census tract racial/ethnic composition (%)", expand=c(0,0), limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  scale_x_discrete(expand=c(0,0)) 

#Overall average for the state across study period
#Across all years
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(year="2006-2020 average")
plot2<-reardon_plot1_group_overall %>% 
  ggplot(aes(x=peak_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  facet_wrap(~year,  nrow=1) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'), expand=c(0,0), breaks=c(25,50,75,100)) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))

#2018
plot4<- reardon_plot1_group %>% filter(year==2018) %>%
  ggplot(aes(x=peak_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),  breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'),expand = c(0,0), breaks=c(25,50,75,100)
                     )

#2016
plot3<- reardon_plot1_group %>% filter(year==2016) %>%
  ggplot(aes(x=peak_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'),expand = c(0,0),  breaks=c(25,50,75,100), limits=c(1,115)
                     )

#2020
plot5<- reardon_plot1_group %>% filter(year==2020) %>%
  ggplot(aes(x=peak_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'),expand = c(0,0), breaks=c(25,50,75,100)
                     )

plot1 + plot2 + plot3  + plot4 + plot5 + guide_area()   + plot_layout(ncol=6, guides = "collect")
ggsave("longterm-pm/analysis/reardon/peak_week_v2.png", dpi=300, height=5, width=20, units="in" )

#####Non-zero days
wf %>% ggplot() + geom_histogram(aes(non_zero_days)) + facet_wrap(~year)
summary(wf$non_zero_days)
quantile(wf$non_zero_days, probs = .999)
zero_seq = seq(0, 128, 8)

reardon_plot1 <- wf %>%
  dplyr::select(geoid, year, non_zero_days, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(year) %>% 
  mutate(zero_cat = cut(non_zero_days, zero_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(zero_cat2 = fct_explicit_na(zero_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(zero_cat2, year, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(year,zero_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(zero_cat2, year, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$zero_cat2 <- factor(reardon_plot1_group$zero_cat2)
glimpse(reardon_plot1_group)
# reardon_plot1_group <- reardon_plot1_group %>% group_by(pm_cat, year) %>%
#   mutate(mean_prop=mean_p/sum(mean_p))

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(zero_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))
# reardon_plot1_group_overall <- reardon_plot1_group_overall %>% group_by(pm_cat) %>%
#   mutate(mean_prop=mean_p/sum(mean_p))

#Multiplying by 4 because we grouped by three above when generating categories
reardon_plot1_group$zero_cat2 <- unclass(reardon_plot1_group$zero_cat2)*8
reardon_plot1_group_overall$zero_cat2 <- unclass(reardon_plot1_group_overall$zero_cat2)*8
#Replace high category as 0 because it's really 0
reardon_plot1_group <- reardon_plot1_group %>% mutate(zero_cat2 = ifelse(zero_cat2==136,0,zero_cat2))
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(zero_cat2 = ifelse(zero_cat2==136,0,zero_cat2))

#Distribution of California population plot
overall_percents <- overall_percents %>% mutate(prop_race = as.numeric(sprintf("%0.2f", round(mean_race, digits = 3))),
                                                year="Overall composition")
plot1<-overall_percents %>% 
  ggplot(aes(y=prop_race, x="", fill = factor(race) )) +
  geom_col() +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  xlab(NULL) +
  scale_y_continuous("Census tract racial/ethnic composition (%)", expand=c(0,0), limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  scale_x_discrete(expand=c(0,0)) 

#Overall average for the state across study period
#Across all years
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(year="2006-2020 average")
plot2<-reardon_plot1_group_overall %>% 
  ggplot(aes(x=zero_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  facet_wrap(~year,  nrow=1) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous("Days (n)", expand=c(0,0)) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))

#2018
plot4<- reardon_plot1_group %>% filter(year==2018) %>%
  ggplot(aes(x=zero_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),  breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Days (n)",expand = c(0,0)
                     )

#2016
plot3<- reardon_plot1_group %>% filter(year==2016) %>%
  ggplot(aes(x=zero_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Days (n)",expand = c(0,0), breaks=c(0,25,50,75,100,125), limits=c(0,128)
                     )

#2020
plot5<- reardon_plot1_group %>% filter(year==2020) %>%
  ggplot(aes(x=zero_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Days (n)",expand = c(0,0)
                     )

plot1 + plot2 + plot3  + plot4 + plot5 + guide_area()   + plot_layout(ncol=6, guides = "collect")
ggsave("longterm-pm/analysis/reardon/non_zero_v2.png", dpi=300, height=5, width=20, units="in" )

#####Frequency of weeks over 5
wf %>% ggplot() + geom_histogram(aes(pm_freq)) + facet_wrap(~year)
summary(wf$pm_freq)
quantile(wf$pm_freq, probs = .999)
freq_seq = seq(0, 12, 1)

reardon_plot1 <- wf %>%
  dplyr::select(geoid, year, pm_freq, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(year) %>% 
  mutate(zero_cat = cut(pm_freq, freq_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(zero_cat2 = fct_explicit_na(zero_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(zero_cat2, year, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(year,zero_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(zero_cat2, year, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$zero_cat2 <- factor(reardon_plot1_group$zero_cat2, levels=c("0" ,"(0,1]" ,  "(1,2]",  "(2,3] " , "(3,4]", "(4,5]", "(5,6]",   "(6,7]",   "(7,8]",  "(8,9]",  "(9,10]", "(10,11]", "(11,12]"))
glimpse(reardon_plot1_group)

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(zero_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))
# reardon_plot1_group_overall <- reardon_plot1_group_overall %>% group_by(pm_cat) %>%
#   mutate(mean_prop=mean_p/sum(mean_p))

#Unclass and subtract one so we have a 0 category
reardon_plot1_group$zero_cat2 <- unclass(reardon_plot1_group$zero_cat2)
table(reardon_plot1_group$zero_cat2)
reardon_plot1_group <- reardon_plot1_group %>% mutate(zero_cat2 = zero_cat2 - 1)
reardon_plot1_group_overall$zero_cat2 <- unclass(reardon_plot1_group_overall$zero_cat2)
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(zero_cat2 = zero_cat2 - 1)


#Distribution of California population plot
overall_percents <- overall_percents %>% mutate(prop_race = as.numeric(sprintf("%0.2f", round(mean_race, digits = 3))),
                                                year="Overall composition")
plot1<-overall_percents %>% 
  ggplot(aes(y=prop_race, x="", fill = factor(race) )) +
  geom_col() +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  xlab(NULL) +
  scale_y_continuous("Census tract racial/ethnic composition (%)", expand=c(0,0), limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  scale_x_discrete(expand=c(0,0)) 

#Overall average for the state across study period
#Across all years
reardon_plot1_group_overall <- reardon_plot1_group_overall %>% mutate(year="2006-2020 average")
plot2<-reardon_plot1_group_overall %>% 
  ggplot(aes(x=zero_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  facet_wrap(~year,  nrow=1) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous("Weeks (n)", expand=c(0,0), breaks=c(3,6,9,12)) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))

#2018
plot4<- reardon_plot1_group %>% filter(year==2018) %>%
  ggplot(aes(x=zero_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),  breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Weeks (n)",expand = c(0,0), breaks=c(3,6,9,12)
                     )

#2016
plot3<- reardon_plot1_group %>% filter(year==2016) %>%
  ggplot(aes(x=zero_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Weeks (n)",expand = c(0,0), breaks=c(3,6,9,12), limits=c(0,12)
                     )

#2020
plot5<- reardon_plot1_group %>% filter(year==2020) %>%
  ggplot(aes(x=zero_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  facet_wrap(~year,  nrow=1) +
  theme_classic(base_size = 13)  +
  theme(axis.text = element_text(size = 12),legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), strip.text = element_text(size=14)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_y_continuous("", expand = c(0,0),breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100)) + scale_x_continuous("Weeks (n)",expand = c(0,0), breaks=c(3,6,9,12)
                     )

plot1 + plot2 + plot3  + plot4 + plot5 + guide_area()   + plot_layout(ncol=6, guides = "collect")
ggsave("longterm-pm/analysis/reardon/weeks_over5_v2.png", dpi=300, height=5, width=20, units="in" )

####Supplementary Figure 11
##Secondary analysis monthly reardon plots for 2020
wf_month <- readRDS("longterm-pm/long-term-pm-exposure/wfpm_metrics_yearmonth_aug2023.RData")

#Need to run code above to load up the ces dataframe 
wf_month <- wf_month %>% mutate(geoid=as.numeric(geoid))
wf_month <- left_join(wf_month, ces, by=c("geoid"="Tract"))
glimpse(wf_month)

#Remove missing tracts
n_distinct(wf_month$geoid) #>8000, need to subset to those from wf
wf_tracts <- wf %>% dplyr::select(geoid) %>% distinct() #7919 is correct
wf_tracts <- wf_tracts %>% mutate(keep=1)
wf_tracts <- wf_tracts %>% mutate(geoid=as.numeric(geoid))
wf_month <- left_join(wf_month, wf_tracts)
wf_month <- wf_month %>% filter(keep==1)
n_distinct(wf_month$geoid) #7919 is correct

####SHOW non-zero for  2020
wf_month <- wf_month %>% dplyr::select(geoid, year, month, non_zero_days, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr)
wf_month <- wf_month %>% filter(year==2020)
summary(wf_month$non_zero_days)
wf_month %>% group_by(month) %>% summarize(quant=quantile(non_zero_days, probs = .999))
wave_seq = seq(0, 31, 1)

reardon_plot1 <- wf_month %>%
  dplyr::select(geoid, month, non_zero_days, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(month) %>% 
  mutate(smoke_cat = cut(non_zero_days, wave_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

reardon_plot1 <- reardon_plot1 %>% mutate(other_race=ifelse(other_race<0,0,other_race))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(smoke_cat2 = fct_explicit_na(smoke_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(smoke_cat2, month, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(month,smoke_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(smoke_cat2, month, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat <- factor(reardon_plot1_group$smoke_cat2)
glimpse(reardon_plot1_group)

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(smoke_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat2 <- unclass(reardon_plot1_group$smoke_cat2)
reardon_plot1_group_overall$smoke_cat2 <- unclass(reardon_plot1_group_overall$smoke_cat2)
reardon_plot1_group <- reardon_plot1_group %>% mutate(smoke_cat2 = ifelse(smoke_cat==0, 0,smoke_cat2))

#Overall average for the state across study period
#Across all years
plot2<-reardon_plot1_group %>% 
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 12)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  facet_wrap(~month,  nrow=3) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous("Days (n)", expand=c(0,0)) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))
plot2

ggsave("longterm-pm/analysis/reardon/supp_fig11_monthly_non0days_2020.png", dpi=300, height=5, width=10, units="in" )

##Weeks over 5ug/m3
wf_month <- wf_month %>% dplyr::select(geoid, year, month, pm_freq, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr)
wf_month <- wf_month %>% filter(year==2020)
wf_month %>% group_by(month) %>% summarize(quant=quantile(pm_freq, probs = .999))
wave_seq = seq(0, 5, 1)

reardon_plot1 <- wf_month %>%
  dplyr::select(geoid, month, pm_freq, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(month) %>% 
  mutate(smoke_cat = cut(pm_freq, wave_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

reardon_plot1 <- reardon_plot1 %>% mutate(other_race=ifelse(other_race<0,0,other_race))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(smoke_cat2 = fct_explicit_na(smoke_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(smoke_cat2, month, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(month,smoke_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(smoke_cat2, month, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat <- factor(reardon_plot1_group$smoke_cat2)

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(smoke_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat2 <- unclass(reardon_plot1_group$smoke_cat2)
reardon_plot1_group_overall$smoke_cat2 <- unclass(reardon_plot1_group_overall$smoke_cat2)
reardon_plot1_group <- reardon_plot1_group %>% mutate(smoke_cat2 = ifelse(smoke_cat==0, 0,smoke_cat2))

#Overall average for the state across study period
#Across all years
plot2<-reardon_plot1_group %>% 
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 12)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  facet_wrap(~month,  nrow=3) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous("Weeks (n)", expand=c(0,0)) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))
plot2
ggsave("longterm-pm/analysis/reardon/supp_fig11_monthly_weeksgt5_2020.png", dpi=300, height=5, width=10, units="in" )

#Keep only needed columns
wf_month <- wf_month %>% dplyr::select(geoid, year, month, peak_pm, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr)
wf_month <- wf_month %>% filter(year==2020)
wf_month %>% group_by(month) %>% summarize(quant=quantile(peak_pm, probs = .999))

wave_seq = seq(0, 132, 6)

reardon_plot1 <- wf_month %>%
  dplyr::select(geoid, month, peak_pm, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(month) %>% 
  mutate(smoke_cat = cut(peak_pm, wave_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

reardon_plot1 <- reardon_plot1 %>% mutate(other_race=ifelse(other_race<0,0,other_race))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(smoke_cat2 = fct_explicit_na(smoke_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(smoke_cat2, month, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(month,smoke_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(smoke_cat2, month, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat <- factor(reardon_plot1_group$smoke_cat2)

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(smoke_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))


reardon_plot1_group$smoke_cat2 <- unclass(reardon_plot1_group$smoke_cat2)*6
reardon_plot1_group_overall$smoke_cat2 <- unclass(reardon_plot1_group_overall$smoke_cat2)*6
reardon_plot1_group <- reardon_plot1_group %>% mutate(smoke_cat2 = ifelse(smoke_cat==0, 0,smoke_cat2))

#Distribution of California population plot
overall_percents <- overall_percents %>% mutate(prop_race = as.numeric(sprintf("%0.2f", round(mean_race, digits = 3))),
                                                year="Overall CA racial/ethnic composition")
plot1<-overall_percents %>% 
  na.omit() %>%  ggplot(aes(y=mean_race, x="", fill = factor(race) )) +
  geom_col() +
  theme_classic(base_size = 12
  )  +
  facet_wrap(~year)+
  theme(axis.text = element_text(size = 10)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"),  legend.title=element_text(size=13), 
        legend.text=element_text(size=11)
  ) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  xlab(NULL) +
  scale_y_continuous("Census tract racial/ethnic composition (%)", expand=c(0,0), limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  scale_x_discrete(expand=c(0,0)) 

#Overall average for the state across study period
#Across all years
plot2<-reardon_plot1_group %>% 
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 12)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  facet_wrap(~month,  nrow=3) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'), expand=c(0,0) 
  ) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))
plot2
ggsave("longterm-pm/analysis/reardon/supp_fig11_monthly_peakweek_2020.png", dpi=300, height=5, width=10, units="in" )

wf_month <- wf_month %>% dplyr::select(geoid, year, month, mon_wfpm_avg, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr)
wf_month <- wf_month %>% filter(year==2020)
wf_month %>% group_by(month) %>% summarize(quant=quantile(mon_wfpm_avg, probs = .99))
wave_seq = seq(0, 38, 2)

reardon_plot1 <- wf_month %>%
  dplyr::select(geoid, month, mon_wfpm_avg, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(month) %>% 
  mutate(smoke_cat = cut(mon_wfpm_avg, wave_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))

reardon_plot1 <- reardon_plot1 %>% mutate(other_race=ifelse(other_race<0,0,other_race))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(smoke_cat2 = fct_explicit_na(smoke_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(smoke_cat2, month, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(month,smoke_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(smoke_cat2, month, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat <- factor(reardon_plot1_group$smoke_cat2)
glimpse(reardon_plot1_group)

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(smoke_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))


reardon_plot1_group$smoke_cat2 <- unclass(reardon_plot1_group$smoke_cat2)*2
reardon_plot1_group_overall$smoke_cat2 <- unclass(reardon_plot1_group_overall$smoke_cat2)*2
reardon_plot1_group <- reardon_plot1_group %>% mutate(smoke_cat2 = ifelse(smoke_cat==0, 0,smoke_cat2))

#Distribution of California population plot
overall_percents <- overall_percents %>% mutate(prop_race = as.numeric(sprintf("%0.2f", round(mean_race, digits = 3))),
                                                year="Overall CA racial/ethnic composition")
plot1<-overall_percents %>% 
  na.omit() %>%  ggplot(aes(y=mean_race, x="", fill = factor(race) )) +
  geom_col() +
  theme_classic(base_size = 12
  )  +
  facet_wrap(~year)+
  theme(axis.text = element_text(size = 10)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"),  legend.title=element_text(size=13), 
        legend.text=element_text(size=11)
  ) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  xlab(NULL) +
  scale_y_continuous("Census tract racial/ethnic composition (%)", expand=c(0,0), limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  scale_x_discrete(expand=c(0,0)) 

#Overall average for the state across study period
#Across all years
plot2<-reardon_plot1_group %>% 
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 12)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  facet_wrap(~month,  nrow=3) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous(bquote('Wildfire PM'[2.5]~'('*µg~m^-3*')'), expand=c(0,0), breaks=c(10,20,30)) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))
plot2
ggsave("longterm-pm/analysis/reardon/supp_fig11_monthly_average.png", dpi=300, height=5, width=10, units="in" )


#smoke waves
wf_month <- wf_month %>% dplyr::select(geoid, year, month, smoke_waves, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr)
quantile(wf_month$smoke_waves, probs = .999)
wave_seq = seq(0, 3, 1)

reardon_plot1 <- wf_month %>%
  dplyr::select(geoid, month, smoke_waves, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr) %>%
  group_by(month) %>% 
  mutate(smoke_cat = cut(smoke_waves, wave_seq),
         other_race=100-(hisp_p + nhw_p + nhb_p + nha_p + nhain_p+  nh2ormr))
reardon_plot1 <- reardon_plot1 %>% mutate(other_race = ifelse(other_race<0, 0, other_race))

#Recode NA as 0 (which they are)
reardon_plot1 <- reardon_plot1 %>% mutate(smoke_cat2 = fct_explicit_na(smoke_cat, na_level="0"))

reardon_plot1 <- reardon_plot1 %>% dplyr::select(smoke_cat2, month, hisp_p, nhw_p, nhb_p, nha_p, nhain_p, nh2ormr, other_race) %>%
  pivot_longer(-c(month,smoke_cat2), names_to = "race", values_to="percent")

reardon_plot1$race <- factor(reardon_plot1$race , levels=c("nhw_p", "hisp_p", "nha_p", "nhb_p", "nh2ormr", "nhain_p", "other_race"),
                             labels=c("NH white", "Hispanic", "NH Asian", "NH Black", "NH 2+ races", "NH American Indian", "Other"))

#Overall percents
overall_percents <- reardon_plot1 %>% group_by(race) %>%
  summarize(mean_race=mean(percent, na.rm=T))

reardon_plot1_group <- reardon_plot1 %>% group_by(smoke_cat2, month, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat <- factor(reardon_plot1_group$smoke_cat2)
reardon_plot1_group_overall$smoke_cat <- factor(reardon_plot1_group_overall$smoke_cat2)

glimpse(reardon_plot1_group)

#Across all years
reardon_plot1_group_overall <- reardon_plot1 %>% group_by(smoke_cat2, race) %>%
  summarize(mean_p=mean(percent, na.rm=T))

reardon_plot1_group$smoke_cat2 <- unclass(reardon_plot1_group$smoke_cat2)
reardon_plot1_group_overall$smoke_cat2 <- unclass(reardon_plot1_group_overall$smoke_cat2)
table(reardon_plot1_group$smoke_cat2)
reardon_plot1_group <- reardon_plot1_group %>% mutate(smoke_cat2 = ifelse(smoke_cat==0, 0,smoke_cat2))

#Distribution of California population plot
overall_percents <- overall_percents %>% mutate(prop_race = as.numeric(sprintf("%0.2f", round(mean_race, digits = 3))),
                                                year="Overall CA racial/ethnic composition")
plot1<-overall_percents %>% 
  na.omit() %>%  ggplot(aes(y=mean_race, x="", fill = factor(race) )) +
  geom_col() +
  theme_classic(base_size = 12
  )  +
  facet_wrap(~year)+
  theme(axis.text = element_text(size = 10)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"),  legend.title=element_text(size=13), 
        legend.text=element_text(size=11)
  ) +
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  xlab(NULL) +
  scale_y_continuous("Census tract racial/ethnic composition (%)", expand=c(0,0), breaks=c(0,20,40,60,80,100)) +
  scale_x_discrete(expand=c(0,0)) 

#Overall average for the state across study period
#Across all years
plot2<-reardon_plot1_group %>% 
  ggplot(aes(x=smoke_cat2, y=mean_p, fill = factor(race))) +
  geom_area(stat="identity", position = "fill") +
  theme_classic(base_size = 12)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  facet_wrap(~month,  nrow=3) +
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"))+
  scale_fill_manual("Race/ethnicity", values=met.brewer("Archambault", 7)) +
  scale_x_continuous("Smoke waves (n)", expand=c(0,0), breaks=c(1,2,3,4)) +
  scale_y_continuous("", expand=c(0,0), breaks=c(0,.20,.40,.60,.80,1),
                     labels=c(0,20,40,60,80,100))
plot2
ggsave("longterm-pm/analysis/reardon/supp_fig11_monthly_smoke_waves.png", dpi=300, height=5, width=10, units="in" )

####Supplementary Figures 12--16 are in the regression code

######Supplementary Table 1
tribal <- read_csv("longterm-pm/tl_2019_us_aiannh copy/tribal_tract.csv")
tribal <- tribal %>% dplyr::select(Tract, tract_cat, pop_den_sq_mi, rural)


#Summary across the study period###
wf_2 <- wf_reg %>%
  group_by(geoid) %>%
  mutate(
    total_5ug = sum(pm_freq),
    total_non_zero = sum(non_zero_days),
    mean_peak_week = mean(peak_pm),
    mean_wf_pm = mean(ann_wfpm_avg),
    total_smoke_waves = sum(smoke_waves))
total_exposure <- wf_2 %>% group_by(geoid) %>% slice(1)

total_exposure <-
  total_exposure %>% dplyr::select(total_5ug,
                                   total_non_zero,
                                   mean_peak_week,
                                   mean_wf_pm,
                                   total_smoke_waves,
                                   geoid)

total_exposure <- left_join(total_exposure, tribal, by = c("geoid"="Tract") )

total_exposure %>%
  dplyr::group_by(tract_cat) %>%
  skim()

####Supplementary Table 2
#Let's look at percentiles of exposure
wf %>% group_by(year) %>% summarize(quant=quantile(ann_wfpm_avg, probs = .99))

summary(wf$peak_pm)
summary(wf$smoke_waves)

##RR for wildfire stuff
wf <- wf %>% mutate(ann_avg_exposed = ifelse(ann_wfpm_avg>=25,1,0))
rr_calc <- wf %>% dplyr::select(ann_avg_exposed, year, hisp, nhw, nhb, nha, nhain_p, nh2ormr_c, total_pop)
rr_calc <- rr_calc %>% mutate(other_race = total_pop- (hisp+ nhw+ nhb+ nha+ nhain_p+ nh2ormr_c))
rr_calc_long <- rr_calc %>% dplyr::select(ann_avg_exposed, year, hisp, nhw, nhb, nha, nhain_p, nh2ormr_c, total_pop, other_race) %>%
  pivot_longer(-c(year, ann_avg_exposed, total_pop), names_to = "race", values_to="count")

rr_hisp<-rr_calc_long %>% filter(race=="hisp") %>% group_by(year) %>%
  summarize(rr_hisp  =
              (sum(count[ann_avg_exposed == 1], na.rm=T) /
                 (sum(count, na.rm=T))) /
              (sum(total_pop[ann_avg_exposed==1], na.rm=T) /
                 sum(total_pop[ann_avg_exposed==0], na.rm=T)))

rr_nhw<-rr_calc_long %>% filter(race=="nhw") %>% group_by(year) %>%
  summarize(rr_nhw  =
              (sum(count[ann_avg_exposed == 1], na.rm=T) /
                 (sum(count, na.rm=T))) /
              (sum(total_pop[ann_avg_exposed==1], na.rm=T) /
                 sum(total_pop[ann_avg_exposed==0], na.rm=T)))

rr_ai<-rr_calc_long %>% filter(race=="nhain_p") %>% group_by(year) %>%
  summarize(rr_ai  =
              (sum(count[ann_avg_exposed == 1], na.rm=T) /
                 (sum(count, na.rm=T))) /
              (sum(total_pop[ann_avg_exposed==1], na.rm=T) /
                 sum(total_pop[ann_avg_exposed==0], na.rm=T)))        

rr_nhb<-rr_calc_long %>% filter(race=="nhb") %>% group_by(year) %>%
  summarize(rr_nhb  =
              (sum(count[ann_avg_exposed == 1], na.rm=T) /
                 (sum(count, na.rm=T))) /
              (sum(total_pop[ann_avg_exposed==1], na.rm=T) /
                 sum(total_pop[ann_avg_exposed==0], na.rm=T)))    

rr_nha<-rr_calc_long %>% filter(race=="nha") %>% group_by(year) %>%
  summarize(rr_nha  =
              (sum(count[ann_avg_exposed == 1], na.rm=T) /
                 (sum(count, na.rm=T))) /
              (sum(total_pop[ann_avg_exposed==1], na.rm=T) /
                 sum(total_pop[ann_avg_exposed==0], na.rm=T)))  

rr_nh2ormr_c<-rr_calc_long %>% filter(race=="nh2ormr_c") %>% group_by(year) %>%
  summarize(rr_nh2or  =
              (sum(count[ann_avg_exposed == 1], na.rm=T) /
                 (sum(count, na.rm=T))) /
              (sum(total_pop[ann_avg_exposed==1], na.rm=T) /
                 sum(total_pop[ann_avg_exposed==0], na.rm=T)))  

test<-rr_calc_long %>% filter((race=="nh2ormr_c" | race=="nhw") & year==2016)

all_rr_annual_avg <- left_join(rr_hisp,rr_nha )
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_ai)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nhb)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nh2ormr_c)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nhw)
write_csv(all_rr_annual_avg, "longterm-pm/analysis/rr_wfpm25_gt25ugm3.csv")

##Try continuous relative risk
##RR for wildfire stuff
rr_calc <- wf %>% dplyr::select(ann_wfpm_avg, year, hisp, nhw, nhb, nha, nhain_p, nh2ormr_c, total_pop)
rr_calc <- rr_calc %>% mutate(other_race = total_pop- (hisp+ nhw+ nhb+ nha+ nhain_p+ nh2ormr_c))
rr_calc_long <- rr_calc %>% dplyr::select(ann_wfpm_avg, year, hisp, nhw, nhb, nha, nhain_p, nh2ormr_c, total_pop, other_race) %>%
  pivot_longer(-c(year, ann_wfpm_avg, total_pop), names_to = "race", values_to="count")

rr_hisp <- rr_calc_long %>% filter(race=="hisp") %>% 
  group_by(year) %>%
  summarize(rr_hisp = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))

rr_nhw<-rr_calc_long %>% filter(race=="nhw") %>% group_by(year) %>%
  summarize(rr_nhw  =
              (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))

rr_ai<-rr_calc_long %>% filter(race=="nhain_p") %>% group_by(year) %>%
  summarize(rr_ai  = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))      

rr_nhb<-rr_calc_long %>% filter(race=="nhb") %>% group_by(year) %>%
  summarize(rr_nhb = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))

rr_nha<-rr_calc_long %>% filter(race=="nha") %>% group_by(year) %>%
  summarize(rr_nha  = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))

rr_nh2ormr_c<-rr_calc_long %>% filter(race=="nh2ormr_c") %>% group_by(year) %>%
  summarize(rr_nh2or = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))


all_rr_annual_avg <- left_join(rr_hisp,rr_nha )
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_ai)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nhb)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nh2ormr_c)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nhw)
write_csv(all_rr_annual_avg, "longterm-pm/analysis/rr_annualavgpm25.csv")

####Supplementary Table 3 -- is in the regression code 
