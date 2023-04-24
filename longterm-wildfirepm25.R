#Long-term wildfire PM2.5 paper script
#Edited by Joan Casey on 24 April 2023

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
options(tigris_use_cache = TRUE)

##CES data from https://oehha.ca.gov/calenviroscreen/report/calenviroscreen-40 and https://oehha.ca.gov/calenviroscreen/report/calenviroscreen-30
#CES 3--uses 2010 census tract boundaries
ces3 <- st_read("calenviroscreen-3.0-shapefile/CES3Results.shp")

#CES 4--uses 2010 census tract boundaries
ces4 <- st_read("calenviroscreen40shpf2021shp/CES4 Final Shapefile.shp")

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

#Bring in race/ethnicity data from NHGIS.org
race <- read_csv("race_eth_2010_ct_ca_counts.csv")
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
table(ces3$missing_ces) #106 = -999
table(ces4$missing_ces) #103 = -999

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
ces <- ces %>% dplyr::select(-missing_race, -missing_ces)
summary(ces$CIscore)
summary(ces$CIscoreP)

#replace -999 with NA
ces <- ces %>% replace_with_na(replace = list(CIscore = -999))
ces <- ces %>% replace_with_na(replace = list(CIscoreP = -999))

#Bring in annual wildfire data
wf <- readRDS("/Users/joancasey/Documents/Columbia/Grants/Wildfire R01/manuscripts/longterm-pm/wfpm_metrics_yearly.RData")
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
glimpse(wf)

#Maps of wildfire without missing data


wf <- wf %>% mutate(geoid=as.numeric(geoid))
wf <- left_join(wf, ces, by=c("geoid"="Tract", "ces_year"="ces_year"))
glimpse(wf)
table(wf$missing)


#Remove missing tracts
wf <- wf %>% filter(missing!=1)
n_distinct(wf$geoid) #7919 is correct

#Let's look at percentiles of exposure
wf %>% group_by(year) %>% summarize(quant=quantile(ann_wfpm_avg, probs = .99))


####Relative risk for exposure to wf smoke compared to race and ethnicity overall representation
##Continuous relative risk
rr_calc <- wf %>% dplyr::select(ann_wfpm_avg, year, hisp, nhw, nhb, nha, nhaian, nh2ormore_c, total_pop)
rr_calc <- rr_calc %>% mutate(other_race = total_pop- (hisp+ nhw+ nhb+ nha+ nhaian+ nh2ormore_c))
rr_calc_long <- rr_calc %>% dplyr::select(ann_wfpm_avg, year, hisp, nhw, nhb, nha, nhaian, nh2ormore_c, total_pop, other_race) %>%
  pivot_longer(-c(year, ann_wfpm_avg, total_pop), names_to = "race", values_to="count")

rr_hisp <- rr_calc_long %>% filter(race=="hisp") %>% 
  group_by(year) %>%
  summarize(rr_hisp = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))

rr_nhw<-rr_calc_long %>% filter(race=="nhw") %>% group_by(year) %>%
  summarize(rr_nhw  =
              (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))

rr_ai<-rr_calc_long %>% filter(race=="nhaian") %>% group_by(year) %>%
  summarize(rr_ai  = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))      

rr_nhb<-rr_calc_long %>% filter(race=="nhb") %>% group_by(year) %>%
  summarize(rr_nhb = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))

rr_nha<-rr_calc_long %>% filter(race=="nha") %>% group_by(year) %>%
  summarize(rr_nha  = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))

rr_nh2ormore_c<-rr_calc_long %>% filter(race=="nh2ormore_c") %>% group_by(year) %>%
  summarize(rr_nh2or = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop)))

all_rr_annual_avg <- left_join(rr_hisp,rr_nha )
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_ai)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nhb)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nh2ormore_c)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nhw)

##RR for annual exposure binary >=25ug/m3
wf <- wf %>% mutate(ann_avg_exposed = ifelse(ann_wfpm_avg>=25,1,0))
rr_calc <- wf %>% dplyr::select(ann_avg_exposed, year, hisp, nhw, nhb, nha, nhaian, nh2ormore_c, total_pop)
rr_calc <- rr_calc %>% mutate(other_race = total_pop- (hisp+ nhw+ nhb+ nha+ nhaian+ nh2ormore_c))
rr_calc_long <- rr_calc %>% dplyr::select(ann_avg_exposed, year, hisp, nhw, nhb, nha, nhaian, nh2ormore_c, total_pop, other_race) %>%
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

rr_ai<-rr_calc_long %>% filter(race=="nhaian") %>% group_by(year) %>%
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

rr_nh2ormore_c<-rr_calc_long %>% filter(race=="nh2ormore_c") %>% group_by(year) %>%
  summarize(rr_nh2or  =
              (sum(count[ann_avg_exposed == 1], na.rm=T) /
                 (sum(count, na.rm=T))) /
              (sum(total_pop[ann_avg_exposed==1], na.rm=T) /
                 sum(total_pop[ann_avg_exposed==0], na.rm=T)))  

test<-rr_calc_long %>% filter((race=="nh2ormore_c" | race=="nhw") & year==2016)

all_rr_annual_avg <- left_join(rr_hisp,rr_nha )
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_ai)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nhb)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nh2ormore_c)
all_rr_annual_avg <- left_join(all_rr_annual_avg, rr_nhw)


