
library(dplyr)
library(tidyr)
library(readr)


#Load data
load("wfpm25_CT_2006to2020_updated_Aug2023.RData")

wfpm <- wfpm25.updated
rm(wfpm25.updated)

names(wfpm)[3:4] <- c("ml_pm25", "wf_pm25_imp")
# summary(wfpm)
wfpm$year <- as.numeric(format(wfpm$date, "%Y"))

popn.dta <- readRDS("rr_est.RDS")
popn.dta <- popn.dta[,c("geoid", "year", "hisp", "nhw", "nhb", "nha", "nhaian", "nh2ormore", "total_pop")]
popn.dta <- popn.dta %>% mutate(other_race = total_pop - (hisp + nhw + nhb + nha + nhaian + nh2ormore))

names(popn.dta)

all.races <- c("hisp", "nhw", "nhb", "nha", "nhaian", "nh2ormore", "other_race")
all.years <- unique(wfpm$year)

### FIRST ESTIMATE POINT ESTIMATES IN FULL POPULATION
wf <- wfpm %>% 
  # mutate(nonwfpm = (ml_pm25-wf_pm25_imp)) %>%
  group_by(geoid, year) %>% 
  summarize(ann_wfpm_avg = mean(wf_pm25_imp)) %>% #,
            # ann_nonwfpm = mean(nonwfpm)) %>% 
  ungroup() %>%
  select(geoid,year,ann_wfpm_avg)

wf <- merge(wf, popn.dta, by=c("geoid", "year"), all.y = TRUE)

rr_calc      <- wf %>% dplyr::select(ann_wfpm_avg, year, hisp, nhw, nhb, nha, nhaian, nh2ormore, other_race, total_pop)
rr_calc_long <- rr_calc %>% 
  dplyr::select(ann_wfpm_avg, year, hisp, nhw, nhb, nha, nhaian, nh2ormore, total_pop, other_race) %>%
  pivot_longer(-c(year, ann_wfpm_avg, total_pop), names_to = "race", values_to="count")

main.res     <- as.data.frame(matrix(NA, ncol=(length(all.races)+1), nrow=length(all.years)))
main.res[,1] <- all.years

for (r in 1:length(all.races)){

main.res[,(r+1)] <- rr_calc_long %>% filter(race==all.races[r]) %>% 
  group_by(year) %>%
  summarize(rr = (sum(count*ann_wfpm_avg)/sum(count))/
              (sum(total_pop*ann_wfpm_avg)/sum(total_pop))) %>%
  dplyr::select(rr)
}

names(main.res) <- c("year", paste0("rr_", all.races))
write_csv(main.res, "main_results.csv")

## NOW TIME TO BOOT! 

wfpm <- wfpm %>% dplyr::select(geoid, date, wf_pm25_imp, year)

n.boot <- 100

allBoot.res <- array(NA, dim=c(n.boot, length(all.years), length(all.races)))

set.seed(212)

for (b in 1:n.boot){ ## b = 1
print(b)
  
sample.b <- sample(row.names(wfpm), replace=TRUE)
wfpm.boot <- wfpm[sample.b, ]  
  
wf.boot <- wfpm.boot %>% 
  group_by(geoid, year) %>% 
  summarize(ann_wfpm_avg = mean(wf_pm25_imp)) %>% 
  ungroup() %>%
  select(geoid,year,ann_wfpm_avg)

wf.boot <- merge(wf.boot, popn.dta, by=c("geoid", "year"), all.y = TRUE)

rr_calc.b      <- wf.boot %>% 
  dplyr::select(ann_wfpm_avg, year, hisp, nhw, nhb, nha, nhaian, nh2ormore, other_race, total_pop)
rr_calc_long.b <- rr_calc.b %>% 
  dplyr::select(ann_wfpm_avg, year, hisp, nhw, nhb, nha, nhaian, nh2ormore, total_pop, other_race) %>%
  pivot_longer(-c(year, ann_wfpm_avg, total_pop), names_to = "race", values_to="count")

for (r in 1:length(all.races)){
  
  rr <- rr_calc_long.b %>% filter(race==all.races[r]) %>% 
    group_by(year) %>%
    summarize(rr = (sum(count*ann_wfpm_avg)/sum(count))/
                (sum(total_pop*ann_wfpm_avg)/sum(total_pop))) %>%
    dplyr::select(rr)
  
  allBoot.res[b,,r] <- rr$rr
  rm(rr)
  }

rm(sample.b, wf.boot, wfpm.boot, rr_calc_long.b, rr_calc.b)

}

bootRes.mean <- apply(allBoot.res, c(2,3), mean, na.rm=TRUE)
bootRes.med  <- apply(allBoot.res, c(2,3), median, na.rm=TRUE)
bootRes.lci  <- apply(allBoot.res, c(2,3), quantile, 0.025, na.rm=TRUE)
bootRes.uci  <- apply(allBoot.res, c(2,3), quantile, 0.975, na.rm=TRUE)

colnames(bootRes.mean) <- paste0("mean_", all.races)
colnames(bootRes.med)  <- paste0("median_", all.races)
colnames(bootRes.lci)  <- paste0("lci_", all.races)
colnames(bootRes.uci)  <- paste0("uci_", all.races)

bootRes      <- data.frame(bootRes.mean, bootRes.med, bootRes.lci, bootRes.uci)
bootRes$year <- all.years

write_csv(bootRes, "allBootResults.csv")

save(allBoot.res, file="FullBootRes.RData")
