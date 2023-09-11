#Code to re-run regression for long-term wildfire paper
#Regression models
#Updated 09 Sept 2023

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

#Data
wf_reg <- readRDS(file = "regression_data_v2.rds")

#Paragraph with 90th%tile
wf_reg %>% group_by(high_ces) %>% filter(year==2020) %>%
  summarise(peak = quantile(peak_pm, c(0.25, 0.5, 0.90)), q = c(0.25, 0.5, 0.90)) 
wf_reg %>% group_by(high_ces) %>%filter(year==2020) %>%
  summarise(fiveweek = quantile(pm_freq, c(0.25, 0.5, 0.90)), q = c(0.25, 0.5, 0.90)) 


#Correlation plot
wf_corr <- wf_reg %>% dplyr::select(pm_freq,non_zero_days,peak_pm,  smoke_waves,   ann_wfpm_avg, ann_nonwfpm,
                                    PM2_5, CIscoreP, high_ces, hisp_p, nhw_p, nhb_p, nha_p, nhaian_p, nh2ormore,popden,year)

ggcorr(wf_corr, method=c("pairwise", "spearman"),
       nbreaks = 10,digits = 2, label=TRUE, hjust = 0.8, vjust =.5, size = 3, 
       label_size = 3,
       color = "grey50")

#Get lat/long centroid of census tracts
wf_reg <- st_as_sf(wf_reg)
class(wf_reg$high_ces)
wf_reg$high_ces <- factor(wf_reg$high_ces)
wf_reg$year <- factor(wf_reg$year)
centroids <- st_centroid(wf_reg)
coordinates <- st_coordinates(centroids)
wf_reg$latit <- coordinates[1:118774,1]
wf_reg$lon <- coordinates[1:118774,2]

#Run models for annual average first
model1 <- gam(pm_freq ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2 <- gam(non_zero_days ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3 <- gam(peak_pm ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4 <- gam(smoke_waves ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5 <- gam(ann_wfpm_avg ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)
model1_marginal<-avg_comparisons(model1, variables = c("high_ces"))
model2_marginal<-avg_comparisons(model2, variables = c("high_ces"))
model3_marginal<-avg_comparisons(model3, variables = c("high_ces"))
model4_marginal<-avg_comparisons(model4, variables = c("high_ces"))
model5_marginal<-avg_comparisons(model5, variables = c("high_ces"))

#Save estimates and add year variable
model1_marginal<-model1_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal<-model1_marginal %>% mutate(year="2006-2020")
model1_marginal$year <- factor(model1_marginal$year)
model2_marginal<-model2_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal<-model2_marginal %>% mutate(year="2006-2020")
model2_marginal$year <- factor(model2_marginal$year)
model3_marginal<-model3_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal<-model3_marginal %>% mutate(year="2006-2020")
model3_marginal$year <- factor(model3_marginal$year)
model4_marginal<-model4_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal<-model4_marginal %>% mutate(year="2006-2020")
model4_marginal$year <- factor(model4_marginal$year)
model5_marginal<-model5_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal<-model5_marginal %>% mutate(year="2006-2020")
model5_marginal$year <- factor(model5_marginal$year)

model1_int<- gam(pm_freq ~ high_ces*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_int<- gam(non_zero_days ~ high_ces*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_int <- gam(peak_pm ~ high_ces*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_int<- gam(smoke_waves ~ high_ces*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_int <- gam(ann_wfpm_avg ~ high_ces*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = wf_reg)


#Check residuals for annual average wildfire pm
resid_mod5 <- resid(model5_int)

# Also adding deciles of residuals for plotting
wf_reg$resid <- resid_mod5

#Mapping continuous residuals for annual mean
residual_ct_map <- wf_reg %>%
  ggplot() +
  geom_sf(data = wf_reg,
          color="white", aes(geometry = geometry, fill = resid
          ), lwd = NA, alpha = 0.8) +
  scale_fill_viridis("Residuals") +
  theme_void(base_size=14) +
  facet_wrap(~factor(year)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect=element_blank())
residual_ct_map #Supplemental Figure 23

ggsave("supp_fig23_spatial_residuals_annual.png", dpi=300, height=6, width=6, units="in" )

#Marginal means wow this is amazing
model1_int_marginal <- avg_comparisons(model1_int,  variables = "high_ces", by = "year")
model2_int_marginal <- avg_comparisons(model2_int,  variables = "high_ces", by = "year")
model3_int_marginal <- avg_comparisons(model3_int,  variables = "high_ces", by = "year")
model4_int_marginal <- avg_comparisons(model4_int,  variables = "high_ces", by = "year")
model5_int_marginal <- avg_comparisons(model5_int, variables = "high_ces", by = "year")


model1_int_marginal <- model1_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal <- model2_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal <- model3_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal <- model4_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal <- model5_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
#Add 2006-2020 estimates
model1_int_marginal$year <- factor(model1_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot <- rbind(model1_int_marginal,model1_marginal)
model1_int_marginal_plot <- model1_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model2_int_marginal$year <- factor(model2_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot <- rbind(model2_int_marginal,model2_marginal)
model2_int_marginal_plot <- model2_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model3_int_marginal$year <- factor(model3_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot <- rbind(model3_int_marginal,model3_marginal)
model3_int_marginal_plot <- model3_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model4_int_marginal$year <- factor(model4_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot <- rbind(model4_int_marginal,model4_marginal)
model4_int_marginal_plot <- model4_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model5_int_marginal$year <- factor(model5_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot <- rbind(model5_int_marginal,model5_marginal)
model5_int_marginal_plot <- model5_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

colnames(model5_int_marginal)
pd <- position_dodge(0.1)


#Weeks over 5 difference
reg1_plot <- ggplot(model1_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod)))+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model1_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous( "Mean difference (weeks)"  ) +
  scale_color_manual("", values=c("black","#D25625")) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3)) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1))

#Non-zero days difference
reg2_plot <-ggplot(model2_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model2_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (days)",
                     breaks=c(-3,0,3,6)) +
  ggtitle(bquote(Wildfire~PM[2.5]~">0")) + 
  scale_color_manual("", values=c("black","#D25625")) + 
  
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1))


#Peak wildfire difference
reg3_plot <-ggplot(model3_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model3_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) + scale_color_manual("", values=c("black","#D25625")) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5]))

#Smoke waves
reg4_plot <-ggplot(model4_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model4_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (smoke waves)", breaks=c(-0.5, -0.25, 0, 0.25)) + 
  scale_color_manual("", values=c("black","#D25625")) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle("Smoke waves")

#5 Annual average wildfire PM
reg5_plot <-ggplot(model5_int_marginal_plot, aes(x=year, y=estimate,color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model5_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_color_manual("", values=c("black","#D25625")) + 
  scale_y_continuous(bquote(Mean~difference~(µg/m^3)), breaks=c(-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5)) +
  scale_x_discrete("") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5]))

reg1_plot + reg2_plot + reg3_plot + reg4_plot + reg5_plot +  plot_annotation(tag_levels = 'A') +  plot_layout(ncol = 2)
ggsave("longterm-pm/analysis/regression/figure4_ces_difference_regression_v2.png", dpi=300, height=11, width=8, units="in" )

#Save regression results
model5_int_marginal_plot <- model5_int_marginal_plot %>% mutate(exposure=5)
model4_int_marginal_plot <- model4_int_marginal_plot %>% mutate(exposure=4)
model3_int_marginal_plot <- model3_int_marginal_plot %>% mutate(exposure=3)
model2_int_marginal_plot <- model2_int_marginal_plot %>% mutate(exposure=2)
model1_int_marginal_plot <- model1_int_marginal_plot %>% mutate(exposure=1)
model_int_marginal_plot_data <- rbind(model1_int_marginal_plot, model2_int_marginal_plot, model3_int_marginal_plot, model4_int_marginal_plot, model5_int_marginal_plot)

combined <- with(model_int_marginal_plot_data, sprintf('%.2f (%.2f, %.2f)', estimate, conf.low, conf.high))

model_int_marginal_plot_data_2 <- cbind(model_int_marginal_plot_data,combined) 
##Supplementary Table 3
write_csv(model_int_marginal_plot_data_2, "longterm-pm/analysis/regression/supp_table3_ces_regression_resultsv2.csv")

####RACE ETHNICITY REGRESSIONS###
#Run models for annual average first
#I think I will do this for a SD increase in each
summary(wf_reg$hisp_p)
summary(wf_reg$nhaian_p)
summary(wf_reg$nha_p)
summary(wf_reg$nhb_p)
sd(wf_reg$hisp_p)
sd(wf_reg$nhaian_p)
sd(wf_reg$nha_p)
sd(wf_reg$nhb_p)
sd(wf_reg$nh2ormore)
sd(wf_reg$nhw_p)

#Overall Hispanic models
model1_hisp <- gam(pm_freq ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_hisp <- gam(non_zero_days ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_hisp <- gam(peak_pm ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_hisp <- gam(smoke_waves ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_hisp <- gam(ann_wfpm_avg ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
summary(model1_hisp)
summary(model2_hisp)
summary(model3_hisp)
summary(model4_hisp)
summary(model5_hisp)
model1_marginal_hisp<-avg_comparisons(model1_hisp, variables = list(hisp_p="sd"))
model2_marginal_hisp<-avg_comparisons(model2_hisp, variables = list(hisp_p="sd"))
model3_marginal_hisp<-avg_comparisons(model3_hisp, variables = list(hisp_p="sd"))
model4_marginal_hisp<-avg_comparisons(model4_hisp, variables = list(hisp_p="sd"))
model5_marginal_hisp<-avg_comparisons(model5_hisp, variables = list(hisp_p="sd"))

#Overall AIAN models
model1_aian <- gam(pm_freq ~ nhaian_p +   s(popden, fx=TRUE, k=9)  + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_aian <- gam(non_zero_days ~ nhaian_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_aian <- gam(peak_pm ~ nhaian_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = subset(wf_reg))
model4_aian <- gam(smoke_waves ~ nhaian_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_aian <- gam(ann_wfpm_avg ~ nhaian_p + s(popden) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_aian<-avg_comparisons(model1_aian, variables = list(nhaian_p="sd"))
model2_marginal_aian<-avg_comparisons(model2_aian, variables = list(nhaian_p="sd"))
model3_marginal_aian<-avg_comparisons(model3_aian, variables = list(nhaian_p="sd"))
model4_marginal_aian<-avg_comparisons(model4_aian, variables = list(nhaian_p="sd"))
model5_marginal_aian<-avg_comparisons(model5_aian, variables = list(nhaian_p="sd"))

#Overall Asian models
model1_asian <- gam(pm_freq ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_asian <- gam(non_zero_days ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_asian <- gam(peak_pm ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_asian <- gam(smoke_waves ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_asian <- gam(ann_wfpm_avg ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_asian<-avg_comparisons(model1_asian, variables = list(nha_p="sd"))
model2_marginal_asian<-avg_comparisons(model2_asian, variables = list(nha_p="sd"))
model3_marginal_asian<-avg_comparisons(model3_asian, variables = list(nha_p="sd"))
model4_marginal_asian<-avg_comparisons(model4_asian, variables = list(nha_p="sd"))
model5_marginal_asian<-avg_comparisons(model5_asian, variables = list(nha_p="sd"))

#Overall Black models
model1_black <- gam(pm_freq ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_black <- gam(non_zero_days ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_black <- gam(peak_pm ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_black <- gam(smoke_waves ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_black <- gam(ann_wfpm_avg ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_black<-avg_comparisons(model1_black, variables = list(nhb_p="sd"))
model2_marginal_black<-avg_comparisons(model2_black, variables = list(nhb_p="sd"))
model3_marginal_black<-avg_comparisons(model3_black, variables = list(nhb_p="sd"))
model4_marginal_black<-avg_comparisons(model4_black, variables = list(nhb_p="sd"))
model5_marginal_black<-avg_comparisons(model5_black, variables = list(nhb_p="sd"))

#Overall nh2_
model1_2more <- gam(pm_freq ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_2more <- gam(non_zero_days ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_2more <- gam(peak_pm ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_2more <- gam(smoke_waves ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_2more <- gam(ann_wfpm_avg ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_2more<-avg_comparisons(model1_2more, variables = list(nh2ormore="sd"))
model2_marginal_2more<-avg_comparisons(model2_2more, variables = list(nh2ormore="sd"))
model3_marginal_2more<-avg_comparisons(model3_2more, variables = list(nh2ormore="sd"))
model4_marginal_2more<-avg_comparisons(model4_2more, variables = list(nh2ormore="sd"))
model5_marginal_2more<-avg_comparisons(model5_2more, variables = list(nh2ormore="sd"))

#Overall white
model1_white <- gam(pm_freq ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_white <- gam(non_zero_days ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_white <- gam(peak_pm ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_white <- gam(smoke_waves ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_white <- gam(ann_wfpm_avg ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_white<-avg_comparisons(model1_white, variables = list(nhw_p="sd"))
model2_marginal_white<-avg_comparisons(model2_white, variables = list(nhw_p="sd"))
model3_marginal_white<-avg_comparisons(model3_white, variables = list(nhw_p="sd"))
model4_marginal_white<-avg_comparisons(model4_white, variables = list(nhw_p="sd"))
model5_marginal_white<-avg_comparisons(model5_white, variables = list(nhw_p="sd"))


#Make 5 plots with labels from figure 1 (text/code) for each exposure and use colors from figure 3 for race/ethnicity designations
race <- c("hisp_p", "nha_p", "nhaian_p", "nhb_p", "nh2ormore", "nhw_p")
race <- factor(race,  levels=c("hisp_p", "nha_p", "nhaian_p", "nhb_p", "nh2ormore", "nhw_p"),
               labels=c("Hispanic", "NH Asian", "NH American Indian", "NH Black", "NH 2+", "NH white"))

model1_marginal_plot_race <- rbind(model1_marginal_hisp, model1_marginal_aian, model1_marginal_asian, model1_marginal_black,
                                   model1_marginal_2more, model1_marginal_white)
model1_marginal_plot_race$race <- race
model2_marginal_plot_race <- rbind(model2_marginal_hisp, model2_marginal_aian, model2_marginal_asian, model2_marginal_black,
                                   model2_marginal_2more, model2_marginal_white)
model2_marginal_plot_race$race <- race

model3_marginal_plot_race <- rbind(model3_marginal_hisp, model3_marginal_aian, model3_marginal_asian, model3_marginal_black,
                                   model3_marginal_2more, model3_marginal_white)
model3_marginal_plot_race$race <- race

model4_marginal_plot_race <- rbind(model4_marginal_hisp, model4_marginal_aian, model4_marginal_asian, model4_marginal_black,
                                   model4_marginal_2more, model4_marginal_white)
model4_marginal_plot_race$race <- race

model5_marginal_plot_race <- rbind(model5_marginal_hisp, model5_marginal_aian, model5_marginal_asian, model5_marginal_black,
                                   model5_marginal_2more, model5_marginal_white)
model5_marginal_plot_race$race <- race

#Plotting across the study period the expected change in the exposures for a 1-SD increase in race/ethnicity percent
exp_1_race_plot <- ggplot(model1_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean differences (weeks)", breaks=c(-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5)) +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3))


exp_2_race_plot <- ggplot(model2_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (days)", breaks=c(-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5)) +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">0"))

exp_3_race_plot <-ggplot(model3_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5]))

exp_4_race_plot <- ggplot(model4_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (smoke waves)") +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle("Smoke waves")

exp_5_race_plot <-  ggplot(model5_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5]))

exp_1_race_plot + exp_2_race_plot + exp_3_race_plot + exp_4_race_plot + exp_5_race_plot + plot_annotation(tag_levels = 'A') +  plot_layout(nrow = 1)

#####YEAR SPECIFIC BY RACE#####
#Overall Hispanic models
model1_hisp_int <- gam(pm_freq ~ hisp_p*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_hisp_int <- gam(non_zero_days ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_hisp_int <- gam(peak_pm ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_hisp_int <- gam(smoke_waves ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_hisp_int <- gam(ann_wfpm_avg ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)

model1_int_marginal_hisp<-avg_comparisons(model1_hisp_int, variables = list(hisp_p="sd"), by="year")
model2_int_marginal_hisp<-avg_comparisons(model2_hisp_int, variables = list(hisp_p="sd"), by="year")
model3_int_marginal_hisp<-avg_comparisons(model3_hisp_int, variables = list(hisp_p="sd"), by="year")
model4_int_marginal_hisp<-avg_comparisons(model4_hisp_int, variables = list(hisp_p="sd"), by="year")
model5_int_marginal_hisp<-avg_comparisons(model5_hisp_int, variables = list(hisp_p="sd"), by="year")

#Overall AIAN models
model1_aian_int <- gam(pm_freq ~ nhaian_p*year  + s(popden, fx=TRUE, k=4) + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_aian_int <- gam(non_zero_days ~ nhaian_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_aian_int <- gam(peak_pm ~ nhaian_p*year  + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_aian_int <- gam(smoke_waves ~ nhaian_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_aian_int <- gam(ann_wfpm_avg ~ nhaian_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_aian_int<-avg_comparisons(model1_aian_int, variables = list(nhaian_p="sd"), by="year")
model2_marginal_aian_int<-avg_comparisons(model2_aian_int, variables = list(nhaian_p="sd"), by="year")
model3_marginal_aian_int<-avg_comparisons(model3_aian_int, variables = list(nhaian_p="sd"), by="year")
model4_marginal_aian_int<-avg_comparisons(model4_aian_int, variables = list(nhaian_p="sd"), by="year")
model5_marginal_aian_int<-avg_comparisons(model5_aian_int, variables = list(nhaian_p="sd"), by="year")

#Overall Asian models
model1_asian_int <- gam(pm_freq ~ nha_p*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_asian_int <- gam(non_zero_days ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_asian_int <- gam(peak_pm ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_asian_int <- gam(smoke_waves ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_asian_int <- gam(ann_wfpm_avg ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_asian_int<-avg_comparisons(model1_asian_int, variables = list(nha_p="sd"), by="year")
model2_marginal_asian_int<-avg_comparisons(model2_asian_int, variables = list(nha_p="sd"), by="year")
model3_marginal_asian_int<-avg_comparisons(model3_asian_int, variables = list(nha_p="sd"), by="year")
model4_marginal_asian_int<-avg_comparisons(model4_asian_int, variables = list(nha_p="sd"), by="year")
model5_marginal_asian_int<-avg_comparisons(model5_asian_int, variables = list(nha_p="sd"), by="year")

#Overall Black models
model1_black_int <- gam(pm_freq ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_black_int <- gam(non_zero_days ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_black_int <- gam(peak_pm ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_black_int <- gam(smoke_waves ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_black_int <- gam(ann_wfpm_avg ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_black_int<-avg_comparisons(model1_black_int, variables = list(nhb_p="sd"), by="year")
model2_marginal_black_int<-avg_comparisons(model2_black_int, variables = list(nhb_p="sd"), by="year")
model3_marginal_black_int<-avg_comparisons(model3_black_int, variables = list(nhb_p="sd"), by="year")
model4_marginal_black_int<-avg_comparisons(model4_black_int, variables = list(nhb_p="sd"), by="year")
model5_marginal_black_int<-avg_comparisons(model5_black_int, variables = list(nhb_p="sd"), by="year")

#Overall nh2_
model1_2more_int <- gam(pm_freq ~ nh2ormore*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_2more_int <- gam(non_zero_days ~ nh2ormore*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_2more_int <- gam(peak_pm ~ nh2ormore*year  + s(popden, fx=TRUE, k=9) +  s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_2more_int <- gam(smoke_waves ~ nh2ormore*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_2more_int <- gam(ann_wfpm_avg ~ nh2ormore*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_2more_int<-avg_comparisons(model1_2more_int, variables = list(nh2ormore="sd"), by="year")
model2_marginal_2more_int<-avg_comparisons(model2_2more_int, variables = list(nh2ormore="sd"), by="year")
model3_marginal_2more_int<-avg_comparisons(model3_2more_int, variables = list(nh2ormore="sd"), by="year")
model4_marginal_2more_int<-avg_comparisons(model4_2more_int, variables = list(nh2ormore="sd"), by="year")
model5_marginal_2more_int<-avg_comparisons(model5_2more_int, variables = list(nh2ormore="sd"), by="year")

#Overall white
model1_white_int <- gam(pm_freq ~ nhw_p*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model2_white_int <- gam(non_zero_days ~ nhw_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg, family = nb())
model3_white_int <- gam(peak_pm ~ nhw_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model4_white_int <- gam(smoke_waves ~ nhw_p*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model5_white_int <- gam(ann_wfpm_avg ~ nhw_p*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = wf_reg)
model1_marginal_white_int<-avg_comparisons(model1_white_int, variables = list(nhw_p="sd"), by="year")
model2_marginal_white_int<-avg_comparisons(model2_white_int, variables = list(nhw_p="sd"), by="year")
model3_marginal_white_int<-avg_comparisons(model3_white_int, variables = list(nhw_p="sd"), by="year")
model4_marginal_white_int<-avg_comparisons(model4_white_int, variables = list(nhw_p="sd"), by="year")
model5_marginal_white_int<-avg_comparisons(model5_white_int, variables = list(nhw_p="sd"), by="year")

#Model 1
model1_marginal_hisp<-model1_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_hisp<-model1_marginal_hisp %>% mutate(year="2006-2020")
model1_marginal_asian<-model1_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_asian<-model1_marginal_asian %>% mutate(year="2006-2020")
model1_marginal_aian<-model1_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_aian<-model1_marginal_aian %>% mutate(year="2006-2020")
model1_marginal_black<-model1_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_black<-model1_marginal_black %>% mutate(year="2006-2020")
model1_marginal_2more<-model1_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_2more<-model1_marginal_2more %>% mutate(year="2006-2020")
model1_marginal_white<-model1_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_white<-model1_marginal_white %>% mutate(year="2006-2020")

model1_int_marginal_hisp <- model1_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_asian <- model1_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_aian <- model1_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_black <- model1_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_2more <- model1_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_white <- model1_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model1_int_marginal_hisp$year <- factor(model1_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_hisp <- rbind(model1_int_marginal_hisp,model1_marginal_hisp)
model1_int_marginal_plot_hisp <- model1_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model1_int_marginal_plot_hisp <- model1_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_asian$year <- factor(model1_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_asian <- rbind(model1_int_marginal_asian,model1_marginal_asian)
model1_int_marginal_plot_asian <- model1_int_marginal_plot_asian %>% mutate(race="NH Asian")
model1_int_marginal_plot_asian <- model1_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_aian$year <- factor(model1_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_aian <- rbind(model1_int_marginal_aian,model1_marginal_aian)
model1_int_marginal_plot_aian <- model1_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model1_int_marginal_plot_aian <- model1_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_black$year <- factor(model1_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_black <- rbind(model1_int_marginal_black,model1_marginal_black)
model1_int_marginal_plot_black <- model1_int_marginal_plot_black %>% mutate(race="NH Black")
model1_int_marginal_plot_black <- model1_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_2more$year <- factor(model1_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_2more <- rbind(model1_int_marginal_2more,model1_marginal_2more)
model1_int_marginal_plot_2more <- model1_int_marginal_plot_2more %>% mutate(race="NH 2+")
model1_int_marginal_plot_2more <- model1_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_white$year <- factor(model1_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_white <- rbind(model1_int_marginal_white,model1_marginal_white)
model1_int_marginal_plot_white <- model1_int_marginal_plot_white %>% mutate(race="NH white")
model1_int_marginal_plot_white <- model1_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model1_int_marginal_plot_all <- rbind(model1_int_marginal_plot_hisp,model1_int_marginal_plot_asian,
                                      model1_int_marginal_plot_aian,model1_int_marginal_plot_black,
                                      model1_int_marginal_plot_2more,model1_int_marginal_plot_white)

model1_int_marginal_plot_all_df <- data.frame(model1_int_marginal_plot_all)
model1_int_marginal_plot_all_df <- model1_int_marginal_plot_all_df %>% mutate(exposure=1)
#Race colors
race_values <- c("#341C5D", "#754E71" , "#DA8940",  "#E19A8F", "#9C372B", "#8D9FD7" )

exp_1_race_plot_all <-  ggplot(model1_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous("Mean difference (weeks)") +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3)) +theme(axis.title.x=element_blank(),
                                                         axis.text.x=element_blank(),
                                                         axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_1_race_plot_all
ggsave("longterm-pm/analysis/regression/supp_fig14_race_exp1_meandiff_v2.png", dpi=300, height=11, width=8, units="in" )


#Model 2
model2_marginal_hisp<-model2_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_hisp<-model2_marginal_hisp %>% mutate(year="2006-2020")
model2_marginal_asian<-model2_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_asian<-model2_marginal_asian %>% mutate(year="2006-2020")
model2_marginal_aian<-model2_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_aian<-model2_marginal_aian %>% mutate(year="2006-2020")
model2_marginal_black<-model2_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_black<-model2_marginal_black %>% mutate(year="2006-2020")
model2_marginal_2more<-model2_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_2more<-model2_marginal_2more %>% mutate(year="2006-2020")
model2_marginal_white<-model2_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_white<-model2_marginal_white %>% mutate(year="2006-2020")

model2_int_marginal_hisp <- model2_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_asian <- model2_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_aian <- model2_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_black <- model2_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_2more <- model2_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_white <- model2_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model2_int_marginal_hisp$year <- factor(model2_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_hisp <- rbind(model2_int_marginal_hisp,model2_marginal_hisp)
model2_int_marginal_plot_hisp <- model2_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model2_int_marginal_plot_hisp <- model2_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_asian$year <- factor(model2_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_asian <- rbind(model2_int_marginal_asian,model2_marginal_asian)
model2_int_marginal_plot_asian <- model2_int_marginal_plot_asian %>% mutate(race="NH Asian")
model2_int_marginal_plot_asian <- model2_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_aian$year <- factor(model2_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_aian <- rbind(model2_int_marginal_aian,model2_marginal_aian)
model2_int_marginal_plot_aian <- model2_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model2_int_marginal_plot_aian <- model2_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_black$year <- factor(model2_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_black <- rbind(model2_int_marginal_black,model2_marginal_black)
model2_int_marginal_plot_black <- model2_int_marginal_plot_black %>% mutate(race="NH Black")
model2_int_marginal_plot_black <- model2_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_2more$year <- factor(model2_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_2more <- rbind(model2_int_marginal_2more,model2_marginal_2more)
model2_int_marginal_plot_2more <- model2_int_marginal_plot_2more %>% mutate(race="NH 2+")
model2_int_marginal_plot_2more <- model2_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_white$year <- factor(model2_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_white <- rbind(model2_int_marginal_white,model2_marginal_white)
model2_int_marginal_plot_white <- model2_int_marginal_plot_white %>% mutate(race="NH white")
model2_int_marginal_plot_white <- model2_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model2_int_marginal_plot_all <- rbind(model2_int_marginal_plot_hisp,model2_int_marginal_plot_asian,
                                      model2_int_marginal_plot_aian,model2_int_marginal_plot_black,
                                      model2_int_marginal_plot_2more,model2_int_marginal_plot_white)

model2_int_marginal_plot_all_df <- data.frame(model2_int_marginal_plot_all)
model2_int_marginal_plot_all_df <- model2_int_marginal_plot_all_df %>% mutate(exposure=2)

exp_2_race_plot_all <-  ggplot(model2_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous("Mean difference (days)") +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">0"~"µg/m"^3)) +theme(axis.title.x=element_blank(),
                                                         axis.text.x=element_blank(),
                                                         axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_2_race_plot_all
ggsave("longterm-pm/analysis/regression/supp_fig14_race_exp2_meandiff.png", dpi=300, height=11, width=8, units="in" )

#Model 3 -- peak PM 
model3_marginal_hisp<-model3_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_hisp<-model3_marginal_hisp %>% mutate(year="2006-2020")
model3_marginal_asian<-model3_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_asian<-model3_marginal_asian %>% mutate(year="2006-2020")
model3_marginal_aian<-model3_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_aian<-model3_marginal_aian %>% mutate(year="2006-2020")
model3_marginal_black<-model3_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_black<-model3_marginal_black %>% mutate(year="2006-2020")
model3_marginal_2more<-model3_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_2more<-model3_marginal_2more %>% mutate(year="2006-2020")
model3_marginal_white<-model3_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_white<-model3_marginal_white %>% mutate(year="2006-2020")

model3_int_marginal_hisp <- model3_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_asian <- model3_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_aian <- model3_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_black <- model3_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_2more <- model3_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_white <- model3_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model3_int_marginal_hisp$year <- factor(model3_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_hisp <- rbind(model3_int_marginal_hisp,model3_marginal_hisp)
model3_int_marginal_plot_hisp <- model3_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model3_int_marginal_plot_hisp <- model3_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_asian$year <- factor(model3_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_asian <- rbind(model3_int_marginal_asian,model3_marginal_asian)
model3_int_marginal_plot_asian <- model3_int_marginal_plot_asian %>% mutate(race="NH Asian")
model3_int_marginal_plot_asian <- model3_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_aian$year <- factor(model3_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_aian <- rbind(model3_int_marginal_aian,model3_marginal_aian)
model3_int_marginal_plot_aian <- model3_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model3_int_marginal_plot_aian <- model3_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_black$year <- factor(model3_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_black <- rbind(model3_int_marginal_black,model3_marginal_black)
model3_int_marginal_plot_black <- model3_int_marginal_plot_black %>% mutate(race="NH Black")
model3_int_marginal_plot_black <- model3_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_2more$year <- factor(model3_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_2more <- rbind(model3_int_marginal_2more,model3_marginal_2more)
model3_int_marginal_plot_2more <- model3_int_marginal_plot_2more %>% mutate(race="NH 2+")
model3_int_marginal_plot_2more <- model3_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_white$year <- factor(model3_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_white <- rbind(model3_int_marginal_white,model3_marginal_white)
model3_int_marginal_plot_white <- model3_int_marginal_plot_white %>% mutate(race="NH white")
model3_int_marginal_plot_white <- model3_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model3_int_marginal_plot_all <- rbind(model3_int_marginal_plot_hisp,model3_int_marginal_plot_asian,
                                      model3_int_marginal_plot_aian,model3_int_marginal_plot_black,
                                      model3_int_marginal_plot_2more,model3_int_marginal_plot_white)

model3_int_marginal_plot_all_df <- data.frame(model3_int_marginal_plot_all)
model3_int_marginal_plot_all_df <- model3_int_marginal_plot_all_df %>% mutate(exposure=3)

exp_3_race_plot_all <-  ggplot(model3_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5])) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_3_race_plot_all
ggsave("longterm-pm/analysis/regression/supp_fig14_race_exp3_meandiff.png", dpi=300, height=11, width=8, units="in" )

#Model 4 -- smoke waves
model4_marginal_hisp<-model4_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_hisp<-model4_marginal_hisp %>% mutate(year="2006-2020")
model4_marginal_asian<-model4_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_asian<-model4_marginal_asian %>% mutate(year="2006-2020")
model4_marginal_aian<-model4_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_aian<-model4_marginal_aian %>% mutate(year="2006-2020")
model4_marginal_black<-model4_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_black<-model4_marginal_black %>% mutate(year="2006-2020")
model4_marginal_2more<-model4_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_2more<-model4_marginal_2more %>% mutate(year="2006-2020")
model4_marginal_white<-model4_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_white<-model4_marginal_white %>% mutate(year="2006-2020")

model4_int_marginal_hisp <- model4_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_asian <- model4_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_aian <- model4_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_black <- model4_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_2more <- model4_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_white <- model4_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model4_int_marginal_hisp$year <- factor(model4_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_hisp <- rbind(model4_int_marginal_hisp,model4_marginal_hisp)
model4_int_marginal_plot_hisp <- model4_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model4_int_marginal_plot_hisp <- model4_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_asian$year <- factor(model4_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_asian <- rbind(model4_int_marginal_asian,model4_marginal_asian)
model4_int_marginal_plot_asian <- model4_int_marginal_plot_asian %>% mutate(race="NH Asian")
model4_int_marginal_plot_asian <- model4_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_aian$year <- factor(model4_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_aian <- rbind(model4_int_marginal_aian,model4_marginal_aian)
model4_int_marginal_plot_aian <- model4_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model4_int_marginal_plot_aian <- model4_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_black$year <- factor(model4_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_black <- rbind(model4_int_marginal_black,model4_marginal_black)
model4_int_marginal_plot_black <- model4_int_marginal_plot_black %>% mutate(race="NH Black")
model4_int_marginal_plot_black <- model4_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_2more$year <- factor(model4_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_2more <- rbind(model4_int_marginal_2more,model4_marginal_2more)
model4_int_marginal_plot_2more <- model4_int_marginal_plot_2more %>% mutate(race="NH 2+")
model4_int_marginal_plot_2more <- model4_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_white$year <- factor(model4_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_white <- rbind(model4_int_marginal_white,model4_marginal_white)
model4_int_marginal_plot_white <- model4_int_marginal_plot_white %>% mutate(race="NH white")
model4_int_marginal_plot_white <- model4_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model4_int_marginal_plot_all <- rbind(model4_int_marginal_plot_hisp,model4_int_marginal_plot_asian,
                                      model4_int_marginal_plot_aian,model4_int_marginal_plot_black,
                                      model4_int_marginal_plot_2more,model4_int_marginal_plot_white)

model4_int_marginal_plot_all_df <- data.frame(model4_int_marginal_plot_all)
model4_int_marginal_plot_all_df <- model4_int_marginal_plot_all_df %>% mutate(exposure=4)

exp_4_race_plot_all <-  ggplot(model4_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous("Mean difference (smoke waves)") +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle("Smoke waves") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_4_race_plot_all
ggsave("longterm-pm/analysis/regression/supp_fig14_race_exp4_meandiff.png", dpi=300, height=11, width=8, units="in" )


#Grab main effect from 2006-2020 for model 5
model5_marginal_hisp<-model5_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_hisp<-model5_marginal_hisp %>% mutate(year="2006-2020")
model5_marginal_asian<-model5_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_asian<-model5_marginal_asian %>% mutate(year="2006-2020")
model5_marginal_aian<-model5_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_aian<-model5_marginal_aian %>% mutate(year="2006-2020")
model5_marginal_black<-model5_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_black<-model5_marginal_black %>% mutate(year="2006-2020")
model5_marginal_2more<-model5_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_2more<-model5_marginal_2more %>% mutate(year="2006-2020")
model5_marginal_white<-model5_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_white<-model5_marginal_white %>% mutate(year="2006-2020")

model5_int_marginal_hisp <- model5_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_asian <- model5_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_aian <- model5_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_black <- model5_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_2more <- model5_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_white <- model5_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model5_int_marginal_hisp$year <- factor(model5_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_hisp <- rbind(model5_int_marginal_hisp,model5_marginal_hisp)
model5_int_marginal_plot_hisp <- model5_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model5_int_marginal_plot_hisp <- model5_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_asian$year <- factor(model5_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_asian <- rbind(model5_int_marginal_asian,model5_marginal_asian)
model5_int_marginal_plot_asian <- model5_int_marginal_plot_asian %>% mutate(race="NH Asian")
model5_int_marginal_plot_asian <- model5_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_aian$year <- factor(model5_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_aian <- rbind(model5_int_marginal_aian,model5_marginal_aian)
model5_int_marginal_plot_aian <- model5_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model5_int_marginal_plot_aian <- model5_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_black$year <- factor(model5_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_black <- rbind(model5_int_marginal_black,model5_marginal_black)
model5_int_marginal_plot_black <- model5_int_marginal_plot_black %>% mutate(race="NH Black")
model5_int_marginal_plot_black <- model5_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_2more$year <- factor(model5_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_2more <- rbind(model5_int_marginal_2more,model5_marginal_2more)
model5_int_marginal_plot_2more <- model5_int_marginal_plot_2more %>% mutate(race="NH 2+")
model5_int_marginal_plot_2more <- model5_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_white$year <- factor(model5_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_white <- rbind(model5_int_marginal_white,model5_marginal_white)
model5_int_marginal_plot_white <- model5_int_marginal_plot_white %>% mutate(race="NH white")
model5_int_marginal_plot_white <- model5_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model5_int_marginal_plot_all <- rbind(model5_int_marginal_plot_hisp,model5_int_marginal_plot_asian,
                                      model5_int_marginal_plot_aian,model5_int_marginal_plot_black,
                                      model5_int_marginal_plot_2more,model5_int_marginal_plot_white)

model5_int_marginal_plot_all_df <- data.frame(model5_int_marginal_plot_all)
model5_int_marginal_plot_all_df <- model5_int_marginal_plot_all_df %>% mutate(exposure=5)

exp_5_race_plot_all <-  ggplot(model5_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5])) +theme(axis.title.x=element_blank(),
                                                       axis.text.x=element_blank(),
                                                       axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_5_race_plot_all
ggsave("longterm-pm/analysis/regression/supp_fig14_race_exp5_meandiff.png", dpi=300, height=11, width=8, units="in" )

#Save estimates from regression for race:
models_int_marginal_plot_all_df <-  rbind(model1_int_marginal_plot_all_df,model2_int_marginal_plot_all_df,model3_int_marginal_plot_all_df,model4_int_marginal_plot_all_df,model5_int_marginal_plot_all_df)
write_csv(models_int_marginal_plot_all_df, "longterm-pm/analysis/regression/race_regression_results.csv")

###URBAN/RURAL REGRESSION
#Create specific cuts for urban/rural
wf_reg <- wf_reg %>% mutate(high_ces=ifelse(((CIscore>=39.5 & ces_year=="ces3" & rural==0) | (CIscore>=27.3 & ces_year=="ces3" & rural==1) |
                                               (CIscore>=40.1 & ces_year=="ces4" & rural==0)|(CIscore>=28.8 & ces_year=="ces4" & rural==1)),1,0))
table(wf$high_ces,wf$rural)

#Get data frames for urban and rural alone
urban <- wf_reg %>% filter(rural==0)
rural <- wf_reg %>% filter(rural==1)

#Run models for annual average first
#RURAL#
model1 <- gam(pm_freq ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data =  rural)
model2 <- gam(non_zero_days ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3 <- gam(peak_pm ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model4 <- gam(smoke_waves ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model5 <- gam(ann_wfpm_avg ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)
model1_marginal<-avg_comparisons(model1, variables = c("high_ces"))
model2_marginal<-avg_comparisons(model2, variables = c("high_ces"))
model3_marginal<-avg_comparisons(model3, variables = c("high_ces"))
model4_marginal<-avg_comparisons(model4, variables = c("high_ces"))
model5_marginal<-avg_comparisons(model5, variables = c("high_ces"))

#Save estimates and add year variable
model1_marginal<-model1_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal<-model1_marginal %>% mutate(year="2006-2020")
model1_marginal$year <- factor(model1_marginal$year)
model2_marginal<-model2_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal<-model2_marginal %>% mutate(year="2006-2020")
model2_marginal$year <- factor(model2_marginal$year)
model3_marginal<-model3_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal<-model3_marginal %>% mutate(year="2006-2020")
model3_marginal$year <- factor(model3_marginal$year)
model4_marginal<-model4_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal<-model4_marginal %>% mutate(year="2006-2020")
model4_marginal$year <- factor(model4_marginal$year)
model5_marginal<-model5_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal<-model5_marginal %>% mutate(year="2006-2020")
model5_marginal$year <- factor(model5_marginal$year)

model1_int<- gam(pm_freq ~ high_ces*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_int<- gam(non_zero_days ~ high_ces*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_int <- gam(peak_pm ~ high_ces*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_int<- gam(smoke_waves ~ high_ces*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_int <- gam(ann_wfpm_avg ~ high_ces*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = rural)

summary(model5_int)

#Marginal means wow this is amazing
model1_int_marginal <- avg_comparisons(model1_int,  variables = "high_ces", by = "year")
model2_int_marginal <- avg_comparisons(model2_int,  variables = "high_ces", by = "year")
model3_int_marginal <- avg_comparisons(model3_int,  variables = "high_ces", by = "year")
model4_int_marginal <- avg_comparisons(model4_int,  variables = "high_ces", by = "year")
model5_int_marginal <- avg_comparisons(model5_int, variables = "high_ces", by = "year")

model1_int_marginal <- model1_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal <- model2_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal <- model3_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal <- model4_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal <- model5_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
#Add 2006-2020 estimates
model1_int_marginal$year <- factor(model1_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot <- rbind(model1_int_marginal,model1_marginal)
model1_int_marginal_plot <- model1_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model2_int_marginal$year <- factor(model2_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot <- rbind(model2_int_marginal,model2_marginal)
model2_int_marginal_plot <- model2_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model3_int_marginal$year <- factor(model3_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot <- rbind(model3_int_marginal,model3_marginal)
model3_int_marginal_plot <- model3_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model4_int_marginal$year <- factor(model4_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot <- rbind(model4_int_marginal,model4_marginal)
model4_int_marginal_plot <- model4_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model5_int_marginal$year <- factor(model5_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot <- rbind(model5_int_marginal,model5_marginal)
model5_int_marginal_plot <- model5_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

colnames(model5_int_marginal)
pd <- position_dodge(0.1)


#Weeks over 5 difference
reg1_plot <- ggplot(model1_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod)))+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model1_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous( "Mean difference (weeks)"  ) +
  scale_color_manual("", values=c("black","#D25625")) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3)) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3)) 


#Non-zero days difference
reg2_plot <-ggplot(model2_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model2_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (days)") +
  ggtitle(bquote(Wildfire~PM[2.5]~">0")) + 
  scale_color_manual("", values=c("black","#D25625")) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle(bquote(Wildfire~PM[2.5]~">0"))


#Peak wildfire difference
reg3_plot <-ggplot(model3_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model3_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) + scale_color_manual("", values=c("black","#D25625")) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5]))

#Smoke waves
reg4_plot <-ggplot(model4_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model4_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (smoke waves)") + 
  scale_color_manual("", values=c("black","#D25625")) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle("Smoke waves")

#5 Annual average wildfire PM
reg5_plot <-ggplot(model5_int_marginal_plot, aes(x=year, y=estimate,color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model5_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_color_manual("", values=c("black","#D25625")) + 
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_x_discrete("") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5]))

reg1_plot + reg2_plot + reg3_plot + reg4_plot + reg5_plot +  plot_annotation(tag_levels = 'A') +  plot_layout(ncol = 2)
ggsave("/Users/joancasey/Documents/Columbia/Grants/Wildfire R01/manuscripts/longterm-pm/analysis/regression/supp_fig13.png", dpi=300, height=11, width=8, units="in" )

#Save regression results
model5_int_marginal_plot <- model5_int_marginal_plot %>% mutate(exposure=5)
model4_int_marginal_plot <- model4_int_marginal_plot %>% mutate(exposure=4)
model3_int_marginal_plot <- model3_int_marginal_plot %>% mutate(exposure=3)
model2_int_marginal_plot <- model2_int_marginal_plot %>% mutate(exposure=2)
model1_int_marginal_plot <- model1_int_marginal_plot %>% mutate(exposure=1)
model_int_marginal_plot_data <- rbind(model1_int_marginal_plot, model2_int_marginal_plot, model3_int_marginal_plot, model4_int_marginal_plot, model5_int_marginal_plot)

combined <- with(model_int_marginal_plot_data, sprintf('%.2f (%.2f, %.2f)', estimate, conf.low, conf.high))

model_int_marginal_plot_data_2 <- cbind(model_int_marginal_plot_data,combined) 
write_csv(model_int_marginal_plot_data_2, "/Users/joancasey/Documents/Columbia/Grants/Wildfire R01/manuscripts/longterm-pm/analysis/regression/ces_regression_results_rural.csv")

##URBAN ONLY##
model1 <- gam(pm_freq ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model2 <- gam(non_zero_days ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3 <- gam(peak_pm ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model4 <- gam(smoke_waves ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model5 <- gam(ann_wfpm_avg ~ high_ces + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)
model1_marginal<-avg_comparisons(model1, variables = c("high_ces"))
model2_marginal<-avg_comparisons(model2, variables = c("high_ces"))
model3_marginal<-avg_comparisons(model3, variables = c("high_ces"))
model4_marginal<-avg_comparisons(model4, variables = c("high_ces"))
model5_marginal<-avg_comparisons(model5, variables = c("high_ces"))

#Save estimates and add year variable
model1_marginal<-model1_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal<-model1_marginal %>% mutate(year="2006-2020")
model1_marginal$year <- factor(model1_marginal$year)
model2_marginal<-model2_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal<-model2_marginal %>% mutate(year="2006-2020")
model2_marginal$year <- factor(model2_marginal$year)
model3_marginal<-model3_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal<-model3_marginal %>% mutate(year="2006-2020")
model3_marginal$year <- factor(model3_marginal$year)
model4_marginal<-model4_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal<-model4_marginal %>% mutate(year="2006-2020")
model4_marginal$year <- factor(model4_marginal$year)
model5_marginal<-model5_marginal %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal<-model5_marginal %>% mutate(year="2006-2020")
model5_marginal$year <- factor(model5_marginal$year)

model1_int<- gam(pm_freq ~ high_ces*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_int<- gam(non_zero_days ~ high_ces*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_int <- gam(peak_pm ~ high_ces*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_int<- gam(smoke_waves ~ high_ces*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_int <- gam(ann_wfpm_avg ~ high_ces*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = urban)

summary(model5_int)

#Marginal means wow this is amazing
model1_int_marginal <- avg_comparisons(model1_int,  variables = "high_ces", by = "year")
model2_int_marginal <- avg_comparisons(model2_int,  variables = "high_ces", by = "year")
model3_int_marginal <- avg_comparisons(model3_int,  variables = "high_ces", by = "year")
model4_int_marginal <- avg_comparisons(model4_int,  variables = "high_ces", by = "year")
model5_int_marginal <- avg_comparisons(model5_int, variables = "high_ces", by = "year")

model1_int_marginal <- model1_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal <- model2_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal <- model3_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal <- model4_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal <- model5_int_marginal %>% dplyr::select(year, estimate, conf.low, conf.high)
#Add 2006-2020 estimates
model1_int_marginal$year <- factor(model1_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot <- rbind(model1_int_marginal,model1_marginal)
model1_int_marginal_plot <- model1_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model2_int_marginal$year <- factor(model2_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot <- rbind(model2_int_marginal,model2_marginal)
model2_int_marginal_plot <- model2_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model3_int_marginal$year <- factor(model3_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot <- rbind(model3_int_marginal,model3_marginal)
model3_int_marginal_plot <- model3_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model4_int_marginal$year <- factor(model4_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot <- rbind(model4_int_marginal,model4_marginal)
model4_int_marginal_plot <- model4_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))
model5_int_marginal$year <- factor(model5_int_marginal$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot <- rbind(model5_int_marginal,model5_marginal)
model5_int_marginal_plot <- model5_int_marginal_plot %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

colnames(model5_int_marginal)
pd <- position_dodge(0.1)

#Weeks over 5 difference
reg1_plot <- ggplot(model1_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod)))+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model1_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous( "Mean difference (weeks)"  ) +
  scale_color_manual("", values=c("black","#D25625")) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3)) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1))

#Non-zero days difference
reg2_plot <-ggplot(model2_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model2_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (days)") +
  ggtitle(bquote(Wildfire~PM[2.5]~">0")) + 
  scale_color_manual("", values=c("black","#D25625")) + 
  
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1))


#Peak wildfire difference
reg3_plot <-ggplot(model3_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model3_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) + scale_color_manual("", values=c("black","#D25625")) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5]))

#Smoke waves
reg4_plot <-ggplot(model4_int_marginal_plot, aes(x=year, y=estimate, color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model4_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (smoke waves)") + 
  scale_color_manual("", values=c("black","#D25625")) + 
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle("Smoke waves")

#5 Annual average wildfire PM
reg5_plot <-ggplot(model5_int_marginal_plot, aes(x=year, y=estimate,color=factor(fullperiod))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=factor(fullperiod)), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.75, linetype="dotted") +
  geom_hline(aes(yintercept=model5_int_marginal_plot[16,2]), linewidth=.75,linetype="dashed", color="#D25625") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_color_manual("", values=c("black","#D25625")) + 
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_x_discrete("") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) +
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5]))

reg1_plot + reg2_plot + reg3_plot + reg4_plot + reg5_plot +  plot_annotation(tag_levels = 'A') +  plot_layout(ncol = 2)
ggsave("longterm-pm/analysis/regression/supp_fig12.png", dpi=300, height=11, width=8, units="in" )

#Save regression results
model5_int_marginal_plot <- model5_int_marginal_plot %>% mutate(exposure=5)
model4_int_marginal_plot <- model4_int_marginal_plot %>% mutate(exposure=4)
model3_int_marginal_plot <- model3_int_marginal_plot %>% mutate(exposure=3)
model2_int_marginal_plot <- model2_int_marginal_plot %>% mutate(exposure=2)
model1_int_marginal_plot <- model1_int_marginal_plot %>% mutate(exposure=1)
model_int_marginal_plot_data <- rbind(model1_int_marginal_plot, model2_int_marginal_plot, model3_int_marginal_plot, model4_int_marginal_plot, model5_int_marginal_plot)

combined <- with(model_int_marginal_plot_data, sprintf('%.2f (%.2f, %.2f)', estimate, conf.low, conf.high))

model_int_marginal_plot_data_2 <- cbind(model_int_marginal_plot_data,combined) 
write_csv(model_int_marginal_plot_data_2, "longterm-pm/analysis/regression/ces_regression_results_urban.csv")

####RACE ANALYSES WITH URBAN/RURAL DESIGNATION####
rural <- wf_reg %>% dplyr::filter(rural==1)
urban <- wf_reg %>% dplyr::filter(rural==0)

# ###Summary across the study period### Supplemental Figure 2
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

#Overall Hispanic models
model1_hisp <- gam(pm_freq ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_hisp <- gam(non_zero_days ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_hisp <- gam(peak_pm ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_hisp <- gam(smoke_waves ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_hisp <- gam(ann_wfpm_avg ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
summary(model1_hisp)
summary(model2_hisp)
summary(model3_hisp)
summary(model4_hisp)
summary(model5_hisp)
model1_marginal_hisp<-avg_comparisons(model1_hisp, variables = list(hisp_p="sd"))
model2_marginal_hisp<-avg_comparisons(model2_hisp, variables = list(hisp_p="sd"))
model3_marginal_hisp<-avg_comparisons(model3_hisp, variables = list(hisp_p="sd"))
model4_marginal_hisp<-avg_comparisons(model4_hisp, variables = list(hisp_p="sd"))
model5_marginal_hisp<-avg_comparisons(model5_hisp, variables = list(hisp_p="sd"))

#Overall AIAN models
model1_aian <- gam(pm_freq ~ nhaian_p +   s(popden, fx=TRUE, k=9)  + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_aian <- gam(non_zero_days ~ nhaian_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_aian <- gam(peak_pm ~ nhaian_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_aian <- gam(smoke_waves ~ nhaian_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_aian <- gam(ann_wfpm_avg ~ nhaian_p + s(popden) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_aian<-avg_comparisons(model1_aian, variables = list(nhaian_p="sd"))
model2_marginal_aian<-avg_comparisons(model2_aian, variables = list(nhaian_p="sd"))
model3_marginal_aian<-avg_comparisons(model3_aian, variables = list(nhaian_p="sd"))
model4_marginal_aian<-avg_comparisons(model4_aian, variables = list(nhaian_p="sd"))
model5_marginal_aian<-avg_comparisons(model5_aian, variables = list(nhaian_p="sd"))

#Overall Asian models
model1_asian <- gam(pm_freq ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_asian <- gam(non_zero_days ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_asian <- gam(peak_pm ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_asian <- gam(smoke_waves ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_asian <- gam(ann_wfpm_avg ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_asian<-avg_comparisons(model1_asian, variables = list(nha_p="sd"))
model2_marginal_asian<-avg_comparisons(model2_asian, variables = list(nha_p="sd"))
model3_marginal_asian<-avg_comparisons(model3_asian, variables = list(nha_p="sd"))
model4_marginal_asian<-avg_comparisons(model4_asian, variables = list(nha_p="sd"))
model5_marginal_asian<-avg_comparisons(model5_asian, variables = list(nha_p="sd"))

#Overall Black models
model1_black <- gam(pm_freq ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_black <- gam(non_zero_days ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_black <- gam(peak_pm ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_black <- gam(smoke_waves ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_black <- gam(ann_wfpm_avg ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_black<-avg_comparisons(model1_black, variables = list(nhb_p="sd"))
model2_marginal_black<-avg_comparisons(model2_black, variables = list(nhb_p="sd"))
model3_marginal_black<-avg_comparisons(model3_black, variables = list(nhb_p="sd"))
model4_marginal_black<-avg_comparisons(model4_black, variables = list(nhb_p="sd"))
model5_marginal_black<-avg_comparisons(model5_black, variables = list(nhb_p="sd"))

#Overall nh2_
model1_2more <- gam(pm_freq ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_2more <- gam(non_zero_days ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_2more <- gam(peak_pm ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_2more <- gam(smoke_waves ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_2more <- gam(ann_wfpm_avg ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_2more<-avg_comparisons(model1_2more, variables = list(nh2ormore="sd"))
model2_marginal_2more<-avg_comparisons(model2_2more, variables = list(nh2ormore="sd"))
model3_marginal_2more<-avg_comparisons(model3_2more, variables = list(nh2ormore="sd"))
model4_marginal_2more<-avg_comparisons(model4_2more, variables = list(nh2ormore="sd"))
model5_marginal_2more<-avg_comparisons(model5_2more, variables = list(nh2ormore="sd"))

#Overall white
model1_white <- gam(pm_freq ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_white <- gam(non_zero_days ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_white <- gam(peak_pm ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_white <- gam(smoke_waves ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_white <- gam(ann_wfpm_avg ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_white<-avg_comparisons(model1_white, variables = list(nhw_p="sd"))
model2_marginal_white<-avg_comparisons(model2_white, variables = list(nhw_p="sd"))
model3_marginal_white<-avg_comparisons(model3_white, variables = list(nhw_p="sd"))
model4_marginal_white<-avg_comparisons(model4_white, variables = list(nhw_p="sd"))
model5_marginal_white<-avg_comparisons(model5_white, variables = list(nhw_p="sd"))


#Make 5 plots with labels from figure 1 (text/code) for each exposure and use colors from figure 3 for race/ethnicity designations
race <- c("hisp_p", "nha_p", "nhaian_p", "nhb_p", "nh2ormore", "nhw_p")
race <- factor(race,  levels=c("hisp_p", "nha_p", "nhaian_p", "nhb_p", "nh2ormore", "nhw_p"),
               labels=c("Hispanic", "NH Asian", "NH American Indian", "NH Black", "NH 2+", "NH white"))

model1_marginal_plot_race <- rbind(model1_marginal_hisp, model1_marginal_asian,  model1_marginal_aian,model1_marginal_black,
                                   model1_marginal_2more, model1_marginal_white)
model1_marginal_plot_race$race <- race
model2_marginal_plot_race <- rbind(model2_marginal_hisp, model2_marginal_asian, model2_marginal_aian, model2_marginal_black,
                                   model2_marginal_2more, model2_marginal_white)
model2_marginal_plot_race$race <- race

model3_marginal_plot_race <- rbind(model3_marginal_hisp, model3_marginal_asian,model3_marginal_aian, model3_marginal_black,
                                   model3_marginal_2more, model3_marginal_white)
model3_marginal_plot_race$race <- race

model4_marginal_plot_race <- rbind(model4_marginal_hisp, model4_marginal_asian, model4_marginal_aian, model4_marginal_black,
                                   model4_marginal_2more, model4_marginal_white)
model4_marginal_plot_race$race <- race

model5_marginal_plot_race <- rbind(model5_marginal_hisp, model5_marginal_asian,model5_marginal_aian, model5_marginal_black,
                                   model5_marginal_2more, model5_marginal_white)
model5_marginal_plot_race$race <- race

pd <- position_dodge(0.1)

#Plotting across the study period the expected change in the exposures for a 1-SD increase in race/ethnicity percent
exp_1_race_plot <- ggplot(model1_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean differences (weeks)") +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3))


exp_2_race_plot <- ggplot(model2_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (days)") +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">0"))

exp_3_race_plot <-ggplot(model3_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5]))

exp_4_race_plot <- ggplot(model4_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (smoke waves)") +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle("Smoke waves")

exp_5_race_plot <-  ggplot(model5_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5]))

exp_1_race_plot + exp_2_race_plot + exp_3_race_plot + exp_4_race_plot + exp_5_race_plot + plot_annotation(tag_levels = 'A') +  plot_layout(nrow = 1)

#####YEAR SPECIFIC BY RACE#####
#Overall Hispanic models
model1_hisp_int <- gam(pm_freq ~ hisp_p*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_hisp_int <- gam(non_zero_days ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_hisp_int <- gam(peak_pm ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_hisp_int <- gam(smoke_waves ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_hisp_int <- gam(ann_wfpm_avg ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)

model1_int_marginal_hisp<-avg_comparisons(model1_hisp_int, variables = list(hisp_p="sd"), by="year")
model2_int_marginal_hisp<-avg_comparisons(model2_hisp_int, variables = list(hisp_p="sd"), by="year")
model3_int_marginal_hisp<-avg_comparisons(model3_hisp_int, variables = list(hisp_p="sd"), by="year")
model4_int_marginal_hisp<-avg_comparisons(model4_hisp_int, variables = list(hisp_p="sd"), by="year")
model5_int_marginal_hisp<-avg_comparisons(model5_hisp_int, variables = list(hisp_p="sd"), by="year")

#Overall AIAN models
model1_aian_int <- gam(pm_freq ~ nhaian_p*year  + s(popden, fx=TRUE, k=4) + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_aian_int <- gam(non_zero_days ~ nhaian_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_aian_int <- gam(peak_pm ~ nhaian_p*year  + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_aian_int <- gam(smoke_waves ~ nhaian_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_aian_int <- gam(ann_wfpm_avg ~ nhaian_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_aian_int<-avg_comparisons(model1_aian_int, variables = list(nhaian_p="sd"), by="year")
model2_marginal_aian_int<-avg_comparisons(model2_aian_int, variables = list(nhaian_p="sd"), by="year")
model3_marginal_aian_int<-avg_comparisons(model3_aian_int, variables = list(nhaian_p="sd"), by="year")
model4_marginal_aian_int<-avg_comparisons(model4_aian_int, variables = list(nhaian_p="sd"), by="year")
model5_marginal_aian_int<-avg_comparisons(model5_aian_int, variables = list(nhaian_p="sd"), by="year")

#Overall Asian models
model1_asian_int <- gam(pm_freq ~ nha_p*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_asian_int <- gam(non_zero_days ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_asian_int <- gam(peak_pm ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_asian_int <- gam(smoke_waves ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_asian_int <- gam(ann_wfpm_avg ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_asian_int<-avg_comparisons(model1_asian_int, variables = list(nha_p="sd"), by="year")
model2_marginal_asian_int<-avg_comparisons(model2_asian_int, variables = list(nha_p="sd"), by="year")
model3_marginal_asian_int<-avg_comparisons(model3_asian_int, variables = list(nha_p="sd"), by="year")
model4_marginal_asian_int<-avg_comparisons(model4_asian_int, variables = list(nha_p="sd"), by="year")
model5_marginal_asian_int<-avg_comparisons(model5_asian_int, variables = list(nha_p="sd"), by="year")

#Overall Black models
model1_black_int <- gam(pm_freq ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_black_int <- gam(non_zero_days ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_black_int <- gam(peak_pm ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_black_int <- gam(smoke_waves ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_black_int <- gam(ann_wfpm_avg ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_black_int<-avg_comparisons(model1_black_int, variables = list(nhb_p="sd"), by="year")
model2_marginal_black_int<-avg_comparisons(model2_black_int, variables = list(nhb_p="sd"), by="year")
model3_marginal_black_int<-avg_comparisons(model3_black_int, variables = list(nhb_p="sd"), by="year")
model4_marginal_black_int<-avg_comparisons(model4_black_int, variables = list(nhb_p="sd"), by="year")
model5_marginal_black_int<-avg_comparisons(model5_black_int, variables = list(nhb_p="sd"), by="year")

#Overall nh2_
model1_2more_int <- gam(pm_freq ~ nh2ormore*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_2more_int <- gam(non_zero_days ~ nh2ormore*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_2more_int <- gam(peak_pm ~ nh2ormore*year  + s(popden, fx=TRUE, k=9) +  s(lon, latit, fx=TRUE, k=21), data = rural)
model4_2more_int <- gam(smoke_waves ~ nh2ormore*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_2more_int <- gam(ann_wfpm_avg ~ nh2ormore*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_2more_int<-avg_comparisons(model1_2more_int, variables = list(nh2ormore="sd"), by="year")
model2_marginal_2more_int<-avg_comparisons(model2_2more_int, variables = list(nh2ormore="sd"), by="year")
model3_marginal_2more_int<-avg_comparisons(model3_2more_int, variables = list(nh2ormore="sd"), by="year")
model4_marginal_2more_int<-avg_comparisons(model4_2more_int, variables = list(nh2ormore="sd"), by="year")
model5_marginal_2more_int<-avg_comparisons(model5_2more_int, variables = list(nh2ormore="sd"), by="year")

#Overall white
model1_white_int <- gam(pm_freq ~ nhw_p*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = rural)
model2_white_int <- gam(non_zero_days ~ nhw_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural, family = nb())
model3_white_int <- gam(peak_pm ~ nhw_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = rural)
model4_white_int <- gam(smoke_waves ~ nhw_p*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = rural)
model5_white_int <- gam(ann_wfpm_avg ~ nhw_p*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = rural)
model1_marginal_white_int<-avg_comparisons(model1_white_int, variables = list(nhw_p="sd"), by="year")
model2_marginal_white_int<-avg_comparisons(model2_white_int, variables = list(nhw_p="sd"), by="year")
model3_marginal_white_int<-avg_comparisons(model3_white_int, variables = list(nhw_p="sd"), by="year")
model4_marginal_white_int<-avg_comparisons(model4_white_int, variables = list(nhw_p="sd"), by="year")
model5_marginal_white_int<-avg_comparisons(model5_white_int, variables = list(nhw_p="sd"), by="year")

#Model 1
model1_marginal_hisp<-model1_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_hisp<-model1_marginal_hisp %>% mutate(year="2006-2020")
model1_marginal_asian<-model1_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_asian<-model1_marginal_asian %>% mutate(year="2006-2020")
model1_marginal_aian<-model1_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_aian<-model1_marginal_aian %>% mutate(year="2006-2020")
model1_marginal_black<-model1_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_black<-model1_marginal_black %>% mutate(year="2006-2020")
model1_marginal_2more<-model1_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_2more<-model1_marginal_2more %>% mutate(year="2006-2020")
model1_marginal_white<-model1_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_white<-model1_marginal_white %>% mutate(year="2006-2020")

model1_int_marginal_hisp <- model1_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_asian <- model1_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_aian <- model1_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_black <- model1_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_2more <- model1_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_white <- model1_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model1_int_marginal_hisp$year <- factor(model1_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_hisp <- rbind(model1_int_marginal_hisp,model1_marginal_hisp)
model1_int_marginal_plot_hisp <- model1_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model1_int_marginal_plot_hisp <- model1_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_asian$year <- factor(model1_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_asian <- rbind(model1_int_marginal_asian,model1_marginal_asian)
model1_int_marginal_plot_asian <- model1_int_marginal_plot_asian %>% mutate(race="NH Asian")
model1_int_marginal_plot_asian <- model1_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_aian$year <- factor(model1_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_aian <- rbind(model1_int_marginal_aian,model1_marginal_aian)
model1_int_marginal_plot_aian <- model1_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model1_int_marginal_plot_aian <- model1_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_black$year <- factor(model1_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_black <- rbind(model1_int_marginal_black,model1_marginal_black)
model1_int_marginal_plot_black <- model1_int_marginal_plot_black %>% mutate(race="NH Black")
model1_int_marginal_plot_black <- model1_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_2more$year <- factor(model1_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_2more <- rbind(model1_int_marginal_2more,model1_marginal_2more)
model1_int_marginal_plot_2more <- model1_int_marginal_plot_2more %>% mutate(race="NH 2+")
model1_int_marginal_plot_2more <- model1_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_white$year <- factor(model1_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_white <- rbind(model1_int_marginal_white,model1_marginal_white)
model1_int_marginal_plot_white <- model1_int_marginal_plot_white %>% mutate(race="NH white")
model1_int_marginal_plot_white <- model1_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model1_int_marginal_plot_all <- rbind(model1_int_marginal_plot_hisp,model1_int_marginal_plot_asian,
                                      model1_int_marginal_plot_aian,model1_int_marginal_plot_black,
                                      model1_int_marginal_plot_2more,model1_int_marginal_plot_white)

model1_int_marginal_plot_all_df <- data.frame(model1_int_marginal_plot_all)
model1_int_marginal_plot_all_df <- model1_int_marginal_plot_all_df %>% mutate(exposure=1)
#Race colors
race_values <- c("#341C5D", "#754E71" , "#DA8940",  "#E19A8F", "#9C372B", "#8D9FD7" )

exp_1_race_plot_all <-  ggplot(model1_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous("Mean difference (weeks)") +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3)) +theme(axis.title.x=element_blank(),
                                                         axis.text.x=element_blank(),
                                                         axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_1_race_plot_all
ggsave("/supp_fig16_regression/race_exp1_meandiff_rural_v2.png", dpi=300, height=11, width=8, units="in" )


#Model 2
model2_marginal_hisp<-model2_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_hisp<-model2_marginal_hisp %>% mutate(year="2006-2020")
model2_marginal_asian<-model2_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_asian<-model2_marginal_asian %>% mutate(year="2006-2020")
model2_marginal_aian<-model2_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_aian<-model2_marginal_aian %>% mutate(year="2006-2020")
model2_marginal_black<-model2_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_black<-model2_marginal_black %>% mutate(year="2006-2020")
model2_marginal_2more<-model2_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_2more<-model2_marginal_2more %>% mutate(year="2006-2020")
model2_marginal_white<-model2_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_white<-model2_marginal_white %>% mutate(year="2006-2020")

model2_int_marginal_hisp <- model2_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_asian <- model2_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_aian <- model2_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_black <- model2_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_2more <- model2_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_white <- model2_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model2_int_marginal_hisp$year <- factor(model2_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_hisp <- rbind(model2_int_marginal_hisp,model2_marginal_hisp)
model2_int_marginal_plot_hisp <- model2_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model2_int_marginal_plot_hisp <- model2_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_asian$year <- factor(model2_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_asian <- rbind(model2_int_marginal_asian,model2_marginal_asian)
model2_int_marginal_plot_asian <- model2_int_marginal_plot_asian %>% mutate(race="NH Asian")
model2_int_marginal_plot_asian <- model2_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_aian$year <- factor(model2_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_aian <- rbind(model2_int_marginal_aian,model2_marginal_aian)
model2_int_marginal_plot_aian <- model2_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model2_int_marginal_plot_aian <- model2_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_black$year <- factor(model2_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_black <- rbind(model2_int_marginal_black,model2_marginal_black)
model2_int_marginal_plot_black <- model2_int_marginal_plot_black %>% mutate(race="NH Black")
model2_int_marginal_plot_black <- model2_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_2more$year <- factor(model2_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_2more <- rbind(model2_int_marginal_2more,model2_marginal_2more)
model2_int_marginal_plot_2more <- model2_int_marginal_plot_2more %>% mutate(race="NH 2+")
model2_int_marginal_plot_2more <- model2_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_white$year <- factor(model2_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_white <- rbind(model2_int_marginal_white,model2_marginal_white)
model2_int_marginal_plot_white <- model2_int_marginal_plot_white %>% mutate(race="NH white")
model2_int_marginal_plot_white <- model2_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model2_int_marginal_plot_all <- rbind(model2_int_marginal_plot_hisp,model2_int_marginal_plot_asian,
                                      model2_int_marginal_plot_aian,model2_int_marginal_plot_black,
                                      model2_int_marginal_plot_2more,model2_int_marginal_plot_white)

model2_int_marginal_plot_all_df <- data.frame(model2_int_marginal_plot_all)
model2_int_marginal_plot_all_df <- model2_int_marginal_plot_all_df %>% mutate(exposure=2)

exp_2_race_plot_all <-  ggplot(model2_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous("Mean difference (days)") +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">0"~"µg/m"^3)) +theme(axis.title.x=element_blank(),
                                                         axis.text.x=element_blank(),
                                                         axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_2_race_plot_all
ggsave("supp_fig16_race_exp2_meandiff_rural_v2.png", dpi=300, height=11, width=8, units="in" )

#Model 3 -- peak PM 
model3_marginal_hisp<-model3_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_hisp<-model3_marginal_hisp %>% mutate(year="2006-2020")
model3_marginal_asian<-model3_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_asian<-model3_marginal_asian %>% mutate(year="2006-2020")
model3_marginal_aian<-model3_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_aian<-model3_marginal_aian %>% mutate(year="2006-2020")
model3_marginal_black<-model3_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_black<-model3_marginal_black %>% mutate(year="2006-2020")
model3_marginal_2more<-model3_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_2more<-model3_marginal_2more %>% mutate(year="2006-2020")
model3_marginal_white<-model3_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_white<-model3_marginal_white %>% mutate(year="2006-2020")

model3_int_marginal_hisp <- model3_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_asian <- model3_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_aian <- model3_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_black <- model3_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_2more <- model3_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_white <- model3_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model3_int_marginal_hisp$year <- factor(model3_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_hisp <- rbind(model3_int_marginal_hisp,model3_marginal_hisp)
model3_int_marginal_plot_hisp <- model3_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model3_int_marginal_plot_hisp <- model3_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_asian$year <- factor(model3_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_asian <- rbind(model3_int_marginal_asian,model3_marginal_asian)
model3_int_marginal_plot_asian <- model3_int_marginal_plot_asian %>% mutate(race="NH Asian")
model3_int_marginal_plot_asian <- model3_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_aian$year <- factor(model3_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_aian <- rbind(model3_int_marginal_aian,model3_marginal_aian)
model3_int_marginal_plot_aian <- model3_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model3_int_marginal_plot_aian <- model3_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_black$year <- factor(model3_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_black <- rbind(model3_int_marginal_black,model3_marginal_black)
model3_int_marginal_plot_black <- model3_int_marginal_plot_black %>% mutate(race="NH Black")
model3_int_marginal_plot_black <- model3_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_2more$year <- factor(model3_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_2more <- rbind(model3_int_marginal_2more,model3_marginal_2more)
model3_int_marginal_plot_2more <- model3_int_marginal_plot_2more %>% mutate(race="NH 2+")
model3_int_marginal_plot_2more <- model3_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_white$year <- factor(model3_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_white <- rbind(model3_int_marginal_white,model3_marginal_white)
model3_int_marginal_plot_white <- model3_int_marginal_plot_white %>% mutate(race="NH white")
model3_int_marginal_plot_white <- model3_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model3_int_marginal_plot_all <- rbind(model3_int_marginal_plot_hisp,model3_int_marginal_plot_asian,
                                      model3_int_marginal_plot_aian,model3_int_marginal_plot_black,
                                      model3_int_marginal_plot_2more,model3_int_marginal_plot_white)

model3_int_marginal_plot_all_df <- data.frame(model3_int_marginal_plot_all)
model3_int_marginal_plot_all_df <- model3_int_marginal_plot_all_df %>% mutate(exposure=3)

exp_3_race_plot_all <-  ggplot(model3_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5])) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_3_race_plot_all
ggsave("supp_fig16_race_exp3_meandiff_rural_v2.png", dpi=300, height=11, width=8, units="in" )

#Model 4 -- smoke waves
model4_marginal_hisp<-model4_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_hisp<-model4_marginal_hisp %>% mutate(year="2006-2020")
model4_marginal_asian<-model4_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_asian<-model4_marginal_asian %>% mutate(year="2006-2020")
model4_marginal_aian<-model4_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_aian<-model4_marginal_aian %>% mutate(year="2006-2020")
model4_marginal_black<-model4_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_black<-model4_marginal_black %>% mutate(year="2006-2020")
model4_marginal_2more<-model4_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_2more<-model4_marginal_2more %>% mutate(year="2006-2020")
model4_marginal_white<-model4_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_white<-model4_marginal_white %>% mutate(year="2006-2020")

model4_int_marginal_hisp <- model4_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_asian <- model4_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_aian <- model4_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_black <- model4_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_2more <- model4_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_white <- model4_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model4_int_marginal_hisp$year <- factor(model4_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_hisp <- rbind(model4_int_marginal_hisp,model4_marginal_hisp)
model4_int_marginal_plot_hisp <- model4_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model4_int_marginal_plot_hisp <- model4_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_asian$year <- factor(model4_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_asian <- rbind(model4_int_marginal_asian,model4_marginal_asian)
model4_int_marginal_plot_asian <- model4_int_marginal_plot_asian %>% mutate(race="NH Asian")
model4_int_marginal_plot_asian <- model4_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_aian$year <- factor(model4_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_aian <- rbind(model4_int_marginal_aian,model4_marginal_aian)
model4_int_marginal_plot_aian <- model4_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model4_int_marginal_plot_aian <- model4_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_black$year <- factor(model4_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_black <- rbind(model4_int_marginal_black,model4_marginal_black)
model4_int_marginal_plot_black <- model4_int_marginal_plot_black %>% mutate(race="NH Black")
model4_int_marginal_plot_black <- model4_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_2more$year <- factor(model4_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_2more <- rbind(model4_int_marginal_2more,model4_marginal_2more)
model4_int_marginal_plot_2more <- model4_int_marginal_plot_2more %>% mutate(race="NH 2+")
model4_int_marginal_plot_2more <- model4_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_white$year <- factor(model4_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_white <- rbind(model4_int_marginal_white,model4_marginal_white)
model4_int_marginal_plot_white <- model4_int_marginal_plot_white %>% mutate(race="NH white")
model4_int_marginal_plot_white <- model4_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model4_int_marginal_plot_all <- rbind(model4_int_marginal_plot_hisp,model4_int_marginal_plot_asian,
                                      model4_int_marginal_plot_aian,model4_int_marginal_plot_black,
                                      model4_int_marginal_plot_2more,model4_int_marginal_plot_white)

model4_int_marginal_plot_all_df <- data.frame(model4_int_marginal_plot_all)
model4_int_marginal_plot_all_df <- model4_int_marginal_plot_all_df %>% mutate(exposure=4)

exp_4_race_plot_all <-  ggplot(model4_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous("Mean difference (smoke waves)") +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle("Smoke waves") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_4_race_plot_all
ggsave("supp_fig16_race_exp4_meandiff_rural_v2.png", dpi=300, height=11, width=8, units="in" )


#Grab main effect from 2006-2020 for model 5
model5_marginal_hisp<-model5_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_hisp<-model5_marginal_hisp %>% mutate(year="2006-2020")
model5_marginal_asian<-model5_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_asian<-model5_marginal_asian %>% mutate(year="2006-2020")
model5_marginal_aian<-model5_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_aian<-model5_marginal_aian %>% mutate(year="2006-2020")
model5_marginal_black<-model5_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_black<-model5_marginal_black %>% mutate(year="2006-2020")
model5_marginal_2more<-model5_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_2more<-model5_marginal_2more %>% mutate(year="2006-2020")
model5_marginal_white<-model5_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_white<-model5_marginal_white %>% mutate(year="2006-2020")

model5_int_marginal_hisp <- model5_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_asian <- model5_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_aian <- model5_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_black <- model5_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_2more <- model5_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_white <- model5_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model5_int_marginal_hisp$year <- factor(model5_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_hisp <- rbind(model5_int_marginal_hisp,model5_marginal_hisp)
model5_int_marginal_plot_hisp <- model5_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model5_int_marginal_plot_hisp <- model5_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_asian$year <- factor(model5_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_asian <- rbind(model5_int_marginal_asian,model5_marginal_asian)
model5_int_marginal_plot_asian <- model5_int_marginal_plot_asian %>% mutate(race="NH Asian")
model5_int_marginal_plot_asian <- model5_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_aian$year <- factor(model5_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_aian <- rbind(model5_int_marginal_aian,model5_marginal_aian)
model5_int_marginal_plot_aian <- model5_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model5_int_marginal_plot_aian <- model5_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_black$year <- factor(model5_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_black <- rbind(model5_int_marginal_black,model5_marginal_black)
model5_int_marginal_plot_black <- model5_int_marginal_plot_black %>% mutate(race="NH Black")
model5_int_marginal_plot_black <- model5_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_2more$year <- factor(model5_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_2more <- rbind(model5_int_marginal_2more,model5_marginal_2more)
model5_int_marginal_plot_2more <- model5_int_marginal_plot_2more %>% mutate(race="NH 2+")
model5_int_marginal_plot_2more <- model5_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_white$year <- factor(model5_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_white <- rbind(model5_int_marginal_white,model5_marginal_white)
model5_int_marginal_plot_white <- model5_int_marginal_plot_white %>% mutate(race="NH white")
model5_int_marginal_plot_white <- model5_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model5_int_marginal_plot_all <- rbind(model5_int_marginal_plot_hisp,model5_int_marginal_plot_asian,
                                      model5_int_marginal_plot_aian,model5_int_marginal_plot_black,
                                      model5_int_marginal_plot_2more,model5_int_marginal_plot_white)

model5_int_marginal_plot_all_df <- data.frame(model5_int_marginal_plot_all)
model5_int_marginal_plot_all_df <- model5_int_marginal_plot_all_df %>% mutate(exposure=5)

exp_5_race_plot_all <-  ggplot(model5_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5])) +theme(axis.title.x=element_blank(),
                                                       axis.text.x=element_blank(),
                                                       axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_5_race_plot_all
ggsave("supp_fig16_race_exp5_meandiff_rural_v2.png", dpi=300, height=11, width=8, units="in" )

#Save estimates from regression for race:
models_int_marginal_plot_all_df <-  rbind(model1_int_marginal_plot_all_df,model2_int_marginal_plot_all_df,model3_int_marginal_plot_all_df,model4_int_marginal_plot_all_df,model5_int_marginal_plot_all_df)
write_csv(models_int_marginal_plot_all_df, "longterm-pm/analysis/regression/race_regression_results_rural.csv")

##SAME BUT FOR URBAN##
#Overall Hispanic models
model1_hisp <- gam(pm_freq ~ hisp_p + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_hisp <- gam(non_zero_days ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_hisp <- gam(peak_pm ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_hisp <- gam(smoke_waves ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_hisp <- gam(ann_wfpm_avg ~ hisp_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
summary(model1_hisp)
summary(model2_hisp)
summary(model3_hisp)
summary(model4_hisp)
summary(model5_hisp)
model1_marginal_hisp<-avg_comparisons(model1_hisp, variables = list(hisp_p="sd"))
model2_marginal_hisp<-avg_comparisons(model2_hisp, variables = list(hisp_p="sd"))
model3_marginal_hisp<-avg_comparisons(model3_hisp, variables = list(hisp_p="sd"))
model4_marginal_hisp<-avg_comparisons(model4_hisp, variables = list(hisp_p="sd"))
model5_marginal_hisp<-avg_comparisons(model5_hisp, variables = list(hisp_p="sd"))

#Overall AIAN models
model1_aian <- gam(pm_freq ~ nhaian_p +   s(popden, fx=TRUE, k=9)  + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_aian <- gam(non_zero_days ~ nhaian_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_aian <- gam(peak_pm ~ nhaian_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_aian <- gam(smoke_waves ~ nhaian_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_aian <- gam(ann_wfpm_avg ~ nhaian_p + s(popden) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model1_marginal_aian<-avg_comparisons(model1_aian, variables = list(nhaian_p="sd"))
model2_marginal_aian<-avg_comparisons(model2_aian, variables = list(nhaian_p="sd"))
model3_marginal_aian<-avg_comparisons(model3_aian, variables = list(nhaian_p="sd"))
model4_marginal_aian<-avg_comparisons(model4_aian, variables = list(nhaian_p="sd"))
model5_marginal_aian<-avg_comparisons(model5_aian, variables = list(nhaian_p="sd"))

#Overall Asian models
model1_asian <- gam(pm_freq ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_asian <- gam(non_zero_days ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_asian <- gam(peak_pm ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_asian <- gam(smoke_waves ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_asian <- gam(ann_wfpm_avg ~ nha_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
summary(model1_asian)
model1_marginal_asian<-avg_comparisons(model1_asian, variables = list(nha_p="sd"))
model2_marginal_asian<-avg_comparisons(model2_asian, variables = list(nha_p="sd"))
model3_marginal_asian<-avg_comparisons(model3_asian, variables = list(nha_p="sd"))
model4_marginal_asian<-avg_comparisons(model4_asian, variables = list(nha_p="sd"))
model5_marginal_asian<-avg_comparisons(model5_asian, variables = list(nha_p="sd"))

#Overall Black models
model1_black <- gam(pm_freq ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_black <- gam(non_zero_days ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_black <- gam(peak_pm ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_black <- gam(smoke_waves ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_black <- gam(ann_wfpm_avg ~ nhb_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model1_marginal_black<-avg_comparisons(model1_black, variables = list(nhb_p="sd"))
model2_marginal_black<-avg_comparisons(model2_black, variables = list(nhb_p="sd"))
model3_marginal_black<-avg_comparisons(model3_black, variables = list(nhb_p="sd"))
model4_marginal_black<-avg_comparisons(model4_black, variables = list(nhb_p="sd"))
model5_marginal_black<-avg_comparisons(model5_black, variables = list(nhb_p="sd"))

#Overall nh2_
model1_2more <- gam(pm_freq ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_2more <- gam(non_zero_days ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_2more <- gam(peak_pm ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_2more <- gam(smoke_waves ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_2more <- gam(ann_wfpm_avg ~ nh2ormore + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model1_marginal_2more<-avg_comparisons(model1_2more, variables = list(nh2ormore="sd"))
model2_marginal_2more<-avg_comparisons(model2_2more, variables = list(nh2ormore="sd"))
model3_marginal_2more<-avg_comparisons(model3_2more, variables = list(nh2ormore="sd"))
model4_marginal_2more<-avg_comparisons(model4_2more, variables = list(nh2ormore="sd"))
model5_marginal_2more<-avg_comparisons(model5_2more, variables = list(nh2ormore="sd"))

#Overall white
model1_white <- gam(pm_freq ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_white <- gam(non_zero_days ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_white <- gam(peak_pm ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_white <- gam(smoke_waves ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_white <- gam(ann_wfpm_avg ~ nhw_p + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model1_marginal_white<-avg_comparisons(model1_white, variables = list(nhw_p="sd"))
model2_marginal_white<-avg_comparisons(model2_white, variables = list(nhw_p="sd"))
model3_marginal_white<-avg_comparisons(model3_white, variables = list(nhw_p="sd"))
model4_marginal_white<-avg_comparisons(model4_white, variables = list(nhw_p="sd"))
model5_marginal_white<-avg_comparisons(model5_white, variables = list(nhw_p="sd"))

#Make 5 plots with labels from figure 1 (text/code) for each exposure and use colors from figure 3 for race/ethnicity designations
race <- c("hisp_p",  "nhaian_p", "nha_p","nhb_p", "nh2ormore", "nhw_p")
race <- factor(race,  levels=c("hisp_p", "nhaian_p", "nha_p", "nhb_p", "nh2ormore", "nhw_p"),
               labels=c("Hispanic",  "NH American Indian", "NH Asian","NH Black", "NH 2+", "NH white"))

model1_marginal_plot_race <- rbind(model1_marginal_hisp, model1_marginal_aian, model1_marginal_asian, model1_marginal_black,
                                   model1_marginal_2more, model1_marginal_white)
model1_marginal_plot_race$race <- race
model2_marginal_plot_race <- rbind(model2_marginal_hisp, model2_marginal_aian, model2_marginal_asian, model2_marginal_black,
                                   model2_marginal_2more, model2_marginal_white)
model2_marginal_plot_race$race <- race

model3_marginal_plot_race <- rbind(model3_marginal_hisp, model3_marginal_aian, model3_marginal_asian, model3_marginal_black,
                                   model3_marginal_2more, model3_marginal_white)
model3_marginal_plot_race$race <- race

model4_marginal_plot_race <- rbind(model4_marginal_hisp, model4_marginal_aian, model4_marginal_asian, model4_marginal_black,
                                   model4_marginal_2more, model4_marginal_white)
model4_marginal_plot_race$race <- race

model5_marginal_plot_race <- rbind(model5_marginal_hisp, model5_marginal_aian, model5_marginal_asian, model5_marginal_black,
                                   model5_marginal_2more, model5_marginal_white)
model5_marginal_plot_race$race <- race

#Plotting across the study period the expected change in the exposures for a 1-SD increase in race/ethnicity percent
exp_1_race_plot <- ggplot(model1_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean differences (weeks)") +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3))


exp_2_race_plot <- ggplot(model2_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (days)") +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">0"))

exp_3_race_plot <-ggplot(model3_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5]))

exp_4_race_plot <- ggplot(model4_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous("Mean difference (smoke waves)") +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle("Smoke waves")

exp_5_race_plot <-  ggplot(model5_marginal_plot_race, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.75) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver"), legend.position="none") +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=met.brewer("Hokusai1", 6)) +
  scale_x_discrete("")+ theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5]))

exp_1_race_plot + exp_2_race_plot + exp_3_race_plot + exp_4_race_plot + exp_5_race_plot + plot_annotation(tag_levels = 'A') +  plot_layout(nrow = 1)

#####YEAR SPECIFIC BY RACE#####
#Overall Hispanic models
model1_hisp_int <- gam(pm_freq ~ hisp_p*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_hisp_int <- gam(non_zero_days ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_hisp_int <- gam(peak_pm ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_hisp_int <- gam(smoke_waves ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_hisp_int <- gam(ann_wfpm_avg ~ hisp_p*year + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)

model1_int_marginal_hisp<-avg_comparisons(model1_hisp_int, variables = list(hisp_p="sd"), by="year")
model2_int_marginal_hisp<-avg_comparisons(model2_hisp_int, variables = list(hisp_p="sd"), by="year")
model3_int_marginal_hisp<-avg_comparisons(model3_hisp_int, variables = list(hisp_p="sd"), by="year")
model4_int_marginal_hisp<-avg_comparisons(model4_hisp_int, variables = list(hisp_p="sd"), by="year")
model5_int_marginal_hisp<-avg_comparisons(model5_hisp_int, variables = list(hisp_p="sd"), by="year")

#Overall AIAN models
model1_aian_int <- gam(pm_freq ~ nhaian_p*year  + s(popden, fx=TRUE, k=4) + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_aian_int <- gam(non_zero_days ~ nhaian_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_aian_int <- gam(peak_pm ~ nhaian_p*year  + s(popden, fx=TRUE, k=9) + year + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_aian_int <- gam(smoke_waves ~ nhaian_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_aian_int <- gam(ann_wfpm_avg ~ nhaian_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model1_marginal_aian_int<-avg_comparisons(model1_aian_int, variables = list(nhaian_p="sd"), by="year")
model2_marginal_aian_int<-avg_comparisons(model2_aian_int, variables = list(nhaian_p="sd"), by="year")
model3_marginal_aian_int<-avg_comparisons(model3_aian_int, variables = list(nhaian_p="sd"), by="year")
model4_marginal_aian_int<-avg_comparisons(model4_aian_int, variables = list(nhaian_p="sd"), by="year")
model5_marginal_aian_int<-avg_comparisons(model5_aian_int, variables = list(nhaian_p="sd"), by="year")

#Overall Asian models
model1_asian_int <- gam(pm_freq ~ nha_p*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_asian_int <- gam(non_zero_days ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_asian_int <- gam(peak_pm ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_asian_int <- gam(smoke_waves ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_asian_int <- gam(ann_wfpm_avg ~ nha_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model1_marginal_asian_int<-avg_comparisons(model1_asian_int, variables = list(nha_p="sd"), by="year")
model2_marginal_asian_int<-avg_comparisons(model2_asian_int, variables = list(nha_p="sd"), by="year")
model3_marginal_asian_int<-avg_comparisons(model3_asian_int, variables = list(nha_p="sd"), by="year")
model4_marginal_asian_int<-avg_comparisons(model4_asian_int, variables = list(nha_p="sd"), by="year")
model5_marginal_asian_int<-avg_comparisons(model5_asian_int, variables = list(nha_p="sd"), by="year")

#Overall Black models
model1_black_int <- gam(pm_freq ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_black_int <- gam(non_zero_days ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_black_int <- gam(peak_pm ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_black_int <- gam(smoke_waves ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_black_int <- gam(ann_wfpm_avg ~ nhb_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model1_marginal_black_int<-avg_comparisons(model1_black_int, variables = list(nhb_p="sd"), by="year")
model2_marginal_black_int<-avg_comparisons(model2_black_int, variables = list(nhb_p="sd"), by="year")
model3_marginal_black_int<-avg_comparisons(model3_black_int, variables = list(nhb_p="sd"), by="year")
model4_marginal_black_int<-avg_comparisons(model4_black_int, variables = list(nhb_p="sd"), by="year")
model5_marginal_black_int<-avg_comparisons(model5_black_int, variables = list(nhb_p="sd"), by="year")

#Overall nh2_
model1_2more_int <- gam(pm_freq ~ nh2ormore*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_2more_int <- gam(non_zero_days ~ nh2ormore*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_2more_int <- gam(peak_pm ~ nh2ormore*year  + s(popden, fx=TRUE, k=9) +  s(lon, latit, fx=TRUE, k=21), data = urban)
model4_2more_int <- gam(smoke_waves ~ nh2ormore*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_2more_int <- gam(ann_wfpm_avg ~ nh2ormore*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model1_marginal_2more_int<-avg_comparisons(model1_2more_int, variables = list(nh2ormore="sd"), by="year")
model2_marginal_2more_int<-avg_comparisons(model2_2more_int, variables = list(nh2ormore="sd"), by="year")
model3_marginal_2more_int<-avg_comparisons(model3_2more_int, variables = list(nh2ormore="sd"), by="year")
model4_marginal_2more_int<-avg_comparisons(model4_2more_int, variables = list(nh2ormore="sd"), by="year")
model5_marginal_2more_int<-avg_comparisons(model5_2more_int, variables = list(nh2ormore="sd"), by="year")

#Overall white
model1_white_int <- gam(pm_freq ~ nhw_p*year + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = urban)
model2_white_int <- gam(non_zero_days ~ nhw_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban, family = nb())
model3_white_int <- gam(peak_pm ~ nhw_p*year  + s(popden, fx=TRUE, k=9)  + s(lon, latit, fx=TRUE, k=21), data = urban)
model4_white_int <- gam(smoke_waves ~ nhw_p*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = urban)
model5_white_int <- gam(ann_wfpm_avg ~ nhw_p*year  + s(popden, fx=TRUE, k=9) + s(lon, latit, fx=TRUE, k=21), data = urban)
model1_marginal_white_int<-avg_comparisons(model1_white_int, variables = list(nhw_p="sd"), by="year")
model2_marginal_white_int<-avg_comparisons(model2_white_int, variables = list(nhw_p="sd"), by="year")
model3_marginal_white_int<-avg_comparisons(model3_white_int, variables = list(nhw_p="sd"), by="year")
model4_marginal_white_int<-avg_comparisons(model4_white_int, variables = list(nhw_p="sd"), by="year")
model5_marginal_white_int<-avg_comparisons(model5_white_int, variables = list(nhw_p="sd"), by="year")

#Model 1
model1_marginal_hisp<-model1_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_hisp<-model1_marginal_hisp %>% mutate(year="2006-2020")
model1_marginal_asian<-model1_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_asian<-model1_marginal_asian %>% mutate(year="2006-2020")
model1_marginal_aian<-model1_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_aian<-model1_marginal_aian %>% mutate(year="2006-2020")
model1_marginal_black<-model1_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_black<-model1_marginal_black %>% mutate(year="2006-2020")
model1_marginal_2more<-model1_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_2more<-model1_marginal_2more %>% mutate(year="2006-2020")
model1_marginal_white<-model1_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model1_marginal_white<-model1_marginal_white %>% mutate(year="2006-2020")

model1_int_marginal_hisp <- model1_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_asian <- model1_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_aian <- model1_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_black <- model1_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_2more <- model1_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model1_int_marginal_white <- model1_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model1_int_marginal_hisp$year <- factor(model1_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_hisp <- rbind(model1_int_marginal_hisp,model1_marginal_hisp)
model1_int_marginal_plot_hisp <- model1_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model1_int_marginal_plot_hisp <- model1_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_asian$year <- factor(model1_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_asian <- rbind(model1_int_marginal_asian,model1_marginal_asian)
model1_int_marginal_plot_asian <- model1_int_marginal_plot_asian %>% mutate(race="NH Asian")
model1_int_marginal_plot_asian <- model1_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_aian$year <- factor(model1_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_aian <- rbind(model1_int_marginal_aian,model1_marginal_aian)
model1_int_marginal_plot_aian <- model1_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model1_int_marginal_plot_aian <- model1_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_black$year <- factor(model1_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_black <- rbind(model1_int_marginal_black,model1_marginal_black)
model1_int_marginal_plot_black <- model1_int_marginal_plot_black %>% mutate(race="NH Black")
model1_int_marginal_plot_black <- model1_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_2more$year <- factor(model1_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_2more <- rbind(model1_int_marginal_2more,model1_marginal_2more)
model1_int_marginal_plot_2more <- model1_int_marginal_plot_2more %>% mutate(race="NH 2+")
model1_int_marginal_plot_2more <- model1_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model1_int_marginal_white$year <- factor(model1_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model1_int_marginal_plot_white <- rbind(model1_int_marginal_white,model1_marginal_white)
model1_int_marginal_plot_white <- model1_int_marginal_plot_white %>% mutate(race="NH white")
model1_int_marginal_plot_white <- model1_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model1_int_marginal_plot_all <- rbind(model1_int_marginal_plot_hisp,model1_int_marginal_plot_asian,
                                      model1_int_marginal_plot_aian,model1_int_marginal_plot_black,
                                      model1_int_marginal_plot_2more,model1_int_marginal_plot_white)

model1_int_marginal_plot_all_df <- data.frame(model1_int_marginal_plot_all)
model1_int_marginal_plot_all_df <- model1_int_marginal_plot_all_df %>% mutate(exposure=1)
#Race colors
race_values <- c("#341C5D", "#754E71" , "#DA8940",  "#E19A8F", "#9C372B", "#8D9FD7" )

exp_1_race_plot_all <-  ggplot(model1_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous("Mean difference (weeks)") +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">5"~"µg/m"^3)) +theme(axis.title.x=element_blank(),
                                                         axis.text.x=element_blank(),
                                                         axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_1_race_plot_all
ggsave("supp_fig15_race_exp1_meandiff_urban_v2.png", dpi=300, height=11, width=8, units="in" )


#Model 2
model2_marginal_hisp<-model2_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_hisp<-model2_marginal_hisp %>% mutate(year="2006-2020")
model2_marginal_asian<-model2_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_asian<-model2_marginal_asian %>% mutate(year="2006-2020")
model2_marginal_aian<-model2_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_aian<-model2_marginal_aian %>% mutate(year="2006-2020")
model2_marginal_black<-model2_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_black<-model2_marginal_black %>% mutate(year="2006-2020")
model2_marginal_2more<-model2_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_2more<-model2_marginal_2more %>% mutate(year="2006-2020")
model2_marginal_white<-model2_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model2_marginal_white<-model2_marginal_white %>% mutate(year="2006-2020")

model2_int_marginal_hisp <- model2_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_asian <- model2_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_aian <- model2_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_black <- model2_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_2more <- model2_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model2_int_marginal_white <- model2_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model2_int_marginal_hisp$year <- factor(model2_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_hisp <- rbind(model2_int_marginal_hisp,model2_marginal_hisp)
model2_int_marginal_plot_hisp <- model2_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model2_int_marginal_plot_hisp <- model2_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_asian$year <- factor(model2_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_asian <- rbind(model2_int_marginal_asian,model2_marginal_asian)
model2_int_marginal_plot_asian <- model2_int_marginal_plot_asian %>% mutate(race="NH Asian")
model2_int_marginal_plot_asian <- model2_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_aian$year <- factor(model2_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_aian <- rbind(model2_int_marginal_aian,model2_marginal_aian)
model2_int_marginal_plot_aian <- model2_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model2_int_marginal_plot_aian <- model2_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_black$year <- factor(model2_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_black <- rbind(model2_int_marginal_black,model2_marginal_black)
model2_int_marginal_plot_black <- model2_int_marginal_plot_black %>% mutate(race="NH Black")
model2_int_marginal_plot_black <- model2_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_2more$year <- factor(model2_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_2more <- rbind(model2_int_marginal_2more,model2_marginal_2more)
model2_int_marginal_plot_2more <- model2_int_marginal_plot_2more %>% mutate(race="NH 2+")
model2_int_marginal_plot_2more <- model2_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model2_int_marginal_white$year <- factor(model2_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model2_int_marginal_plot_white <- rbind(model2_int_marginal_white,model2_marginal_white)
model2_int_marginal_plot_white <- model2_int_marginal_plot_white %>% mutate(race="NH white")
model2_int_marginal_plot_white <- model2_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model2_int_marginal_plot_all <- rbind(model2_int_marginal_plot_hisp,model2_int_marginal_plot_asian,
                                      model2_int_marginal_plot_aian,model2_int_marginal_plot_black,
                                      model2_int_marginal_plot_2more,model2_int_marginal_plot_white)

model2_int_marginal_plot_all_df <- data.frame(model2_int_marginal_plot_all)
model2_int_marginal_plot_all_df <- model2_int_marginal_plot_all_df %>% mutate(exposure=2)

exp_2_race_plot_all <-  ggplot(model2_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous("Mean difference (days)") +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Wildfire~PM[2.5]~">0"~"µg/m"^3)) +theme(axis.title.x=element_blank(),
                                                         axis.text.x=element_blank(),
                                                         axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_2_race_plot_all
ggsave("supp_fig15_race_exp2_meandiff_urban_v2.png", dpi=300, height=11, width=8, units="in" )

#Model 3 -- peak PM 
model3_marginal_hisp<-model3_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_hisp<-model3_marginal_hisp %>% mutate(year="2006-2020")
model3_marginal_asian<-model3_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_asian<-model3_marginal_asian %>% mutate(year="2006-2020")
model3_marginal_aian<-model3_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_aian<-model3_marginal_aian %>% mutate(year="2006-2020")
model3_marginal_black<-model3_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_black<-model3_marginal_black %>% mutate(year="2006-2020")
model3_marginal_2more<-model3_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_2more<-model3_marginal_2more %>% mutate(year="2006-2020")
model3_marginal_white<-model3_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model3_marginal_white<-model3_marginal_white %>% mutate(year="2006-2020")

model3_int_marginal_hisp <- model3_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_asian <- model3_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_aian <- model3_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_black <- model3_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_2more <- model3_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model3_int_marginal_white <- model3_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model3_int_marginal_hisp$year <- factor(model3_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_hisp <- rbind(model3_int_marginal_hisp,model3_marginal_hisp)
model3_int_marginal_plot_hisp <- model3_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model3_int_marginal_plot_hisp <- model3_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_asian$year <- factor(model3_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_asian <- rbind(model3_int_marginal_asian,model3_marginal_asian)
model3_int_marginal_plot_asian <- model3_int_marginal_plot_asian %>% mutate(race="NH Asian")
model3_int_marginal_plot_asian <- model3_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_aian$year <- factor(model3_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_aian <- rbind(model3_int_marginal_aian,model3_marginal_aian)
model3_int_marginal_plot_aian <- model3_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model3_int_marginal_plot_aian <- model3_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_black$year <- factor(model3_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_black <- rbind(model3_int_marginal_black,model3_marginal_black)
model3_int_marginal_plot_black <- model3_int_marginal_plot_black %>% mutate(race="NH Black")
model3_int_marginal_plot_black <- model3_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_2more$year <- factor(model3_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_2more <- rbind(model3_int_marginal_2more,model3_marginal_2more)
model3_int_marginal_plot_2more <- model3_int_marginal_plot_2more %>% mutate(race="NH 2+")
model3_int_marginal_plot_2more <- model3_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model3_int_marginal_white$year <- factor(model3_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model3_int_marginal_plot_white <- rbind(model3_int_marginal_white,model3_marginal_white)
model3_int_marginal_plot_white <- model3_int_marginal_plot_white %>% mutate(race="NH white")
model3_int_marginal_plot_white <- model3_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model3_int_marginal_plot_all <- rbind(model3_int_marginal_plot_hisp,model3_int_marginal_plot_asian,
                                      model3_int_marginal_plot_aian,model3_int_marginal_plot_black,
                                      model3_int_marginal_plot_2more,model3_int_marginal_plot_white)

model3_int_marginal_plot_all_df <- data.frame(model3_int_marginal_plot_all)
model3_int_marginal_plot_all_df <- model3_int_marginal_plot_all_df %>% mutate(exposure=3)

exp_3_race_plot_all <-  ggplot(model3_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Peak~week~mean~wildfire~PM[2.5])) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_3_race_plot_all
ggsave("supp_fig15_race_exp3_meandiff_urban_v2.png", dpi=300, height=11, width=8, units="in" )

#Model 4 -- smoke waves
model4_marginal_hisp<-model4_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_hisp<-model4_marginal_hisp %>% mutate(year="2006-2020")
model4_marginal_asian<-model4_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_asian<-model4_marginal_asian %>% mutate(year="2006-2020")
model4_marginal_aian<-model4_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_aian<-model4_marginal_aian %>% mutate(year="2006-2020")
model4_marginal_black<-model4_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_black<-model4_marginal_black %>% mutate(year="2006-2020")
model4_marginal_2more<-model4_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_2more<-model4_marginal_2more %>% mutate(year="2006-2020")
model4_marginal_white<-model4_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model4_marginal_white<-model4_marginal_white %>% mutate(year="2006-2020")

model4_int_marginal_hisp <- model4_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_asian <- model4_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_aian <- model4_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_black <- model4_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_2more <- model4_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model4_int_marginal_white <- model4_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model4_int_marginal_hisp$year <- factor(model4_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_hisp <- rbind(model4_int_marginal_hisp,model4_marginal_hisp)
model4_int_marginal_plot_hisp <- model4_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model4_int_marginal_plot_hisp <- model4_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_asian$year <- factor(model4_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_asian <- rbind(model4_int_marginal_asian,model4_marginal_asian)
model4_int_marginal_plot_asian <- model4_int_marginal_plot_asian %>% mutate(race="NH Asian")
model4_int_marginal_plot_asian <- model4_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_aian$year <- factor(model4_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_aian <- rbind(model4_int_marginal_aian,model4_marginal_aian)
model4_int_marginal_plot_aian <- model4_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model4_int_marginal_plot_aian <- model4_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_black$year <- factor(model4_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_black <- rbind(model4_int_marginal_black,model4_marginal_black)
model4_int_marginal_plot_black <- model4_int_marginal_plot_black %>% mutate(race="NH Black")
model4_int_marginal_plot_black <- model4_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_2more$year <- factor(model4_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_2more <- rbind(model4_int_marginal_2more,model4_marginal_2more)
model4_int_marginal_plot_2more <- model4_int_marginal_plot_2more %>% mutate(race="NH 2+")
model4_int_marginal_plot_2more <- model4_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model4_int_marginal_white$year <- factor(model4_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model4_int_marginal_plot_white <- rbind(model4_int_marginal_white,model4_marginal_white)
model4_int_marginal_plot_white <- model4_int_marginal_plot_white %>% mutate(race="NH white")
model4_int_marginal_plot_white <- model4_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model4_int_marginal_plot_all <- rbind(model4_int_marginal_plot_hisp,model4_int_marginal_plot_asian,
                                      model4_int_marginal_plot_aian,model4_int_marginal_plot_black,
                                      model4_int_marginal_plot_2more,model4_int_marginal_plot_white)

model4_int_marginal_plot_all_df <- data.frame(model4_int_marginal_plot_all)
model4_int_marginal_plot_all_df <- model4_int_marginal_plot_all_df %>% mutate(exposure=4)

exp_4_race_plot_all <-  ggplot(model4_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous("Mean difference (smoke waves)") +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle("Smoke waves") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_4_race_plot_all
ggsave("supp_fig15_race_exp4_meandiff_urban_v2.png", dpi=300, height=11, width=8, units="in" )


#Grab main effect from 2006-2020 for model 5
model5_marginal_hisp<-model5_marginal_hisp %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_hisp<-model5_marginal_hisp %>% mutate(year="2006-2020")
model5_marginal_asian<-model5_marginal_asian %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_asian<-model5_marginal_asian %>% mutate(year="2006-2020")
model5_marginal_aian<-model5_marginal_aian %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_aian<-model5_marginal_aian %>% mutate(year="2006-2020")
model5_marginal_black<-model5_marginal_black %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_black<-model5_marginal_black %>% mutate(year="2006-2020")
model5_marginal_2more<-model5_marginal_2more %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_2more<-model5_marginal_2more %>% mutate(year="2006-2020")
model5_marginal_white<-model5_marginal_white %>% dplyr::select(estimate, conf.low, conf.high)
model5_marginal_white<-model5_marginal_white %>% mutate(year="2006-2020")

model5_int_marginal_hisp <- model5_int_marginal_hisp %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_asian <- model5_marginal_asian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_aian <- model5_marginal_aian_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_black <- model5_marginal_black_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_2more <- model5_marginal_2more_int %>% dplyr::select(year, estimate, conf.low, conf.high)
model5_int_marginal_white <- model5_marginal_white_int %>% dplyr::select(year, estimate, conf.low, conf.high)

#Add 2006-2020 estimates
model5_int_marginal_hisp$year <- factor(model5_int_marginal_hisp$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_hisp <- rbind(model5_int_marginal_hisp,model5_marginal_hisp)
model5_int_marginal_plot_hisp <- model5_int_marginal_plot_hisp %>% mutate(race="Hispanic")
model5_int_marginal_plot_hisp <- model5_int_marginal_plot_hisp %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_asian$year <- factor(model5_int_marginal_asian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_asian <- rbind(model5_int_marginal_asian,model5_marginal_asian)
model5_int_marginal_plot_asian <- model5_int_marginal_plot_asian %>% mutate(race="NH Asian")
model5_int_marginal_plot_asian <- model5_int_marginal_plot_asian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_aian$year <- factor(model5_int_marginal_aian$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_aian <- rbind(model5_int_marginal_aian,model5_marginal_aian)
model5_int_marginal_plot_aian <- model5_int_marginal_plot_aian %>% mutate(race="NH American Indian")
model5_int_marginal_plot_aian <- model5_int_marginal_plot_aian %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_black$year <- factor(model5_int_marginal_black$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_black <- rbind(model5_int_marginal_black,model5_marginal_black)
model5_int_marginal_plot_black <- model5_int_marginal_plot_black %>% mutate(race="NH Black")
model5_int_marginal_plot_black <- model5_int_marginal_plot_black %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_2more$year <- factor(model5_int_marginal_2more$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_2more <- rbind(model5_int_marginal_2more,model5_marginal_2more)
model5_int_marginal_plot_2more <- model5_int_marginal_plot_2more %>% mutate(race="NH 2+")
model5_int_marginal_plot_2more <- model5_int_marginal_plot_2more %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))

model5_int_marginal_white$year <- factor(model5_int_marginal_white$year, levels=c("2006-2020", 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))
model5_int_marginal_plot_white <- rbind(model5_int_marginal_white,model5_marginal_white)
model5_int_marginal_plot_white <- model5_int_marginal_plot_white %>% mutate(race="NH white")
model5_int_marginal_plot_white <- model5_int_marginal_plot_white %>% mutate(fullperiod=ifelse(year=="2006-2020",1,0))


model5_int_marginal_plot_all <- rbind(model5_int_marginal_plot_hisp,model5_int_marginal_plot_asian,
                                      model5_int_marginal_plot_aian,model5_int_marginal_plot_black,
                                      model5_int_marginal_plot_2more,model5_int_marginal_plot_white)

model5_int_marginal_plot_all_df <- data.frame(model5_int_marginal_plot_all)
model5_int_marginal_plot_all_df <- model5_int_marginal_plot_all_df %>% mutate(exposure=5)

exp_5_race_plot_all <-  ggplot(model5_int_marginal_plot_all_df, aes(x=race, y=estimate, color=race)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour=race), width=0, position=pd) +
  geom_line(position=pd, aes(group=race, color=race)) +
  geom_point(position=pd, size=1.75) +
  facet_wrap(~year) +
  geom_hline(aes(yintercept=0), linewidth=.25, linetype="dotted") +
  theme_classic(base_size = 10)  +
  theme(axis.text = element_text(size = 10),legend.title=element_text(size=13), 
        legend.text=element_text(size=11)) + 
  theme(strip.text.y = element_blank(),
        text=element_text(family="Gulliver")) +
  scale_y_continuous(bquote(Mean~difference~(µg/m^3))) +
  scale_color_manual("Race/ethnicity", values=race_values) +
  labs(x="Race/ethnicty") + theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust=1)) + 
  ggtitle(bquote(Annual~mean~wildfire~PM[2.5])) +theme(axis.title.x=element_blank(),
                                                       axis.text.x=element_blank(),
                                                       axis.ticks.x=element_blank())+ theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = .5))
exp_5_race_plot_all
ggsave("supp_fig15_race_exp5_meandiff_urban_v2.png", dpi=300, height=11, width=8, units="in" )

#Save estimates from regression for race:
models_int_marginal_plot_all_df <-  rbind(model1_int_marginal_plot_all_df,model2_int_marginal_plot_all_df,model3_int_marginal_plot_all_df,model4_int_marginal_plot_all_df,model5_int_marginal_plot_all_df)
write_csv(models_int_marginal_plot_all_df, "longterm-pm/analysis/regression/race_regression_results_urban.csv")


