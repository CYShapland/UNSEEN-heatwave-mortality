################################################################################
# Create figures for manuscript
# Last updated: 05/09/2025
################################################################################

#load libraries
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(dplyr)
library(reshape2)
library(patchwork)

# set directories
Dir<-"YourDir"
setwd(Dir)
output_dir<-paste0(Dir, "YourDir")

# LOAD functions
source(paste0(Dir,"scripts/Functions.R"))

#####################################
# Preparing data for Figure 4 and 5
#####################################

# Input: Proportion: DLNM exposure-mortality model

# lead time
leadtime_range<-c("2021-06-26", "2021-06-22", "2021-06-18")
leadtime_name<-c(3, 7, 11)
leadtime_lab<-c("three", "seven", "eleven")
total_ensemble<-c(51,51,251)

# LOAD OBSERVED DATA - DAILY MORTALITY & TEMPERATURE-SERIES BETWEEN 1981-2015 in regions in Canada
regions <- c("abbotsford", "vancouver")
reg <- regions[1]
reg1 <- gsub(" ", "", reg, fixed=TRUE)

forecast<-c("Pi", "ENS","Incr")
forecast_lab<-c("pre-indust.","current","future")

plot_list<-list()

for (l in 1:3){
  #l<-3
  forcast_all<-list()
  forcast_temp_all<-list()
  
  for(i in 1:3){
    #i<-2
    #load results from DLNM exposure-mortality model
    load(paste0(dtaDir, forecast[i],"/daily_attributable_summer3dayslag_MMT50to100_TAVG_scalingfactor_",
                leadtime_range[l],"_nsim1000",reg1,".RData"))
    ansim_allens[[1]]$temp
    
    #Convert negative values to zero
    varlist<-paste0("sim", 1:100)
    allvar<-c("est", varlist, "expdeaths")
    convert_neg<-function(x){ifelse(x<0,0,x)}
    ansim_allens_noNeg<-lapply(ansim_allens, function(x) x %>% mutate_at(allvar, convert_neg))
    
    #extract heatwave days
    HWPeriod<-lapply(ansim_allens_noNeg, function(x) x[x$date >= "2021-06-26" & x$date <= "2021-06-30",])
    
    # estimate sum of deaths and average temp during heatwave period
    sum_deaths <- do.call("cbind",lapply(HWPeriod,function(x) sum(x[, 'est'])))
    avg_temp <- do.call("cbind",lapply(HWPeriod,function(x) mean(x[, 2])))
    
    # sd of estimate
    varlist<-paste0("sim", 1:1000)
    sum_deaths_allsim <- lapply(HWPeriod,function(x) apply(x[varlist],2,sum))
    sd_deaths<-as.data.frame(do.call("cbind",lapply(sum_deaths_allsim, function(x) sd(x))))
    
    #expected all-cause deaths
    exp_deaths<-mean(do.call("cbind",lapply(ansim_allens, function(x) sum(x[x$date >= "2021-06-26" & x$date <= "2021-06-30"
                                                                            ,"expdeaths"]))))
    
    #rank computes return time
    forcast_all[[i]] <- data.frame(deaths=as.numeric(sum_deaths), 
                                   stdev=as.numeric(sd_deaths), 
                                   rank=total_ensemble[l]/rank(-sum_deaths), 
                                   prop=as.numeric(sum_deaths)/exp_deaths,
                                   condition=rep(forecast_lab[i], length(sum_deaths)))
    forcast_temp_all[[i]] <- data.frame(avgtemp=as.numeric(avg_temp), 
                                        rank=total_ensemble[l]/rank(-avg_temp),
                                        condition=rep(forecast_lab[i], length(sum_deaths)))
    
  }
  
  #show predicted heat-related number of deaths from ERA5
  load(paste0(dtaDir, "ERA5_daily_attributable_deaths_EDVancouver2008to2015_summer3dayslag_MMT50to100_TAVG_nsim100",
              "_",reg1,".RData"))
  ansim_noNeg<-ansim %>% mutate_at(allvar, convert_neg)
  ERA5_deaths<-sum(ansim_noNeg[ansim_noNeg$date >= "2021-06-26" & ansim_noNeg$date <= "2021-06-30","est"])
  avg_obtemp<-mean(ansim[ansim$date >= "2021-06-26" & ansim$date <= "2021-06-30","temp"])
  
  dta_plot<-do.call("rbind", forcast_all)
  dta_plot['leadtime']<-rep(leadtime_lab[l], nrow(dta_plot))
  tbl_hw_obs<-dta_plot %>% dplyr:::group_by(condition) %>% dplyr:::filter(abs(prop - ERA5_deaths/exp_deaths) == min(abs(prop - ERA5_deaths/exp_deaths)))
  
  test<-dta_plot %>% dplyr:::filter(condition=="pre-indust.")
  
  dta_temp_plot<-do.call("rbind", forcast_temp_all)
  tbl_temp_obs<-dta_temp_plot %>% dplyr:::group_by(condition) %>% dplyr:::filter(abs(avgtemp - avg_obtemp) == min(abs(avgtemp - avg_obtemp)))
  
  print("proportion of Heat-related deaths to be the same as observed event")
  print(tbl_hw_obs)
  print("temperature change to be the same as observed event")
  print(tbl_temp_obs)
  
  plot_list[[l]]<-dta_plot
  
  rm(forcast_all, forcast_temp_all, tbl_hw_obs, tbl_temp_obs)
}

plot_all<-do.call("rbind", plot_list)

bar_plt<-data.frame(conditions=c("preInd", "future"), 
                    days3=rep(0,2), days7=rep(0,2), 
                    days11=rep(0,2)) 

unseen_plt<-data.frame(conditions="current", 
                       days3=rep(0,1), days7=rep(0,1), 
                       days11=rep(0,1)) 

for (i in 1:3) {
  #i<-3
  # formating data to plot % attribution
  obsAga<-plot_list[[i]] %>% filter(condition=="current") %>% dplyr:::filter(abs(deaths - ERA5_deaths) == min(abs(deaths - ERA5_deaths)))
  obsAga_future<-plot_list[[i]] %>% filter(condition=="future") %>% dplyr:::filter(abs(rank - obsAga$rank) == min(abs(rank - obsAga$rank)))
  obsAga_preInd<-plot_list[[i]] %>% filter(condition=="pre-indust.") %>% dplyr:::filter(abs(rank - obsAga$rank) == min(abs(rank - obsAga$rank)))
  
  bar_plt[2,i+1]<-((obsAga_future$death-obsAga$death)/obsAga$death)*100
  bar_plt[1,i+1]<-((obsAga_preInd$death-obsAga$death)/obsAga$death)*100
  
  print(((obsAga$death-obsAga_preInd$death)/exp_deaths)*100)
  
  # formating data to plot unseen deaths
  data_sub<-plot_list[[i]] %>% filter(condition=="current") %>% arrange(desc(rank))
  unseen_plt[1, i+1]<-((mean(data_sub[1:3,"deaths"])-obsAga$death)/obsAga$death)*100
  
}
#####################################
# Plotting Figure 4
#####################################

ylim_reg<-ifelse(reg1=="vancouver", 30,15)
subplot_label<-c("3-day", "7-day","11-day")
plot_predDeath<-list()
for (i in 1:3){
  plot_predDeath[[i]]<- ggplot(subset(plot_all, leadtime==leadtime_lab[i]), 
                               aes(rank, deaths, color = condition)) + 
    geom_point(size = 1.2) +
    scale_x_log10() +
    geom_hline(yintercept=ERA5_deaths, linetype="dashed") +
    xlab("conditional exceedance probability") +
    ylab("heat-related deaths") +
    labs(title=subplot_label[i]) +
    ylim(0,ylim_reg) +
    theme(plot.title = element_text(size=12, hjust = 0.5), 
          axis.title=element_text(size=12),
          legend.title = element_text(size = 12),
          legend.position="bottom") +
    scale_color_manual(values=c("grey", "red", "blue"))
}

ggsave(
  paste0("/DLNM/PredHeatdeaths_scalingfactor_", reg1, ".pdf"),
  plot_predDeath[[1]]+ plot_predDeath[[2]]+ plot_predDeath[[3]]+
    plot_annotation(tag_levels = "A")  +
    plot_layout(axis_titles = "collect", guides = 'collect')  &
    theme(legend.position = "bottom", plot.tag = element_text(face = 'bold')), 
  width=11, height=8.5)

#####################################
# Plotting Figure 5
#####################################

gfg<-melt(bar_plt, id.vars="conditions")
unseen<-melt(unseen_plt, id.vars="conditions")

ylim_unseen<-ifelse(reg1=="vancouver", 40,60)

plot_unseen<-ggplot(unseen,aes(x = variable, y =value, fill = conditions)) + 
  geom_bar(stat = "identity", position = "dodge", fill="grey")+ 
  scale_x_discrete(labels=c('3-day', '7-day', '11-day')) +
  xlab("lead time") +
  ylab("% unseen heat mortality") +
  theme_light() +
  theme(plot.title = element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ylim(0,ylim_unseen) 

plot_attribution<-ggplot(gfg,aes(x = variable, y =value, fill = conditions)) + 
  geom_bar(stat = "identity", position = "dodge")+ 
  scale_x_discrete(labels=c('3-day', '7-day', '11-day')) +
  xlab("lead time") +
  ylab("% Attribution") +
  ylim(-30,30) +
  theme_light() +
  theme(plot.title = element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size = 11), 
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position="bottom") +
  scale_fill_manual(name="condition",labels = c("future", "pre-indust."),values=c( "red","blue"))

ggsave(
  paste0("/DLNM/Unseen_PercAttr_", reg1, ".pdf"),
  plot_unseen + plot_attribution + plot_annotation(tag_levels = "A")  &
    theme(plot.tag = element_text(face = 'bold')), 
  width=11, height=8.5) 



