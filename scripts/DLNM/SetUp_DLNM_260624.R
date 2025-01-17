################################################################################
# Title: Setup data for DLNM
#
# Last updated: 24/11/2023
# Author: Chin Yang Shapland
################################################################################

#####################################################
##                                   ################
## Load library, directory and data  ################
##                                   ################
#####################################################

##rm(list=ls())
library(splines)
library(Hmisc)
library(BMA)
library(coda)
library(RColorBrewer)
library(dplyr)
library(tidyr)

wkDir<-"M:/working/data/Health/Canada/Mortality/cleaned"
setwd(wkDir)

### load data ###
load("data_Canada_HWcities_151223.RData")
lapply(dat, dim)

#####################################################
##                                   ################
## Reformatting                      ################
##                                   ################
#####################################################

### variables ###
cities <- c("abbotsford", "victoria", "vancouver","lytton")
ncities <- length(cities)

### subset to all cause mortality data with age and without separating gender ###
for(i in 1:ncities) {
  #i<-1
  #head(dat[[i]])
  dat[[i]]<-dat[[i]] %>% select(!starts_with("non_acc") & 
                                !starts_with("digestive") & 
                                !starts_with("external") & 
                                !starts_with("nervous")) %>%
                          mutate(death=all_mortality, 
                                death65plus=rowSums(across(c(all_mortality_65_74, all_mortality_75))),
                                death65plus_female=rowSums(across(c(all_mortality_65_74_female, 
                                                                    all_mortality_75_female))),
                                death65plus_male=rowSums(across(c(all_mortality_65_74_male, 
                                                                  all_mortality_75_male))))
}
health_dat <- dat
rm(dat)

save(health_dat,file="M:/working/data/Health/Canada/Mortality/cleaned/DLNM/Canadian_HWcity_DLNMdata_260624.RData")
