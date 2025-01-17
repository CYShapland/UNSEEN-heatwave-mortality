################################################################################
# SUPPLEMENTAL MATERIAL of the article:
#   "Health impact projections under climate change scenarios: a tutorial."
#   Ana M. Vicedo-Cabrera, Francesco Sera, Antonio Gasparrini.
#
# This code reproduces the analysis described in the tutorial.
#
# [06/11/2024]
# * an updated version of this code compatible with future
#   versions of the software is available at the personal website of the
#   last author (https://github.com/gasparrini/....)
################################################################################

# This script is adapted by C.Y. Shapland for Regions in Pacific West 
# Using Canadian data provided by Eric Lavigne 1981-2015
# To retrospectively estimate daily temperature-related deaths 
# For years 1981-2015
# With ER curves from 2005-the previous year
# i.e. 2015 mortality estimated from 2005-2014 ER curve
#      2016 mortality estimated from 2005-2015 ER curve
#      ...
#      2020 mortality estimated from 2005-2018 ER curve (ONS mortality data currently until 2018)
# Without climate model projections
# Adapted from Ana M. Vicedo-Cabrera's:
#  - 01EstimationERassociation.r & 
#  - 02ProjTempMortalitySeries.r & 
#  - 04_05_06ExtCurveProjUncert.r 
# Updated 25/01/2021
# Updated 04/03/2021

# LOAD THE PACKAGES
library(MASS) ; library(dplyr)
library(plyr) ; library(PCICt); library(tidyverse)
library(clock) ; library(lubridate)

# set directories
Dir<-"M:/projects/ieu3/p2/001/working/"
setwd(Dir)
output_dir<-paste0(Dir, "results/DLNM/")

# data dir
climate_dir<-paste0(Dir,"data/Health/Canada/Mortality/cleaned/nleach/")

# LOAD functions
source(paste0(Dir,"scripts/Functions.R"))

################################################################################
# 01 ESTIMATION OF THE EXPOSURE-RESPONSE ASSOCIATIONS
################################################################################

# Use the observed daily temperature-mortality series to estimate
#     the coefficients defining the exposure-response association.

# LOAD OBSERVED DATA - DAILY MORTALITY & TEMPERATURE-SERIES BETWEEN 1981-2015 in regions in Canada
regions <- c("abbotsford", "victoria", "vancouver", "lytton")

# mortality and temperature data 1981-2015, used for establishing ER curves
load(paste0(Dir,"data/Health/Canada/Mortality/cleaned/DLNM/Canadian_HWcity_DLNMdata_151223.RData"))

#names of HW attribution data
forecast<-c("Pi","ENS","Incr")
forecast_cond<-c("b2mg","1","b2mh")

# lead time
leadtime_range<-c("2021-06-26", "2021-06-22")

for (j in 1:length(forecast)){
  
  #j<-1
  print(forecast[j])
  
  #for (k in 1:length(regions)){
    k<-3
    ## city
    reg1<-gsub(" ", "", regions[k], fixed=TRUE)
    print(reg1)
    
    ## load health data for specific city
    obs<-health_dat[[reg1]]

    ### 00 Cleaning and reformat ###
    obs_sub <- obs %>% drop_na(death) %>% # removing days with missing mortality data
                subset(month %in% c(6:9)) # subset to summer
    
    ### 01 ESTIMATION OF THE EXPOSURE-RESPONSE ASSOCIATIONS ###
    argvar<-list(fun="ns", knots=quantile(obs_sub$tavg,c(50,90)/100, na.rm=T), Bound=range(obs_sub$tavg,na.rm=T))
  
    ERass_res<-EstERass_HW_01(data=obs_sub, exposure= obs_sub$tavg, outcome=obs_sub$death, 
                              argVAR=argvar, maxlag_sp=3, obs_sub$year)

    for (l in 1:length(leadtime_range)){
      #j=2
      #l=1
      ### 02 CLIMATOLOGICALLY EXPECTED MORTALITY SERIES ###
      obs_ext<-reformat_ncf_to_DLNM(ncf_dat=paste0(climate_dir,forecast[j], "/EricCoord_", reg1,"_", 
                                                    forecast_cond[j],"_", leadtime_range[l],".nc", sep=""),
                                      var = 't2m', is_ERA5=F)
      
      ### 03 SCALING FACTOR ###
      #results extracted using scalingfactor_v2.py from JASMIN
      if (forecast[j]=="Pi"){
        obs_ext_adj<-adjust_scalingfactor(obs_ext,
                                        paste0(climate_dir,"ScalingFactor/scalingfactor_ENS_Pi_",
                                               reg1, "_", leadtime_range[l],".nc.txt"), "subtract")
      } else if (forecast[j]=="Incr"){
        obs_ext_adj<-adjust_scalingfactor(obs_ext,
                                          paste0(climate_dir,"ScalingFactor/scalingfactor_ENS_Incr_",
                                                 reg1, "_", leadtime_range[l],".nc.txt"), "add")
      } else {
        obs_ext_adj<-obs_ext
      }
      print(obs_ext_adj$V1)
      
      # years of interest for estimating daily temp-related deaths
      iyrs <- 2021
      
      # loop years
      # ONLY ONE YEAR
      y<-1
      
      # Year for which daily attributable deaths will be computed
      if (forecast[j]=="ENS"){
        yr <- obs_ext_adj 
        } else {
          yr <- obs_ext_adj %>% dplyr:::select(-V1)
        }

      # EXPECTED MORTALITY SERIES:
      # It is computed as the average mortality for each day of the year 
      #   from daily observed deaths, then repeated along the same projection period
      #   of the modelled temperature series.
      if (reg1=="vancouver"){
        deathexp<-ExpMort_02(data=obs, rm_bound_change=T, data_to_compute=yr)
      } else {
        deathexp<-ExpMort_02(data=obs, rm_bound_change=F, data_to_compute=yr)
      }
      print(leadtime_range[l])
      print(paste0("ED=", sum(deathexp[which(yr$date=="2021-06-26"):which(yr$date=="2021-06-30")])))
      
      ### 04 EXTRAPOLATION OF THE EXPOSURE-RESPONSE CURVE ###
      ### 05 PROJECTION & QUANTIFICATION OF THE IMPACT ###
      ### 06 ENSEMBLE ESTIMATES & QUANTIFICATION OF THE UNCERTAINTY ###

      # The three last steps of the analysis (extrapolation of the curve, 
      #   impact projections and quantification of the uncertainty) can be performed 
      #   sequentially using the following code.

      # In brief, once we extrapolate the curve, we estimate the daily number of attributable 
      #   deaths (AN). 
      # By dividing between the total mortality, we estimate the corresponding 
      #   attributable fractions (AFs).
      # Uncertainty of the estimated impacts is expressed in terms of empirical 
      #   confidence intervals, defined as the 2.5th and 97.5th percentiles of the 
      #   empirical distribution of the impacts across coefficients samples. 
      #   The distribution is obtained through Monte Carlo simulation of the coefficients.
      nsim<-1000
    
      ansim_allens = list()
      for (e in 1:51) { 
        dat_ens <- yr %>% dplyr:::select(date,paste0("temp.",e)) %>% dplyr:::rename(temp=paste0("temp.",e)) 
        ansim_allens[[e]]<-ExtrapoPred_040506(Nsim=nsim, clim_data=dat_ens, argVAR=argvar, CEN=ERass_res$cen, RED=ERass_res, expectDeath=deathexp, seed=13041975)
      }
    
      ### save output ###
      fname <- paste0(output_dir, forecast[j],"/daily_attributable_deathsVancouver2008to2015_summer3dayslag_MMT50to100_TAVG_scalingfactor_nsim",
                    nsim,"_",leadtime_range[l],"_", reg1,".RData")
      save(ansim_allens, file=fname)
    }
  
  rm(argvar, ERass_res, ansim_allens)
  #}
}

