################################################################################
# SUPPLEMENTAL MATERIAL of the article:
#   "Health impact projections under climate change scenarios: a tutorial."
#   Ana M. Vicedo-Cabrera, Francesco Sera, Antonio Gasparrini.
#
# This code reproduces the analysis described in the tutorial.
#
# [date]
# * an updated version of this code compatible with future
#   versions of the software is available at the personal website of the
#   last author (https://github.com/gasparrini/....)
################################################################################

# This script is adapted by Chin Yang Shapland for Regions in Pacific Northwest of Canada
# Using Canadian data provided by Eric Lavigne 1981-2015
# Without additional linear fitting tests
# Adapted from Ana M. Vicedo-Cabrera's 01EstimationERassociation.r

# LOAD THE PACKAGES
library(MASS); library(tidyverse)

# set directories
datpath <- "YourDir"

# LOAD functions
source(paste0(datpath,"scripts/Functions.R"))

################################################################################
# 01 ESTIMATION OF THE EXPOSURE-RESPONSE ASSOCIATIONS
################################################################################

# Use the observed daily temperature-mortality series to estimate
#     the coefficients defining the exposure-response association.
# In this example, we estimate the temperature-related mortality association 
#     using data from Canada between 1981 and 2015.

load(paste0(datpath, "/DLNM/Canadian_HWcity_DLNMdata_260624.RData"))
lapply(health_dat, dim)

# outcome of interest
outcome_var<-names(health_dat[[1]] %>% dplyr::select(contains("all_mortality")))

# LOAD OBSERVED DATA - DAILY TEMPERATURE-SERIES BETWEEN 1981-2015 in Regions in Pacific Northwest of Canada
cities <- c("abbotsford", "victoria", "vancouver", "lytton")
city_lab <- c("Fraser East", "South Vancouver Island", "Vancouver", "Thompson/Cariboo")

ncities <- length(cities)

for(i in 1:ncities) {
  
  #i<-3
  
  ## city
  city <- cities[[i]]
  print(city)
  
  for(j in 1:length(outcome_var)) {
    #j<-1
    
    # reformatting for DLNM
    obs <- health_dat[[city]] %>% rename(outcome=outcome_var[j]) %>% #isolate outcome of interest
      subset(month %in% c(6:9)) %>%      #isolate summer months
      drop_na(outcome)                   #removing days with missing mortality data
    
    # DEFINITION OF THE CROSS-BASIS FOR TEMPERATURE
    # - SPECIFICATION PARAMETERS OF THE EXPOSURE-RESPONSE DIMENSION OF THE CROSS-BASIS
    
    # argvar: main model, cubic natural spline with two internal knots in 
    #   the 50th, 90th percentiles of the temperature distribution
    argvar <- list(fun="ns", knots = quantile(obs$tavg,c(50,90)/100, na.rm=T), Bound=range(obs$tavg,na.rm=T))
    
    # - SPECIFICATION PARAMETERS OF THE LAG-ASSOCIATION DIMENSION OF THE CROSS-BASIS
    # Definition of the maximum lag, that is, 3 days
    maxlag <- 3
    # arglag: main model, unconstrained 
    #   includes intercept, i.e. minimum lag=0 by default
    arglag <- list(fun="integer")
    
    # - CREATE CROSSBASIS OBJECTS
    cb <- crossbasis(obs$tavg,maxlag,argvar,arglag,group=obs$year)  
    
    # FIT THE MODEL
    # Include in the model the crossbasis term of temperature, along with the 
    #    indicator for day of day of the week (dow) and 
    #    natural cubic spline of day of year with 4 dfs, 
    #    and an interaction of this spline with year 
    m <- glm(outcome ~ cb + dow + ns(doy,df=4)*factor(year), data=obs, family=quasipoisson)
    
    # GET PREDICTIONS & PLOT
    # - DEFINE PROVVISIONAL CENTERING POINT TO HAVE THE INITIAL PREDICTION
    varcen <- 19 
    
    # constrain centering point to certain percentiles
    cenlim <- round(quantile(obs$tavg, c(0.5,0.98), na.rm=T),1) #changed from c(0.5,1) as we don't want the highest temp.
    
    # - ESTIMATE MMT FROM THE PREDICTED EXPOSURE-RESPONSE ASSOCIATION 
    # MMT corresponds to the temperature of minimum mortality, which will be used as
    #    as reference to estimate relative risks and as temperature threshold 
    #    to differentiate the contribution of heat and cold to the total mortality 
    #    attributable to non-optimal temperatures.
    cp <- crosspred(cb,m,cen=varcen,by=0.1,from=cenlim[1],to=cenlim[2])
    cen <- cp$predvar[which.min(cp$allRRfit)] 
    
    # - RE-CENTER & GET PREDICTIONS FOR EACH MODEL CENTERING ON THE MMT 
    pred <- crosspred(cb, m, cen=cen, by=1)   
    
    # print MMT and its percentile in summer T distribution
    mmp <- round(mean(obs$tavg <= cen)*100, digits=3)
    cat(paste0("MMT = ",cen," deg C\n"))
    cat(paste0("MMP = ",mmp," th p\n"))
    
    # PLOT - FIGURE 1
    xlab <- expression(paste("Temperature (",degree,"C)"))
    
    pdf(paste0(datpath, "/DLNM/3D_2D_summer",maxlag, "dayslag_MMT50to100_mdl2_", 
               outcome_var[j], "_", city, ".pdf"),height=4.5,width=8)
    layout(matrix(c(1,1,2,3,2,3,2,3),ncol=2,byrow=T))
    
    par(mar=c(0,0,0,0))
    plot(0:5,axes=F,ylab="",xlab="",type="n")
    text(3.5,3.5,paste0("Mortality in ", city_lab[i]),cex=2, font=2)
    text(3.5,2,"1981-2015",cex=1.4)
    
    # PLOT - 3D
    par(mar=c(2,1,1,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
    plot(pred,"3d",ltheta=150,xlab="Temperature (\u00B0C)",ylab="Lag",zlab="RR", col=gray(0.9), main="Exposure-lag-response")
    fig_label("A", cex=1.4, font=2)
    
    ### 04 EXTRAPOLATION OF THE EXPOSURE-RESPONSE CURVE ###
    red <- crossreduce(cb, m, cen=cen)
    coef <- coef(red)
    vcov <- vcov(red)
    
    # load climate data
    obs_ext<-reformat_ncf_to_DLNM(ncf_dat=paste0(datpath, 
                                                 "/Incr/250Ensemble_EricCoord_", 
                                                 city, "_b2mq_2021-06-18.nc", sep=""),
                                  var = 't2m', is_ERA5=F)
    
    bvar <- do.call(onebasis,c(list(x=as.vector(unlist(obs_ext %>% dplyr::select(contains("temp"))))),argvar))
    pred_exp <- crosspred(bvar,coef=coef,vcov=vcov,model.link="log", cen=cen,by=0.1)
    
    # OVERALL with exploration
    # The plots show the cumulative exposure-response association, in terms of 
    #    relative risks (RR) and centered in the MMT, across the 21 days of lag.
    par(mar=c(5,3,1,2),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
    plot(pred_exp,type="n",ylim=c(0.5,4),lab=c(6,5,7),xlab="",ylab="RR", main=paste0("Overall, MMT=",round(cen,digits=1),"\u00B0C"))
    ind1 <- pred_exp$predvar<=min(obs$tavg)
    ind2 <- pred_exp$predvar>=min(obs$tavg) & pred_exp$predvar<=cen
    ind3 <- pred_exp$predvar>=cen & pred_exp$predvar<=max(obs$tavg)
    ind4 <- pred_exp$predvar>=max(obs$tavg)
    if(any(ind1)) lines(pred_exp$predvar[ind1],pred_exp$allRRfit[ind1],col=4,lwd=1.5,lty=2)
    lines(pred_exp$predvar[ind2],pred_exp$allRRfit[ind2],col=4,lwd=1.5)
    lines(pred_exp$predvar[ind3],pred_exp$allRRfit[ind3],col=2,lwd=1.5)
    lines(pred_exp$predvar[ind4],pred_exp$allRRfit[ind4],col=2,lwd=1.5,lty=2)
    abline(v=cen,lty=3)
    abline(v=max(obs$tavg),lty=2)
    fig_label("B", cex=1.4, font=2)
    
    layout(1)
    dev.off()
    
    save(pred, file=paste0(datpath,"results/DLNM/DLNM_Res_summer",maxlag, "dayslag_MMT50to100_mdl2_", 
                           outcome_var[j], "_", city, ".RData"))
    rm(pred, obs, m)
  }
}
