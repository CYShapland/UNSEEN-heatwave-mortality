################################################################################
# Functions for DLNM 
# Author: Chin Yang Shapland (Adapted from Vicedo-Cabrera et al.)
# Last updated: 05/11/2024
################################################################################

require(ncdf4)
require(dlnm)
require(splines)

#########################################
###            DLNM functions         ###
#########################################

EstERass_HW_01<- function(data, exposure, outcome, argVAR, maxlag_sp, grp_sp){
  
  #exposure = data$tavg
  #knots_sp = quantile(exposure,c(50,90)/100, na.rm=T)
  #grp_sp = data$year
  #maxlag_sp=3
  
  # DEFINITION OF THE CROSS-BASIS FOR TEMPERATURE
  # - SPECIFICATION PARAMETERS OF THE EXPOSURE-RESPONSE DIMENSION OF THE CROSS-BASIS

  # argvar: main model, cubic natural spline with two internal knots in 
  #   the 50th, 90th percentiles of the summer temperature distribution
  #argvar <- argVAR
  
  # - SPECIFICATION PARAMETERS OF THE LAG-ASSOCIATION DIMENSION OF THE CROSS-BASIS
  # Definition of the maximum lag, that is, 3 days
  maxlag <- maxlag_sp
  # arglag: main model, unconstrained
  arglag <- list(fun="integer")
  
  # - CREATE CROSSBASIS OBJECTS
  cb <- crossbasis(exposure,lag=maxlag,argVAR,arglag=arglag,group=grp_sp)  
  
  # FIT THE MODEL
  # Include in the model the crossbasis term of temperature, along with the 
  #    indicator for day of the week (dow) and natural cubic spline of day of year (doy) 
  #    with 4 degrees of freedom, and an interaction of this spline with year.
  m <- glm(outcome ~ cb + dow + ns(date, df=round(4*length(date)/122)), data=data, family=quasipoisson)  
  
  # - DEFINE PROVVISIONAL CENTERING POINT TO HAVE THE INITIAL PREDICTION
  varcen <- 19 
  
  # constrain centering point to certain percentiles
  cenlim <- round(quantile(exposure, c(0.5,0.98), na.rm=T),1)
  
  # - ESTIMATE MMT FROM THE PREDICTED EXPOSURE-RESPONSE ASSOCIATION 
  # MMT corresponds to the temperature of minimum mortality, which will be used as
  #    as reference to estimate relative risks and as temperature threshold 
  #    to differentiate the contribution of heat and cold to the total mortality 
  #    attributable to non-optimal temperatures.
  cp <- crosspred(cb,m,cen=varcen,by=0.1,from=cenlim[1],to=cenlim[2])
  cen <- cp$predvar[which.min(cp$allRRfit)] #mmt is the lowest RR.
  
  # EXTRACT COEFFICIENTS AND VCOV FROM THE MODEL IN STAGE 1
  # With the crossreduce function we can reduce the fit of the bidimensional DLNM
  #   (of the stage 1 using the observed temperature-mortality series)
  #   to summaries defined in the exposure-response dimension. We can then 
  #   extract the coefficients and the covariance matrix defining the overall 
  #   exposure-response association.
  red <- crossreduce(cb, m, cen=cen)
  
  return(red=red)
  
}

ExpMort_02<- function(data, rm_bound_change, data_to_compute){
  
  # use year-round data because some April days are included in certain summers
  obs_nl <- data[substr(data$date,6,10)!="02-29",]
  
  if (rm_bound_change==T){
    obs_nl <- subset(obs_nl, year>=2008)
  }
  
  # first, year-round average deaths per day
  deathdoy <- tapply(obs_nl$death,as.numeric(format(obs_nl$date,"%j")),mean,na.rm=T)[seq(365)]
  while(any(isna <- is.na(deathdoy)))
    deathdoy[isna] <- rowMeans(Lag(deathdoy,c(-1,1)),na.rm=T)[isna]
  # then, either repeat them or subset them
  if (nrow(data_to_compute)>=365) {
    # has to start from Jan 1st
    expect_deaths <- rep(deathdoy,length=nrow(data_to_compute))
  } else {
    # start from any date
    dinx1 <- yday(data_to_compute$date[1])
    dinx2 <- yday(data_to_compute$date[nrow(data_to_compute)])
    expect_deaths <- as.numeric(deathdoy[dinx1:dinx2])
  }
  
  return(expect_deaths)
  
}

ExtrapoPred_040506<-function(Nsim, clim_data, argVAR, CEN, RED, expectDeath, seed) {
  
  coef <- coef(RED)
  vcov <- vcov(RED)
  
  # Loop through ensemble members of tmax, if any
  # Monitor mode only has one real-world realisation of Ts
  
  # DEFINE THE DATAFRAME TO STORE THE ESTIMATED AN
  ANSIM <- data.frame(matrix(NA, nrow=length(clim_data$date), ncol=Nsim+4,
                             dimnames=list(NULL,c("date","temp","mmt","est",paste0("sim",seq(Nsim))))))
  
  # STORE DATES, TMean AND MMT 
  ANSIM[,"date"] <- clim_data$date
  ANSIM[,"temp"] <- clim_data$temp
  ANSIM[,"mmt"] <- CEN
  
  # (4) EXTRAPOLATION OF THE CURVE: 
  # - DERIVE THE CENTERED BASIS USING THE PROJECTED TEMPERATURE SERIES
  #   AND EXTRACT PARAMETERS
  bvar <- do.call(onebasis,c(list(x=clim_data$temp),argVAR))
  cenvec <- do.call(onebasis,c(list(x=CEN),argVAR))
  bvarcen <- scale(bvar,center=cenvec,scale=F)
  
  # (5) IMPACT PROJECTIONS:
  # - COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
  an <- (1-exp(-bvarcen%*%coef(RED)))*expectDeath
  
  # - STORE AN IN ARRAY BEFORE THE ITERATIONS
  ANSIM[,"est"] <- an
  
  # (6) ESTIMATE UNCERTAINTY OF THE PROJECTED AN:
  # - SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
  set.seed(seed)
  coefsim <- mvrnorm(Nsim,coef,vcov)
  
  # - LOOP ACROSS ITERATIONS
  for(s in seq(Nsim)) {
    
    # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
    an <- (1-exp(-bvarcen%*%coefsim[s,]))*expectDeath
    
    # STORE THE ATTRIBUTABLE MORTALITY
    ANSIM[,s+4] <- an
  }
  ANSIM["expdeaths"]<-expectDeath

  return(ANSIM)
  
}


#########################################
###     Reformat and cleaning         ###
#########################################

# reformat nleach's data to fit with DLNM
reformat_ncf_to_DLNM <- function(ncf_dat, var, is_ERA5){ 
  
  #ncf_dat <- paste0(climate_dir,forecast[j], "/250Ensemble_EricCoord_", reg1,"_", forecast_cond[j],"_2021-06-18.nc", sep="")
  #var <- 't2m'
  
  # LOAD FORECAST DATA
  ncin <- nc_open(ncf_dat)
  #print(ncin)
  
  ### get time ###
  time <- ncvar_get(ncin,"time")
  tunits <- ncatt_get(ncin,"time","units")
  nt <- dim(time)
  # use tunits to find the start time
  datetime <- as.POSIXct(gsub('days since ', '', tunits$value), "Canada/Pacific")+days(time)
  
  ### get variables ###
  var_array <- ncvar_get(ncin,var)
  nvar_array <- dim(var_array)
  
  nc_close(ncin)
  
  # create dataframe and combine with time variables
  if (is_ERA5==T){
    nleach_dat <- data.frame(temp=var_array, date=as.Date(datetime, format="%Y-%m-%d"))
  } else{
    nleach_dat <- data.frame(temp=t(var_array), date=as.Date(datetime, format="%Y-%m-%d"))
  }
  
  return(nleach_dat)
  
}  


######################
### Scaling factor ###
######################

adjust_scalingfactor <- function(original_data, scalingfactor_data, direction){
  scalingfactor_city <- read.table(scalingfactor_data, quote="\"", comment.char="")
  GlobalMean_time <- read.table("ScalingFactor/GlobalMean_time.txt", quote="\"", comment.char="", stringsAsFactors=F)
  GlobalMean_time["date"]<-gsub("\\T.*","",GlobalMean_time$V1)
  scalingfactor_city["date"]<-tail(GlobalMean_time$date,n=nrow(scalingfactor_city))
  
  if (direction=="add"){
    obs_adjust<-merge(original_data, scalingfactor_city, by = "date") %>%
      dplyr:::mutate(across(starts_with("temp."), ~ .x + V1))
    } else if(direction=="subtract"){
    obs_adjust<-merge(original_data, scalingfactor_city, by = "date") %>%
      dplyr:::mutate(across(starts_with("temp."), ~ .x - V1))
    } else {
      obs_adjust<-merge(original_data, scalingfactor_city, by = "date")
  }
  
  return(obs_adjust)
}










