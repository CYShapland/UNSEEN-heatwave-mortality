import xarray as xr
import pandas as pd
import fiona
import iris
import os
import numpy as np

import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

import cartopy.crs as ccrs
import matplotlib.cm as mpl_cm

import iris.plot as iplt
import iris.quickplot as qplt

### Functions n###

def preproc_ds(ds):
    
    ds = ds.copy().squeeze()
    
    fname = ds.encoding['source'].split('/')[-1].split('.')[0]
    
    #expver = fname.split('_')[0]
    #ds = ds.expand_dims({'experiment':[expver]})
    
    # set up aux data
    inidate = pd.to_datetime(ds.time[0].values)
    
    # expand dimensions to include extra info
    #if not 'hDate' in ds:
    #    ds = ds.expand_dims({'inidate':[inidate]})
        
    if not 'number' in ds:
    	ds = ds.expand_dims({'number':[0]})
       
    # put time dimension at front
    ds = ds.transpose('time',...)
            
    return ds

# load convert kelvin To Celsius
def kelvinToCelsius(kelvin):
    return kelvin - 273.15

# anomaliy mean

def AnomAveTempProc(directory1, directory2, condition1, condition2, time):
 
	# load control data
	filename_cf= os.path.join(directory1, 'cf', f'{condition1}_{time}')
	ds_cf = xr.open_mfdataset(format(filename_cf),preprocess=preproc_ds)['t2m']
	
	# load control data
	filename_pf= os.path.join(directory1, 'pf', f'{condition1}_{time}')
	ds_pf = xr.open_mfdataset(format(filename_pf),preprocess=preproc_ds)['t2m']

	# merged data
	ds_con1 = xr.concat([ds_cf, ds_pf], dim="number")

	# convert to celsius
	ds_con1 = kelvinToCelsius(ds_con1)
	ds_con1.attrs["units"] = "celsius"

	#making ensemble mean
	ds_con1 = ds_con1.resample(time='D').mean()

	### load future global mean data ###
	# load control data
	filename_cf_con2= os.path.join(directory2, 'cf', f'{condition2}_{time}')
	ds_con2_cf = xr.open_mfdataset(format(filename_cf_con2),preprocess=preproc_ds)['t2m']

	# load perturbed data
	filename_pf_con2 = os.path.join(directory2, 'pf', f'{condition2}_{time}')
	ds_con2_pf = xr.open_mfdataset(format(filename_pf_con2),preprocess=preproc_ds)['t2m']

	# merged data
	ds_con2 = xr.concat([ds_con2_cf, ds_con2_pf], dim="number")

	# convert to celsius
	ds_con2 = kelvinToCelsius(ds_con2)
	ds_con2.attrs["units"] = "celsius"

	#making daily temp.
	ds_con2 = ds_con2.resample(time='D').mean()

	#calculate anomalies
	ds=ds_con2 - ds_con1

	# making ensemble mean
	ds=ds.reduce(np.mean,dim="number")

	return ds

#####################################################
# Finding regression coefficient
#####################################################

### Directories ###
ENSDir='/ENS/pnw025/sfc/'
PiDir='/pIC/pi/pnw025/sfc/'
IncrDir='/pIC/incr/pnw025/sfc/'

# global means
GloM_ENSDir='/postproc/glomean/land/ENS/sfc/'
GloM_PiDir='/postproc/glomean/land/pIC/pi/sfc/'
GloM_IncrDir='/postproc/glomean/land/pIC/incr/sfc/'

inidate=('2021-06-26.nc', '2021-06-22.nc', '2021-06-18.nc')

### load land sea mask ###
ds_mask = xr.open_mfdataset('/ERA5/land_sea_mask_PNW.nc')
ds_mask = ds_mask.lsm.isel(time=-1)
mask = (ds_mask > 0.5)

#x='2021-06-26.nc'

day1_gloM=np.empty(3, dtype=float)
day2_gloM=np.empty(3, dtype=float)

### load current global mean data ###

for d, b in zip(inidate, range(3)):
	if d=='2021-06-18.nc':
		ds_gloM=AnomAveTempProc(GloM_PiDir, GloM_IncrDir, 'b2mp','b2mq',d)
	else:
		ds_gloM=AnomAveTempProc(GloM_PiDir, GloM_IncrDir, 'b2mg','b2mh',d) 
	
	day1_gloM[b]= float(ds_gloM.loc["2021-06-26"].values)
	day2_gloM[b]= float(ds_gloM.loc["2021-06-27"].values)

### BC temperatures ###

day1_sub_loc=np.empty(3, dtype=float)
day2_sub_loc=np.empty(3, dtype=float)

for x, c in zip(inidate, range(3)):
	if x=='2021-06-18.nc':
		ds_locTemp=AnomAveTempProc(PiDir, IncrDir, 'b2mp','b2mq',x)  
	else:
		ds_locTemp=AnomAveTempProc(PiDir, IncrDir, 'b2mg','b2mh',x)
	
	ds_LandOnly=ds_locTemp.where(mask)
	ds_locM=ds_LandOnly.mean(('longitude', 'latitude'))
	day1_sub_loc[c]= float(ds_locM.loc["2021-06-26"].values)
	day2_sub_loc[c]= float(ds_locM.loc["2021-06-27"].values)
	
# estimate the coefficient between BC and global mean temperature 
model_day1 = LinearRegression(fit_intercept=False).fit(day1_gloM.reshape((-1,1)), day1_sub_loc)
model_day2 = LinearRegression(fit_intercept=False).fit(day2_gloM.reshape((-1,1)), day2_sub_loc) 
print(sum(model_day1.coef_, model_day2.coef_)/2)


