import xarray as xr
import pandas as pd
import fiona
import iris
import os
import numpy as np

import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
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

### Directories ###
ENSDir='/ENS/pnw025/sfc/'
PiDir='/pIC/pi/pnw025/sfc/'
IncrDir='/pIC/incr/pnw025/sfc/'

# global means
GloM_ENSDir='/postproc/glomean/land/ENS/sfc/'
GloM_PiDir='/postproc/glomean/land/pIC/pi/sfc/'
GloM_IncrDir='/postproc/glomean/land/pIC/incr/sfc/'

inidate=('2021-06-26.nc', '2021-06-22.nc', '2021-06-18.nc')

#x='2021-06-26.nc'

day1_gloM=np.empty(3, dtype=float)
day2_gloM=np.empty(3, dtype=float)

### load current global mean data ###

for d, b in zip(inidate, range(3)):
	ds_gloM=AnomAveTempProc(GloM_PiDir, GloM_IncrDir, 'b2mg','b2mh',d) 
	day1_gloM[b]= float(ds_gloM.loc["2021-06-26"].values)
	day2_gloM[b]= float(ds_gloM.loc["2021-06-27"].values)

#remove last day as it is only one time point
#ds_gloM=ds_gloM[:-1]

### Local temperatures ###

# coordinates for abbotsford, victoria, vancouver
lats =  [49.06, 48.43, 49.24, 50.23]
lons =  [-122.25, -123.37, -123.11, -121.58]
city_name=['abbotsford', 'victoria', 'vancouver', 'lytton']


for f, b, d in zip(lons, lats, city_name):
	day1_sub_loc=np.empty(3, dtype=float)
	day2_sub_loc=np.empty(3, dtype=float)

	for x, c in zip(inidate, range(3)):
		ds_locM=AnomAveTempProc(PiDir, IncrDir, 'b2mg','b2mh',x)
		ds_sub_loc = ds_locM.sel(longitude=f, latitude=b, method='nearest')
		day1_sub_loc[c]= float(ds_sub_loc.loc["2021-06-26"].values)
		day2_sub_loc[c]= float(ds_sub_loc.loc["2021-06-27"].values)
	
	# separate cities and save output
	model_day1 = LinearRegression(fit_intercept=False).fit(day1_gloM.reshape((-1,1)), day1_sub_loc)
	model_day2 = LinearRegression(fit_intercept=False).fit(day2_gloM.reshape((-1,1)), day2_sub_loc) 
	print(d)
	print(sum(model_day1.coef_, model_day2.coef_)/2)
		

#create the amount that needs to be increase
#results is from above loop
slope_city=[0.65601952, 1.0339223, 0.65749857, 0.47766658]

for f, d in zip(slope_city, city_name):
	for x in inidate:
		ds_gloM=AnomAveTempProc(GloM_PiDir, GloM_ENSDir, 'b2mg','1', x) 
		ScaleFactor=(1.6-ds_gloM.values)*f
		np.savetxt(f'scalingfactor_ENS_Pi_{d}_{x}.txt', ScaleFactor)

for f, d in zip(slope_city, city_name):
	for x in inidate:
		ds_gloM=AnomAveTempProc(GloM_IncrDir, GloM_ENSDir, 'b2mh','1', x)
		ScaleFactor=(1.6-ds_gloM.values)*f
		np.savetxt(f'scalingfactor_ENS_Incr_{d}_{x}.txt', ScaleFactor)

x='2021-06-26.nc'
ds_gloM=AnomAveTempProc(GloM_ENSDir, GloM_PiDir, '1','b2mg', x)
ds_gloM=AnomAveTempProc(GloM_ENSDir, GloM_IncrDir, '1','b2mh', x)

np.savetxt(f'GlobalMean_time.txt', ds_gloM.time.values,fmt="%s")
