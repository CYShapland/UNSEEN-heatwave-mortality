import xarray as xr
import pandas as pd
import fiona
import iris
import os
import numpy as np

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm

import iris.plot as iplt
import iris.quickplot as qplt

from sklearn.linear_model import LinearRegression
from matplotlib.transforms import ScaledTranslation

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
def Comb_cf_pf(directory, condition, time):
	
	# load control data
	filename_cf= os.path.join(directory, 'cf', f'{condition}_{time}')
	ds_cf = xr.open_mfdataset(format(filename_cf),preprocess=preproc_ds)['t2m']
    
	# load control data
	filename_pf= os.path.join(directory, 'pf', f'{condition}_{time}')
	ds_pf = xr.open_mfdataset(format(filename_pf),preprocess=preproc_ds)['t2m']
    
	# merged data
	ds_con = xr.concat([ds_cf, ds_pf], dim="number")
    
	# convert to celsius
	ds_con = kelvinToCelsius(ds_con)
	ds_con.attrs["units"] = "celsius"
    
	#making ensemble mean
	ds_con = ds_con.resample(time='D').mean()
	
	return ds_con


def AnomAveTempProc(directory1, directory2, condition1, condition2, time):
 
	### process first data ###
	ds_con1=Comb_cf_pf(directory1, condition1, time)	
    
	### process second data ###
	ds_con2=Comb_cf_pf(directory2, condition2, time)
    
	#calculate anomalies
	ds=ds_con2 - ds_con1
    
	# making ensemble mean
	ds=ds.reduce(np.mean,dim="number")
    
	return ds

def scalingfactor(x):
	return (1.6-x)*0.77272048

def letter_annotation(ax, xoffset, yoffset, letter):
    ax.text(xoffset, yoffset, letter, transform=ax.transAxes, size=12, weight='bold')


### Directories ###

# global means
GloM_ENSDir='/postproc/glomean/land/ENS/sfc/'
GloM_PiDir='/postproc/glomean/land/pIC/pi/sfc/'
GloM_IncrDir='/postproc/glomean/land/pIC/incr/sfc/'

### load land sea mask ###
ds_mask = xr.open_mfdataset('/ERA5/land_sea_mask_PNW.nc')
ds_mask = ds_mask.lsm.isel(time=-1)
mask = (ds_mask > 0.5)

###########################################################
# Current - pre-industrial with 250 ensembles
###########################################################

### loop over for ENS data ###
inidate=('2021-06-26.nc', '2021-06-22.nc')
inidate_lab=('3', '7')
ENSPi_lab=('A', 'B')
#inidate=('2021-06-26.nc', '2021-06-22.nc', '2021-06-18.nc')
ENSDir='/ENS/pnw025/sfc/'
PiDir='/pIC/pi/pnw025/sfc/'
IncrDir='/pIC/incr/pnw025/sfc/'
OutputDir='YourDir'

# subset heatwave period
day_range=list(range(26,31))

fig, axes = plt.subplot_mosaic([['A', 'B','C'], ['D', 'E', 'F']], figsize=(15, 12))
fig = plt.figure(figsize=(20, 10))
fig.tight_layout()
fig.subplots_adjust(right=0.8)
plt.rcParams.update({'font.size': 20})

t=1

for x, f, z in zip(inidate, inidate_lab, ENSPi_lab):
	
	#current scenario
	ds_ENS=Comb_cf_pf(ENSDir, '1', x)
	LandOnly_ENS=ds_ENS.where(mask)
	ds_ENS=LandOnly_ENS.sel(time=LandOnly_ENS.time.dt.day.isin(day_range))
	
	#pre-industrial scenario
	ds_Pi=Comb_cf_pf(PiDir, 'b2mg',x)
	ds_Pi=ds_Pi.sel(time=ds_Pi.time.dt.day.isin(day_range))
	
	# Scaling factor 
	ds_gloM_Pi=AnomAveTempProc(GloM_PiDir, GloM_ENSDir, 'b2mg','1', x) 
	
	#ds_gloM_Incr=AnomAveTempProc(GloM_IncrDir, GloM_ENSDir, 'b2mh','1', x)
	#ScaleFactor_Incr=scalingfactor(ds_gloM_Incr.values)
	
	results_Pi=[]
	#results_Incr=[]
	for SF, group in ds_Pi.groupby('time'):
		results_Pi.append(group-scalingfactor(ds_gloM_Pi.sel(time=SF).values))	
		print(scalingfactor(ds_gloM_Pi.sel(time=SF).values))
		#results_Incr.append(group-scalingfactor(ds_gloM_Incr.sel(time=SF).values))
	
	results_Pi=xr.combine_by_coords(results_Pi)
	LandOnly_Pi=results_Pi.where(mask)
	
	#results_Incr=xr.combine_by_coords(results_Incr)
	
	### Temperature anomalies ###
	
	# subset heatwave period
	ds_ENS = ds_ENS.sel(latitude=slice(70,30), longitude=slice(-150,-100))
	ds_ENS = ds_ENS.reduce(np.mean,dim="time")
	ds_ENS = ds_ENS.reduce(np.mean,dim="number")
	
	ds_Pi = LandOnly_Pi.sel(latitude=slice(70,30), longitude=slice(-150,-100))
	ds_Pi = ds_Pi.reduce(np.mean,dim="time")
	ds_Pi = ds_Pi.reduce(np.mean,dim="number")
	
	ds_anom = ds_ENS - ds_Pi
	av_anom=ds_anom['t2m'].mean()
	print(av_anom.values)
	
	### plot ###
	ds_anom.to_netcdf('temp.nc')
	in_cube = iris.load_cube('temp.nc')
	levels =  np.linspace(-2,2,9)
	plt.subplot(2,3,t)
	contour_result = iplt.contourf(in_cube, levels=levels,  cmap="RdBu_r", extend='both')
	plt.gca().coastlines()
	plt.title(f'{f}-day')
	plt.title(z, loc='left', weight='bold', x=-0.1, y=1.1)
	t=t+1

### 250 ensembles ###
 
t=3

x = '2021-06-18.nc'
f = '11' 

#current scenario
ds_ENS=Comb_cf_pf(ENSDir, 'b2mr', x)
LandOnly_ENS=ds_ENS.where(mask)
ds_ENS=LandOnly_ENS.sel(time=LandOnly_ENS.time.dt.day.isin(day_range))

#pre-industrial scenario
ds_Pi=Comb_cf_pf(PiDir, 'b2mp', x)
ds_Pi=ds_Pi.sel(time=ds_Pi.time.dt.day.isin(day_range))

# Scaling factor
ds_gloM_Pi=AnomAveTempProc(GloM_PiDir, GloM_ENSDir, 'b2mg','1', x)
ds_gloM_Pi=ds_gloM_Pi.sel(time=ds_gloM_Pi.time.dt.day.isin(day_range))

results_Pi=[]
for SF, group in ds_Pi.groupby('time'):
	results_Pi.append(group-scalingfactor(ds_gloM_Pi.sel(time=SF).values))

results_Pi=xr.combine_by_coords(results_Pi)
LandOnly_Pi=results_Pi.where(mask)

### Temperature anomalies ###
	
# subset heatwave period
ds_ENS = ds_ENS.sel(latitude=slice(70,30), longitude=slice(-150,-100))
ds_ENS = ds_ENS.reduce(np.mean,dim="time")
ds_ENS = ds_ENS.reduce(np.mean,dim="number")
	
ds_Pi = LandOnly_Pi.sel(latitude=slice(70,30), longitude=slice(-150,-100))
ds_Pi = ds_Pi.reduce(np.mean,dim="time")
ds_Pi = ds_Pi.reduce(np.mean,dim="number")
	
ds_anom = ds_ENS - ds_Pi
av_anom=ds_anom['t2m'].mean()
print(av_anom.values)

ds_anom.to_netcdf('temp.nc')
in_cube250Pi = iris.load_cube('temp.nc')
levels =  np.linspace(-2,2,9)
plt.subplot(2,3,t)
contour_result = iplt.contourf(in_cube250Pi, levels=levels,  cmap="RdBu_r", extend='both')
plt.gca().coastlines()
plt.title(f'{f}-day')
plt.title('C', loc='left', weight='bold', x=-0.1,y=1.1)
#letter_annotation(axes[0,0], -.2, 1.1, 'A')

# Future - current with 250 ensembles
ENSIncr_lab=('D', 'E')

t=4

for x, f, z in zip(inidate, inidate_lab, ENSIncr_lab):
	
	print(x,f)
	
	#current scenario
	ds_ENS=Comb_cf_pf(ENSDir, '1', x)
	LandOnly_ENS=ds_ENS.where(mask)
	ds_ENS=LandOnly_ENS.sel(time=LandOnly_ENS.time.dt.day.isin(day_range))
	
	#pre-industrial scenario
	ds_Incr=Comb_cf_pf(IncrDir, 'b2mh',x)
	ds_Incr=ds_Incr.sel(time=ds_Incr.time.dt.day.isin(day_range))
	
	# Scaling factor 
	ds_gloM_Incr=AnomAveTempProc(GloM_IncrDir, GloM_ENSDir, 'b2mh','1', x)
	
	results_Incr=[]
	
	for SF, group in ds_Incr.groupby('time'):	
		print(scalingfactor(ds_gloM_Incr.sel(time=SF).values))
		results_Incr.append(group+scalingfactor(ds_gloM_Incr.sel(time=SF).values))
	
	results_Incr=xr.combine_by_coords(results_Incr)
	LandOnly_Incr=results_Incr.where(mask)
	
	### Temperature anomalies ###
	
	# subset heatwave period
	ds_ENS = ds_ENS.sel(latitude=slice(70,30), longitude=slice(-150,-100))
	ds_ENS = ds_ENS.reduce(np.mean,dim="time")
	ds_ENS = ds_ENS.reduce(np.mean,dim="number")
	
	ds_Incr = LandOnly_Incr.sel(latitude=slice(70,30), longitude=slice(-150,-100))
	ds_Incr = ds_Incr.reduce(np.mean,dim="time")
	ds_Incr = ds_Incr.reduce(np.mean,dim="number")
	
	ds_anom = ds_Incr - ds_ENS
	av_anom=ds_anom['t2m'].mean()
	print(av_anom.values)
	
	### plot ###
	ds_anom.to_netcdf('temp.nc')
	in_cube = iris.load_cube('temp.nc')
	levels =  np.linspace(-2,2,9)
	plt.subplot(2,3,t)
	contour_result = iplt.contourf(in_cube, levels=levels,  cmap="RdBu_r", extend='both')
	plt.gca().coastlines()
	plt.title(f'{f}-day')
	plt.title(z, loc='left', weight='bold', x=-0.1, y=1.1)
	t=t+1

### 250 ensembles ###
 
t=6

x = '2021-06-18.nc'
f = '11' 

#current scenario
ds_ENS=Comb_cf_pf(ENSDir, 'b2mr', x)
LandOnly_ENS=ds_ENS.where(mask)
ds_ENS=LandOnly_ENS.sel(time=LandOnly_ENS.time.dt.day.isin(day_range))

#pre-industrial scenario
ds_Incr=Comb_cf_pf(IncrDir, 'b2mq', x)
ds_Incr=ds_Incr.sel(time=ds_Incr.time.dt.day.isin(day_range))

# Scaling factor
ds_gloM_Incr=AnomAveTempProc(GloM_IncrDir, GloM_ENSDir, 'b2mh','1', x)
ds_gloM_Incr=ds_gloM_Incr.sel(time=ds_gloM_Incr.time.dt.day.isin(day_range))

results_Incr=[]
for SF, group in ds_Incr.groupby('time'):
	results_Incr.append(group+scalingfactor(ds_gloM_Incr.sel(time=SF).values))

results_Incr=xr.combine_by_coords(results_Incr)
LandOnly_Incr=results_Incr.where(mask)

### Temperature anomalies ###
	
# subset heatwave period
ds_ENS = ds_ENS.sel(latitude=slice(70,30), longitude=slice(-150,-100))
ds_ENS = ds_ENS.reduce(np.mean,dim="time")
ds_ENS = ds_ENS.reduce(np.mean,dim="number")
	
ds_Incr = LandOnly_Incr.sel(latitude=slice(70,30), longitude=slice(-150,-100))
ds_Incr = ds_Incr.reduce(np.mean,dim="time")
ds_Incr = ds_Incr.reduce(np.mean,dim="number")
	
ds_anom = ds_Incr - ds_ENS
av_anom=ds_anom['t2m'].mean()
print(av_anom.values)

ds_anom.to_netcdf('temp.nc')
in_cube250Incr = iris.load_cube('temp.nc')
levels =  np.linspace(-2,2,9)
plt.subplot(2,3,t)
contour_result = iplt.contourf(in_cube250Incr, levels=levels,  cmap="RdBu_r", extend='both')
plt.gca().coastlines()
plt.title(f'{f}-day')
plt.title('F', loc='left', weight='bold', x=-0.1,y=1.1)

cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.5])
fig.colorbar(contour_result,  cax=cbar_ax, ticks=range(-2,2,1))
#cbar = plt.colorbar(contour_result, ticks=range(-2,2,1))
filename_output= os.path.join(OutputDir, f'Quickplot_AvgTemp_PNW_29Jun21_250Ensembles_scalingfactor_ALL.png')
plt.savefig(format(filename_output))


