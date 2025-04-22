

# This code creates hysteresis based on simulations with full ocean module MOM
#%%
#%%
# ###  LOAD PACKAGES
#-------------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1 import make_axes_locatable

savedir = "C:\\Users\\anjaka\\Nextcloud\\PhD\\03_MonsoonPlanet_Hysteresis\\Figures"

#%%###  LOAD DATA
#-------------------

# monthly data of 10 years
# adapt directory to the folder fullocean
data_dir = 'C:/Users/anjaka/Nextcloud/PhD/03_MonsoonPlanet_Hysteresis/data/fullocean/full_ocean_atmos_month_10years_sel.nc'

#%%
# Extract the variables

data = xr.open_dataset(data_dir)
sw = data['swdn_toa']
tsurf = data['t_ref']-273.15
pr = data['precip']*86400
wvp= data['WVP']

#%%
# Extract India 

lat_low_il2 = 17 
lat_up_il2 = 26
lon_low_il2 = 76
lon_up_il2 = 82

sw_lat= sw.sel(lat=slice(lat_low_il2,lat_up_il2))
weights = np.cos(np.deg2rad(sw_lat.lat))
sw_lat_mean = sw_lat.weighted(weights).mean(dim=['lat'])
sw_lat_mean_lon = sw_lat_mean.sel(lon = slice(lon_low_il2,lon_up_il2))
sw_india = sw_lat_mean_lon.mean('lon')

pr_lat= pr.sel(lat=slice(lat_low_il2,lat_up_il2))
weights = np.cos(np.deg2rad(pr_lat.lat))
pr_lat_mean = pr_lat.weighted(weights).mean(dim=['lat'])
pr_lat_mean_lon = pr_lat_mean.sel(lon = slice(lon_low_il2,lon_up_il2))
pr_india = pr_lat_mean_lon.mean('lon')

#%%
# calcuating monthly means of 10 years
sw_india_mean = sw_india.groupby('time.month').mean(dim='time')
pr_india_mean = pr_india.groupby('time.month').mean(dim='time')


#%%
# Plot 

plt.figure()
plt.plot(sw_india, pr_india, color = "gray", alpha = 0.3)
plt.plot(sw_india_mean, pr_india_mean, color = "black",  marker = '.')
plt.plot([sw_india_mean[0], sw_india_mean[-1]], [pr_india_mean[0], pr_india_mean[-1]], color='black')
plt.plot(sw_india_mean[0], pr_india_mean[0], marker = '>', color = 'red', markersize = 10)
plt.plot(sw_india_mean[6], pr_india_mean[6], marker = 'o', color = 'red', markersize = 10)
plt.xlabel('Incoming Solar Radiation ($W/m^2$)')
plt.ylabel('Monsoon Rainfall (mm/day)')
plt.savefig(savedir + '/hysteresis_fullocean.pdf',bbox_inches = 'tight')

# %%
