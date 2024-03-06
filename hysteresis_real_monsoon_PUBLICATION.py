#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 13:59:49 2023

@author: anjakatzenberger
"""
# This code creates Hysteresis plots based on Observations
# Plots are provided for different monsoon regions, 
# including India, Bay of Bengal, East Asia, Australia
# Hysteresis is shown depending on solar radiation as well as land surface temperature

#%%
### LOAD PACKAGES
#-------------------

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# Where to save the figures
savedir = r'C:\Users\anjaka\Nextcloud\PhD\03_MonsoonPlanet_Hysteresis\Figures\hysteresis_real_world\\'

# observed shortwave down (top of atmosphere) data from CERES (monthly 2001-2019)
# https://doi.org/10.1175/JCLI-D-17-0208.1
dir_sw = "C:/Users/anjaka/Nextcloud/PhD/03_MonsoonPlanet_Hysteresis/data/CERES_EBAF-TOA_Ed4.2_Subset_200101-201912_ymonmean2.nc"

# observed precipitation data from GPCC (monthly 2001-2019, ymonmean)
# 10.5676/DWD_GPCC/FD_D_V2020_100
dir_pr = "C:/Users/anjaka/Nextcloud/PhD/03_MonsoonPlanet_Hysteresis/data/GPCC_daily_2001_2019_ymonmean.nc"
dir_pr_20 = "C:/Users/anjaka/Nextcloud/PhD/03_MonsoonPlanet_Hysteresis/data/GPCC_daily_monmean.nc"

# land surface temperatures from TERRA-MODIS (monthly 2001-2018, ymonmean)
# https://dx.doi.org/10.5285/32d7bc64c7b740e9ad7a43589ab91592
dir_tsurf = "C:/Users/anjaka/Nextcloud/PhD/03_MonsoonPlanet_Hysteresis/data/terramodis_ymonmean.nc"
dir_tsurf_18 = "C:/Users/anjaka/Nextcloud/PhD/03_MonsoonPlanet_Hysteresis/data/terramodis.nc"



#%%
### LOAD DATA
#-------------------

data_sw = xr.open_dataset(dir_sw)
swdn = data_sw["solar_mon"]

data_pr = xr.open_dataset(dir_pr)
precip = data_pr["precip"]

data_tsurf = xr.open_dataset(dir_tsurf)
tsurf = data_tsurf["lst"]-273.15


#%%
#  PROCESS DATA
#-------------------

# India (upper and lower thresholds for region)

lat_low = 17 
lat_up = 26
lon_low = 76
lon_up = 82

sw_lat= swdn.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(sw_lat.lat))
sw_lat_mean = sw_lat.weighted(weights).mean(dim=['lat'])
sw_lat_mean_lon = sw_lat_mean.sel(lon = slice(lon_low,lon_up))
sw_il2 = sw_lat_mean_lon.mean('lon')

pr_lat= precip.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(pr_lat.lat))
pr_lat_mean = pr_lat.weighted(weights).mean(dim=['lat'])
pr_lat_mean_lon = pr_lat_mean.sel(lon = slice(lon_low,lon_up))
pr_il2 = pr_lat_mean_lon.mean('lon')


# Bay_of_Bengal

lat_low = 15
lat_up = 30
lon_low = 80
lon_up = 100

sw_lat= swdn.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(sw_lat.lat))
sw_lat_mean = sw_lat.weighted(weights).mean(dim=['lat'])
sw_lat_mean_lon = sw_lat_mean.sel(lon = slice(lon_low,lon_up))
sw_b = sw_lat_mean_lon.mean('lon')

pr_lat= precip.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(pr_lat.lat))
pr_lat_mean = pr_lat.weighted(weights).mean(dim=['lat'])
pr_lat_mean_lon = pr_lat_mean.sel(lon = slice(lon_low,lon_up))
pr_b = pr_lat_mean_lon.mean('lon')


# China

lat_low = 25
lat_up = 30
lon_low = 100
lon_up = 110

sw_lat= swdn.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(sw_lat.lat))
sw_lat_mean = sw_lat.weighted(weights).mean(dim=['lat'])
sw_lat_mean_lon = sw_lat_mean.sel(lon = slice(lon_low,lon_up))
sw_c = sw_lat_mean_lon.mean('lon')

pr_lat= precip.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(pr_lat.lat))
pr_lat_mean = pr_lat.weighted(weights).mean(dim=['lat'])
pr_lat_mean_lon = pr_lat_mean.sel(lon = slice(lon_low,lon_up))
pr_c = pr_lat_mean_lon.mean('lon')


# Australia

lat_low = -18
lat_up = -10
lon_low = 120
lon_up = 150

sw_lat= swdn.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(sw_lat.lat))
sw_lat_mean = sw_lat.weighted(weights).mean(dim=['lat'])
sw_lat_mean_lon = sw_lat_mean.sel(lon = slice(lon_low,lon_up))
sw_a = sw_lat_mean_lon.mean('lon')

pr_lat= precip.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(pr_lat.lat))
pr_lat_mean = pr_lat.weighted(weights).mean(dim=['lat'])
pr_lat_mean_lon = pr_lat_mean.sel(lon = slice(lon_low,lon_up))
pr_a = pr_lat_mean_lon.mean('lon')

#%%
# PROCESSING 18/19 YEARS DATA
#-------------------

### TSURF
data_tsurf_18 = xr.open_dataset(dir_tsurf_18)
tsurf_18 = data_tsurf_18["lst"]-273.15

# adapt format of tsurf
lst_values = tsurf.values
lst_list = lst_values.tolist()
lst_flat_list = [item for sublist1 in lst_list for sublist2 in sublist1 for item in sublist2]
tsurf = lst_flat_list

tsurf_18 = tsurf_18.isel(lat=0, lon=0)



# PRECIPITATION

data_pr = xr.open_dataset(dir_pr_20)
precip = data_pr["precip"]

lat_low = 17 
lat_up = 26
lon_low = 76
lon_up = 82

pr_lat= precip.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(pr_lat.lat))
pr_lat_mean = pr_lat.weighted(weights).mean(dim=['lat'])
pr_lat_mean_lon = pr_lat_mean.sel(lon = slice(lon_low,lon_up))
pr20_il2 = pr_lat_mean_lon.mean('lon')

pr20_il2_18 = pr20_il2.sel(time = slice("2001-01", "2018-12"))


# Bay_of_Bengal

lat_low = 15
lat_up = 30
lon_low = 80
lon_up = 100

sw_lat= swdn.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(sw_lat.lat))
sw_lat_mean = sw_lat.weighted(weights).mean(dim=['lat'])
sw_lat_mean_lon = sw_lat_mean.sel(lon = slice(lon_low,lon_up))
sw_b = sw_lat_mean_lon.mean('lon')

pr_lat= precip.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(pr_lat.lat))
pr_lat_mean = pr_lat.weighted(weights).mean(dim=['lat'])
pr_lat_mean_lon = pr_lat_mean.sel(lon = slice(lon_low,lon_up))
pr20_b = pr_lat_mean_lon.mean('lon')


# China

lat_low = 25
lat_up = 30
lon_low = 100
lon_up = 110

sw_lat= swdn.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(sw_lat.lat))
sw_lat_mean = sw_lat.weighted(weights).mean(dim=['lat'])
sw_lat_mean_lon = sw_lat_mean.sel(lon = slice(lon_low,lon_up))
sw_c = sw_lat_mean_lon.mean('lon')

pr_lat= precip.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(pr_lat.lat))
pr_lat_mean = pr_lat.weighted(weights).mean(dim=['lat'])
pr_lat_mean_lon = pr_lat_mean.sel(lon = slice(lon_low,lon_up))
pr20_c = pr_lat_mean_lon.mean('lon')



# Australia

lat_low = -18
lat_up = -10
lon_low = 120
lon_up = 150

sw_lat= swdn.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(sw_lat.lat))
sw_lat_mean = sw_lat.weighted(weights).mean(dim=['lat'])
sw_lat_mean_lon = sw_lat_mean.sel(lon = slice(lon_low,lon_up))
sw_a = sw_lat_mean_lon.mean('lon')

pr_lat= precip.sel(lat=slice(lat_low,lat_up))
weights = np.cos(np.deg2rad(pr_lat.lat))
pr_lat_mean = pr_lat.weighted(weights).mean(dim=['lat'])
pr_lat_mean_lon = pr_lat_mean.sel(lon = slice(lon_low,lon_up))
pr20_a = pr_lat_mean_lon.mean('lon')


#%%
### PLOT: Hysteresis sw panel 
#-------------------

# Create a 2x2 panel plot
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
l=16

# Plot on the first subplot (top-left)
axes[0, 0].plot(sw_il2, pr_il2, color = "black", marker = ".")
axes[0, 0].plot(np.tile(sw_il2.data,19),pr20_il2, color = "gray", alpha = 0.3)
axes[0, 0].plot(sw_il2[0], pr_il2[0],  marker = '.', color = "red")
axes[0, 0].plot(sw_il2[11], pr_il2[11],  marker = '.', color = "blue")
axes[0, 0].set_title('India')
axes[0, 0].set_ylabel('Monsoon Rainfall (mm/day)')
axes[0, 0].set_ylim([0,l])
axes[0, 0].set_xlim([240,475])
#axes[0, 0].set_xticks([])  # Hide yticks for Bay of Bengal panel
axes[0, 0].plot([sw_il2[0], sw_il2[-1]], [pr_il2[0], pr_il2[-1]], color='black')
axes[0, 0].set_xlabel('Solar Insolation ($W/m^2$)')
axes[0, 0].set_xlim([240,475])

# Plot on the second subplot (top-right)
axes[0, 1].plot(sw_b, pr_b, color='black')
axes[0, 1].plot(np.tile(sw_b.data,19),pr20_b, color = "gray", alpha = 0.3)
axes[0, 1].plot(sw_b[0], pr_b[0],  marker = '.', color = "red")
axes[0, 1].plot(sw_b[11], pr_b[11],  marker = '.', color = "blue")
axes[0, 1].set_title('Bay of Bengal')
axes[0, 1].set_ylim([0,l])
axes[0, 1].set_xlim([240,475])
axes[0, 1].set_yticks([])  # Hide yticks for Bay of Bengal panel
axes[0, 1].set_xticks([])  # Hide yticks for Bay of Bengal panel
axes[0, 1].plot([sw_b[0], sw_b[-1]], [pr_b[0], pr_b[-1]], color='black')


# Plot on the third subplot (bottom-left)
axes[1, 0].plot(sw_c, pr_c, color='black')
axes[1, 0].plot(np.tile(sw_c.data,19),pr20_c, color = "gray", alpha = 0.3)
axes[1, 0].plot(sw_c[0], pr_c[0],  marker = '.', color = "red")
axes[1, 0].plot(sw_c[11], pr_c[11],  marker = '.', color = "blue")
axes[1, 0].set_title('East Asia')
axes[1, 0].set_ylabel('Monsoon Rainfall  (mm/day)')
axes[1, 0].set_ylim([0,l])
axes[1, 0].set_xlabel('Solar Insolation ($W/m^2$)')
axes[1, 0].set_xlim([240,475])
axes[1, 0].plot([sw_c[0], sw_c[-1]], [pr_c[0], pr_c[-1]], color='black')
axes[1, 1].plot([sw_a[0], sw_a[-1]], [pr_a[0], pr_a[-1]], color='black')


# Plot on the fourth subplot (bottom-right)
axes[1, 1].plot(sw_a, pr_a, color='black')
axes[1, 1].plot(np.tile(sw_a.data,19),pr20_a, color = "gray", alpha = 0.3)
axes[1, 1].plot(sw_a[0], pr_a[0],  marker = '.', color = "red")
axes[1, 1].plot(sw_a[11], pr_a[11],  marker = '.', color = "blue")
axes[1, 1].set_title('Australia')
axes[1, 1].set_ylim([0,l])
axes[1, 1].set_yticks([])  # Hide yticks for Bay of Bengal panel
axes[1, 1].set_xlabel('Solar Insolation ($W/m^2$)')
axes[1, 1].set_xlim([240,475])


# Adjust layout to prevent overlapping
plt.tight_layout()
plt.savefig(savedir + r"/hysteresis_observations.pdf")

#%%
### PLOT: Hysteresis tsurf
#-------------------

plt.figure()
plt.plot(tsurf, pr_il2, color = "black", marker = ".")
plt.plot(tsurf_18, pr20_il2_18, color = "gray", alpha = 0.3)
plt.plot(tsurf[0], pr_il2[0],  marker = '.', color = "red")
plt.plot(tsurf[11], pr_il2[11],  marker = '.', color = "blue")
plt.title('India')
plt.ylabel('Monsoon Rainfall (mm/day)')
plt.xlim([25,55])
plt.plot([tsurf[0], tsurf[-1]], [pr_il2[0], pr_il2[-1]], color='black')
plt.xlabel('Surface temperature ($°C$)')
plt.tight_layout()

plt.savefig(savedir + r'hysteresis_observations_tsurf.pdf')


#%%
### PLOT: Hysteresis 3 regions
#-------------------

# Create a 1x3 panel plot with quadratic panels
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(14, 5))
l = 16

# Plot on the second subplot (top-left)
axes[0].plot(sw_b, pr_b, color='black', marker='.')
axes[0].plot(np.tile(sw_b.data, 19), pr20_b, color="gray", alpha=0.3)
axes[0].plot(sw_b[0], pr_b[0], marker='.', color="red")
axes[0].plot(sw_b[11], pr_b[11], marker='.', color="blue")
axes[0].set_title('Bay of Bengal', fontsize=14)
axes[0].set_ylim([0, l])
axes[0].set_xlim([240, 475])
axes[0].plot([sw_b[0], sw_b[-1]], [pr_b[0], pr_b[-1]], color='black')
axes[0].set_xlabel('Solar Insolation ($W/m^2$)', fontsize=14)
axes[0].set_ylabel('Monsoon Rainfall ($kg/m^2$)', fontsize=14)
axes[0].set_xlim([240, 475])
axes[0].tick_params(axis='both', which='both', labelsize=14)

# Plot on the third subplot (top-right)
axes[1].plot(sw_c, pr_c, color='black', marker='.')
axes[1].plot(np.tile(sw_c.data, 19), pr20_c, color="gray", alpha=0.3)
axes[1].plot(sw_c[0], pr_c[0], marker='.', color="red")
axes[1].plot(sw_c[11], pr_c[11], marker='.', color="blue")
axes[1].set_title('East Asia', fontsize=14)
axes[1].set_ylim([0, l])
axes[1].set_xlabel('Solar Insolation ($W/m^2$)', fontsize=14)
axes[1].set_xlim([240, 475])
axes[1].plot([sw_c[0], sw_c[-1]], [pr_c[0], pr_c[-1]], color='black')
axes[1].tick_params(axis='both', which='both', labelsize=14)
axes[1].set_yticks([]) 

# Plot on the fourth subplot (bottom-left)
axes[2].plot(sw_a, pr_a, color='black', marker='.')
axes[2].plot(np.tile(sw_a.data, 19), pr20_a, color="gray", alpha=0.3)
axes[2].plot(sw_a[0], pr_a[0], marker='.', color="red")
axes[2].plot(sw_a[11], pr_a[11], marker='.', color="blue")
axes[2].set_title('Australia', fontsize=14)
axes[2].set_ylim([0, l])
axes[2].set_yticks([])  # Hide yticks for Bay of Bengal panel
axes[2].set_xlabel('Solar Insolation ($W/m^2$)', fontsize=14)
axes[2].set_xlim([240, 475])
axes[2].plot([sw_a[0], sw_a[-1]], [pr_a[0], pr_a[-1]], color='black')
axes[2].tick_params(axis='both', which='both', labelsize=14)


plt.tight_layout()
plt.subplots_adjust(wspace=0.2)  # Default white space between all subplots

plt.savefig(savedir + r"/hysteresis_observations_updated.pdf")

#%%
### PLOT: Hysteresis 1 region (India)
#-----
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))
l = 16

# Plot on the first subplot
axes.plot(sw_il2, pr_il2, color='black', marker='.')
axes.plot(np.tile(sw_il2.data, 19), pr20_il2, color="gray", alpha=0.3)
axes.plot(sw_il2[0], pr_il2[0], marker='.', color="red")
axes.plot(sw_il2[11], pr_il2[11], marker='.', color="blue")
axes.set_title('India', fontsize=14)
axes.set_ylim([0, l])
axes.set_xlim([240, 475])
axes.plot([sw_il2[0], sw_il2[-1]], [pr_il2[0], pr_il2[-1]], color='black')
axes.set_xlabel('Solar Insolation ($W/m^2$)', fontsize=14)
axes.set_ylabel('Monsoon Rainfall ($kg/m^2$)', fontsize=14)
axes.set_xlim([240, 475])
axes.tick_params(axis='both', which='both', labelsize=14)
axes.set_yticks([0, 5, 10, 15])
axes.set_xticks([300, 400])

plt.tight_layout()
plt.savefig(savedir + r"/hysteresis_india.pdf")



#%%
### PLOT: Hysteresis 1 region (Bay of Bengal)
#-----

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))
l = 16

# Plot on the first subplot
axes.plot(sw_b, pr_b, color='black', marker='.')
axes.plot(np.tile(sw_b.data, 19), pr20_b, color="gray", alpha=0.3)
axes.plot(sw_b[0], pr_b[0], marker='.', color="red")
axes.plot(sw_b[11], pr_b[11], marker='.', color="blue")
axes.set_title('Bay of Bengal', fontsize=14)
axes.set_ylim([0, l])
axes.set_xlim([240, 475])
axes.plot([sw_b[0], sw_b[-1]], [pr_b[0], pr_b[-1]], color='black')
axes.set_xlabel('Solar Insolation ($W/m^2$)', fontsize=14)
axes.set_ylabel('Monsoon Rainfall ($kg/m^2$)', fontsize=14)
axes.set_xlim([240, 475])
axes.tick_params(axis='both', which='both', labelsize=14)
axes.set_yticks([0, 5, 10, 15])
axes.set_xticks([300, 400])

plt.tight_layout()
plt.savefig(savedir + r"/hysteresis_bayofbengal.pdf")



#%%
### PLOT: Hysteresis 1 region (East Asia)
#-----
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))
l = 16

# Plot on the first subplot
axes.plot(sw_c, pr_c, color='black', marker='.')
axes.plot(np.tile(sw_c.data, 19), pr20_c, color="gray", alpha=0.3)
axes.plot(sw_c[0], pr_c[0], marker='.', color="red")
axes.plot(sw_c[11], pr_c[11], marker='.', color="blue")
axes.set_title('East Asia', fontsize=14)
axes.set_ylim([0, l])
axes.set_xlim([240, 475])
axes.plot([sw_c[0], sw_c[-1]], [pr_c[0], pr_c[-1]], color='black')
axes.set_xlabel('Solar Insolation ($W/m^2$)', fontsize=14)
axes.set_ylabel('Monsoon Rainfall ($kg/m^2$)', fontsize=14)
axes.set_xlim([240, 475])
axes.tick_params(axis='both', which='both', labelsize=14)
axes.set_yticks([0, 5, 10, 15])
axes.set_xticks([300, 400])

plt.tight_layout()
plt.savefig(savedir + r"/hysteresis_eastasia.pdf")


#%%
### PLOT: Hysteresis 1 region (Australia)
#-----

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))
l = 16

# Plot on the first subplot
axes.plot(sw_a, pr_a, color='black', marker='.')
axes.plot(np.tile(sw_a.data, 19), pr20_a, color="gray", alpha=0.3)
axes.plot(sw_a[0], pr_a[0], marker='.', color="red")
axes.plot(sw_a[11], pr_a[11], marker='.', color="blue")
axes.set_title('Australia', fontsize=14)
axes.set_ylim([0, l])
axes.set_xlim([240, 475])
axes.plot([sw_a[0], sw_a[-1]], [pr_a[0], pr_a[-1]], color='black')
axes.set_xlabel('Solar Insolation ($W/m^2$)', fontsize=14)
axes.set_ylabel('Monsoon Rainfall ($kg/m^2$)', fontsize=14)
axes.set_xlim([240, 475])
axes.tick_params(axis='both', which='both', labelsize=14)
axes.set_yticks([0, 5, 10, 15])
axes.set_xticks([300, 400])

plt.tight_layout()
plt.savefig(savedir + r"/hysteresis_australia.pdf")

#%%
# PLOT: Precipitation, Solar Radiation, Land Surface Temperature throughout the year

fig, ax1 = plt.subplots(figsize=(10, 6)) 

# SW plot
ax1.plot(range(1,13), sw_il2, marker='.', color='darkorange')
ax1.set_ylabel('SW ($W/m^2$)', color='darkorange')
ax1.tick_params(axis='y', labelcolor='darkorange')
ax1.set_xlabel('Months')

# TSURF plot
ax2 = ax1.twinx()
ax2.spines['left'].set_position(('outward', 80))
ax2.yaxis.set_ticks_position('left')
ax2.yaxis.set_label_position('left')
ax2.plot(range(1,13), tsurf, marker='.', color='darkred')
ax2.set_ylabel('TSURF (°C)', color='darkred')
ax2.tick_params(axis='y', labelcolor='darkred')

# PRECIP plot
ax3 = ax1.twinx()
ax3.spines['right'].set_position(('outward', 0))
ax3.plot(range(1,13), pr_il2, marker='.', color='darkblue')
ax3.set_ylabel('Monsoon Rainfall ($mm/day$)', color='darkblue')
ax3.set_ylim([0, 11])
ax3.tick_params(axis='y', labelcolor='darkblue')


fig.tight_layout()
plt.xlabel('Months')
plt.title('India')
plt.savefig(savedir + '/observations_year_india.pdf', bbox_inches='tight')
# %%
