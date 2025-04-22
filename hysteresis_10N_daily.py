#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# AUTHOR: Anja Katzenberger

# This code creates the hysteresis plots for the Monsoon Planet
# Inclung precipitation depending on solar radiation and surface temperature
# and water vapour content depending on solar radiation and surface temperature

# The code is based on daily Monsoon Planet data that is averaged over 10 years

#%%
### LOAD PACKAGES
#-------------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pyts.decomposition import SingularSpectrumAnalysis

save_dir = 'C:/Users/anjaka/Nextcloud/PhD/03_MonsoonPlanet_Hysteresis/Figures'

#%%
###  LOAD DATA
#-------------------
data_dir = 'C:/Users/anjaka/Nextcloud/PhD/03_MonsoonPlanet_Hysteresis/data/slab_10years_ydaymean'

slab_list = [50,100,200,500]  # available slabs
slab_list_sel = [50,100,200,500] # selected slabs

# Initialize dictionaries to store the data for each slab
ds_dict = {}

pr = {}
sw = {}
tsurf = {}
wvp = {}

for slab in slab_list:
    # Load the dataset
    ds_dict[slab] = xr.open_dataset(f'{data_dir}/slab{slab}m.nc')
    
    # Extract the variables
    sw[slab] = ds_dict[slab]['swdn_toa']
    tsurf[slab] = ds_dict[slab]['t_surf']-273.15
    pr[slab] = ds_dict[slab]['precip']*86400
    wvp[slab] = ds_dict[slab]['WVP']



#%%
#  PROCESS DATA
#-------------------

l_m = 10
u_m = 20


# Initialize dictionaries to store the results
pr_m = {}
sw_m = {}
tsurf_m = {}
wvp_m = {}

pr_m_z = {}
sw_m_z = {}
tsurf_m_z = {}
wvp_m_z = {}

pr_m_z_mean = {}
sw_m_z_mean = {}
tsurf_m_z_mean = {}
wvp_m_z_mean = {}


for slab in slab_list:
    # Select the latitude slice and calculate the mean over longitude
    pr_m[slab] = pr[slab].sel(lat=slice(l_m,u_m)).mean('lon')
    sw_m[slab] = sw[slab].sel(lat=slice(l_m,u_m)).mean('lon')
    tsurf_m[slab] = tsurf[slab].sel(lat=slice(l_m,u_m)).mean('lon')
    wvp_m[slab] = wvp[slab].sel(lat=slice(l_m,u_m)).mean('lon')

    # Calculate the weighted mean over latitude
    weights = np.cos(np.deg2rad(pr_m[slab].lat))
    pr_m_z_mean[slab] = pr_m[slab].weighted(weights).mean(dim=['lat'])

    weights = np.cos(np.deg2rad(sw_m[slab].lat))
    sw_m_z_mean[slab] = sw_m[slab].weighted(weights).mean(dim=['lat'])

    weights = np.cos(np.deg2rad(tsurf_m[slab].lat))
    tsurf_m_z_mean[slab] = tsurf_m[slab].weighted(weights).mean(dim=['lat'])

    weights = np.cos(np.deg2rad(wvp_m[slab].lat))
    wvp_m_z_mean[slab] = wvp_m[slab].weighted(weights).mean(dim=['lat'])

#%%
# MERIDIONAL RADIATION PLOT
#-------------------

# select timesteps with same radiation (427W/m^2)
sw_march = sw[50].isel(time=84-1) # 25th march
sw_sept = sw[50].isel(time=252-1) # 9th september

# create meridional profiles
sw_march_lat = sw_march.mean('lon')
sw_sept_lat = sw_sept.mean('lon')

# plot
plt.figure()
plt.axvspan(10, 20, color='red', alpha=0.2, label = "Monsoon Region")
plt.plot(sw_march_lat.lat, sw_march_lat, label='25th March', color = "darkred")
plt.plot(sw_sept_lat.lat, sw_sept_lat, label='09th September', color = "darkorange")
plt.axvline(x=0, color='black', linestyle='--')
plt.axvline(x=10, color='black', linestyle='-')
plt.axvline(x=60, color='black', linestyle='-')
plt.xlabel('Latitude')
plt.ylabel('Solar Radiation ($W/m^2$)')
plt.legend(loc = 'upper left', frameon = False,bbox_to_anchor=(1, 0.5))
plt.xlim(-90,90)
xticks = [-90, -60, -30, 0, 30, 60, 90]
xtick_labels = ['-90', '-60', '-30', '0', '30', '60', '90']
plt.xticks(ticks=xticks, labels=xtick_labels)
plt.savefig(save_dir + '/MonsoonPlanet_Radiation.pdf', bbox_inches='tight')

#%%
# Singular Spectrum Analysis
#-------------------

# Singular Spectrum Analysis
L = 20 # window_size

pr_m_z_mean_ssa = {}
tsurf_m_z_mean_ssa = {}
sw_m_z_mean_ssa = {}
wvp_m_z_mean_ssa = {}

# Perform the Singular Spectrum Analysis for 'tsurf'
for i in range(len(slab_list)):
    F = tsurf_m_z_mean[slab_list[i]]
    F_arr = np.array(F)
    F_in = np.array([F_arr])
    ssa = SingularSpectrumAnalysis(window_size = L)
    X_ssa = ssa.transform(F_in)
    tsurf_m_z_mean_ssa[slab_list[i]] = X_ssa[0, 0, :]

# Perform the Singular Spectrum Analysis for 'sw'
for i in range(len(slab_list)):
    F = sw_m_z_mean[slab_list[i]]
    F_arr = np.array(F)
    F_in = np.array([F_arr])
    ssa = SingularSpectrumAnalysis(window_size = L)
    X_ssa = ssa.transform(F_in)
    sw_m_z_mean_ssa[slab_list[i]] = X_ssa[0, 0, :]

# Perform the Singular Spectrum Analysis for 'pr'
for i in range(len(slab_list)):
    F = pr_m_z_mean[slab_list[i]]
    F_arr = np.array(F)
    F_in = np.array([F_arr])
    ssa = SingularSpectrumAnalysis(window_size = L)
    X_ssa = ssa.transform(F_in)
    pr_m_z_mean_ssa[slab_list[i]] = X_ssa[0, 0, :]

# Perform the Singular Spectrum Analysis for 'wvp'
for i in range(len(slab_list)):
    F = wvp_m_z_mean[slab_list[i]]
    F_arr = np.array(F)
    F_in = np.array([F_arr])
    ssa = SingularSpectrumAnalysis(window_size = L)
    X_ssa = ssa.transform(F_in)
    wvp_m_z_mean_ssa[slab_list[i]] = X_ssa[0, 0, :]

#%%
# PLOTTING DAILY DATA WITHOUT SSA 
#-------------------   

# Create a 2x2 panel
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

# Define slab depths
slab_list = [50, 100, 200, 500]

# Define colors for each slab depth
colors = ["#DC143C", "#FF8000", "#33A1C9", "#191970"]

# Plot 1
for i, slab in enumerate(slab_list):
    axs[1, 0].plot(sw_m_z_mean[slab], wvp_m_z_mean[slab], label=f'{slab}m', color=colors[i])
axs[1, 0].set_xlabel('Solar Radiation ($W/m^2$)')
axs[1, 0].set_ylabel('Water Vapour Content ($kg/m^2$)')

# Plot 2
for i, slab in enumerate(slab_list):
    axs[0, 0].plot(sw_m_z_mean[slab], pr_m_z_mean[slab], label=f'{slab}m', color=colors[i])
axs[0, 0].set_xlabel('Solar Radiation ($W/m^2$)')
axs[0, 0].set_ylabel('Monsoon Rainfall ($kg/m^2$)')

# Plot 3
for i, slab in enumerate(slab_list):
    axs[0, 1].plot(tsurf_m_z_mean[slab], pr_m_z_mean[slab], label=f'{slab}m', color=colors[i])
axs[0, 1].set_xlabel('Surface Temperature ($°C$)')
axs[0, 1].set_ylabel('Monsoon Rainfall ($kg/m^2$)')

# Plot 4
for i, slab in enumerate(slab_list):
    axs[1, 1].plot(tsurf_m_z_mean[slab], wvp_m_z_mean[slab], label=f'{slab}m', color=colors[i])
axs[1, 1].set_xlabel('Surface Temperature ($°C$)')
axs[1, 1].set_ylabel('Water Vapour Content ($kg/m^2$)')

# Adjust layout to prevent clipping of titles
plt.tight_layout()
fig.subplots_adjust(bottom=0.1)

plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol=4, title = "Slab Depth")
plt.tight_layout()

plt.savefig(save_dir + '/MonsoonPlanet_Hysteresis_4panels_daily.pdf', bbox_inches='tight')


#%%
# PLOTTING DAILY DATA WITH SSA
#-------------------   

# Create a 2x2 panel
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

# Define colors for each slab depth
colors = ["#DC143C", "#FF8000", "#33A1C9", "#191970"]

# Plot 1
for i, slab in enumerate(slab_list):
    sw = sw_m_z_mean_ssa[slab]
    wvp = wvp_m_z_mean_ssa[slab]
    axs[1, 0].plot(sw, wvp, label=f'{slab}m', color=colors[i])
    
    # Connect the first and last data points
    axs[1, 0].plot(sw[[0, -1]], wvp[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[1, 0].plot(sw[0], wvp[0], marker=">", color=colors[i])

axs[1, 0].set_xlabel('Solar Radiation ($W/m^2$)')
axs[1, 0].set_ylabel('Water Vapour Content ($kg/m^2$)')

# Plot 2
for i, slab in enumerate(slab_list):
    sw = sw_m_z_mean_ssa[slab]
    pr = pr_m_z_mean_ssa[slab]
    axs[0, 0].plot(sw, pr, label=f'{slab}m', color=colors[i])
    
    # Connect the first and last data points
    axs[0, 0].plot(sw[[0, -1]], pr[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[0, 0].plot(sw[0], pr[0], marker=">", color=colors[i])

axs[0, 0].set_xlabel('Solar Radiation ($W/m^2$)')
axs[0, 0].set_ylabel('Monsoon Rainfall ($kg/m^2$)')

# Plot 3
for i, slab in enumerate(slab_list):
    tsurf = tsurf_m_z_mean_ssa[slab]
    pr = pr_m_z_mean_ssa[slab]
    axs[0, 1].plot(tsurf, pr, label=f'{slab}m', color=colors[i])
    
    # Connect the first and last data points
    axs[0, 1].plot(tsurf[[0, -1]], pr[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[0, 1].plot(tsurf[0], pr[0], marker=">", color=colors[i])

axs[0, 1].set_xlabel('Surface Temperature ($°C$)')
axs[0, 1].set_ylabel('Monsoon Rainfall ($kg/m^2$)')

# Plot 4
for i, slab in enumerate(slab_list):
    tsurf = tsurf_m_z_mean_ssa[slab]
    wvp = wvp_m_z_mean_ssa[slab]
    axs[1, 1].plot(tsurf, wvp, label=f'{slab}m', color=colors[i])
    
    # Connect the first and last data points
    axs[1, 1].plot(tsurf[[0, -1]], wvp[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[1, 1].plot(tsurf[0], wvp[0], marker=">", color=colors[i])

axs[1, 1].set_xlabel('Surface Temperature ($°C$)')
axs[1, 1].set_ylabel('Water Vapour Content ($kg/m^2$)')


# Adjust layout to prevent clipping of titles
plt.tight_layout()
fig.subplots_adjust(bottom=0.1)

# Set x-tick and y-tick labels font size for all subplots
for ax in axs.flatten():
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

# Move the legend to the right side
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 0.5), fancybox=True, shadow=True, ncol=1, title="Slab Depth",
           frameon=False, fontsize=14, title_fontsize=14)

# Adjust layout to prevent clipping of titles
plt.tight_layout()
fig.subplots_adjust(bottom=0.1, right=0.8)

plt.savefig(save_dir + '/MonsoonPlanet_Hysteresis_4panels.pdf', bbox_inches='tight')


# %%
### Plot only main hysteresis 
#-------------------

for i, slab in enumerate(slab_list):
    sw = sw_m_z_mean_ssa[slab]
    pr = pr_m_z_mean_ssa[slab]
    plt.plot(sw, pr, label=f'{slab}m', color=colors[i])
    
    # Connect the first and last data points
    plt.plot(sw[[0, -1]], pr[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    plt.plot(sw[0], pr[0], marker=">", color=colors[i])
    plt.plot(sw[6], pr[6], marker=".", color=colors[i])


plt.xlabel('Solar Radiation ($W/m^2$)')
plt.ylabel('Monsoon Rainfall ($kg/m^2$)')

plt.savefig(save_dir + '/MonsoonPlanet_Hysteresis_1panel.pdf', bbox_inches='tight')

#%%

# Plot only hysteresis depending on solar radiation 
#-------------------

fig, axs = plt.subplots(1, 2, figsize=(11, 5))  # Change the values as needed

# Define colors for each slab depth
colors = ["#DC143C", "#FF8000", "#33A1C9", "#191970"]

# Plot 1
for i, slab in enumerate(slab_list):
    sw = sw_m_z_mean_ssa[slab]
    wvp = wvp_m_z_mean_ssa[slab]
    axs[1].plot(sw, wvp, label=f'{slab}m', color=colors[i])
    
    # Connect the first and last data points
    axs[1].plot(sw[[0, -1]], wvp[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[1].plot(sw[0], wvp[0], marker=">", color=colors[i], markersize = 14)
    axs[1].plot(sw[183], wvp[183], marker=".", color=colors[i], markersize = 14)


axs[1].set_xlabel('Solar Radiation ($W/m^2$)')
axs[1].set_ylabel('Water Vapour Content ($kg/m^2$)')

# Plot 2
for i, slab in enumerate(slab_list):
    sw = sw_m_z_mean_ssa[slab]
    pr = pr_m_z_mean_ssa[slab]
    axs[0].plot(sw, pr, label=f'{slab}m', color=colors[i])
    
    # Connect the first and last data points
    axs[0].plot(sw[[0, -1]], pr[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[0].plot(sw[0], pr[0], marker=">", color=colors[i], markersize = 14)
    axs[0].plot(sw[183], pr[183], marker=".", color=colors[i], markersize = 14)


#axs[0].legend(loc='best', title = "Slab Depth", frameon = False, bbox_to_anchor=(0.35, 0.85)) 
axs[0].set_xlabel('Solar Radiation ($W/m^2$)')
axs[0].set_ylabel('Monsoon Rainfall ($mm/day$)')

# Or use text function to add labels at specific positions
axs[0].text(0.15, 0.95, 'A.', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
axs[1].text(0.15, 0.95, 'B.', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')


plt.savefig(save_dir + '/MonsoonPlanet_Hysteresis_2panel_sw.pdf', bbox_inches='tight')

#%%

# Plot only hysteresis depending on tsurf
#-------------------

fig, axs = plt.subplots(1, 2, figsize=(10, 5))  # Change the values as needed

# Define colors for each slab depth
colors = ["#DC143C", "#FF8000", "#33A1C9", "#191970"]

# Plot 1
for i, slab in enumerate(slab_list):
    tsurf = tsurf_m_z_mean_ssa[slab]
    wvp = wvp_m_z_mean_ssa[slab]
    axs[1].plot(tsurf, wvp, label=f'{slab}m', color=colors[i])
    
    # Connect the first and last data points
    axs[1].plot(tsurf[[0, -1]], wvp[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[1].plot(tsurf[0], wvp[0], marker=">", color=colors[i])
    axs[1].plot(tsurf[183], wvp[183], marker="o", color=colors[i])

axs[1].set_xlabel('Surface Temperature ($°C$)')
axs[1].set_ylabel('Water Vapour Content ($kg/m^2$)')

# Plot 2
for i, slab in enumerate(slab_list):
    tsurf = tsurf_m_z_mean_ssa[slab]
    pr = pr_m_z_mean_ssa[slab]
    axs[0].plot(tsurf, pr, label=f'{slab}m', color=colors[i])
    
    # Connect the first and last data points
    axs[0].plot(tsurf[[0, -1]], pr[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[0].plot(tsurf[0], pr[0], marker=">", color=colors[i])
    axs[0].plot(tsurf[183], pr[183], marker="o", color=colors[i])

#axs[1].legend(loc='best', title = "Slab Depth", frameon = False)
axs[0].set_xlabel('Surface temperature ($°C$)')
axs[0].set_ylabel('Monsoon Rainfall ($mm/day$)')

plt.savefig(save_dir + '/MonsoonPlanet_Hysteresis_2panel_tsurf.pdf', bbox_inches='tight')

# %%
# HYSTERESIS WITH STATES
#-------------------

sw = sw_m_z_mean_ssa[slab]
wvp = wvp_m_z_mean_ssa[slab]
plt.plot(sw, wvp, label=f'{slab}m', color=colors[0])
    
# Connect the first and last data points
plt.plot(sw[[0, -1]], wvp[[0, -1]], color=colors[0], linestyle='-')
    
# Mark only the first data point with a marker "."
plt.plot(sw[0], wvp[0], marker=">", color=colors[0])

plt.xlabel('Solar Radiation ($W/m^2$)')
plt.ylabel('Water Vapour Content ($kg/m^2$)')

# %%
