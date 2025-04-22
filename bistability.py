# Author: Anja Katzenberger, anja.katzenberger@pik-potsdam.de

# This code creates figures of bistability
# Two stable rainfall states are possible for the same insolation
# Timeseries of insolation and precipitation are presented 

#%%
### LOAD PACKAGES
#-------------------

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
from pyts.decomposition import SingularSpectrumAnalysis


# Where to save the figures
savedir = r'C:\Users\anjaka\Nextcloud\PhD\03_MonsoonPlanet_Hysteresis\Figures\bistability\\'

# Directories 
dir = r"C:\Users\anjaka\Nextcloud\PhD\03_MonsoonPlanet_Hysteresis\Codes_and_Data_for_Monsoon_Hysteresis_reveals_Atmospheric_Memory\data_reviewprocess\bistability"
files = os.listdir(dir)

#%%
### LOAD DATA
#-------------------

# Get list of keys for down simulations
keys_down = [file.replace("output_", "").replace("down_sel.nc", "") for file in files if file.startswith("output_") and file.endswith("down_sel.nc")]
keys_up =  [file.replace("output_", "").replace("up_sel.nc", "") for file in files if file.startswith("output_") and file.endswith("up_sel.nc")]


sw_dict_down = {}
pr_dict_down = {}
for key in keys_down:
    data = xr.open_dataset(dir + r"\output_" + key + "down_sel.nc")
    sw_dict_down[key] = data["swdn_toa"]
    pr_dict_down[key] = data["precip"]*86400

sw_dict_up = {}
pr_dict_up = {}
for key in keys_up:
    data = xr.open_dataset(dir + r"\output_" + key + "up_sel.nc")
    sw_dict_up[key] = data["swdn_toa"]
    pr_dict_up[key] = data["precip"]*86400

    
    

#%%
### PROCESSING DATA
#-------------------
    
# Monsoon region
lat_low = 10
lat_up = 20

pr_down = {}
sw_down = {}
for key in keys_down:
    pr_lat = pr_dict_down[key].mean('lon')
    pr_lat_mr = pr_lat.sel(lat=slice(lat_low,lat_up))
    weights = np.cos(np.deg2rad(pr_lat_mr.lat))
    pr_down[key] = pr_lat_mr.weighted(weights).mean(dim=['lat'])

    sw_lat = sw_dict_down[key].mean('lon')
    sw_lat_mr = sw_lat.sel(lat=slice(lat_low,lat_up))
    weights = np.cos(np.deg2rad(sw_lat_mr.lat))
    sw_down[key] = sw_lat_mr.weighted(weights).mean(dim=['lat'])

pr_up = {}
sw_up = {}
for key in keys_up:
    pr_lat = pr_dict_up[key].mean('lon')
    pr_lat_mr = pr_lat.sel(lat=slice(lat_low,lat_up))
    weights = np.cos(np.deg2rad(pr_lat_mr.lat))
    pr_up[key] = pr_lat_mr.weighted(weights).mean(dim=['lat'])

    sw_lat = sw_dict_up[key].mean('lon')
    sw_lat_mr = sw_lat.sel(lat=slice(lat_low,lat_up))
    weights = np.cos(np.deg2rad(sw_lat_mr.lat))
    sw_up[key] = sw_lat_mr.weighted(weights).mean(dim=['lat'])


#%%
### SINGULAR SPECTRUM ANALYSIS
#-------------------

# Smoothening time series by using Singular Spectrum Analysis    
L = 10 # window size

pr_down_ssa = {}
for key in keys_down:
    F = pr_down[key]
    F_arr = np.array(F)
    F_in = np.array([F_arr])
    ssa = SingularSpectrumAnalysis(window_size = L)
    X_ssa = ssa.transform(F_in)
    pr_down_ssa[key] = X_ssa[0, 0, :]

#%%
pr_up_ssa = {}
for key in keys_up:
    F = pr_up[key]
    F_arr = np.array(F)
    F_in = np.array([F_arr])
    ssa = SingularSpectrumAnalysis(window_size = L)
    X_ssa = ssa.transform(F_in)
    pr_up_ssa[key] = X_ssa[0, 0, :]



# %%
### OVERVIEW: ALL RUNS 2x2
#----------------------------

fig, axs = plt.subplots(2, 2, figsize=(10, 10))  # Create a 2x2 grid of subplots

# PLOT 1: SW DOWN
for key in keys_down:
    if key != '0' and key.endswith('firstyeardaily_'):
        axs[0, 1].plot(sw_down[key], label=key)
axs[0, 1].plot(sw_down['0'], label="Standard", linestyle = '--', color = "black")
#axs[0, 1].set_ylabel("Insolation ($W/m^2$)")
#axs[0, 1].set_xlabel("Days of the Year")
axs[0, 1].set_xticks([])
axs[0, 1].set_yticks([])


# PLOT 2: PR DOWN
for key in keys_down:
    if key != '0'  and key.endswith('firstyeardaily_'):
        axs[1, 1].plot(pr_down_ssa[key], label=key)
axs[1, 1].plot(pr_down_ssa['0'], label="Standard", linestyle = '--', color = "black")
axs[1, 1].legend()
#axs[1, 1].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[1, 1].set_yticks([])
axs[1, 1].set_ylim(0, 13)
axs[1, 1].set_xlabel("Days of the Year")

# PLOT 3: SW UP
for key in keys_up:
    if key != '0' and key.endswith('firstyeardaily_'):
        axs[0, 0].plot(sw_up[key], label=key)
axs[0, 0].plot(sw_up['0'], label="Standard", linestyle = '--', color = "black")
axs[0, 0].set_ylabel("Insolation ($W/m^2$)")
#axs[0, 0].set_xlabel("Days of the Year")
axs[0, 0].set_xticks([])


# PLOT 4: PR UP
for key in keys_up:
    if key != '0' and key.endswith('firstyeardaily_'):
        axs[1, 0].plot(pr_up_ssa[key], label=key)
axs[1, 0].plot(pr_up_ssa['0'], label="Standard", linestyle = '--', color = "black")
axs[1, 0].legend()
axs[1, 0].set_ylim(0, 13)
axs[1, 0].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[1, 0].set_xlabel("Days of the Year")

plt.tight_layout()  # Adjust the layout to prevent overlap
plt.savefig(savedir + "bistability_2x2_panel.pdf", bbox_inches = 'tight')  # Save the figure


# %%
### PLOT: OVERVIEW PROCEDURE
#----------------------------

key_sel = '4276firstyeardaily_'
lim = 365

fig, axs = plt.subplots(2, 2, figsize=(10, 10))  # Create a 2x2 grid of subplots

# PLOT 1: SW DOWN
axs[0, 1].plot(sw_down[key_sel], color = "red")
axs[0, 1].plot(sw_down['0'], linestyle = '--', color = "black")
#axs[0, 1].set_ylabel("Insolation ($W/m^2$)")
#axs[0, 1].set_xlabel("Days of the Year")
axs[0, 1].set_xticks([])
axs[0, 1].set_yticks([])
axs[0, 1].set_xlim(0, lim)
axs[0, 1].axvline(252, color = "black", linestyle = "--")


# PLOT 2: PR DOWN
axs[1, 1].plot(pr_down_ssa[key_sel],  color = "red")
axs[1, 1].plot(pr_down_ssa['0'],  linestyle = '--', color = "black")
#axs[1, 1].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[1, 1].set_yticks([])
axs[1, 1].set_ylim(0, 13)
axs[1, 1].set_xlabel("Days of the Year")
axs[1, 1].set_xlim(0, lim)
axs[1, 1].axvline(252, color = "black", linestyle = "--")


# PLOT 3: SW UP
axs[0, 0].plot(sw_up[key_sel], color = "blue")
axs[0, 0].plot(sw_up['0'], linestyle = '--', color = "black")
axs[0, 0].set_ylabel("Insolation ($W/m^2$)")
#axs[0, 0].set_yticks([])
#axs[0, 0].set_xlabel("Days of the Year")
axs[0, 0].set_xticks([])
axs[0, 0].set_xlim(0, lim)
axs[0, 0].axvline(84, color = "black", linestyle = "--")

# PLOT 4: PR UP
axs[1, 0].plot(pr_up_ssa[key_sel], label=key_sel[0:3] + "$W/m^2$", color = "blue")
axs[1, 0].plot(pr_up_ssa['0'], label="Standard", linestyle = '--', color = "black")
axs[1, 0].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[1, 0].set_ylim(0, 13)
#axs[1, 0].set_yticks([])
axs[1, 0].set_xlabel("Days of the Year")
axs[1, 0].set_xlim(0, lim)
axs[1, 0].axvline(84, color = "black", linestyle = "--")


# Create a single legend for all subplots
fig.legend(loc='lower center', ncol=2,frameon = False)

plt.tight_layout()  # Adjust the layout to prevent overlap
plt.subplots_adjust(bottom=0.1)  # Make space for the legend
plt.savefig(savedir + "bistability_2x2_panel_selected.pdf", bbox_inches = 'tight')  # Save the figure

# %%
# INDIVIDUAL RUN - FIRST YEAR
#------------------

key_sel = '4272firstyeardaily_'
plt.figure(figsize = (5,5))
plt.plot(pr_down[key_sel], label=key_sel, color = "red")
plt.plot(pr_up[key_sel], label=key_sel, color = "blue")
plt.plot(pr_up_ssa['0'], label="Standard", linestyle = '--', color = "black")
plt.axvline(84, color = "black", linestyle = "--")
plt.axvline(252, color = "black", linestyle = "--")
plt.ylabel("Monsoon rainfall ($mm/day$)")
plt.ylim(0, 16)
plt.xlabel("Days of the Year")
plt.xlim(0,365)
plt.savefig(savedir + "bistability_1_panel_" + key_sel + ".pdf", bbox_inches = 'tight')  # Save the figure

# %%
# INDIVIDUAL RUN - 60 MONTHS LONGTERM
#------------------

key_sel = '4272monthly_'
key_firstyear = '4272firstyeardailymonmean_'

plt.figure(figsize = (5,5))
plt.plot(pr_down[key_sel], color = "red")
plt.plot(pr_up[key_sel], color = "blue")
plt.axvline(12, color = "green", linestyle = "--")
plt.ylabel("Monsoon rainfall ($mm/day$)")
plt.ylim(0, 11.75)
plt.xlabel("Months")
plt.xlim(0,60)
plt.savefig(savedir + "bistability_1_panel_longterm" + key_sel + "_60.pdf", bbox_inches = 'tight')  # Save the figure

# %%
# INDIVIDUAL RUN - 720 MONTHS LONGTERM
#------------------

key_sel = '4276monthly_'
key_firstyear = '4276firstyeardailymonmean_'

plt.figure(figsize = (5,5))
plt.plot(pr_down[key_sel], color = "red")
plt.plot(pr_up[key_sel], color = "blue")
plt.axvline(60, color = "green", linestyle = "--")
plt.ylabel("Monsoon rainfall ($mm/day$)")
plt.ylim(0, 11.75)
plt.xlabel("Months")
plt.xlim(0,900)
plt.savefig(savedir + "bistability_1_panel_longterm" + key_sel + ".pdf", bbox_inches = 'tight')  # Save the figure


# %%# %%
# PLOT: MULTI-RUN
#-------------------

#key_sel = '427'
#plt.figure(figsize = (5,5))
#plt.plot(pr_down_ssa[key_sel], label=key_sel, color = "red")
#plt.plot(pr_down_ssa[key_sel + '2'], label=key_sel, color = "red",alpha = 0.4)
#plt.plot(pr_down_ssa[key_sel + '3'], label=key_sel, color = "red", alpha = 0.5)
#plt.plot(pr_down_ssa[key_sel + '4'], label=key_sel, color = "red", alpha = 0.6)
#plt.plot(pr_down_ssa[key_sel + '6'], label=key_sel, color = "red",alpha = 0.7)
#plt.plot(pr_up_ssa[key_sel + '2'], color = "blue",alpha = 0.4)
#plt.plot(pr_up_ssa[key_sel + '3'], color = "blue",alpha = 0.5)
#plt.plot(pr_up_ssa[key_sel + '4'], color = "blue",alpha = 0.6)
#plt.plot(pr_up_ssa[key_sel + '5'], color = "blue",alpha = 0.7)
#plt.plot(pr_up_ssa[key_sel + '6'], color = "blue",alpha = 0.8)

# Put all the time series into a list
#time_series_list = [pr_down_ssa[key_sel][:617], pr_down_ssa[key_sel + '2'], pr_down_ssa[key_sel + '3'], pr_down_ssa[key_sel + '4'], pr_down_ssa[key_sel + '5']]
#time_series_array = np.array(time_series_list)
#mean_time_series = np.mean(time_series_array, axis=0)
#std_time_series = np.std(time_series_array, axis=0)
#upper_time_series = mean_time_series + std_time_series
#lower_time_series = mean_time_series - std_time_series
#plt.plot(mean_time_series, color = "darkred")
#plt.fill_between(range(len(mean_time_series)), lower_time_series, upper_time_series, color='darkred', alpha=0.3)

# Put all the time series into a list
#time_series_list = [pr_up_ssa[key_sel + '2'][:617], pr_up_ssa[key_sel + '3'][:617], pr_up_ssa[key_sel + '4'][:617], pr_up_ssa[key_sel + '5'][:617], pr_up_ssa[key_sel + '6'][:617]]
#time_series_array = np.array(time_series_list)
#mean_time_series_up = np.mean(time_series_array, axis=0)
#std_time_series_up = np.std(time_series_array, axis=0)
#upper_time_series_up = mean_time_series_up + std_time_series_up
#lower_time_series_up = mean_time_series_up - std_time_series_up
#plt.plot(mean_time_series_up, color = "darkblue")
#plt.fill_between(range(len(mean_time_series_up)), lower_time_series_up, upper_time_series_up, color='darkblue', alpha=0.3)


#plt.plot(pr_up_ssa['0'], label="Standard", linestyle = '--', color = "black")
#plt.ylabel("Monsoon rainfall ($mm/day$)")
#plt.ylim(0, 13)
#plt.xlim(0,400)
#plt.xlabel("Days of the Year")
#plt.savefig(savedir + "bistability_1_panel_" + key_sel + ".pdf", bbox_inches = 'tight')  # Save the figure

#%%
### One Run - 3 panels
#-------------------

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15,5))

# ENSEMBLE - 1 year
key_sel = '4276'


axs[0].text(0.04, 0.96, '5A.', horizontalalignment='left', verticalalignment='top', transform=axs[0].transAxes, fontsize=20, fontweight='bold')
axs[1].text(0.04, 0.96, '5B.', horizontalalignment='left', verticalalignment='top', transform=axs[1].transAxes, fontsize=20, fontweight='bold')
axs[2].text(0.04, 0.96, '5C.', horizontalalignment='left', verticalalignment='top', transform=axs[2].transAxes, fontsize=20, fontweight='bold')

axs[0].plot(pr_down[key_sel + 'firstyeardaily_'], color = "red")
axs[0].plot(pr_up[key_sel + 'firstyeardaily_'], color = "blue")
axs[0].axvline(252, color = "red", linestyle = "-", linewidth = 2)
axs[0].axvline(84, color = "blue", linestyle = "-")
axs[0].plot(pr_up_ssa['0'], label="Standard", linestyle = '--', color = "black")
axs[0].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[0].set_ylim(0, 13)
axs[0].set_xlim(0,365)
axs[0].set_xlabel("Days of the Year")

# ENSEMBLE - 60 months
axs[1].plot(pr_down[key_sel + 'monthly_'], label=key_sel, color = "red")
axs[1].plot(pr_up[key_sel + 'monthly_'], color = "blue")
axs[1].axvline(12, color = "green", linestyle = "--")
axs[1].set_ylim(0, 13)
axs[1].set_xlim(0,60)
axs[1].set_xlabel("Months")
axs[1].set_yticklabels([])

# ENSEMBLE - 360 months
axs[2].plot(pr_down[key_sel + 'monthly_'], label=key_sel, color = "red")
axs[2].plot(pr_up[key_sel  + 'monthly_'], color = "blue")
axs[2].axvline(60, color = "green", linestyle = "--")
axs[2].set_ylim(0, 13)
axs[2].set_xlim(0,360)
axs[2].set_xlabel("Months")
axs[2].set_yticklabels([])

plt.subplots_adjust(wspace=0.1)  # Adjust this value to bring the subplots closer together

plt.savefig(savedir + "bistability_singlerun_" + key_sel + "_combined.pdf", bbox_inches = 'tight')  # Save the figure

#%%
### One Run - 3 panels (SSA)
#-------------------

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15,5))

# ENSEMBLE - 1 year
key_sel = '4272'


axs[0].text(0.04, 0.96, 'A.', horizontalalignment='left', verticalalignment='top', transform=axs[0].transAxes, fontsize=20)
axs[1].text(0.04, 0.96, 'B.', horizontalalignment='left', verticalalignment='top', transform=axs[1].transAxes, fontsize=20)
axs[2].text(0.04, 0.96, 'C.', horizontalalignment='left', verticalalignment='top', transform=axs[2].transAxes, fontsize=20)

axs[0].plot(pr_down_ssa[key_sel + 'firstyeardaily_'], color = "red")
axs[0].plot(pr_up_ssa[key_sel + 'firstyeardaily_'], color = "blue")
axs[0].axvline(252, color = "black", linestyle = "--")
axs[0].axvline(84, color = "black", linestyle = "--")
axs[0].plot(pr_up_ssa['0'], label="Standard", linestyle = '--', color = "black")
axs[0].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[0].set_ylim(0, 13)
axs[0].set_xlim(0,365)
axs[0].set_xlabel("Days of the Year")

# ENSEMBLE - 60 months
axs[1].plot(pr_down_ssa[key_sel + 'monthly_'], label=key_sel, color = "red")
axs[1].plot(pr_up_ssa[key_sel + 'monthly_'], color = "blue")
axs[1].axvline(12, color = "green", linestyle = "--")
axs[1].set_ylim(0, 12)
axs[1].set_xlim(0,60)
axs[1].set_xlabel("Months")
axs[1].set_yticklabels([])

# ENSEMBLE - 360 months
axs[2].plot(pr_down_ssa[key_sel + 'monthly_'], label=key_sel, color = "red")
axs[2].plot(pr_up_ssa[key_sel  + 'monthly_'], color = "blue")
axs[2].axvline(60, color = "green", linestyle = "--")
axs[2].set_ylim(0, 12)
axs[2].set_xlim(0,360)
axs[2].set_xlabel("Months")
axs[2].set_yticklabels([])

plt.subplots_adjust(wspace=0.1)  # Adjust this value to bring the subplots closer together

plt.savefig(savedir + "bistability_singlerun_" + key_sel + "_combined_ssa.pdf", bbox_inches = 'tight')  # Save the figure


# %%# %%
# ENSEMBLE - 1 year
#-------------------

key_sel = '427'
plt.figure(figsize = (5,5))
plt.plot(pr_down[key_sel + '1' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.4)
plt.plot(pr_down[key_sel + '2' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.4)
#plt.plot(pr_down[key_sel + '3' + 'firstyeardaily_'], label=key_sel, color = "red", alpha = 0.5)
#plt.plot(pr_down[key_sel + '4' + 'firstyeardaily_'], label=key_sel, color = "red", alpha = 0.6)
plt.plot(pr_down[key_sel + '6' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.7)
#plt.plot(pr_up[key_sel + '1' + 'firstyeardaily_'], color = "blue",alpha = 0.4)
plt.plot(pr_up[key_sel + '2' + 'firstyeardaily_'], color = "blue",alpha = 0.4)
#plt.plot(pr_up[key_sel + '3' + 'firstyeardaily_'], color = "blue",alpha = 0.5)
#plt.plot(pr_up[key_sel + '4' + 'firstyeardaily_'], color = "blue",alpha = 0.6)
plt.plot(pr_up[key_sel + '6' + 'firstyeardaily_'], color = "blue",alpha = 0.8)
plt.axvline(252, color = "black", linestyle = "--")
plt.axvline(84, color = "black", linestyle = "--")

plt.plot(pr_up_ssa['0'], label="Standard", linestyle = '--', color = "black")
plt.ylabel("Monsoon rainfall ($mm/day$)")
plt.ylim(0, 13)
plt.xlim(0,365)
plt.xlabel("Days of the Year")
plt.savefig(savedir + "bistability_ensemble_" + key_sel + "_1year.pdf", bbox_inches = 'tight')  # Save the figure

# %%# %%
# ENSEMBLE - 60 months
#-------------------

key_sel = '427'
plt.figure(figsize = (5,5))
plt.plot(pr_down[key_sel + '1' + 'monthly_'], label=key_sel, color = "red",alpha = 0.4)
plt.plot(pr_down[key_sel + '2' + 'monthly_'], label=key_sel, color = "red",alpha = 0.4)
#plt.plot(pr_down[key_sel + '3' + 'monthly_'], label=key_sel, color = "red", alpha = 0.5)
#plt.plot(pr_down[key_sel + '4' + 'monthly_'], label=key_sel, color = "red", alpha = 0.6)
plt.plot(pr_down[key_sel + '6' + 'monthly_'], label=key_sel, color = "red",alpha = 0.7)
#plt.plot(pr_up[key_sel + '1' + 'monthly_'], color = "blue",alpha = 0.4)
plt.plot(pr_up[key_sel + '2' + 'monthly_'], color = "blue",alpha = 0.4)
#plt.plot(pr_up[key_sel + '3' + 'monthly_'], color = "blue",alpha = 0.5)
#plt.plot(pr_up[key_sel + '4' + 'monthly_'], color = "blue",alpha = 0.6)
plt.plot(pr_up[key_sel + '6' + 'monthly_'], color = "blue",alpha = 0.8)
plt.axvline(12, color = "black", linestyle = "--")

plt.ylabel("Monsoon rainfall ($mm/day$)")
plt.ylim(0, 13)
plt.xlim(0,60)
plt.xlabel("Months")
plt.savefig(savedir + "bistability_ensemble_" + key_sel + "_60months.pdf", bbox_inches = 'tight')  # Save the figure


# %%# %%
# ENSEMBLE - 360 months
#-------------------

key_sel = '427'
plt.figure(figsize = (5,5))
plt.plot(pr_down[key_sel + '1' + 'monthly_'], label=key_sel, color = "red",alpha = 0.4)
plt.plot(pr_down[key_sel + '2' + 'monthly_'], label=key_sel, color = "red",alpha = 0.4)
#plt.plot(pr_down[key_sel + '3' + 'monthly_'], label=key_sel, color = "red", alpha = 0.5)
#plt.plot(pr_down[key_sel + '4' + 'monthly_'], label=key_sel, color = "red", alpha = 0.6)
plt.plot(pr_down[key_sel + '6' + 'monthly_'], label=key_sel, color = "red",alpha = 0.7)
#plt.plot(pr_up[key_sel + '1' + 'monthly_'], color = "blue",alpha = 0.4)
plt.plot(pr_up[key_sel + '2' + 'monthly_'], color = "blue",alpha = 0.4)
#plt.plot(pr_up[key_sel + '3' + 'monthly_'], color = "blue",alpha = 0.5)
#plt.plot(pr_up[key_sel + '4' + 'monthly_'], color = "blue",alpha = 0.6)
plt.plot(pr_up[key_sel + '6' + 'monthly_'], color = "blue",alpha = 0.8)
plt.axvline(60, color = "green", linestyle = "--")

plt.ylabel("Monsoon rainfall ($mm/day$)")
plt.ylim(0, 13)
plt.xlim(0,360)
plt.xlabel("Months")
plt.savefig(savedir + "bistability_ensemble_" + key_sel + "_60months.pdf", bbox_inches = 'tight')  # Save the figure


#%%
### ENSEMBLE 3 PANEL
#-------------------

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15,5))

# ENSEMBLE - 1 year
key_sel = '427'


axs[0].text(0.04, 0.96, 'A.', horizontalalignment='left', verticalalignment='top', transform=axs[0].transAxes, fontsize=20)
axs[1].text(0.04, 0.96, 'B.', horizontalalignment='left', verticalalignment='top', transform=axs[1].transAxes, fontsize=20)
axs[2].text(0.04, 0.96, 'C.', horizontalalignment='left', verticalalignment='top', transform=axs[2].transAxes, fontsize=20)

axs[0].plot(pr_down[key_sel + '1' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.4)
axs[0].plot(pr_down[key_sel + '2' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.5)
axs[0].plot(pr_down[key_sel + '3' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.6)
axs[0].plot(pr_down[key_sel + '4' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.7)
axs[0].plot(pr_down[key_sel + '6' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.8)
axs[0].plot(pr_up[key_sel + '1' + 'firstyeardaily_'], color = "blue",alpha = 0.4)
axs[0].plot(pr_up[key_sel + '2' + 'firstyeardaily_'], color = "blue",alpha = 0.5)
axs[0].plot(pr_up[key_sel + '3' + 'firstyeardaily_'], color = "blue",alpha = 0.6)
axs[0].plot(pr_up[key_sel + '4' + 'firstyeardaily_'], color = "blue",alpha = 0.7)
axs[0].plot(pr_up[key_sel + '6' + 'firstyeardaily_'], color = "blue",alpha = 0.8)
axs[0].axvline(252, color = "black", linestyle = "--")
axs[0].axvline(84, color = "black", linestyle = "--")
axs[0].plot(pr_up_ssa['0'], label="Standard", linestyle = '--', color = "black")
axs[0].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[0].set_ylim(0, 13)
axs[0].set_xlim(0,365)
axs[0].set_xlabel("Days of the Year")

# ENSEMBLE - 60 months
axs[1].plot(pr_down[key_sel + '1' + 'monthly_'], label=key_sel, color = "red",alpha = 0.4)
axs[1].plot(pr_down[key_sel + '2' + 'monthly_'], label=key_sel, color = "red",alpha = 0.5)
axs[1].plot(pr_down[key_sel + '3' + 'monthly_'], label=key_sel, color = "red",alpha = 0.6)
axs[1].plot(pr_down[key_sel + '4' + 'monthly_'], label=key_sel, color = "red",alpha = 0.7)
axs[1].plot(pr_down[key_sel + '6' + 'monthly_'], label=key_sel, color = "red",alpha = 0.8)
axs[1].plot(pr_up[key_sel + '1' + 'monthly_'], color = "blue",alpha = 0.4)
axs[1].plot(pr_up[key_sel + '2' + 'monthly_'], color = "blue",alpha = 0.5)
axs[1].plot(pr_up[key_sel + '3' + 'monthly_'], color = "blue",alpha = 0.6)
axs[1].plot(pr_up[key_sel + '4' + 'monthly_'], color = "blue",alpha = 0.7)
axs[1].plot(pr_up[key_sel + '6' + 'monthly_'], color = "blue",alpha = 0.8)
axs[1].axvline(12, color = "green", linestyle = "--")
axs[1].set_ylim(0, 13)
axs[1].set_xlim(0,60)
axs[1].set_xlabel("Months")
axs[1].set_yticklabels([])

# ENSEMBLE - 360 months
axs[2].plot(pr_down[key_sel + '1' + 'monthly_'], label=key_sel, color = "red",alpha = 0.4)
axs[2].plot(pr_down[key_sel + '2' + 'monthly_'], label=key_sel, color = "red",alpha = 0.5)
axs[2].plot(pr_down[key_sel + '3' + 'monthly_'], label=key_sel, color = "red",alpha = 0.6)
axs[2].plot(pr_down[key_sel + '4' + 'monthly_'], label=key_sel, color = "red",alpha = 0.7)
axs[2].plot(pr_down[key_sel + '6' + 'monthly_'], label=key_sel, color = "red",alpha = 0.8)
axs[2].plot(pr_up[key_sel + '1' + 'monthly_'], color = "blue",alpha = 0.4)
axs[2].plot(pr_up[key_sel + '2' + 'monthly_'], color = "blue",alpha = 0.5)
axs[2].plot(pr_up[key_sel + '3' + 'monthly_'], color = "blue",alpha = 0.6)
axs[2].plot(pr_up[key_sel + '4' + 'monthly_'], color = "blue",alpha = 0.7)
axs[2].plot(pr_up[key_sel + '6' + 'monthly_'], color = "blue",alpha = 0.8)
axs[2].axvline(60, color = "green", linestyle = "--")
axs[2].set_ylim(0, 13)
axs[2].set_xlim(0,360)
axs[2].set_xlabel("Months")
axs[2].set_yticklabels([])

plt.subplots_adjust(wspace=0.1)  # Adjust this value to bring the subplots closer together

plt.savefig(savedir + "bistability_ensemble_" + key_sel + "_combined.pdf", bbox_inches = 'tight') # Save the figure


#%%
### ENSEMBLE 3 PANEL (SSA)
#-------------------

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15,5))

# ENSEMBLE - 1 year
key_sel = '427'


axs[0].text(0.04, 0.96, 'A.', horizontalalignment='left', verticalalignment='top', transform=axs[0].transAxes, fontsize=20, fontweight='bold')
axs[1].text(0.04, 0.96, 'B.', horizontalalignment='left', verticalalignment='top', transform=axs[1].transAxes, fontsize=20, fontweight='bold')
axs[2].text(0.04, 0.96, 'C.', horizontalalignment='left', verticalalignment='top', transform=axs[2].transAxes, fontsize=20, fontweight='bold')

axs[0].plot(pr_down_ssa[key_sel + '1' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.4)
axs[0].plot(pr_down_ssa[key_sel + '2' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.5)
axs[0].plot(pr_down_ssa[key_sel + '3' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.6)
axs[0].plot(pr_down_ssa[key_sel + '4' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.7)
axs[0].plot(pr_down_ssa[key_sel + '6' + 'firstyeardaily_'], label=key_sel, color = "red",alpha = 0.8)
axs[0].plot(pr_up_ssa[key_sel + '1' + 'firstyeardaily_'], color = "blue",alpha = 0.4)
axs[0].plot(pr_up_ssa[key_sel + '2' + 'firstyeardaily_'], color = "blue",alpha = 0.5)
axs[0].plot(pr_up_ssa[key_sel + '3' + 'firstyeardaily_'], color = "blue",alpha = 0.6)
axs[0].plot(pr_up_ssa[key_sel + '4' + 'firstyeardaily_'], color = "blue",alpha = 0.7)
axs[0].plot(pr_up_ssa[key_sel + '6' + 'firstyeardaily_'], color = "blue",alpha = 0.8)
axs[0].axvline(252, color = "red", linestyle = "-", linewidth = 1.5)
axs[0].axvline(84, color = "blue", linestyle = "-")
axs[0].plot(pr_up_ssa['0'], label="Standard", linestyle = '--', color = "black")
axs[0].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[0].set_ylim(0, 13)
axs[0].set_xlim(0,365)
axs[0].set_xlabel("Days of the Year")

# ENSEMBLE - 60 months
axs[1].plot(pr_down_ssa[key_sel + '1' + 'monthly_'], label=key_sel, color = "red",alpha = 0.4)
axs[1].plot(pr_down_ssa[key_sel + '2' + 'monthly_'], label=key_sel, color = "red",alpha = 0.5)
axs[1].plot(pr_down_ssa[key_sel + '3' + 'monthly_'], label=key_sel, color = "red",alpha = 0.6)
axs[1].plot(pr_down_ssa[key_sel + '4' + 'monthly_'], label=key_sel, color = "red",alpha = 0.7)
axs[1].plot(pr_down_ssa[key_sel + '6' + 'monthly_'], label=key_sel, color = "red",alpha = 0.8)
axs[1].plot(pr_up_ssa[key_sel + '1' + 'monthly_'], color = "blue",alpha = 0.4)
axs[1].plot(pr_up_ssa[key_sel + '2' + 'monthly_'], color = "blue",alpha = 0.5)
axs[1].plot(pr_up_ssa[key_sel + '3' + 'monthly_'], color = "blue",alpha = 0.6)
axs[1].plot(pr_up_ssa[key_sel + '4' + 'monthly_'], color = "blue",alpha = 0.7)
axs[1].plot(pr_up_ssa[key_sel + '6' + 'monthly_'], color = "blue",alpha = 0.8)
axs[1].axvline(12, color = "green", linestyle = "--")
axs[1].set_ylim(0, 13)
axs[1].set_xlim(0,60)
axs[1].set_xlabel("Months")
axs[1].set_yticklabels([])

# ENSEMBLE - 360 months
axs[2].plot(pr_down_ssa[key_sel + '1' + 'monthly_'], label=key_sel, color = "red",alpha = 0.4)
axs[2].plot(pr_down_ssa[key_sel + '2' + 'monthly_'], label=key_sel, color = "red",alpha = 0.5)
axs[2].plot(pr_down_ssa[key_sel + '3' + 'monthly_'], label=key_sel, color = "red",alpha = 0.6)
axs[2].plot(pr_down_ssa[key_sel + '4' + 'monthly_'], label=key_sel, color = "red",alpha = 0.7)
axs[2].plot(pr_down_ssa[key_sel + '6' + 'monthly_'], label=key_sel, color = "red",alpha = 0.8)
axs[2].plot(pr_up_ssa[key_sel + '1' + 'monthly_'], color = "blue",alpha = 0.4)
axs[2].plot(pr_up_ssa[key_sel + '2' + 'monthly_'], color = "blue",alpha = 0.5)
axs[2].plot(pr_up_ssa[key_sel + '3' + 'monthly_'], color = "blue",alpha = 0.6)
axs[2].plot(pr_up_ssa[key_sel + '4' + 'monthly_'], color = "blue",alpha = 0.7)
axs[2].plot(pr_up_ssa[key_sel + '6' + 'monthly_'], color = "blue",alpha = 0.8)
axs[2].axvline(60, color = "green", linestyle = "--")
axs[2].set_ylim(0, 13)
axs[2].set_xlim(0,360)
axs[2].set_xlabel("Months")
axs[2].set_yticklabels([])

plt.subplots_adjust(wspace=0.1)  # Adjust this value to bring the subplots closer together

plt.savefig(savedir + "bistability_ensemble_" + key_sel + "_combined_ssa.pdf", bbox_inches = 'tight')  # Save the figure

# %%
