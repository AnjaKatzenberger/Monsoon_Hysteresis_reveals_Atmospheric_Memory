# Author: Anja Katzenberger, anja.katzenberger@pik-potsdam.de

# This code creates figures of the memory effect
# It is used to quantify the duration of the memory effect

#%%
### LOAD PACKAGES
#-------------------

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
from pyts.decomposition import SingularSpectrumAnalysis
import matplotlib.patches as mpatches

# Where to save the figures
savedir = r'C:\Users\anjaka\Nextcloud\PhD\03_MonsoonPlanet_Hysteresis\Figures\memory\\'

# Directories 
dir = r"C:\Users\anjaka\Nextcloud\PhD\03_MonsoonPlanet_Hysteresis\data\memory"
files = os.listdir(dir)

#%%
### LOAD DATA
#-------------------

# Get list of keys for memory simulations
keys = [file.replace("output_", "").replace("memory_selection.nc", "") for file in files if file.startswith("output_") and file.endswith("memory_selection.nc")]

sw_dict = {}
pr_dict = {}
wvp_dict = {}
for key in keys:
    data = xr.open_dataset(dir + r"\output_" + key + "memory_selection.nc")
    sw_dict[key] = data["swdn_toa"]
    pr_dict[key] = data["precip"]*86400
    wvp_dict[key] = data["WVP"]


#%%
### PROCESSING DATA
#-------------------
    
# Monsoon region
lat_low = 10
lat_up = 20

pr = {}
sw = {}
wvp = {}
for key in keys:
    pr_lat = pr_dict[key].mean('lon').isel(time=slice(0, 500))
    pr_lat_mr = pr_lat.sel(lat=slice(lat_low,lat_up))
    weights = np.cos(np.deg2rad(pr_lat_mr.lat))
    pr[key] = pr_lat_mr.weighted(weights).mean(dim=['lat'])

    sw_lat = sw_dict[key].mean('lon').isel(time=slice(0, 500))
    sw_lat_mr = sw_lat.sel(lat=slice(lat_low,lat_up))
    weights = np.cos(np.deg2rad(sw_lat_mr.lat))
    sw[key] = sw_lat_mr.weighted(weights).mean(dim=['lat'])

    wvp_lat = wvp_dict[key].mean('lon').isel(time=slice(0, 500))
    wvp_lat_mr = wvp_lat.sel(lat=slice(lat_low,lat_up))
    weights = np.cos(np.deg2rad(wvp_lat_mr.lat))
    wvp[key] = wvp_lat_mr.weighted(weights).mean(dim=['lat'])



#%%
### SINGULAR SPECTRUM ANALYSIS
#-------------------

# Smoothening time series by using Singular Spectrum Analysis    
L = 15 # window size

pr_ssa = {}
for key in keys:
    F = pr[key]
    F_arr = np.array(F)
    F_in = np.array([F_arr])
    ssa = SingularSpectrumAnalysis(window_size = L)
    X_ssa = ssa.transform(F_in)
    pr_ssa[key] = X_ssa[0, 0, :]

wvp_ssa = {}
for key in keys: 
    F = wvp[key]
    F_arr = np.array(F)
    F_in = np.array([F_arr])
    ssa = SingularSpectrumAnalysis(window_size = L)
    X_ssa = ssa.transform(F_in)
    wvp_ssa[key] = X_ssa[0, 0, :]





# %%
### PLOT: MEMORY FOR AUGUST
#-------------------
    
key_sel = '0801_'
threshold = 213 # day of intervention

fig, axs = plt.subplots(3, 1, figsize=(10, 10)) 

# PLOT 1: SW 
axs[0].set_title("Month:" + key_sel[0]+key_sel[1]+ "; Day:" + key_sel[2]+key_sel[3])
axs[0].plot(sw[key_sel], color = "red")
axs[0].plot(sw['0_'], linestyle = '--', color = "black")
axs[0].set_ylabel("Insolation ($W/m^2$)")
axs[0].set_xlim(0, 350)
axs[0].axvline(threshold, color = "black", linestyle = "-")

# PLOT 2: PR D
axs[1].plot(pr[key_sel],  color = "red",label = "Memory = 19 days")
axs[1].plot(pr_ssa['0_'],  linestyle = '--', color = "black")
axs[1].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[1].set_xlim(0, 350)	
axs[1].axvline(threshold, color = "black", linestyle = "-")
#axs[1].axhline(0.2, color = "green", linestyle = "--")
axs[1].legend()

# PLOT 3: WVP
axs[2].plot(wvp[key_sel], color = "red")
axs[2].plot(wvp_ssa['0_'], linestyle = '--', color = "black")
axs[2].set_ylabel("WVP ($kg/m^2$)")
axs[2].set_xlim(0, 350)
axs[2].axvline(threshold, color = "black", linestyle = "-")
axs[2].axhline(35, color = "black", linestyle = "--")
axs[2].axvline(270, color = "blue", linestyle = "--")
axs[2].axvspan(213, 270, alpha=0.1, color='gray')
axs[1].axvspan(213, 232, alpha=0.1, color='gray')

axs[2].axvline(138, color = "gray", linestyle = ":")
axs[1].axvline(138, color = "gray", linestyle = ":")
axs[2].axvline(232, color = "black", linestyle = "--")
axs[1].axvline(232, color = "black", linestyle = "--")
axs[2].set_xlabel("Days of the Year")
axs[2].legend()

# Create a custom legend handle with gray color
gray_patch1 = mpatches.Patch(color='gray', label='Rainfall Memory \n = 19 days', alpha = 0.4)
gray_patch2 = mpatches.Patch(color='gray', label='WVP Memory \n = 57 days', alpha = 0.4)

# Add the custom handle to the legend
axs[2].legend(handles=[gray_patch2], loc = 6, facecolor = 'white',  framealpha = 1, edgecolor = 'white')
axs[1].legend(handles=[gray_patch1], loc = 6, facecolor = 'white',  framealpha = 1, edgecolor = 'white')

# Add text labels to the panels
axs[0].text(0.02, 0.95, 'A.', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
axs[1].text(0.02, 0.95, 'B.', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
axs[2].text(0.02, 0.95, 'C.', transform=axs[2].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


plt.savefig(savedir + "memory_effect_wvp_"+ key_sel + ".pdf")

# %%
### PLOT: MEMORY FOR SEPTEMBER
#-------------------

key_sel = '0901_'
threshold = 244

fig, axs = plt.subplots(3, 1, figsize=(10, 10))  
# PLOT 1: SW 
axs[0].set_title("Month:" + key_sel[0]+key_sel[1]+ "; Day:" + key_sel[2]+key_sel[3])
axs[0].plot(sw[key_sel], color = "red")
axs[0].plot(sw['0_'], linestyle = '--', color = "black")
axs[0].set_ylabel("Insolation ($W/m^2$)")
axs[0].set_xlim(0, 350)
axs[0].axvline(threshold, color = "black", linestyle = "-")

# PLOT 2: PR D
axs[1].plot(pr[key_sel],  color = "red",label = "Memory = 20 days")
axs[1].plot(pr_ssa['0_'],  linestyle = '--', color = "black")
axs[1].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[1].set_xlim(0, 350)	
axs[1].axvline(threshold, color = "black", linestyle = "-")
axs[1].legend()

# PLOT 3: WVP
axs[2].plot(wvp[key_sel], color = "red", label = "Memory = 58 days")
axs[2].plot(wvp_ssa['0_'], linestyle = '--', color = "black")
axs[2].set_ylabel("WVP ($kg/m^2$)")
axs[2].set_xlim(0, 350)
axs[2].axvline(threshold, color = "black", linestyle = "-")
axs[2].set_xlabel("Days of the Year")
axs[2].legend()
axs[2].axhline(35, color = "black", linestyle = "--")
axs[2].axvline(302, color = "blue", linestyle = "--")
axs[2].axvline(138, color = "gray", linestyle = ":")
axs[1].axvline(138, color = "gray", linestyle = ":")

axs[2].axvline(264, color = "black", linestyle = "--")
axs[1].axvline(264, color = "black", linestyle = "--")

axs[1].axvspan(244, 264, alpha=0.1, color='gray')
axs[2].axvspan(244, 302, alpha=0.1, color='gray')

# Create a custom legend handle with gray color
gray_patch1 = mpatches.Patch(color='gray', label='Rainfall Memory \n = 20 days', alpha = 0.4)
gray_patch2 = mpatches.Patch(color='gray', label='WVP Memory \n = 58 days', alpha = 0.4)

# Add the custom handle to the legend
axs[2].legend(handles=[gray_patch2], loc = 6, facecolor = 'white',  framealpha = 1, edgecolor = 'white')
axs[1].legend(handles=[gray_patch1], loc = 6, facecolor = 'white',  framealpha = 1, edgecolor = 'white')

# Add text labels to the panels
axs[0].text(0.02, 0.95, 'A.', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
axs[1].text(0.02, 0.95, 'B.', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
axs[2].text(0.02, 0.95, 'C.', transform=axs[2].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')



plt.savefig(savedir + "memory_effect_wvp_"+ key_sel + ".pdf")

# %%
### PLOT: MEMORY FOR OCTOBER
#-------------------

key_sel = '1001_'
threshold = 274

fig, axs = plt.subplots(3, 1, figsize=(10, 10))  
# PLOT 1: SW 
axs[0].set_title("Month:" + key_sel[0]+key_sel[1]+ "; Day:" + key_sel[2]+key_sel[3])
axs[0].plot(sw[key_sel], color = "red")
axs[0].plot(sw['0_'], linestyle = '--', color = "black")
axs[0].set_ylabel("Insolation ($W/m^2$)")
axs[0].set_xlim(0, 350)
axs[0].axvline(threshold, color = "black", linestyle = "-")

# PLOT 2: PR D
axs[1].plot(pr[key_sel],  color = "red",label = "Memory = 15 days")
axs[1].plot(pr_ssa['0_'],  linestyle = '--', color = "black")
axs[1].set_ylabel("Monsoon rainfall ($mm/day$)")
axs[1].set_xlim(0, 350)	
axs[1].axvline(threshold, color = "black", linestyle = "-")
axs[1].legend()

# PLOT 3: WVP
axs[2].plot(wvp[key_sel], color = "red", label = "Memory = 66 days")
axs[2].plot(wvp_ssa['0_'], linestyle = '--', color = "black")
axs[2].set_ylabel("WVP ($kg/m^2$)")
axs[2].set_xlim(0, 350)
axs[2].axvline(threshold, color = "black", linestyle = "-")
axs[2].axvline(138, color = "gray", linestyle = ":")
axs[1].axvline(138, color = "gray", linestyle = ":")

axs[2].axvline(289, color = "black", linestyle = "--")
axs[1].axvline(289, color = "black", linestyle = "--")

axs[2].axvline(340, color = "blue", linestyle = "--")

axs[1].axvspan(274, 289, alpha=0.1, color='gray')
axs[2].axvspan(274, 340, alpha=0.1, color='gray')

axs[2].set_xlabel("Days of the Year")
axs[2].legend()
axs[2].axhline(35, color = "black", linestyle = "--")

# Create a custom legend handle with gray color
gray_patch1 = mpatches.Patch(color='gray', label='Rainfall Memory \n = 15 days', alpha = 0.4)
gray_patch2 = mpatches.Patch(color='gray', label='WVP Memory \n = 66 days', alpha = 0.4)

# Add the custom handle to the legend
axs[2].legend(handles=[gray_patch2], loc = 6, facecolor = 'white',  framealpha = 1, edgecolor = 'white')
axs[1].legend(handles=[gray_patch1], loc = 6, facecolor = 'white',  framealpha = 1, edgecolor = 'white')

# Add text labels to the panels
axs[0].text(0.02, 0.95, 'A.', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
axs[1].text(0.02, 0.95, 'B.', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
axs[2].text(0.02, 0.95, 'C.', transform=axs[2].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


plt.savefig(savedir + "memory_effect_wvp_"+ key_sel + ".pdf")

# %%
### PLOT: WVP THRESHOLD AUGUST
#-------------------

key_sel = '0801_'

plt.figure(figsize=(5, 5))
plt.title("Month:" + key_sel[0]+key_sel[1]+ "; Day:" + key_sel[2]+key_sel[3])

plt.scatter(wvp[key_sel],pr[key_sel], s= 10)

plt.xlabel("WVP ($kg/m^2$)")
plt.ylabel("Monsoon rainfall ($mm/day$)")
plt.axvline(35, color = "black", linestyle = "--")
plt.savefig(savedir + "pr_wvp_plot_"+ key_sel + ".pdf")


# %%
### PLOT: WVP THRESHOLD SEPTEMBER
#-------------------

key_sel = '0901_'

plt.figure(figsize=(5, 5))
plt.title("Month:" + key_sel[0]+key_sel[1]+ "; Day:" + key_sel[2]+key_sel[3])

plt.scatter(wvp[key_sel],pr[key_sel], s= 10)
plt.xlabel("WVP ($kg/m^2$)")
plt.ylabel("Monsoon rainfall ($mm/day$)")
plt.axvline(35, color = "black", linestyle = "--")
plt.savefig(savedir + "pr_wvp_plot_"+ key_sel + ".pdf")

# %%
### PLOT: WVP THRESHOLD OCTOBER
#-------------------

key_sel = '1001_'


plt.figure(figsize=(5, 5))
plt.title("Month:" + key_sel[0]+key_sel[1]+ "; Day:" + key_sel[2]+key_sel[3])

plt.scatter(wvp[key_sel],pr[key_sel], s= 10)
plt.xlabel("WVP ($kg/m^2$)")
plt.ylabel("Monsoon rainfall ($mm/day$)")
plt.axvline(35, color = "black", linestyle = "--")
plt.savefig(savedir + "pr_wvp_plot_"+ key_sel + ".pdf")



#%%
### PLOT: WVP THRESHOLD (3 PANEL in 1 FIG)
#-------------------

fig, axs = plt.subplots(3, 1, figsize=(5, 15))  # Create a figure with 3 subplots

keys = ['0801_', '0901_', '1001_']

for i, key_sel in enumerate(keys):
    axs[i].scatter(wvp[key_sel],pr[key_sel], s= 10)
    axs[i].set_title("Month:" + key_sel[0]+key_sel[1]+ "; Day:" + key_sel[2]+key_sel[3])
    if i == 2:
        axs[i].set_xlabel("WVP ($kg/m^2$)")
    if i != 2:
        axs[i].set_xticklabels([])
    axs[i].set_ylabel("Monsoon rainfall ($mm/day$)")
    axs[i].set_xlim(0, 80)
    axs[i].set_ylim(0, 12)
    axs[i].axvline(35, color = "black", linestyle = "--")

plt.tight_layout()  # Adjust the layout so that the plots do not overlap
plt.savefig(savedir + "pr_wvp_plot_combined.pdf")
# %%
