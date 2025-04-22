# AUTHOR: Anja Katzenberger, anja.katzenberger@pik-potsdam.de

# This code analyses robustness of hysteresis behaviour
# with regard to varying land surface albedo


#%%
### IMPORT MODULES
#-------------------
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1 import make_axes_locatable

save_dir = "C:\\Users\\anjaka\\Nextcloud\\PhD\\03_MonsoonPlanet_Hysteresis\\Figures"


#%%###  LOAD DATA
#-------------------

albedo_list_sel = [10,20,30,34]  

# monthly average of last 20 years
# adapt the directory to connect to the folder albedo_20years_ymonmean
data_dir = 'C:/Users/anjaka/Nextcloud/PhD/03_MonsoonPlanet_Hysteresis/data/sensitivity/'

ds_dict = {}

pr = {}
sw = {}
tsurf = {}
wvp = {}

for value in albedo_list_sel:
    # Load the dataset
    ds_dict[value] = xr.open_dataset(f'{data_dir}MonsoonPlanet_albedo_{value}.nc')
    
    # Extract the variables
    sw[value] = ds_dict[value]['swdn_toa']
    tsurf[value] = ds_dict[value]['t_surf']-273.15
    pr[value] = ds_dict[value]['precip']*86400
    wvp[value] = ds_dict[value]['WVP']



    
#%%
###  PROCESSING DATA
#---------------------------

# Monsoon region (m)
l_m= 10
u_m = 20

# Define the latitude slices
slice_m = slice(l_m, u_m)


# M
sw_m = {}
tsurf_m = {}
pr_m = {}
evap_m = {}
wvp_m = {}
rh_m = {}

sw_m_z = {}
tsurf_m_z = {}
pr_m_z = {}
evap_m_z = {}
wvp_m_z = {}
rh_m_z = {}

sw_m_z_mean = {}
tsurf_m_z_mean = {}
pr_m_z_mean = {}
evap_m_z_mean = {}
wvp_m_z_mean = {}
rh_m_z_mean = {}


for value in albedo_list_sel:
    sw_m[value] = sw[value].sel(lat=slice_m)
    tsurf_m[value] = tsurf[value].sel(lat=slice_m)
    pr_m[value] = pr[value].sel(lat=slice_m)
    wvp_m[value] = wvp[value].sel(lat=slice_m)

    sw_m_z[value] = sw_m[value].mean('lon')
    tsurf_m_z[value] = tsurf_m[value].mean('lon') 
    pr_m_z[value] = pr_m[value].mean('lon')
    wvp_m_z[value] = wvp_m[value].mean('lon')

    weights = np.cos(np.deg2rad(sw[value].lat))
    sw_m_z_mean[value] = sw_m_z[value].weighted(weights).mean(dim=['lat'])
    weights = np.cos(np.deg2rad(tsurf[value].lat))
    tsurf_m_z_mean[value] = tsurf_m_z[value].weighted(weights).mean(dim=['lat'])
    weights = np.cos(np.deg2rad(pr[value].lat))
    pr_m_z_mean[value] = pr_m_z[value].weighted(weights).mean(dim=['lat'])
    weights = np.cos(np.deg2rad(wvp[value].lat))
    wvp_m_z_mean[value] = wvp_m_z[value].weighted(weights).mean(dim=['lat'])
   

# %%
### PLOT: ROBUSTNESS CHECK albedo
##---------------------------
plt.rcParams.update({'font.size': 14})


# Create a 2x2 panel
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

# Define colors for each albedo depth
colors = ["#DC143C", "#FF8000", "#33A1C9", "#191970"]

# Plot 1
for i, albedo in enumerate(albedo_list_sel):
    sw = sw_m_z_mean[albedo]
    wvp = wvp_m_z_mean[albedo]
    axs[1, 0].plot(sw, wvp, label=f'{albedo}', color=colors[i])
    
    # Connect the first and last data points
    axs[1, 0].plot(sw[[0, -1]], wvp[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[1, 0].plot(sw[0], wvp[0], marker=">", color=colors[i])
    axs[1, 0].plot(sw[6], wvp[6], marker="o", color=colors[i])


axs[1, 0].set_xlabel('Solar Radiation ($W/m^2$)')
axs[1, 0].set_ylabel('Water Vapour Content ($kg/m^2$)')

# Plot 2
for i, albedo in enumerate(albedo_list_sel):
    sw = sw_m_z_mean[albedo]
    pr = pr_m_z_mean[albedo]
    axs[0, 0].plot(sw, pr, label=f'{albedo}', color=colors[i])
    
    # Connect the first and last data points
    axs[0, 0].plot(sw[[0, -1]], pr[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[0, 0].plot(sw[0], pr[0], marker=">", color=colors[i])
    axs[0, 0].plot(sw[6], pr[6], marker="o", color=colors[i])


axs[0, 0].set_xlabel('Solar Radiation ($W/m^2$)')
axs[0, 0].set_ylabel('Monsoon Rainfall ($mm/day$)')

# Plot 3
for i, albedo in enumerate(albedo_list_sel):
    tsurf = tsurf_m_z_mean[albedo]
    pr = pr_m_z_mean[albedo]
    axs[0, 1].plot(tsurf, pr, label=f'{albedo}', color=colors[i])
    
    # Connect the first and last data points
    axs[0, 1].plot(tsurf[[0, -1]], pr[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[0, 1].plot(tsurf[0], pr[0], marker=">", color=colors[i])
    axs[0, 1].plot(tsurf[6], pr[6], marker="o", color=colors[i])


axs[0, 1].set_xlabel('Surface Temperature ($°C$)')
axs[0, 1].set_ylabel('Monsoon Rainfall ($mm/day$)')

# Plot 4
for i, albedo in enumerate(albedo_list_sel):
    tsurf = tsurf_m_z_mean[albedo]
    wvp = wvp_m_z_mean[albedo]
    axs[1, 1].plot(tsurf, wvp, label=f'{albedo}', color=colors[i])
    
    # Connect the first and last data points
    axs[1, 1].plot(tsurf[[0, -1]], wvp[[0, -1]], color=colors[i], linestyle='-')
    
    # Mark only the first data point with a marker "."
    axs[1, 1].plot(tsurf[0], wvp[0], marker=">", color=colors[i])
    axs[1, 1].plot(tsurf[6], wvp[6], marker="o", color=colors[i])


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
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 0.5), fancybox=True, shadow=True, ncol=1, title="Albedo (%)",
           frameon=False, fontsize=14, title_fontsize=14)

# Adjust layout to prevent clipping of titles
plt.tight_layout()
fig.subplots_adjust(bottom=0.1, right=0.8)

# Add text labels to the panels
axs[0, 0].text(0.13, 0.95, 'A.', transform=axs[0, 0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
axs[0, 1].text(0.13, 0.95, 'B.', transform=axs[0, 1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
axs[1, 0].text(0.13, 0.95, 'C.', transform=axs[1, 0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
axs[1, 1].text(0.13, 0.95, 'D.', transform=axs[1, 1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')



plt.savefig(save_dir + '/robustness/MonsoonPlanet_Hysteresis_4panels_albedo.pdf', bbox_inches='tight')

# %%
