# Author: Anja Katzenberger, anja.katzenberger@pik-potsdam.de

# This code produces figures of wind direction and precipitation distribution on the Monsoon Planet

#%%
### IMPORT MODULES
#-------------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1 import make_axes_locatable

savedir = "C:\\Users\\anjaka\\Nextcloud\\PhD\\03_MonsoonPlanet_Hysteresis\\Figures"


#%%
###  LOAD DATA
#-------------------

slab_list = [50,100,200,500]  # available slabs
slab_list_sel = [50,100,200,500] # selected slabs

# monthly average of last 20 years
# adapt the path to the folder slab_20years_ymonmean
data_dir = 'C:/Users/anjaka/Nextcloud/PhD/01_Monsoon_Planet/A_data/'

# Initialize dictionaries to store the data for each slab
ds_dict = {}

pr = {}
sw = {}
tsurf = {}
evap = {}
wvp = {}
rh = {}
u={}
v={}


for slab in slab_list:
    # Load the dataset
    ds_dict[slab] = xr.open_dataset(f'{data_dir}MonsoonPlanet_slab_{slab}.nc')

    # Extract the variables
    u[slab] = ds_dict[slab]['ucomp']
    v[slab] = ds_dict[slab]['vcomp']

    sw[slab] = ds_dict[slab]['swdn_toa']
    tsurf[slab] = ds_dict[slab]['t_surf']-273.15
    pr[slab] = ds_dict[slab]['precip']*86400
    evap[slab] = ds_dict[slab]['evap']*86400
    wvp[slab] = ds_dict[slab]['WVP']
    rh[slab] = ds_dict[slab]['rh']

#%%
# CALCULATING WIND SPEED
#-------------------    

wind_speed = {}
for slab in slab_list:
    wind_speed[slab] = np.sqrt(u[slab].isel(pfull=-1)**2 + v[slab].isel(pfull=-1)**2)

print(ds_dict[50]['pfull'])





#%%
### PLOT WIND DIRECTIONS
#-------------------

plt.rcParams['font.size'] = 16
month_names = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

# Create a 3x4 panel for all 12 months
fig, axs = plt.subplots(4, 3, figsize=(15, 20))

# Select slab depth
slab = 200

# Assuming that your data includes longitude and latitude coordinates
lon = ds_dict[slab].lon
lat = ds_dict[slab].lat

# Calculate wind components at highest pressure level
u_surf = u[slab].isel(pfull=-1)
v_surf = v[slab].isel(pfull=-1)

# Only every xth wind vector is shown
stride_x = 10
stride_y = 3

# Loop over all 12 months
for month in range(12):
    # Calculate row and column index for the subplot
    i, j = divmod(month, 3)

    # Create a filled contour plot of the wind speed
    cs = axs[i, j].contourf(lon, lat, pr[slab][month,:,:], extend='max', levels=range(0, 15, 1), cmap='Blues')

    # Add wind vectors
    Q = axs[i, j].quiver(lon[::stride_x], lat[::stride_y], u_surf[month, ::stride_y, ::stride_x], v_surf[month, ::stride_y, ::stride_x], color='black')

    # Add horizontal lines
    axs[i, j].axhline(0, color='black', linestyle='--')
    axs[i, j].axhline(10, color='black')
    axs[i, j].axhline(60, color='black')
    axs[i, j].set_ylim(-20, 25)


    # Set the title
    axs[i, j].set_title(month_names[month])

        # Add x label to the bottom row
    if i == 3: 
        axs[i, j].set_xlabel('lon')

    # Add y label to the first column
    if j == 0:
        axs[i, j].set_ylabel('lat')

cbar_ax = fig.add_axes([0.1, -0.07, 0.8, 0.03])  # Adjust as needed
cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Precipitation (mm/day)')

plt.tight_layout()
plt.savefig(savedir + '/wind_slab' + str(slab)+'_pr.pdf', bbox_inches='tight')

# %%
