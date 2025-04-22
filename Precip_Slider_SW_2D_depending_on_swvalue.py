#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### author: anja katzenberger, anja.katzenberger@pik-potsdam.de

# This code creates figure with 5 panels:
# 1: Off-State: Meridional precip distribution
# 2: Off-State: 3D Precip (sphere)
# 3: SW Radiation througout the year
# 4: Monsoon-State: 3D Precip (sphere)
# 5: Monsoon-State: Meridional precip distribution

#%%
### LOAD MODULES
###--------------------------- 
import sys
print(sys.executable)

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, timedelta
import imageio
import os

# coastward boundary of land (°N)
latlow = 10

# Where to save the figures
savedir = r'C:\Users\anjaka\Nextcloud\PhD\03_MonsoonPlanet_Hysteresis\Figures\precip_slider\\'

#%%
### LOAD AND PROCESS DATA
###---------------------------

# Monthly precip data for 1 year (averaged over 20 years)
# adapt to your path to the folder slab_20years_ymonmean
dir_10N = r"C:\Users\anjaka\Nextcloud\PhD\03_MonsoonPlanet_Hysteresis\data_publication\slab_20years_ymonmean\MonsoonPlanet_slab_50.nc"
ds_10N = xr.open_dataset(dir_10N)
precip = ds_10N['precip']
precip_z = precip.mean('lon')

# Daily precip data for one year
# adapt to your path to the folder randomyear_precip_tsurf
dir_day = r"C:\Users\anjaka\Nextcloud\PhD\03_MonsoonPlanet_Hysteresis\data_publication\random_year\01010101.atmos_day_selyear_precip_tsurf.nc"
ds_day = xr.open_dataset(dir_day)
precip_day = ds_day['precip']
precip_day_z = precip_day.mean('lon')

# Daily SW data for one year 
# adapt to your path to the folder randomyear_sw
dir_day_sw = r"C:\Users\anjaka\Nextcloud\PhD\03_MonsoonPlanet_Hysteresis\data_publication\random_year\01010101.atmos_day_selyear_swdntoa.nc"
ds = xr.open_dataset(dir_day_sw)
sw_day = ds['swdn_toa']
sw_z = sw_day.mean('lon')
sw_day_z = sw_day.mean('lon')
sw_day_10_30 = sw_day_z.sel(lat = slice(10,30))
weights = np.cos(np.deg2rad(sw_day_10_30.lat))
sw = sw_day_10_30.weighted(weights).mean(dim=['lat'])


#%%
### TIME FORMAT
###--------------------------

# Function that transforms date to day of year
def date_to_dayofyear(month, day):
    try:
        # Get the current year
        current_year = datetime.now().year
        
        # Create a datetime object for the given date
        date_obj = datetime(current_year, month, day)

        # Get the day of the year (1-based index)
        day_of_year = date_obj.timetuple().tm_yday
        return day_of_year
    except ValueError as e:
        return f"Invalid date format: {e}"

# Example usage for 21. June
month = 6  # 6 represents June
day = 21
base = date_to_dayofyear(6, 21)


# Function that transforms day of year to date
def dayofyear_to_date(day_of_year):
    try:
        # Get the current year
        current_year = datetime.now().year
        
        # Create a datetime object for the given day of the year
        date_obj = datetime(current_year, 1, 1) + timedelta(days=day_of_year - 1)

        # Format the date as "day. Month"
        formatted_date = date_obj.strftime("%d. %B")
        return formatted_date
    except ValueError as e:
        return f"Invalid day of the year: {e}"

# Example usage for day of the year 172
dayofyear_to_date(172)


# Function that transforms day of year to months (e.g. 15 th May = 5.5)
def dayofyear_to_month(day_of_year):
    try:
        # Get the current year
        current_year = datetime.now().year
        
        # Define a list of the number of days in each month
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        # Initialize variables
        month = 1
        day = day_of_year

        # Find the month for the given day_of_year
        for days in days_in_month:
            if day <= days:
                break
            day -= days
            month += 1

        # Calculate the value as months + day/30
        value = month + day / 30.0 - 1

        return value
    except ValueError as e:
        return f"Invalid day of the year: {e}"

# Example usage for 21. June
day_of_year = 172
value = dayofyear_to_month(day_of_year)
print(value)  # Output: 5.7



#%%
### CREATE SLIDER
#---------------------------

sw_max = sw.argmax()

# Function that returns the index of the closest sw_values for a specific radiation
# for each radition value, there are two such values
def index_l_for_rad(rad):
    check_l = []
    for i in range(0,sw_max.values):
        if int(sw.values[i]) > rad: 
            check_l.append(i)
    val_l = check_l[0]
    return val_l

def index_r_for_rad(rad):
    check_r = []
    for i in range(sw_max.values,365):
        if int(sw.values[i]) > rad: 
            check_r.append(i)
    val_r = check_r[len(check_r)-1]
    return val_r
    



def update_plot():
    m_value = int(slider.get())  # Get the value of m from the slider
    plt.clf()  # Clear the previous plot
    plt.plot(sw)
    plt.xlabel("Day of the year")
    plt.ylabel("SW Radiation [$W/m^2$]")
    plt.axvline(base, color="red")
    plt.axvline(base + m_value)
    plt.axvline(base - m_value)
    plt.show()
    



#%%
### PANELS AS INDIVIDUAL FIGURES
###---------------------------

### Panel 1 & 5: Meridional precip distribution
m=4
fig = plt.figure(figsize=(15,5))
ax2 = fig.add_subplot(121)
pos2 = ax2.get_position()
ax2.set_xlim([0,40])
ax2.set_xticks((0,15,30))
ax2.set_yticks([])
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')
ax2.tick_params(axis='x', colors='#015482')
ax2 = plt.plot(precip_day_z.sel(time=precip_day_z['time.dayofyear'] == 266).values.flatten()*86400,precip_day['lat'], color = '#015482')
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')

#%%

### Panel 2 & 4: 3D Precip (sphere)
m = 4
globe_proj = ccrs.Orthographic(central_longitude=80, central_latitude=0)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(2, 1, 1, projection=globe_proj)
im = ax.contourf(precip_day.lon, precip_day.lat, precip_day[m,:,:]*86400, transform=ccrs.PlateCarree(), cmap='Blues', levels = range(0,21), extend = "max")
ax.plot(np.linspace(-180, 180, 1000), np.zeros(1000), color='black', transform=ccrs.PlateCarree(), linewidth = 1.2, linestyle = "--")
ax.plot(np.linspace(-180, 180, 1000), np.full(1000, 10), color='black', transform=ccrs.PlateCarree(), linewidth = 1.2)
ax.plot(np.linspace(-180, 180, 1000), np.full(1000, 60), color='black', transform=ccrs.PlateCarree(), linewidth = 1.2)
ax.set_global()
plt.title(dayofyear_to_date(base-m))

#%%
### Panel 3 - sw
fig, ax = plt.subplots(1,1)
plt.plot(sw, color = "orange")
ax.set_ylabel("Radiation ($W/m^2$")
ax.set_xlabel("Days of the Year")
ax.set_ylim([290,460])
ax.axvline(sw.argmax(),color = "black",linestyle = '--')
ax.set_title("Radiation =" + str(int(np.round(sw.isel(time=base + m).values)))+ "$W/m^2$", color = 'orange',fontsize = 20, y=1.09)

#%%
### Panel 3 - alternative
fig,ax = plt.subplots(figsize = (4,4))
plt.contourf(sw_z['time.month'],sw_z['lat'],sw_z.transpose(),cmap="Oranges",extend="both",levels = range(160,560,10))
plt.ylabel("lat")
plt.xlabel("Months")
plt.ylim([-60,60])
plt.axhline(0, color = "black", linestyle = '--')
plt.axhline(latlow, color = "black")
plt.axhline(latlow + 50, color = "black")
plt.axvline(dayofyear_to_month(base),color = "orange",linestyle = '--')
plt.colorbar(label = "SW down TOA (W$/m^2$)",orientation="horizontal",pad = 0.2)

#%%
### FINAL FIGURE
###---------------------------

globe_proj = ccrs.Orthographic(central_longitude=80, central_latitude=0)
base = date_to_dayofyear(6,21)
size = 19

#for rad in range(int(max(sw).values)-1,310,-1): # until left line reaches january
for rad in range(426,425,-1): # example for one value (line above gives all required values)
    print(rad)
    
    fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1,5,figsize = (30,5))
    
    ax1.remove()
    ax1 = fig.add_subplot(151)
    ax1.set_xlim([0,15])
    ax1.set_xticks((0,5,10,15))
    ax1.set_yticks([])
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top')
    ax1.tick_params(axis='x', colors='#015482')
    ax1.set_xticklabels(ax1.get_xticks(), fontsize=size)
    ax1.plot(precip_day_z.sel(time=precip_day_z['time.dayofyear'] == index_l_for_rad(rad)).values.flatten()*86400,precip_day['lat'], color = '#015482')
    ax1.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    ax1.axhline(latlow, linewidth = 0.7, color = 'black')
    ax1.axhline(latlow + 50, linewidth = 0.7, color = 'black')
    ax1.invert_xaxis()
    
    ax2.remove()
    ax2 = fig.add_subplot(152, projection=globe_proj) # 1x5 grid, subplot nr 2
    ax2fig = ax2.contourf(precip_day.lon, precip_day.lat, precip_day[index_l_for_rad(rad),:,:]*86400, transform=ccrs.PlateCarree(), cmap='Blues', levels = range(0,21), extend = "max")
    ax2.plot(np.linspace(-180, 180, 1000), np.zeros(1000), color='black', transform=ccrs.PlateCarree(), linewidth = 1.2, linestyle = "--")
    ax2.plot(np.linspace(-180, 180, 1000), np.full(1000, latlow), color='black', transform=ccrs.PlateCarree(), linewidth = 1.2)
    ax2.plot(np.linspace(-180, 180, 1000), np.full(1000, latlow+50), color='black', transform=ccrs.PlateCarree(), linewidth = 1.2)
    ax2.set_global()
    plt.title(dayofyear_to_date(index_l_for_rad(rad)),fontsize = size)
    
    ax3fig = ax3.plot(sw, color = "orange")
    ax3.set_ylabel("Radiation 10-30N ($W/m^2$)", fontsize=size)
    ax3.set_xlabel("Days of the Year", fontsize=size)
    ax3.set_ylim([290,460])
    ax3.axvline(index_l_for_rad(rad),color = "black")
    ax3.axvline(index_r_for_rad(rad),color = "black")
    ax3.axvline(sw.argmax(),color = "black",linestyle = '--')
    ax3.set_title("Radiation =" + str(rad)+ "$W/m^2$", color = 'orange',fontsize = 20, y=1.09)
    title_line1 = "Upper atmospheric boundary condition:"
    ax3.text(-0.2, 1.3, title_line1, color='orange', fontsize=20, transform=ax3.transAxes)
    ax3.set_xticks((0,100,200,300))
    ax3.set_yticks((300,350,400,450))
    ax3.set_xticklabels(ax3.get_xticks(), fontsize=size)
    ax3.set_yticklabels(ax3.get_yticks(), fontsize=size)
    
    precip_colorbar_ax = plt.gcf().add_axes([0.21, 0.0, 0.13, 0.03])  # Adjust the position as needed [left, bottom, width, height]

    ax4.remove()
    ax4 = fig.add_subplot(154, projection=globe_proj)
    ax4fig = ax4.contourf(precip_day.lon, precip_day.lat, precip_day[index_r_for_rad(rad),:,:]*86400, transform=ccrs.PlateCarree(), cmap='Blues', levels = range(0,21), extend = "max")
    precip_colorbar = plt.colorbar(ax4fig, cax=precip_colorbar_ax, orientation='horizontal', label='Precipitation (mm/day)')
    precip_colorbar.ax.xaxis.label.set_size(size)  # Adjust the font size (14) as needed
    precip_colorbar.ax.tick_params(axis='x', labelsize=size)
    ax4.plot(np.linspace(-180, 180, 1000), np.zeros(1000), color='black', transform=ccrs.PlateCarree(), linewidth = 1.2, linestyle = "--")
    ax4.plot(np.linspace(-180, 180, 1000), np.full(1000, latlow), color='black', transform=ccrs.PlateCarree(), linewidth = 1.2)
    ax4.plot(np.linspace(-180, 180, 1000), np.full(1000, latlow + 50), color='black', transform=ccrs.PlateCarree(), linewidth = 1.2)
    ax4.set_global()
    plt.title(dayofyear_to_date(index_r_for_rad(rad)),fontsize = size)
    
    ax5.remove()
    ax5 = fig.add_subplot(155)
    ax5.set_xlim([0,15])
    ax5.set_xticks((0,5,10,15))
    ax5.set_yticks([])
    ax5.xaxis.tick_top()
    ax5.xaxis.set_label_position('top')
    ax5.tick_params(axis='x', colors='#015482')
    ax5.set_xticklabels(ax1.get_xticks(), fontsize=size)
    ax5.plot(precip_day_z.sel(time=precip_day_z['time.dayofyear'] == index_r_for_rad(rad)).values.flatten()*86400,precip_day['lat'], color = '#015482')
    ax5.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    ax5.axhline(latlow, linewidth = 0.7, color = 'black')
    ax5.axhline(latlow + 50, linewidth = 0.7, color = 'black')
    
    # Left upper line
    start_x_lu = index_l_for_rad(rad)  # Adjust the x-coordinate of the starting point as needed
    end_x_lu = -300  # Adjust the x-coordinate of the ending point as needed
    start_y_lu = 460  # Adjust the y-coordinate as needed
    end_y_lu = 497
    ax3.annotate('',xy=(end_x_lu, end_y_lu), xytext=(start_x_lu, start_y_lu),arrowprops=dict(arrowstyle='-', color='black', lw=1.5),annotation_clip=False,xycoords='data',textcoords='data')
    plt.draw()
    
    # Right upper line
    start_x = index_r_for_rad(rad)  
    end_x_ru = 650  
    start_y = start_y_lu  
    end_y = end_y_lu
    ax3.annotate('',xy=(end_x_ru, end_y), xytext=(start_x, start_y),arrowprops=dict(arrowstyle='-', color='black', lw=1.5),annotation_clip=False,xycoords='data',textcoords='data')
    plt.draw()
    
    # Right lower line
    start_x = index_r_for_rad(rad) 
    end_x = end_x_ru  
    start_y_rl = 290  
    end_y_rl = 254
    ax3.annotate('',xy=(end_x, end_y_rl), xytext=(start_x, start_y_rl),arrowprops=dict(arrowstyle='-', color='black', lw=1.5),annotation_clip=False,xycoords='data',textcoords='data')
    plt.draw()
    
    # Left lower line
    start_x = index_l_for_rad(rad)  
    end_x = end_x_lu 
    start_y = start_y_rl  
    end_y = end_y_rl
    ax3.annotate('',xy=(end_x, end_y), xytext=(start_x, start_y),arrowprops=dict(arrowstyle='-', color='black', lw=1.5),annotation_clip=False,xycoords='data',textcoords='data')
    plt.draw()
    
    # Adjust the position of ax1 and ax5
    ax1.set_position([0.21, 0.1, 0.06, 0.8])  # Adjust the [left, bottom, width, height]
    ax5.set_position([0.75, 0.1, 0.06, 0.8])  # Adjust the [left, bottom, width, height]
    ax3.set_position([0.45, 0.235,0.13, 0.53])
    plt.savefig(savedir + 'Precip_Flexible_' + str(rad) + '.png',bbox_inches = 'tight')
    
#%%-------------------------
# Create gif
#--------------------------

# Set the directory containing the PNG files
directory = savedir

# Get a list of all the PNG files in the directory
files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.png')]

# Sort the files in ascending order
files = sorted(files, key=lambda x: int(x.split('_')[11].split('.')[0]))
files_r = list(reversed(files))

# Create an imageio writer object
with imageio.get_writer(savedir + 'precip_flexible_sw.gif', mode='I') as writer:
    # Loop through the PNG files and add them to the writer object
    for filename in files_r:
        image = imageio.imread(filename)
        writer.append_data(image)
        
       