# Author: Anja Katzenberger, anja.katzenberger@pik-potsdam.de

# This code applies CDOs for preprocessing the model output


ens = "6" ### ensemble member
updown = "down"


#-------------------
### LOAD PACKAGES
#-------------------

import os
import subprocess
import numpy as np
import sys

#-------------------
### DIRECTORIES
#-------------------
 
dir = os.getcwd()
files = os.listdir(dir)

if updown == "down":
    filename = "009" + ens + "0911.atmos_day.nc"
    print(filename)
    if os.path.exists(filename):
        print("File exists")
    else:
        print("File does not exist!!!!")
        sys.exit()


if updown == "up":
    filename = "009" + ens + "0326.atmos_day.nc"
    print(filename)
    if os.path.exists(filename):
        print("File exists")
    else:
        print("File does not exist!!!")
        sys.exit()


#-------------------
### DAILY FILES
#-------------------
# get a list of daily files 
day_files = [file for file in files if file.endswith('atmos_day.nc')]

print("Daily files: ", day_files)

# CDO1: selname (select variables)
day_files_sel = []
for file in day_files:
    output = file.split('.')[0] + file.split('.')[1] + "_sel.nc"
    cdo_cmd = "cdo selname,swdn_toa,precip " + file + " " + output
    print(cdo_cmd)
    subprocess.check_call(cdo_cmd, shell=True)
    day_files_sel.append(output)

# CDO2: mergetime (merge files)
output = "output_427" + ens + "firstyeardaily_" + updown +  "_sel.nc"
print(cdo_cmd)
cdo_cmd = "cdo mergetime " + " ".join(day_files_sel) + " " + output 
subprocess.check_call(cdo_cmd, shell=True)

# CDO3: monmean (monthly mean)
output_monmean = "output_427" + ens + "firstyeardailymonmean_" + updown +  "_sel.nc"
print(cdo_cmd)
cdo_cmd = "cdo monmean " + output + " " + output_monmean
subprocess.check_call(cdo_cmd, shell=True)


#-------------------
### MONTHLY FILES
#-------------------

# get a list of monthly files (monmean)
month_files = [file for file in files if file.endswith('atmos_month.nc')]

print("Monthly files: ", month_files)

# CDO1: selname (select variables)
month_files_sel = []
for file in month_files:
    output = file.split('.')[0] + file.split('.')[1] + "_sel.nc"
    cdo_cmd = "cdo selname,swdn_toa,precip, " + file + " " + output
    print(cdo_cmd)
    subprocess.check_call(cdo_cmd, shell=True)
    month_files_sel.append(output)

# CDO2: mergetime (merge files)
print(cdo_cmd)
cdo_cmd = "cdo mergetime " + output_monmean + " " + " ".join(month_files_sel) + " " + "output_427" + ens + "monthly_" + updown +  "_sel.nc"
subprocess.check_call(cdo_cmd, shell=True)




