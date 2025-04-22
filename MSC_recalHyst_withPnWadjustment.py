### -----------------------------------------------------
### MINIMALISTIC THEORY OF MOISTURE ADVECTION FEEDBACK
### -----------------------------------------------------
# Author: Anders Levermann, contact: anders.levermann@pik-potsdam.de

# This code implements the minimalistic theory of moisture advection feedback 
# (with delay) 
# as introduced in Katzenberger & Levermann: Seasonal Monsoon Hysteresis reveals Atmospheric Memory


#%%
### LOAD PACKAGES
#------------------

import numpy as np
import matplotlib.pylab as plt
from numpy import genfromtxt
import random as rd

#%%
# Parameters
#-------------------
SecsInDay = 24 * 60 * 60 # s/day
SecsInYear = 365 * 24 * 60 * 60 # s/yr
SecsInMonth = SecsInYear/12 # s/month
DaysInYear = 365
Pfactor = SecsInDay/1000*1000 # kg/m2/s * [SecsInDay s/day / (1000 kg/m3) * (1000 mm/m) ]  = mm/day 

TimeRange = 1 * SecsInYear # s
dt = 10 # s
Tnum = int(TimeRange/dt)
print("TimeRange (s):", TimeRange)
print("Tnum:", Tnum)
print()

# physical variables
Cp = 1.295 # J/m3/K
Height = 15e3 #  m
LatentH = 2.6e6 # J/kg
Depth = 2.725e3 # m
Alpha = 0.6  # m/s/K 
Beta = 0.8 * SecsInDay / 1000 / 1000 # (0.8 mm/day / g/kg ) * 1/1000 m/mm * SecsInDay s/day *1/1000 kg/g = 0.8 * SecsInDay/1000/1000
rhoA = 1.24e-3 # kg / m3
qO = 5.1e-3#
PO = Beta*qO
TO = 5 # K above T_base

# Solar insulation
Rsol_base = 380 # W/m2
Rsol_amp = 70 # W/m2
Rsol_freq = 2*np.pi/SecsInYear
t = SecsInYear/4
Rsol = Rsol_base + Rsol_amp * np.sin(Rsol_freq * t)

# Outgoing radiation
sigma = 5.670374419e-8 # W/m2/K4 
Tbase = 293  # K
Gamma = 4*sigma*Tbase**3 # W/m2/K
Rout_base = sigma*Tbase**4 
Rout = Rout_base + Gamma * TO  # W/m2 for T=5K

# Time constants for slow adjustment
tauW = 2.0 * SecsInDay #3.5 * SecsInDay #3 * SecsInDay # s
tauP = 1/24 * SecsInDay #3 * SecsInDay # s

print("max(Rsol):", Rsol_base+Rsol_amp)
print("Rout(Tval=0):", Rout_base + Gamma*TO)
print()


#%%
#### Example of T=5K and P=7 mm/day
#-----------------------------------

Tval = 5 # K # Temperature difference  Land minus Ocean
Pval = 7 / Pfactor  # kg/m3 mm/day * m/mm * day/s = kg/m3 * m/s = kg/m2/s

print("PO (kg/m2/s):", PO)
print("PO (mm/day):", PO*Pfactor)
print("qO:", qO)
print("qL:", Pval/Beta, "(=Pval/Beta)")
print("Pval:", Pval)
print("Pfactor:", Pfactor)

####### Radiation
Rout = Rout_base + Gamma * (TO+Tval)  # W/m2 for T=5K
Rad = Rsol-Rout

print()
print("Rsol:", Rsol)
print("Rout_base:", Rout_base)
print("Rout Temp:", Gamma * (TO+Tval))
print("Gamma (W/m2/K):", Gamma)
print("Rsol-Rout:", Rsol-Rout)

####### Temperature time evolution with cp*H still on the left-hand side
valT1 = Rad
valT2 = LatentH*Pval
valT3 = -Alpha*Cp*Height/Depth * Tval*Tval

print()
print("Temperature equation in W/m2")
print("Cp*Height*dT/dt=")
print("Radiation: valT1:", valT1)
print("Latent Heat: valT2:", valT2)
print("Convergence valT3:", valT3)
print("sum:", valT1+valT2+valT3)

######## Precipitation 
valP1 = Alpha/Depth*Tval*Beta*qO*Pfactor
valP2 = -Alpha/Depth*Tval*Pval*Pfactor
valP3 = - Beta/rhoA/Height*Pval*Pfactor

print()
print("Precipitation equation in mm/day/s")
print("dP/t=")
print("valP1:", valP1)
print("valP2:", valP2)
print("valP3:", valP3)
print("sum:", valP1+valP2+valP3)

######## Stability
print()
print("Slope with P:", (valP2+valP3)/Pval/Pfactor, " [1/((mm/day)*s)]")
print("Slope with T:", valP1/Tval, " [1/(K*s)]")

######## Prime parameters    
Gamma_prime = Gamma/(Cp*Height)
LatentH_prime = LatentH/(Cp*Height)
Alpha_prime = Alpha/Depth
Beta_prime = Beta/(rhoA*Height)
PO = Beta*qO

Rsol = []
Rout = []
Rnet = []
Temp = []
Prec = []
Time = []

Tval = 4.92401 # K
Pval = 6.87894 / Pfactor  # kg/m3 mm/day * m/mm * day/s = kg/m3 * m/s = kg/m2/s
Tval = 0
qval = 0
Wval = Alpha*Tval
Pval = Beta*qval

#%%
# Timeseries
#-------------------

t0=-SecsInYear/8 # mid winter
t=t0
print("Solar cycle:", Rsol_amp * np.sin(Rsol_freq * t))
for i in range(Tnum):
    t = t + dt
    
    #### Varying solar insulation (starting in winter)
    Rsol_val = Rsol_base + Rsol_amp * np.sin(Rsol_freq * t)
    
    Rout_val = Rout_base + Gamma * (TO+Tval)  # W/m2 for T=5K
    Rnet_val = Rsol_val-Rout_val
    
    R_prime = Rnet_val/(Cp*Height)
 
    valT1 =   R_prime # radiative heat loss
    valT2 =   LatentH_prime*Pval    # Latent heat
    valT3 = - Alpha/Depth * Tval*Tval # Heat advection
    
    valq1 =   Wval/Depth*qO
    valq2 = - Wval/Depth*qval
    valq3 = - Pval/(rhoA*Height)
   
    Tval_new = Tval + dt * (valT1 + valT2 + valT3)
    qval_new = qval + dt * (valq1 + valq2 + valq3)
    
    Wval_new = Wval + dt * (Alpha*Tval - Wval)/tauW
    Pval_new = Pval + dt * (Beta*qval - Pval)/tauP
    
    Wval = Wval_new
    Tval = Tval_new
    Pval = Pval_new
    qval = qval_new
    
    if qval < 0:
        qval = 0
    if Tval < 0:
        Tval = 0
    if Wval < 0:
        Wval = 0
    if Pval < 0:
        Pval = 0

    if t%100000 == 0:
        print(t/SecsInYear, Tval, Wval/Alpha, Wval, Pval*Pfactor,Rsol_val,Rout_val,Rnet_val)
    
        Rsol.append(Rsol_val)
        Rout.append(Rout_val)
        Rnet.append(Rnet_val)
        Temp.append(Tval)
        Prec.append(Pval*Pfactor)
        Time.append((t-t0)/SecsInMonth)

#%%
### VISUALIZATION
#-------------------
        
figwidth = 5

fig, axes = plt.subplots(2,2,figsize=(2*figwidth,2*figwidth))

ax=axes[0,0]
ax.plot(Time,Rsol,'y',label="Solar")
ax.plot(Time,Rnet,'b',label="Net")
ax.plot(Time[0],Rsol[0],'x')
ax.plot(Time[0],Rnet[0],'x')
ax.set_xlabel("Time (Month)")
ax.set_ylabel("Radiation (W/m2)")
ax.legend()

ax=axes[1,0]
ax.plot(Time,Prec,'b')
ax.plot(Time[0],Prec[0],'x')
ax.set_xlabel("Time (Month)")
ax.set_ylabel("Precipitation (mm/day)")

ax=axes[0,1]
ax.plot(Time,Temp,'r',label="T")
ax.plot(Time[0],Temp[0],'x')
ax.plot(Time[-1],Temp[-1],'or')
ax.set_xlabel("Time (Month)")
ax.set_ylabel("(Land) Temperature (°C)")

ax=axes[1,1]
ax.plot(Rsol,Prec,'b')
ax.plot(Rsol[0],Prec[0],'x')
ax.plot(Rsol[-1],Prec[-1],'or')
ax.set_xlabel("Solar insulation (W/m2)")
ax.set_ylabel("Precipitation (mm/day)")

plt.show()


# %%
