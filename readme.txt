README
Data and Codes for “Monsoon Hysteresis reveals Atmospheric Memory” from Katzenberger, Anja and Levermann, Anders

This repository contains code and preprocessed data to reproduce the results for the manuscript “Monsoon Hysteresis reveals Atmospheric Memory”. This README is structured as below:

1.)	Description of codes
2.)	Description of data sets
3.)	List of codes to create individual figures
4.)	List of datasets as required by codes
5.)	Requirements and versions

All results were produced using Python (version 3.11.7). In order to reproduce specific figures, 
select the required code from the list given in (3.) and adapt the paths to datasets following (4.). 
You will also have to adapt the savedir-path by giving a directory where you wish to save the created
figures. The codes where run in the virtual environment new_venv. To simplify the process of installing  
the required Python packages according to this environment, you can install the requirements following  
requirements.txt or the list in (5.). The command to install the virtual environment packages on the target machine 
is  pip install -r requirements.txt 

1.)	DESCRIPTION OF CODES

wind.py: This code create wind direction and precipitation distribution on the Monsoon Planet for every month 

hysteresis_10N_monthly.ipynb: This notebook explores key variables on the Monsoon Planet as precipitation, evaporation, surface temperature, solar radiation, etc. based on monthly data. It also creates the hysteresis plot based on monthly values.

hysteresis_10N_daily.py: This code creates hysteresis plot based on daily data that is smoothened using singular spectrum analysis. 

3D_hysteresis.ipynb: This notebook creates the interactive 3D hysteresis plot depending on upper (solar radiation) and lower (surface temperature) boundary conditions. 

Hysteresis_real_monsoon.py: This code creates hysteresis plot of observed monsoon systems based on satellite data and observed rainfall. 

MonsoonSeasonalCycle_recalibratedHysteresis.py: This code implements the Moisture Advection Feedback and creates the ‘theory hysteresis plot’. 

MSC_recalHyst_withPnWadjustment.py: As in MonsoonSeasonalCycle_recalibratedHysteresis.py but with an implemented delay as described in the equations. 

hysteresis_observations_model_theory.py: This code plots the hysteresis plot of real monsoon systems, the Monsoon Planet hysteresis as well as the theory hystereis plot in one figure to simplify comparison. 

memory.py: This code helps to explore the memory experiments and quantify the memory effect within the atmosphere. 

bistability.py: This code creates the bistability plot and visualizes the experimental setup. 

Precip_Slider_SW_2D_depending_on_swvalue.py, Precip_Slider_tsurf_2D_depending_on_tsurfvalue.py: These two codes create a visualisation of the hysteresis on the Monsoon Planet depending on upper and lower boundary conditions. 

robustness_check_sol.py,
robustness_check_co2.py,
robustness_check_aerosol.py,
robustness_check_albedo.py: These four codes check the robustness of the hysteresis with regard to changes in the baseline setup. 


2.)	DESCRIPTION OF DATASETS
Observational datasets 

Observed shortwave down (top of atmosphere) data from CERES (monthly 2001-2019, ymonmean)
https://doi.org/10.1175/JCLI-D-17-0208.1

Observed precipitation data from GPCC (monthly 2001-2019, ymonmean)
https://dx.doi.org/10.5676/DWD_GPCC/FD_D_V2020_100

Observed Land surface temperatures from TERRA-MODIS (monthly 2001-2018, ymonmean)
https://dx.doi.org/10.5285/32d7bc64c7b740e9ad7a43589ab91592

Monsoon Planet data: 
The Monsoon Planet simulation data was preprocessed using CDOs in order to provide the relevant data in reasonable size.

slab_10years_ydaymean: Daily Monsoon Planet data averaged over 10years using CDO ydaymean for slab depths of 50m, 100m, 200m, 500m

slab_20years_ymonmean: Monthly Monsoon Planet data averaged over 20years using CDO ymonmean for slab depths of 50m, 100m, 200m, 500m

sol_20years_ymonmean: Monthly Monsoon Planet data averaged over 20years using CDO ymonmean for different solar radiation conditions

co2_20years_ymonmean: Monthly Monsoon Planet data averaged over 20years using CDO ymonmean for varying carbon dioxide concentration in the atmosphere

aerosol_20years_ymonmean: Monthly Monsoon Planet data averaged over 20years using CDO ymonmean for different atmospheric sulfate aerosol concentrations

albedo_20years_ymonmean: Monthly Monsoon Planet data averaged over 20years using CDO ymonmean for varying land surface albedo

memory: Daily Monsoon Planet Data with experiments of interrupted solar radiation for August, September and October. The file output_0_memory_selection.nc refers to the baseline simulation with natural solar radiation. 

bistability: Bistability experiments (daily output) with solar radiation kept fixed after reaching specific radiation sthreshold, e.g. 427 W/m^2. 

random_year: Daily Monsoon Planet data for one specific year


3.)	LIST OF CODES TO CREATE INDIVIDUAL FIGURES
Main
Fig 1: Hysteresis_real_monsoon.py
Fig 2 A, B: hysteresis_10N_daily.py
-2C: hysteresis_10N_monthly.ipynb
-2D: hysteresis_oservations_model_theory.py
Fig 3: memory.py
Fig 4: bistability.py



4.)	LIST OF DATASETS AS REQUIRED BY CODES
To reproduce the results, you will need to adapt the directory path at the beginning of the codes according to the following list: 

Hysteresis_real_monsoon.py
-	Observed shortwave down (top of atmosphere) data from CERES (monthly 2001-2019, ymonmean)
https://doi.org/10.1175/JCLI-D-17-0208.1

-	Observed precipitation data from GPCC (monthly 2001-2019, ymonmean)
https://dx.doi.org/10.5676/DWD_GPCC/FD_D_V2020_100

-	Observed Land surface temperatures from TERRA-MODIS (monthly 2001-2018, ymonmean)
https://dx.doi.org/10.5285/32d7bc64c7b740e9ad7a43589ab91592

hysteresis_10N_daily.py
-	slab_10years_ydaymean

hysteresis_10N_monthly.ipynb
-	slab_20years_ymonmean

robustness_check_sol.py
-	sol_20years_ymonmean

robustness_check_co2.py
-	co2_20years_ymonmean

robustness_check_aerosol.py
-	aerosol_20years_ymonmean

robustness_check_albedo.py
-	albedo_20years_ymonmean

memory.py
-	memory 

wind.py
-	slab_20years_ymonmean

3D_hysteresis.ipynb
-	slab_10years_ydaymean

hysteresis_oservations_model_theory.py
-	relevant data snippet plus source is given within the code

bistability.py
-	bistability 

Precip_Slider_SW_2D_depending_on_swvalue.py
-	slab_20years_ymonmean
-	random_year

Precip_Slider_tsurf-2D_depending_on_tsurfvalue.py
-	slab_20years_ymonmean
-	random_year



5.)	REQUIREMENTS AND VERSIONS
The codes were run using Python Version 3.11.7. You can install the required virtual environment packages on the target machine by using the command ‘pip install -r requirements.txt’. Below a list of the integrated version of the modules is given. 

asttokens==2.4.1
attrs==23.2.0
Cartopy==0.22.0
certifi==2023.11.17
cftime==1.6.3
colorama==0.4.6
comm==0.2.1
contourpy==1.2.0
cycler==0.12.1
debugpy==1.8.0
decorator==5.1.1
executing==2.0.1
fastjsonschema==2.19.1
fonttools==4.47.2
imageio==2.34.0
ipykernel==6.29.0
ipython==8.20.0
jedi==0.19.1
joblib==1.3.2
jsonschema==4.21.1
jsonschema-specifications==2023.12.1
jupyter_client==8.6.0
jupyter_core==5.7.1
kaleido==0.2.1
kiwisolver==1.4.5
llvmlite==0.41.1
matplotlib==3.8.2
matplotlib-inline==0.1.6
nbformat==5.9.2
nest-asyncio==1.6.0
netCDF4==1.6.5
numba==0.58.1
numpy==1.26.3
packaging==23.2
pandas==2.2.0
parso==0.8.3
pillow==10.2.0
platformdirs==4.1.0
plotly==5.18.0
prompt-toolkit==3.0.43
psutil==5.9.8
pure-eval==0.2.2
Pygments==2.17.2
pyparsing==3.1.1
pyproj==3.6.1
pyshp==2.3.1
python-dateutil==2.8.2
pyts==0.13.0
pytz==2023.3.post1
pywin32==306
pyzmq==25.1.2
referencing==0.32.1
rpds-py==0.17.1
scikit-learn==1.4.0
scipy==1.12.0
shapely==2.0.3
six==1.16.0
stack-data==0.6.3
tenacity==8.2.3
threadpoolctl==3.2.0
tornado==6.4
traitlets==5.14.1
tzdata==2023.4
wcwidth==0.2.13
xarray==2024.1.0

