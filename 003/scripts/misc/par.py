import numpy as np

###############################################################
# SCALARS
###############################################################

a = 6378.1370e3 # Radius of Earth [m]

sigma = 5.67e-8 # stephan boltzmann constant [W m**-2 K**-4]
g = 9.80665 # gravitational acceleration [m s**-2]

rho_air = 1.225 # density of dry air [kg m**-3]
rho_co2 = 1.98 # density of CO2 [kg m**-3]

Lv = 2.501e6 # latent heat of vaporization [J kg**-1] (AMS glossary)
Lf = 3.337e5 # latent heat of fusion [J kg**-1] (AMS glossary)

###############################################################
# GRIDS
###############################################################

# standard plev grid
p50 = np.logspace(1,5,50)
# p50 = np.linspace(10,1e5,50)

# standard sigma grid
si50 = np.logspace(-1,0,50)
si100 = np.logspace(-1,0,100)
# si50 = np.linspace(0,1,50)
# si100 = np.linspace(0,1,100)

# standard lat/lon grid
lon360 = np.linspace(0,360,3600)
lat180 = np.linspace(-90,90,1800)

###############################################################
# BINS
###############################################################

# R1 bin
r1_bins = np.arange(-55, 150, 10)/100
