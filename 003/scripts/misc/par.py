import numpy as np

# scalars
a = 6378.1370e3 # Radius of Earth [m]

sigma = 5.67e-8 # stephan boltzmann constant [W m**-2 K**-4]
g = 9.80665 # gravitational acceleration [m s**-2]

rho_air = 1.225 # density of dry air [kg m**-3]
rho_co2 = 1.98 # density of CO2 [kg m**-3]

# vectors
si50 = np.linspace(1e-2,1,50)
si100 = np.linspace(1e-2,1,100)

lon360 = np.linspace(0,360,360)
lat180 = np.linspace(-90,90,180)