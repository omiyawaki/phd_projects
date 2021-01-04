from seasonalebm import SeasonalEbm
import climlab
import numpy as np
import os
import pickle
from IPython import get_ipython
get_ipython().magic('reset -sf')

def get_data():
    seasonal = SeasonalEbm(plotpath, datapath, tstep=tstep, tsteps=tsteps, N_lat=N_lat, N_lev=N_lev, mld=mld, Ninit=Ninit, albedo=albedo, D=D, U=U, Cd=Cd, conv_scheme=conv_scheme, transport=transport, forcing=forcing)
    seasonal.run_nyears() # initialize
    seasonal.run_year() # run for one additional year and save monthly data

    # outfile = open('%s/ebm_init.pickle' % (datapath), 'wb')
    # pickle.dump(seasonal, outfile)
    # outfile.close()

    return seasonal

def load_data():
    infile = open('%s/ebm_init.pickle' % (datapath), 'rb')
    seasonal = pickle.load(infile)
    infile.close()

    return seasonal

def make_plots(seasonal):
    seasonal.plot_energy_fluxes('mse_all') # plot all annual mean energy fluxes in MSE equation
    seasonal.plot_mon_lat(0.1, 'r1', 'R1', 'r1_mon_lat') # plot R1
    seasonal.plot_mse('nhmid', 'mse_nhmid') # MSE seasonality
    seasonal.plot_dr1(0.1, 'nhmid', 'dr1_nhmid') # Delta R1 decomposition
    seasonal.plot_temp_sel('temp_sel') # temperature profile at select latitudes

# model parameters
tstep = 24*3600 # long timestep for radiation
tsteps = 3600 # short timestep for everything else
N_lat = 45
N_lev = 30
Cd = 3e-3
phi_h = 25*np.pi/180 # hadley cell edge
Dx = 1.6 # extratropical diffusivity [W/m^2/K]
Dt = 3.2 # tropical diffusivity [W/m^2/K]
lat = climlab.Axis(axis_type='Latitude', num_points=N_lat)
phi = lat.bounds*np.pi/180 # latitude edges
S = 1/2*(1-np.tanh(np.pi*(phi-phi_h)/phi_h)*np.tanh(np.pi*(phi+phi_h)/phi_h))
D = np.tile(Dx + (Dt-Dx)*S, (N_lev, 1))
lev = climlab.Axis(axis_type='Pressure', num_points=N_lev)
D = D*np.expand_dims(lev.points, axis=1) / lev.points[-1] # diffusivity weakens with height
# D = 0.04
U = 0
albedo = 0.25 # surface albedo
conv_scheme = 'MALR' # enter lapse rate in K/km or MALR/DALR
transport = 'dry_diff'
forcing = None # None=default CO2 (380 ppmv), 4xco2=4x default

# mld_list = [40,10,50,80] # mixed layer depth
# mld_list = [50,60,70,80] # mixed layer depth
# mld_list = [10,20,30,40,50,80] # mixed layer depth
mld_list = [20, 40] # mixed layer depth

for mld in mld_list:
    print('d=%g' % (mld))
    # Ninit=1
    Ninit = int(3/4*max(int(mld/2), 10)) # number of years to initialize
    plotpath = './figures/d%g' % (mld)
    datapath = './data/d%g' % (mld)

    # create directories if they doesn't exist
    if forcing is None and albedo == 0.25:
        plotpath = './figures/d%g' % (mld)
    elif albedo != 0.25:
        plotpath = './figures/alb%g/d%g' % (albedo, mld)
    else:
        plotpath = './figures/%s/d%g' % (forcing, mld)
    if not os.path.isdir(plotpath):
        os.makedirs(plotpath)
    datapath = './data/d%g' % (mld)
    if not os.path.isdir(datapath):
        os.makedirs(datapath)

    seasonal = get_data()
    make_plots(seasonal)
