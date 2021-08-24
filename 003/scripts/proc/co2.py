import sys
import os
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.means import global_int
from misc.translate import translate_varname
from misc.dirnames import get_datadir
from misc.filenames import *
from misc import par
from proc.r1 import save_r1
import numpy as np
import pickle
from netCDF4 import Dataset

def save_co2(sim, **kwargs):
    # sim is the name of the simulation, e.g. rcp85

    model = kwargs.get('model') # name of model
    timemean = kwargs.get('timemean', 'yearmean') # do annual mean? (bool)
    yr_span = kwargs.get('yr_span') # considered span of years
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    try_load = kwargs.get('try_load', 1)

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    if sim in ['rcp85', 'historical']:
        # All models follow same CO2 trajectory for RCP8.5 pathway so use MPI-ESM-LR data for all
        co2mass_file = filenames_raw(sim, 'co2mass', model='MPI-ESM-LR', timemean=timemean, yr_span=yr_span)
    else:
        co2mass_file = filenames_raw(sim, 'co2mass', model=model, timemean=timemean, yr_span=yr_span)
    co2mass = co2mass_file.variables['co2mass'][:]
    ps_file = filenames_raw(sim, 'ps', model=model, timemean=timemean, yr_span=yr_span)
    ps = ps_file.variables['ps'][:].filled(fill_value=np.nan)
    
    grid = {}
    grid['lon'] = ps_file.variables['lon'][:].filled(fill_value=np.nan)
    grid['lat'] = ps_file.variables['lat'][:].filled(fill_value=np.nan)

    ps_g = global_int(ps, grid, par.a, dim_lat=1, dim_lon=2)
    co2 = 1e6 * par.g * co2mass/ps_g*par.rho_air/par.rho_co2

    pickle.dump(co2, open(remove_repdots('%s/co2.%s.pickle' % (datadir, timemean)), 'wb'))