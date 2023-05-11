import sys
import os
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.translate import translate_varname
from misc.dirnames import get_datadir
from misc.filenames import *
from proc.r1 import save_r1
import numpy as np
import pickle
from netCDF4 import Dataset

def save_rad(sim, **kwargs):
    # saves various radiation data (e.g., clear vs cloudy sky, lw vs sw) in one file
    # sim is the name of the simulation, e.g. rcp85

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', 'zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', 'yearmean') # do annual mean? (bool)
    yr_span = kwargs.get('yr_span') # considered span of years
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    try_load = kwargs.get('try_load', 1)

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # initialize dictionaries
    file = {}
    grid = {}
    rad = {}

    # counter so grid is only loaded once
    loaded_grid = 0

    # variable names
    if sim == 'echam':
        varnames = ['trad0', 'srad0', 'trads', 'srads', 'traf0', 'sraf0', 'trafs', 'srafs']
    elif sim == 'era5':
        varnames = ['ssr', 'str', 'tsr', 'ttr', 'slhf', 'sshf']
    else:
        # varnames = ['rlut', 'rlutcs', 'rsut', 'rsutcs', 'rsus', 'rsuscs', 'rsds', 'rsdscs', 'rlds', 'rldscs', 'rsdt', 'rlus', 'tmax']
        varnames = ['rlut', 'rlutcs', 'rsut', 'rsutcs', 'rsus', 'rsuscs', 'rsds', 'rsdscs', 'rlds', 'rldscs', 'rsdt', 'rlus']

    # load all variables
    print(varnames)
    for varname in varnames:
        print(varname)
        file[varname] = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

        if loaded_grid == 0:
            grid = {}
            grid['lat'] = file[varname].variables['lat'][:]
            grid['lon'] = file[varname].variables['lon'][:]
            loaded_grid = 1

        rad[translate_varname(varname)] = np.squeeze(file[varname].variables[varname][:])
        if sim == 'era5':
            rad[translate_varname(varname)] = rad[translate_varname(varname)]/86400

    # compute useful derived quantities
    if sim == 'era5' or sim == 'echam':
        rad['ra'] = rad['trad0'] + rad['srad0'] - rad['trads'] - rad['srads'] 
        rad['ra_cs'] = rad['traf0'] + rad['sraf0'] - rad['trafs'] - rad['srafs']
        rad['sw'] = rad['srad0'] - rad['srads'] 
        rad['sw_cs'] = rad['sraf0'] - rad['srafs'] 
        rad['lw'] = rad['trad0'] - rad['trads'] 
        rad['lw_cs'] = rad['traf0'] - rad['trafs'] 
    else:
        rad['ra'] = rad['rsdt'] - rad['rsut'] - rad['rlut'] + rad['rsus'] - rad['rsds'] + rad['rlus'] - rad['rlds']
        rad['ra_cs'] = rad['rsdt'] - rad['rsutcs'] - rad['rlutcs'] + rad['rsuscs'] - rad['rsdscs'] + rad['rlus'] - rad['rldscs']

        rad['sw'] = rad['rsdt'] - rad['rsut'] + rad['rsus'] - rad['rsds']
        rad['sw_cs'] = rad['rsdt'] - rad['rsutcs'] + rad['rsuscs'] - rad['rsdscs']

        rad['lw'] = -rad['rlut'] + rad['rlus'] - rad['rlds']
        rad['lw_cs'] = -rad['rlutcs'] + rad['rlus'] - rad['rldscs']

    # cloud radiative fluxes
    rad['ra_cld'] = rad['ra'] - rad['ra_cs']
    rad['sw_cld'] = rad['sw'] - rad['sw_cs']
    rad['lw_cld'] = rad['lw'] - rad['lw_cs']

    if zonmean:
        for radname in rad:
            print(radname)
            print(model)
            print(rad[radname].shape)
            rad[radname] = np.mean(rad[radname], 2)
    
    pickle.dump([rad, grid], open(remove_repdots('%s/rad.%s.%s.pickle' % (datadir, zonmean, timemean)), 'wb'))
