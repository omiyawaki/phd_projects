import sys
import os
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.translate import translate_varname
from misc.dirnames import get_datadir
from misc.filenames import *
from proc.r1 import save_r1
import numpy as np
from scipy.interpolate import interp1d
import pickle
from netCDF4 import Dataset

def save_hydro(sim, **kwargs):
    # saves various hydroiation data (e.g., clear vs cloudy sky, lw vs sw) in one file
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
    hydro = {}

    # counter so grid is only loaded once
    loaded_grid = 0

    # variable names
    if sim == 'echam':
        print('todo: fix variable names')
        # fix variable names below
        varnames = ['precip', 'aprc']
    elif sim == 'era5':
        varnames = ['cp', 'lsp', 'e']
    else:
        # varnames = ['pr', 'prc', 'prl', 'prfrac', 'evspsbl', 'vhur', 'tas', 't850']
        varnames = ['pr', 'prc', 'prfrac', 'evspsbl']
        # varnames = ['pr', 'prc', 'evspsbl']
        # varnames = ['pr', 'prc', 'clwvi', 'clivi']
        # varnames = ['pr', 'prc', 'evspsbl', 'prw']

    # load all variables
    for varname in varnames:
        file[varname] = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

        if loaded_grid == 0:
            grid = {}
            grid['lat'] = file[varname].variables['lat'][:]
            grid['lon'] = file[varname].variables['lon'][:]
            loaded_grid = 1
        elif varname == 'vhur':
            grid['lat3d'] = file[varname].variables['lat'][:]
            grid['lon3d'] = file[varname].variables['lon'][:]

        if not varname == 't850':
            hydro[translate_varname(varname)] = np.squeeze(file[varname].variables[varname][:])
        else:
            hydro[translate_varname(varname)] = np.squeeze(file[varname].variables['ta'][:])

        if ( ('pr' in translate_varname(varname)) or (translate_varname(varname) == 'evspsbl') ) and not (translate_varname(varname) == 'prw'):
            if sim == 'era5': # convert m accumulated over a day to mm/d
                hydro[translate_varname(varname)] = hydro[translate_varname(varname)]*1e3
            else: # convert kg m**-2 s**-1 to mm/d
                hydro[translate_varname(varname)] = hydro[translate_varname(varname)]*86400

        if varname in ['vhur', 't850'] and not np.any(grid['lat'] == grid['lat3d']):
            vhur_fint = interp1d(grid['lat3d'], hydro[varname].filled(fill_value=np.nan), axis=1, bounds_error=False)
            hydro[varname] = vhur_fint(grid['lat'].filled(fill_value=np.nan))

    # if sim == 'era5':
    #     hydro['pr'] = hydro['prc'] + hydro['prl']
    # else:
    #     hydro['prl'] = hydro['pr'] - hydro['prc']

    # # define convective precipitation fraction
    # hydro['prfrac'] = hydro['prc'] / hydro['pr']
    # if (hydro['prfrac'] > 1).any():
    #     hydro['prfrac'][hydro['prfrac'] > 1] = np.nan

    if zonmean:
        for hydroname in hydro:
            hydro[hydroname] = np.nanmean(hydro[hydroname], 2)

    pickle.dump([hydro, grid], open(remove_repdots('%s/hydro.%s.%s.pickle' % (datadir, zonmean, timemean)), 'wb'))
