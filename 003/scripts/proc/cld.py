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

def save_cld(sim, **kwargs):
    # saves various cldiation data (e.g., clear vs cloudy sky, lw vs sw) in one file
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
    cld = {}

    # counter so grid is only loaded once
    loaded_grid = 0

    # variable names
    if sim == 'echam':
        print('todo: fix variable names')
        # fix variable names below
        varnames = ['clt']
    elif sim == 'era5':
        print('todo: fix variable names')
        varnames = ['clt']
    else:
        varnames = ['clivi','clt']

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
            cld[translate_varname(varname)] = np.squeeze(file[varname].variables[varname][:])
        else:
            cld[translate_varname(varname)] = np.squeeze(file[varname].variables['ta'][:])

    if zonmean:
        for cldname in cld:
            cld[cldname] = np.nanmean(cld[cldname], 2)

    pickle.dump([cld, grid], open(remove_repdots('%s/cld.%s.%s.pickle' % (datadir, zonmean, timemean)), 'wb'))
