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

def save_cl_in(sim, **kwargs):
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
    cl_in = {}

    # counter so grid is only loaded once
    loaded_grid = 0

    # variable names
    if sim == 'echam':
        print('todo: fix variable names')
        # fix variable names below
        # varnames = ['thydro0', 'shydro0', 'thydros', 'shydros', 'ahfl', 'ahfs']
    elif sim == 'era5':
        varnames = ['cp', 'lsp', 'e']
    else:
        varnames = ['ta', 'hur', 'hus']

    # load all variables
    for varname in varnames:
        file[varname] = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

        if loaded_grid == 0:
            grid = {}
            grid['lat'] = file[varname].variables['lat'][:]
            grid['lon'] = file[varname].variables['lon'][:]
            grid['lev'] = file[varname].variables['plev'][:]
            loaded_grid = 1

        cl_in[translate_varname(varname)] = np.squeeze(file[varname].variables[varname][:])

    if zonmean:
        for cl_inname in cl_in:
            cl_in[cl_inname] = np.nanmean(cl_in[cl_inname], 3)
    
    pickle.dump([cl_in, grid], open(remove_repdots('%s/cl_in.%s.%s.pickle' % (datadir, zonmean, timemean)), 'wb'))
