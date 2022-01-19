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

def save_ke(sim, **kwargs):
    # saves various keiation data (e.g., clear vs cloudy sky, lw vs sw) in one file
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
    ke = {}

    # counter so grid is only loaded once
    loaded_grid = 0

    # variable names
    varname = 'veke'
    subvarnames = ['tke', 'mke', 'eke', 'teke', 'seke', 'reke', 'res']

    # load all variables
    file = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

    if loaded_grid == 0:
        grid = {}
        grid['lat'] = file.variables['lat'][:]
        # grid['lon'] = file.variables['lon'][:]
        loaded_grid = 1

    for subvarname in subvarnames:
        # print(subvarname)
        ke[translate_varname(subvarname)] = np.squeeze(file.variables[subvarname][:])

    # variable names
    varname = 'eke'
    mean_subvarnames = ['mtke', 'mmke', 'meke', 'mteke', 'mseke', 'meke', 'mres']

    # load all variables
    file = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

    for mean_subvarname in mean_subvarnames:
        subvarname = mean_subvarname[1:]
        ke[mean_subvarname] = np.squeeze(file.variables[subvarname][:])

    pickle.dump([ke, grid], open(remove_repdots('%s/ke.%s.%s.pickle' % (datadir, zonmean, timemean)), 'wb'))
