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

def save_circ(sim, **kwargs):
    # saves various large scale circulation data (e.g., ua, va, wap, streamfunc) in one file
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
    circ = {}

    # counter so grid is only loaded once
    loaded_grid = 0

    # variable names
    if sim == 'echam':
        varnames = ['tdyn0', 'sdyn0', 'tdyns', 'sdyns', 'ahfl', 'ahfs']
    elif sim == 'era5':
        varnames = ['ssr', 'str', 'tsr', 'ttr', 'slhf', 'sshf']
    else:
        varnames = ['psi']

    # load all variables
    print(varnames)
    for varname in varnames:
        print(varname)
        file[varname] = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

        if loaded_grid == 0:
            grid = {}
            # grid['lon'] = file[varname].variables['lon'][:]
            grid['lat'] = file[varname].variables['lat'][:]
            grid['lev'] = file[varname].variables['plev'][:]
            loaded_grid = 1

            # create subsurface mask
            ps_file = filenames_raw(sim, 'ps', model=model, timemean=timemean, yr_span=yr_span)
            ps = ps_file.variables['ps'][:]
            ps_ext = np.tile(ps, [len(grid['lev']),1,1,1])
            ps_ext = np.transpose(ps_ext, [1,0,2,3])
            pa_ext = np.tile(grid['lev'], [ps.shape[0], ps.shape[1], ps.shape[2], 1])
            pa_ext = np.transpose(pa_ext, [0,3,1,2])
            subsurf = pa_ext > ps_ext

            ps = None; ps_ext = None; pa_ext = None;

        circ[translate_varname(varname)] = np.squeeze(file[varname].variables[varname][:])
        # make subsurface data nan
        circ[translate_varname(varname)] = circ[translate_varname(varname)].filled(fill_value=np.nan)
        if varname is not 'psi':
            circ[translate_varname(varname)][subsurf] = np.nan

    if zonmean:
        for circname in circ:
            if circname is not 'psi':
                print(circname)
                circ[circname] = np.nanmean(circ[circname], 3)
    
    pickle.dump([circ, grid], open(remove_repdots('%s/circ.%s.%s.pickle' % (datadir, zonmean, timemean)), 'wb'))
