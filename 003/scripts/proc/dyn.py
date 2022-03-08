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

def save_dyn(sim, **kwargs):
    # saves various dyniation data (e.g., clear vs cloudy sky, lw vs sw) in one file
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
    dyn = {}

    # counter so grid is only loaded once
    loaded_grid = 0

    # variable names
    if sim == 'echam':
        varnames = ['tdyn0', 'sdyn0', 'tdyns', 'sdyns', 'ahfl', 'ahfs']
    elif sim == 'era5':
        varnames = ['ssr', 'str', 'tsr', 'ttr', 'slhf', 'sshf']
    else:
        varnames = ['aht', 'qaht', 'saht', 'vmmmc', 'vqmmc', 'vsmmc', 'vmse', 'vqse', 'vsse', 'vmte', 'vqte', 'vste', 'daht', 'dqaht', 'dsaht', 'dvmmmc', 'dvqmmc', 'dvsmmc', 'dvmse', 'dvqse', 'dvsse', 'dvmte', 'dvqte', 'dvste','aht_sm', 'qaht_sm', 'saht_sm', 'vmmmc_sm', 'vqmmc_sm', 'vsmmc_sm', 'vmse_sm', 'vqse_sm', 'vsse_sm', 'vmte_sm', 'vqte_sm', 'vste_sm', 'daht_sm', 'dqaht_sm', 'dsaht_sm', 'dvmmmc_sm', 'dvqmmc_sm', 'dvsmmc_sm', 'dvmse_sm', 'dvqse_sm', 'dvsse_sm', 'dvmte_sm', 'dvqte_sm', 'dvste_sm']

    # load all variables
    for varname in varnames:
        file[varname] = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

        if loaded_grid == 0:
            grid = {}
            grid['lat'] = file[varname].variables['lat'][:]
            loaded_grid = 1

        dyn[translate_varname(varname)] = np.squeeze(file[varname].variables[varname][:])
        if sim == 'era5':
            dyn[translate_varname(varname)] = dyn[translate_varname(varname)]/86400

    # if zonmean:
    #     for dynname in dyn:
    #         dyn[dynname] = np.mean(dyn[dynname], 2)
    
    pickle.dump([dyn, grid], open(remove_repdots('%s/dyn.%s.%s.pickle' % (datadir, zonmean, timemean)), 'wb'))
