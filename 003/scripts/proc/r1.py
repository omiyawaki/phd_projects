import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.translate import translate_varname
from misc.dirnames import get_datadir
from misc.filenames import filenames_raw
import numpy as np
import pickle
from netCDF4 import Dataset

def save_r1(sim, **kwargs):
    # computes R1
    # sim is the name of the simulation, e.g. rcp85

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', '.zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '.yearmean') # do annual mean? (bool)
    yr_span = kwargs.get('yr_span') # considered span of years

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # initialize dictionaries
    file = {}
    grid = {}
    flux = {}

    # counter so grid is only loaded once
    loaded_grid = 0

    # variable names
    if sim == 'echam':
        varnames = ['trad0', 'srad0', 'trads', 'srads', 'ahfl', 'ahfs']
    elif sim == 'era5':
        varnames = ['ssr', 'str', 'tsr', 'ttr', 'slhf', 'sshf']
    else:
        varnames = ['rlut', 'rsdt', 'rsut', 'rsus', 'rsds', 'rlds', 'rlus', 'hfls', 'hfss']

    # load all variables required to compute R1
    for varname in varnames:
        file[varname] = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

        if loaded_grid == 0:
            grid = {}
            grid['lat'] = file[varname].variables['lat'][:]
            grid['lon'] = file[varname].variables['lon'][:]
            loaded_grid = 1

        flux[translate_varname(varname)] = file[varname].variables[varname][:]

    if sim == 'era5' or sim == 'echam':
        flux['ra'] = flux['trad0'] + flux['srad0'] - flux['trads'] - flux['srads'] 
        flux['hfls'] = -flux['hfls']
        flux['hfss'] = -flux['hfss']
    else:
        flux['ra'] = flux['rsdt'] - flux['rsut'] - flux['rlut'] + flux['rsus'] - flux['rsds'] + flux['rlus'] - flux['rlds']

    flux['stg_adv'] = flux['ra'] + flux['hfls'] + flux['hfss']

    if not zonmean:
        r1 = flux['stg_adv']/flux['ra']
        stg_adv = flux['stg_adv']
        ra = flux['ra']
    else:
        r1 = np.mean(flux['stg_adv'], 2)/np.mean(flux['ra'], 2)
        stg_adv = np.mean(flux['stg_adv'], 2)
        ra = np.mean(flux['ra'], 2)

    pickle.dump([r1, grid], open('%s/r1%s%s.pickle' % (datadir, zonmean, timemean), 'wb'))
    pickle.dump([stg_adv, grid], open('%s/stg_adv%s%s.pickle' % (datadir, zonmean, timemean), 'wb'))
    pickle.dump([ra, grid], open('%s/ra%s%s.pickle' % (datadir, zonmean, timemean), 'wb'))

    rlut = None; rsdt = None; rsut = None; rsus = None; rsds = None; rlds = None; rlus = None;