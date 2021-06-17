import sys
import os
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.translate import translate_varname, latetrans_grid
from misc.dirnames import get_datadir
from misc.filenames import filenames_raw
from misc.vertconv import pa_to_sigma
from misc import par
import numpy as np
import pickle
from netCDF4 import Dataset

def make_ga_dev(sim, **kwargs):
    # computes vertically integrated lapse rate deviation from a moist adiabat
    # sim is the name of the simulation, e.g. rcp85

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', '.zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    vertcoord = kwargs.get('vertcoord', '.si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    if vertcoord == '.si':
        si_std = np.linspace(1e-2,1,100) # standard sigma grid to convert to

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled lapse rate deviation data
    ga_dev_file = '%s/ga_dev%s%s%s.pickle' % (datadir, vertcoord, zonmean, timemean)

    # load data required to compute reanalysis/model lapse rate
    ga_indata = {}
    file_vertconv = {}
    if not (os.path.isfile(ga_dev_file) and try_load):
        # variable names to load (load surface pressure first because it is required to convert the 3D variables to other vertical coordinate systems)
        if sim == 'echam':
            varnames = ['sp', 't', 'z', 'ts', 'zs']
        elif sim == 'era5':
            varnames = ['sp', 'skt', 't', 'zs', 'z']
        else:
            varnames = ['ps', 'ta', 'zg', 'ts', 'zs']

        for varname in varnames:
            varname_std = translate_varname(varname)
            ga_indata[varname_std] = load_var(sim, varname, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load) # load raw data (in pressure grid)

            if (varname_std in ['ta', 'zg'] and vertcoord == '.si'): # convert 3D data to sigma coordinate
                file_vertconv[varname_std] = '%s/%s%s%s%s.pickle' % (datadir, varname, vertcoord, zonmean, timemean)

                if (os.path.isfile(file_vertconv[varname_std]) and try_load): # load pickled data if available
                    ga_indata[varname_std] = pickle.load(open(file_vertconv, 'rb'))

                else: # if not convert and save
                    if varname_std == 'ta': 
                        ga_indata[varname_std] = pa_to_sigma(varname_std, ga_indata[varname_std], 'ts', ga_indata['ts'], ga_indata['ps'], si_std)
                    elif varname_std == 'zg': 
                        ga_indata[varname_std] = pa_to_sigma(varname_std, ga_indata[varname_std], 'orog', ga_indata['orog'], ga_indata['ps'], si_std)

                    pickle.dump(ga_indata[translate_varname(varname)], open('%s/%s%s%s%s.pickle' % (datadir, varname, vertcoord, zonmean, timemean), 'wb'))

        # compute lapse rate
        make_ga(sim, ga_indata, model=model, vertcoord=vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

    sys.exit()

    [ga_dev, grid] = pickle.load(open(ga_dev_file, 'rb'))



    pickle.dump([ga_dev_vint, grid], open('%s/ga_dev_vint%s%s%s.pickle' % (datadir, vertcoord, zonmean, timemean), 'wb'))

    rlut = None; rsdt = None; rsut = None; rsus = None; rsds = None; rlds = None; rlus = None;

def make_ga_m(sim, **kwargs):
    # Computes the moist adiabatic lapse rate
    placeholder = None

def make_ga(sim, ga_indata, **kwargs):
    # Computes the reanalysis/model lapse rate
    
    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', '.zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    vertcoord = kwargs.get('vertcoord', '.si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled data if available
    file = '%s/ga_%s%s%s.pickle' % (datadir, vertcoord, zonmean, timemean)

    if (os.path.isfile(file) and try_load):
        alldata = pickle.load(open(file, 'rb'))
    else:
        alldata = None

    return alldata

def load_var(sim, varname, **kwargs):

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', '.zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled data if available
    file = '%s/%s%s%s.pickle' % (datadir, varname, zonmean, timemean)

    if (os.path.isfile(file) and try_load):
        alldata = pickle.load(open(file, 'rb'))
    else:
        file_raw = {}
        grid = {}
        alldata = {}

        file_raw[varname] = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

        grid = {}
        grid['lat'] = file_raw[varname].variables[latetrans_grid(sim, 'lat')][:]
        grid['lon'] = file_raw[varname].variables[latetrans_grid(sim, 'lon')][:]
        try: # try to load vertical grid  data
            grid['lev'] = file_raw[varname].variables[latetrans_grid(sim, 'lev')][:]
            if sim == 'era5':
                grid['lev'] = 100 * grid['lev'] # convert hPa to Pa
        except:
            print('%s does not contain vertical level data.' % (varname))

        vardata = file_raw[varname].variables[varname][:]

        if sim == 'era5' and varname in ['z', 'zs']:
            vardata = 1/par.g * vardata # convert geopotential (m**2 s**-2) to geopotential height (m)

        alldata[translate_varname(varname)] = vardata
        alldata['grid'] = grid

        pickle.dump(alldata, open('%s/%s%s%s.pickle' % (datadir, varname, zonmean, timemean), 'wb'))

    return alldata
