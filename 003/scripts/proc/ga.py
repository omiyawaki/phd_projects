import sys
import os
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.translate import *
from misc.dirnames import get_datadir
from misc.filenames import *
from misc.vertconv import pa_to_sigma
from misc import par
import numpy as np
from scipy import interpolate
import pickle
from netCDF4 import Dataset

def varnames_ga_m(sim):
    if sim == 'echam':
        varnames = ['sp', 't']
    elif sim == 'era5':
        varnames = ['sp', 't']
    else:
        varnames = ['ps', 'ta']

    return varnames

def varnames_ga(sim):
    if sim == 'echam':
        varnames = ['sp', 'ts', 't', 'zs', 'z']
    elif sim == 'era5':
        varnames = ['sp', 't2m', 't', 'zs', 'z']
    else:
        varnames = ['ps', 'tas', 'ta', 'orog', 'zg']

    return varnames

def make_ga_dev_vint(sim, vertbnd, **kwargs):
    # computes vertically integrated lapse rate deviation from a moist adiabat

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', 'zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled vertically-integrated lapse rate deviation data
    ga_dev_vint_file = remove_repdots('%s/ga_dev_vint%g.%g.%s.%s.%s.pickle' % (datadir, vertbnd[0], vertbnd[1], vertcoord, zonmean, timemean))

    if not (os.path.isfile(ga_dev_vint_file) and try_load):

        ga_dev = make_ga_dev(sim, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

        f = interpolate.interp1d(ga_dev['grid']['lev'], ga_dev['ga_dev'], axis=1)

        # identify which bound is the lower/upper bound
        if ((vertbnd[0] > vertbnd[1]) and vertcoord in ['si', 'pa'] ):
            vertbnd_lo = vertbnd[0]
            vertbnd_up = vertbnd[1]
        else:
            vertbnd_lo = vertbnd[1]
            vertbnd_up = vertbnd[0]
        # vertically interpolate between the provided bounds
        ga_dev_itp = f(np.linspace(vertbnd_up,vertbnd_lo,ga_dev['grid']['lev'].size))

        ga_dev_vint = {}
        ga_dev_vint['grid'] = {}
        ga_dev_vint['grid']['lon'] = ga_dev['grid']['lon']
        ga_dev_vint['grid']['lat'] = ga_dev['grid']['lat']
        ga_dev_vint['ga_dev_vint'] = np.mean(ga_dev_itp, axis=1)
        ga_dev_itp = None

        # plot_mon_lat_test(ga_dev_vint, 'ga_dev_vint')

        pickle.dump(ga_dev_vint, open(remove_repdots('%s/ga_dev_vint%g.%g.%s.%s.%s.pickle' % (datadir, vertbnd[0], vertbnd[1], vertcoord, zonmean, timemean)), 'wb'))
    else:
        ga_dev_vint = pickle.load(open(ga_dev_vint_file, 'rb'))

    if zonmean:
        ga_dev_vint['ga_dev_vint'] = np.mean(ga_dev_vint['ga_dev_vint'],2)

    return ga_dev_vint


def make_ga_dev(sim, **kwargs):
    # computes lapse rate deviation from a moist adiabat
    # sim is the name of the simulation, e.g. rcp85

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', 'zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    if vertcoord == 'si':
        si_std = par.si50 # standard sigma grid to convert to

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled lapse rate deviation data
    ga_dev_file = remove_repdots('%s/ga_dev.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean))

    if not (os.path.isfile(ga_dev_file) and try_load):

        # location of pickled lapse rate deviation data
        ga_file = remove_repdots('%s/ga.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean))

        if not (os.path.isfile(ga_file) and try_load):

            # load data required to compute reanalysis/model lapse rate
            ga_indata = {}
            file_vertconv = {}

            # variable names to load (load surface pressure first because it is required to convert the 3D variables to other vertical coordinate systems)
            varnames = varnames_ga(sim)

            for varname in varnames:
                varname_std = translate_varname(varname)
                ga_indata[varname_std] = load_var(sim, varname, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load) # load raw data (in pressure grid)

                if (translate_vardim(varname_std) == '3d' and vertcoord == 'si'): # convert 3D data to sigma coordinate
                    file_vertconv[varname_std] = remove_repdots('%s/%s.%s.%s.%s.pickle' % (datadir, varname, vertcoord, zonmean, timemean))

                    if (os.path.isfile(file_vertconv[varname_std]) and try_load): # load pickled data if available
                        ga_indata[varname_std] = (pickle.load(open(file_vertconv[varname_std], 'rb')))

                    else: # if not convert and save
                        ga_indata[varname_std] = pa_to_sigma(varname_std, ga_indata[varname_std], translate_varsfc(varname_std), ga_indata[translate_varsfc(varname_std)], ga_indata['ps'], si_std)

                        pickle.dump(ga_indata[translate_varname(varname)], open(remove_repdots('%s/%s.%s.%s.%s.pickle' % (datadir, varname, vertcoord, zonmean, timemean)), 'wb'))

            # compute lapse rate
            ga = make_ga(sim, ga_indata, model=model, vertcoord=vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            pickle.dump(ga, open(remove_repdots('%s/ga.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean)), 'wb'))

            ga_indata = None

        else:
            ga = pickle.load(open(ga_file, 'rb'))

        # location of pickled MALR data
        ga_m_file = remove_repdots('%s/ga_m.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean))

        if not (os.path.isfile(ga_m_file) and try_load):
            # load data required to compute MALR 
            varnames = varnames_ga_m(sim)

            file_ga_m = {}
            ga_m_indata = {}
            for varname in varnames:
                varname_std = translate_varname(varname)

                if translate_vardim(varname_std) == '3d':
                    file_ga_m[varname_std] = remove_repdots('%s/%s.%s.%s.%s.pickle' % (datadir, varname, vertcoord, zonmean, timemean))
                else:
                    file_ga_m[varname_std] = remove_repdots('%s/%s.%s.%s.pickle' % (datadir, varname, zonmean, timemean))

                ga_m_indata[varname_std] = (pickle.load(open(file_ga_m[varname_std], 'rb')))

            ga_m = make_ga_m(sim, ga_m_indata, model=model, vertcoord=vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            pickle.dump(ga_m, open(remove_repdots('%s/ga_m%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean)), 'wb'))

            ga_m_indata = None

        else:
            ga_m = pickle.load(open(ga_m_file, 'rb'))

        ga_dev = {}
        ga_dev['grid'] = ga['grid']
        ga_dev['ga_dev'] = 1e2 * (ga_m['ga_m'] - ga['ga'])/ga_m['ga_m'] # lapse rate deviation from moist adiabat (percent)

        # plot_lat_si_test(ga_dev, 'ga_dev')

        pickle.dump(ga_dev, open(remove_repdots('%s/ga_dev.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean)), 'wb'))

    else:
        ga_dev = pickle.load(open(ga_dev_file, 'rb'))

    return ga_dev


def make_ga_m(sim, ga_m_indata, **kwargs):
    # Computes the moist adiabatic lapse rate

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', 'zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled data if available
    file = remove_repdots('%s/ga_m.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean))

    ga_m = {}
    ga_m['ga_m'] = np.empty_like(ga_m_indata['ta']['ta'])
    ga_m['grid'] = ga_m_indata['ta']['grid']

    # rename data for easier access
    ta = ga_m_indata['ta']['ta']
    ps = ga_m_indata['ps']['ps']
    si = ga_m_indata['ta']['grid']['lev']

    ps3d = np.tile(ps, (si.size,1,1,1))
    ps3d = np.transpose(ps3d,(1,0,2,3))
    si3d = np.tile(si, (ta.shape[0],ta.shape[2],ta.shape[3],1))
    si3d = np.transpose(si3d,(0,3,1,2))
    pa = si3d*ps3d # compute pressure field in sigma grid
    pa = 1e-2*pa # convert pressure to mb

    si3d = None; ps3d = None

    # calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8 # K/km
    eps = 0.622
    R = 0.287 # J/g/K
    cp = 1.005 # J/g/K
    es0 = 6.11 # mb
    T0 = 273 # K
    L = 2510 - 2.38*(ta - T0)
    es = es0*np.exp(eps*L/R*(1/T0-1/ta))
    desdT = eps*L*es/(R*ta**2)

    ga_m['ga_m'] = dalr * (1+eps*L*es/(pa*R*ta))/(1+(eps*L/(cp*pa))*(desdT))

    return ga_m

def make_ga(sim, ga_indata, **kwargs):
    # Computes the reanalysis/model lapse rate
    
    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', 'zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled data if available
    file = remove_repdots('%s/ga.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean))

    ga = {}
    ga['ga'] = np.empty_like(ga_indata['ta']['ta'])
    ga['grid'] = ga_indata['ta']['grid']

    # compute lapse rate using centered finite difference (forward and backward Euler at the boundaries)
    ga['ga'][:,0,:,:] = (ga_indata['ta']['ta'][:,1,:,:]-ga_indata['ta']['ta'][:,0,:,:])/(ga_indata['zg']['zg'][:,1,:,:]-ga_indata['zg']['zg'][:,0,:,:])

    ga['ga'][:,1:-1,:,:] = (ga_indata['ta']['ta'][:,2:,:,:]-ga_indata['ta']['ta'][:,:-2,:,:])/(ga_indata['zg']['zg'][:,2:,:,:]-ga_indata['zg']['zg'][:,:-2,:,:])

    ga['ga'][:,-1,:,:] = (ga_indata['ta']['ta'][:,-1,:,:]-ga_indata['ta']['ta'][:,-2,:,:])/(ga_indata['zg']['zg'][:,-1,:,:]-ga_indata['zg']['zg'][:,-2,:,:])

    print(np.where(ga['ga']==ga['ga'].max()))
    print(ga_indata['ta']['ta'][6,:,1,304])
    print(ga_indata['zg']['zg'][6,:,1,304])

    ga['ga'] = -1e3 * ga['ga'] # convert K/m to K/km and use the sign convention lapse rate = - dT/dz

    pickle.dump(ga, open(remove_repdots('%s/ga.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean)), 'wb'))

    return ga

def load_var(sim, varname, **kwargs):

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', 'zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled data if available
    file = remove_repdots('%s/%s.%s.%s.pickle' % (datadir, varname, zonmean, timemean))

    if (os.path.isfile(file) and try_load):
        alldata = pickle.load(open(file, 'rb'))
    else:
        file_raw = {}
        grid = {}
        alldata = {}

        file_raw[varname] = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

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

        pickle.dump(alldata, open(remove_repdots('%s/%s.%s.%s.pickle' % (datadir, varname, zonmean, timemean)), 'wb'))

    return alldata

def plot_lat_si_test(ga, varname):
    import matplotlib.pyplot as plt

    [mesh_si, mesh_lat] = np.meshgrid(ga['grid']['lev'], ga['grid']['lat']) # create mesh

    fig, ax = plt.subplots()
    vmin = -50
    vmax = 50
    csf = ax.contourf(mesh_lat, mesh_si, np.transpose(np.squeeze(np.mean(ga[varname],axis=(0,3)))), np.arange(vmin,vmax,1), cmap='RdBu_r', vmin=vmin, vmax=vmax, extend='both')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_xticks(np.arange(-90,91,30))
    ax.set_ylabel('$\sigma$')
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.show()

def plot_mon_lat_test(ga, varname):
    import matplotlib.pyplot as plt

    print(ga[varname].shape[0])

    [mesh_lat, mesh_time] = np.meshgrid(ga['grid']['lat'], range(ga[varname].shape[0])) # create mesh

    fig, ax = plt.subplots()
    vmin = -50
    vmax = 50
    csf = ax.contourf(mesh_time, mesh_lat, np.squeeze(np.mean(ga[varname],axis=2)), np.arange(vmin,vmax,1), cmap='RdBu_r', vmin=vmin, vmax=vmax, extend='both')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_ylabel('Latitude (deg)')
    ax.set_yticks(np.arange(-90,91,30))
    # ax.set_xlabel('$\sigma$')
    # ax.set_xlim(ax.get_ylim()[::-1])
    plt.show()