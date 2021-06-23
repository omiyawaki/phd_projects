import sys
import os
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.translate import *
from misc.dirnames import get_datadir
from misc.filenames import *
from misc.load_data import *
from misc.vertconv import pa_to_sigma
from misc.timeutil import *
from misc import par
import numpy as np
from scipy import interpolate
import pickle
from netCDF4 import Dataset

def varnames_ta(sim):
    if sim == 'echam':
        varnames = ['sp', 'ts', 't', 'zs', 'z']
    elif sim == 'era5':
        varnames = ['sp', 't2m', 't', 'zs', 'z']
    else:
        varnames = ['ps', 'tas', 'ta', 'orog', 'zg']

    return varnames

def timestep(timemean):
    if timemean in ['yearmean', 'djfmean', 'jjamean']:
        dt = 365*86400 # seconds in a year
    elif timemean == '':
        dt = 365*86400/12 # seconds in a month
    return dt

def make_dTdt_lev(sim, vertlev, **kwargs):
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
    dTdt_lev_file = remove_repdots('%s/dTdt_lev.%g.%s.%s.%s.pickle' % (datadir, vertlev, vertcoord, zonmean, timemean))

    if not (os.path.isfile(dTdt_lev_file) and try_load):

        dTdt = make_dTdt(sim, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

        t0 = start_time('Evaluating dTdt at level %g...' % (vertlev))

        f = interpolate.interp1d(dTdt['grid']['lev'], dTdt['dTdt'], axis=1)

        dTdt_lev = {}
        dTdt_lev['grid'] = {}
        dTdt_lev['grid']['lon'] = dTdt['grid']['lon']
        dTdt_lev['grid']['lat'] = dTdt['grid']['lat']
        dTdt_lev['dTdt_lev'] = f(vertlev)

        end_time(t0)

        t0 = start_time('Pickling dTdt at level %g...' % (vertlev))

        pickle.dump(dTdt_lev, open(remove_repdots('%s/dTdt_lev.%g.%s.%s.%s.pickle' % (datadir, vertlev, vertcoord, zonmean, timemean)), 'wb'))

        end_time(t0)
    else:
        t0 = start_time('Reading dTdt at level %g...' % (vertlev))

        dTdt_lev = pickle.load(open(dTdt_lev_file, 'rb'))

        end_time(t0)

    if zonmean:
        t0 = start_time('Computing zonal mean of dTdt at level %g...' % (vertlev))

        dTdt_lev['dTdt_lev'] = np.mean(dTdt_lev['dTdt_lev'],2)

        end_time(t0)

    return dTdt_lev


def make_dTdt(sim, **kwargs):
    # computes 3D temperature tendency
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
    dTdt_file = remove_repdots('%s/dTdt.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean))

    if not (os.path.isfile(dTdt_file) and try_load):

        # location of pickled lapse rate data
        ta_file = remove_repdots('%s/ta.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean))

        if not (os.path.isfile(ta_file) and try_load):

            # load data required to compute reanalysis/model lapse rate
            ta_indata = {}
            file_vertconv = {}

            # variable names to load (load surface pressure first because it is required to convert the 3D variables to other vertical coordinate systems)
            varnames = varnames_ta(sim)

            for varname in varnames:
                varname_std = translate_varname(varname)
                ta_indata[varname_std] = load_var(sim, varname, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load) # load raw data (in pressure grid)

                if (translate_vardim(varname_std) == '3d' and vertcoord == 'si'): # convert 3D data to sigma coordinate
                    file_vertconv[varname_std] = remove_repdots('%s/%s.%s.%s.%s.pickle' % (datadir, varname, vertcoord, zonmean, timemean))

                    if (os.path.isfile(file_vertconv[varname_std]) and try_load): # load pickled data if available
                        ta_indata[varname_std] = (pickle.load(open(file_vertconv[varname_std], 'rb')))

                    else: # if not convert and save
                        ta_indata[varname_std] = pa_to_sigma(varname_std, ta_indata[varname_std], translate_varsfc(varname_std), ta_indata[translate_varsfc(varname_std)], ta_indata['ps'], si_std)

                        pickle.dump(ta_indata[translate_varname(varname)], open(remove_repdots('%s/%s.%s.%s.%s.pickle' % (datadir, varname, vertcoord, zonmean, timemean)), 'wb'), protocol=4)

        else:
            ta = pickle.load(open(ta_file, 'rb'))

        t0 = start_time('Computing temperature tendency...')

        dTdt = {}
        dTdt['grid'] = ta['grid']
        dTdt['dTdt'] = np.empty_like(ta['ta'])

        dt = timestep(timemean)

        # compute lapse rate using centered finite difference (forward and backward Euler at the boundaries)
        dTdt['dTdt'][0,:,:,:] = (ta['ta'][1,:,:,:]-ta['ta'][0,:,:,:])/dt
        dTdt['dTdt'][1:-1,:,:,:] = (ta['ta'][2:,:,:,:]-ta['ta'][:-2,:,:,:])/dt
        dTdt['dTdt'][-1,:,:,:] = (ta['ta'][-1,:,:,:]-ta['ta'][-2,:,:,:])/dt

        # convert K/s to K/decade
        dTdt['dTdt'] = 365*86400*10* dTdt['dTdt']

        end_time(t0)

        # plot_lat_si_test(dTdt, 'dTdt')

        t0 = start_time('Pickling temperature tendency...')

        pickle.dump(dTdt, open(remove_repdots('%s/dTdt.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean)), 'wb'))

        end_time(t0)

    else:
        t0 = start_time('Loading temperature tendency...')

        dTdt = pickle.load(open(dTdt_file, 'rb'))

        end_time(t0)

    return dTdt

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
    vmin = -2e-7
    vmax = 2e-7
    csf = ax.contourf(mesh_time, mesh_lat, np.squeeze(np.mean(ga[varname],axis=2)), np.linspace(vmin,vmax,50), cmap='RdBu_r', vmin=vmin, vmax=vmax, extend='both')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_ylabel('Latitude (deg)')
    ax.set_yticks(np.arange(-90,91,30))
    # ax.set_xlabel('$\sigma$')
    # ax.set_xlim(ax.get_ylim()[::-1])
    plt.show()