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

def make_dta_lev(sim, vertlev, **kwargs):
    # computes temperature response relative to first timestep

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', 'zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled vertically-integrated lapse rate deviation data
    dta_lev_file = remove_repdots('%s/dta_lev.%g.%s.%s.%s.pickle' % (datadir, vertlev, vertcoord, zonmean, timemean))

    if not (os.path.isfile(dta_lev_file) and try_load):

        dta = make_dta(sim, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

        t0 = start_time('Evaluating dta at level %g...' % (vertlev))

        f = interpolate.interp1d(dta['grid']['lev'], dta['dta'], axis=1)

        dta_lev = {}
        dta_lev['grid'] = {}
        dta_lev['grid']['lon'] = dta['grid']['lon']
        dta_lev['grid']['lat'] = dta['grid']['lat']
        dta_lev['dta_lev'] = f(vertlev)

        end_time(t0)

        t0 = start_time('Pickling dta at level %g...' % (vertlev))

        pickle.dump(dta_lev, open(remove_repdots('%s/dta_lev.%g.%s.%s.%s.pickle' % (datadir, vertlev, vertcoord, zonmean, timemean)), 'wb'))

        end_time(t0)
    else:
        t0 = start_time('Reading dta at level %g...' % (vertlev))

        dta_lev = pickle.load(open(dta_lev_file, 'rb'))

        end_time(t0)

    if zonmean:
        t0 = start_time('Computing zonal mean of dta at level %g...' % (vertlev))

        dta_lev['dta_lev'] = np.mean(dta_lev['dta_lev'],2)

        end_time(t0)

    return dta_lev


def make_dta(sim, **kwargs):
    # computes 3D temperature response
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
    dta_file = remove_repdots('%s/dta.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean))

    if not (os.path.isfile(dta_file) and try_load):

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

        t0 = start_time('Computing temperature anomaly...')

        dta = {}
        dta['grid'] = ta['grid']
        dta['dta'] = np.empty_like(ta['ta'])

        # compute temperature anomaly relative to first timestep
        dta['dta'] = ta['ta']-ta['ta'][0,:,:,:]

        end_time(t0)

        # plot_lat_si_test(dta, 'dta')

        t0 = start_time('Pickling temperature anomaly...')

        pickle.dump(dta, open(remove_repdots('%s/dta.%s.%s.%s.pickle' % (datadir, vertcoord, zonmean, timemean)), 'wb'))

        end_time(t0)

    else:
        t0 = start_time('Loading temperature anomaly...')

        dta = pickle.load(open(dta_file, 'rb'))

        end_time(t0)

    return dta