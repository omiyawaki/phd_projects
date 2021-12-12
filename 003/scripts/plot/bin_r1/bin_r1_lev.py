import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import *
from misc.filenames import *
# from misc.translate import translate_varname
from misc.means import lat_mean, global_int
from misc.load_data import *
from misc import par
from proc.r1 import save_r1
from plot.titles import make_title_sim_time_lat
from xaxis_bin_r1_lev import xaxis_bin_r1_lev
import os
import pickle
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def bin_r1_lev(sim, **kwargs):

    categ = 'bin_r1'

    annmean = kwargs.get('annmean', 0)
    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    latbnd = kwargs.get('latbnd', (80,90))
    latstep = kwargs.get('latstep', 0.25) # latitude step size used for interpolation
    plotvar = kwargs.get('plotvar', None) # plot overlay (sic for sea ice, ga_dev for lapse rate deviation, pr for precip, decomp for linear decomposition)?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    legend = kwargs.get('legend', 0)
    if plotvar == 'ga_dev':
        vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
        # vertbnd = kwargs.get('vertbnd', (0.7, 0.3)) # sigma bounds of vertical integral
        vertbnd = kwargs.get('vertbnd', (1.0, 0.9)) # sigma bounds of vertical integral

    lat_int = np.arange(latbnd[0], latbnd[1], latstep)

    if sim == 'longrun':
        model = kwargs.get('model', 'MPIESM12_abrupt4x')
        yr_span = kwargs.get('yr_span', '1000')
        yr_base = 0
    elif sim == 'rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '200601-230012')
        if 'ymonmean' not in timemean:
            yr_base = 2006
        else:
            yr_base = 0
    elif sim == 'historical':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '186001-200512')
        if 'ymonmean' not in timemean:
            yr_base = 1860
        else:
            yr_base = 0
    elif sim == 'echam':
        model = kwargs.get('model', 'rp000140')
        yr_span = kwargs.get('yr_span', '0001_0039')
        yr_base = 0
    elif sim == 'era5':
        model = None
        yr_span = kwargs.get('yr_span', '1979_2019')
        if 'ymonmean' not in timemean:
            yr_base = int(yr_span[0:4])
        else:
            yr_base = 0

    ##########################################
    ## X AXIS SPECIFICATIONS
    ##########################################
    vmin, vmax, xlabel = xaxis_bin_r1_lev(sim, timemean=timemean, latbnd=latbnd, plotvar=plotvar)

    ##################################
    # LOAD DATA
    ##################################
    if isinstance(model, str) or model is None:
        [r1, grid, datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
    else:
        [r1, grid, datadir, plotdir, modelstr, r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 

    # load lapse rate deviation
    if plotvar == 'ga_dev':
        from proc.ga import make_ga_dev

        ga_dev = make_ga_dev(sim, model=model, zonmean=zonmean, timemean=timemean, vertcoord=vertcoord, yr_span=yr_span, try_load=try_load)

        # save vertical grid
        lev = ga_dev['grid']['lev']

        # take zonal mean of the variable to be averaged over r1 bins
        binvar = np.mean(ga_dev['ga_dev'], 3)
        ga_dev = None

        # reshape to (time, lat, lev)
        binvar = np.transpose(binvar, [0,2,1])

    # take annual mean?
    if annmean:
        r1 = np.mean(r1, axis=0, keepdims=True)
        binvar = np.mean(binvar, axis=0, keepdims=True)

    ############################################
    # BIN VARIABLE ACCORDING TO R1
    ############################################

    # categorize r1 data into bins
    id_bins = np.digitize(r1, par.r1_bins)

    # create array to save binned averages
    binvar_avg = np.empty([len(par.r1_bins), binvar.shape[2]])

    # create cos(lat) data for area averaging
    clat = np.cos( np.deg2rad( np.tile(grid['lat'], (r1.shape[0], 1)) ) )

    for bin in range(len(par.r1_bins)):
        idx_bins = np.where( id_bins == bin)
        r1_avg = np.sum( clat[idx_bins]*r1[idx_bins] ) / np.sum(clat[idx_bins])
        binvar_avg[bin,:] = np.sum( np.expand_dims(clat[idx_bins], axis=1) * binvar[idx_bins[0],idx_bins[1],:], axis=0) / np.sum(clat[idx_bins])

    ############################################
    # PLOT
    ############################################

    plotname = '%s/%s.%s' % (plotdir, plotvar, timemean)

    fig, ax = plt.subplots()

    for bin in range(len(par.r1_bins)):
        print(par.r1_bins[bin])
        if par.r1_bins[bin] == 0.05:
            ax.plot(binvar_avg[bin,:], lev, color='tab:red')
        elif par.r1_bins[bin] == 0.85:
            ax.plot(binvar_avg[bin,:], lev, color='tab:blue')
        else:
            ax.plot(binvar_avg[bin,:], lev, color='black')

    # make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    # if 'ymonmean' in timemean:
    #     ax.set_xticks(np.arange(0,12,1))
    #     ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    # else:
    ax.set_xlabel(xlabel)
    #     ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('$\sigma$ (unitless)')
    # ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.set_xlim(vmin,vmax)
    ax.set_ylim([0.3, 1])
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()
