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
from yaxis_r1_mon_hl import yaxis_r1_mon_hl
import os
import pickle
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def stf_mon_hl(sim, **kwargs):

    categ = 'mon_hl'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    latbnd = kwargs.get('latbnd', (80,90))
    latstep = kwargs.get('latstep', 0.25) # latitude step size used for interpolation
    plotover = kwargs.get('plotover', None) # plot overlay (sic for sea ice, ga_dev for lapse rate deviation, pr for precip, decomp for linear decomposition)?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    legend = kwargs.get('legend', 0)
    spread = kwargs.get('spread', 'prc')
    if plotover == 'ga_dev':
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
    ## Y AXIS SPECIFICATIONS
    ##########################################
    vmin, vmax = yaxis_r1_mon_hl(sim, timemean=timemean, latbnd=latbnd, plotover=plotover)

    ##################################
    # LOAD DATA
    ##################################
    if isinstance(model, str) or model is None:
        [flux, grid, datadir, plotdir, modelstr] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [flux, grid, datadir, plotdir, modelstr, flux_mmm] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    ############################################
    # PLOT
    ############################################
    flux_hl = {}
    for varname in flux:
        flux_hl[varname] = lat_mean(flux[varname], grid, lat_int, dim=1)

    if not (isinstance(model, str) or model is None):
        flux_mmm_hl = {}
        for varname in flux_mmm:
            flux_mmm_hl[varname] = {}
            for stat in flux_mmm[varname]:
                flux_mmm_hl[varname][stat] = lat_mean(flux_mmm[varname][stat], grid, lat_int, dim=1)

    time = yr_base + np.arange(flux_hl['dr1'].shape[0]) # create time vector

    plotname = '%s/dr1_dstf_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            ax.fill_between(time, flux_mmm_hl['dr1']['prc25'], flux_mmm_hl['dr1']['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(time, flux_mmm_hl['dr1']['mmm']-flux_mmm_hl['dr1']['std'], flux_mmm_hl['dr1']['mmm']+flux_mmm_hl['dr1']['std'], facecolor='k', alpha=0.2, edgecolor=None)

    lp_r1_dev = ax.plot(time, flux_hl['dr1'], color='black', label='$\Delta R_1$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('$\Delta R_1$ (unitless)')
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.set_xlim(time[0],time[-1])
    # ax.set_ylim(vmin['r1'],vmax['r1'])
    
    ############################################
    # STF
    ###########################################
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(time, flux_mmm_hl['dcstf']['prc25'], flux_mmm_hl['dcstf']['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(time, flux_mmm_hl['dcstf']['mmm']-flux_mmm_hl['dcstf']['std'], flux_mmm_hl['dcstf']['mmm']+flux_mmm_hl['dcstf']['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.plot(time, flux_hl['dcstf'], '-', color='tab:blue', label=r'$\frac{\Delta(LH+SH)}{\overline{R_a}}$')
    # ax.set_ylim(vmin['r1']-r1_hl[0],vmax['r1']-r1_hl[0])
    ax.tick_params(axis='y', labelcolor='black', color='black')
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # add legend
    ax.legend(loc='lower left')

    fig.set_size_inches(4,3)

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

