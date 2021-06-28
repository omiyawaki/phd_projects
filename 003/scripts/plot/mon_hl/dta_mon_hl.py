import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import *
# from misc.translate import translate_varname
from misc.means import lat_mean, global_int
from misc import par
from proc.dta import make_dta_lev
from proc.ga import make_ga_dev_vint
from plot.titles import make_title_sim_time_lat
import os
import pickle
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def dta_mon_hl(sim, **kwargs):

    categ = 'mon_hl'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    latbnd = kwargs.get('latbnd', (80,90))
    latstep = kwargs.get('latstep', 0.25) # latitude step size used for interpolation
    vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    vertlev = kwargs.get('vertlev', 1) # sigma level to evaluate temperature response
    # plotover = kwargs.get('plotover', 'sic') # plot overlay (sic for sea ice, ga_dev for lapse rate deviation)?
    # if plotover == 'ga_dev':
    #     vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    #     vertbnd = kwargs.get('vertbnd', (0.7, 0.3)) # sigma bounds of vertical integral

    lat_int = np.arange(latbnd[0], latbnd[1], latstep)

    if sim == 'longrun':
        model = kwargs.get('model', 'MPIESM12_abrupt4x')
        yr_span = kwargs.get('yr_span', '1000')
        yr_base = 0
    elif sim == 'rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '200601-230012')
        yr_base = 2006
    elif sim == 'echam':
        model = kwargs.get('model', 'rp000140')
        yr_span = kwargs.get('yr_span', '0001_0039')
        yr_base = 0
    elif sim == 'era5':
        model = None
        yr_span = kwargs.get('yr_span', '1980_2005')
        yr_base = 1980

    if latbnd[0] > 0: # NH
        if timemean == 'djfmean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
            vmin_dta_lev = 0.6
            vmax_dta_lev = 1.0
            vmin_sic = -20
            vmax_sic = 100
            vmin_ga_dev = -200
            vmax_ga_dev = 200
        elif timemean == 'jjamean':
            vmin_dta_lev = 0.825
            vmax_dta_lev = 1.0
            vmin_sic = 0
            vmax_sic = 120
        elif timemean == 'yearmean' or timemean == '':
            vmin_dta_lev = 0.7
            vmax_dta_lev = 0.95
            vmin_sic = -40
            vmax_sic = 100
            vmin_ga_dev = -200
            vmax_ga_dev = 175
    else: # SH
        if timemean == 'djfmean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
            vmin_dta_lev = 0.8 
            vmax_dta_lev = 1.4
            vmin_sic = -40
            vmax_sic = 100
            vmin_ga_dev = 50
            vmax_ga_dev = 350
        elif timemean == 'jjamean':
            vmin_dta_lev = 0.8 
            vmax_dta_lev = 1.4
            vmin_sic = -40
            vmax_sic = 100
            vmin_ga_dev = 50
            vmax_ga_dev = 350
        elif timemean == 'yearmean' or timemean == '':
            vmin_dta_lev = 0.8 
            vmax_dta_lev = 1.4
            vmin_sic = -40
            vmax_sic = 100
            vmin_ga_dev = 50
            vmax_ga_dev = 350


    # load data and plot directories
    datadir = get_datadir(sim, model=model, yr_span=yr_span)
    plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

    ############################################
    # LOAD DATA
    ############################################
    dta_lev = make_dta_lev(sim, vertlev, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

    # print(np.reshape(dTdt, (-1,96,12)).shape)
    if timemean == '':
        dta_lev['dta_lev'] = np.mean(np.reshape(dta_lev['dta_lev'], (-1,12,dta_lev['dta_lev'].shape[1])),1)

    ############################################
    # AVERAGE DT ONLY AT HIGH LATITUDES
    ############################################
    dta_lev_hl = lat_mean(dta_lev['dta_lev'], dta_lev['grid'], lat_int, dim=1)

    # rolling_mean = 0; # smooth data using a rolling mean? (units: yr)
    # dta_lev_filt = uniform_filter(dta_lev_hl, [rolling_mean,0]) # apply rolling mean

    time = np.arange(dta_lev_hl.shape[0]) # create time vector

    ############################################
    # PLOT
    ############################################
    plotname = '%s/dta_lev_mon_hl.%g.%g.%s.%g.%s' % (plotdir, latbnd[0], latbnd[1], vertcoord, vertlev, timemean)

    fig, ax = plt.subplots()
    lp_dta_lev = ax.plot(time, dta_lev_hl, color='black')
    make_title_sim_time_lat(ax, sim, model=model, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('$\Delta T_{%0.1f}$ (K)' % (vertlev))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(0,dta_lev_hl.shape[0])
    # ax.set_ylim(vmin_dta_lev,vmax_dta_lev)

    # if plotover == 'sic':
    #     ############################################
    #     # COMPARE WITH SEA ICE CONCENTRATION
    #     ###########################################
    #     plotname = '%s.%s' % (plotname, plotover)

    #     sic_file = filenames_raw(sim, 'sic', model=model, timemean=timemean, yr_span=yr_span)
    #     if not zonmean:
    #         sic = sic_file.variables['sic'][:]
    #     else:
    #         sic = np.mean(np.squeeze(sic_file.variables['sic'][:]),2)

    #     if sim == 'longrun':
    #         sic = 100*sic

    #     sic_hl = lat_mean(sic, grid, lat_int, dim=1)

    #     sax = ax.twinx()
    #     if timemean == 'jjamean':
    #         sax.plot(time, 100-sic_hl, color='tab:blue')
    #         sax.set_ylabel('Ocean fraction (%)', color='tab:blue')
    #     else:
    #         sax.plot(time, sic_hl, color='tab:blue')
    #         sax.set_ylabel('Sea ice fraction (%)', color='tab:blue')
    #     sax.set_ylim(vmin_sic,vmax_sic)
    #     sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
    #     sax.yaxis.set_minor_locator(MultipleLocator(5))
    # elif plotover == 'ga_dev':
    #     ############################################
    #     # COMPARE WITH LAPSE RATE DEVIATION
    #     ###########################################
    #     plotname = '%s.%s.%g.%g' % (plotname, plotover, vertbnd[0], vertbnd[1])

    #     ga_dev_vint = make_ga_dev_vint(sim, vertbnd, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

    #     ga_dev_vint_hl = lat_mean(ga_dev_vint['ga_dev_vint'], ga_dev_vint['grid'], lat_int, dim=1)

    #     sax = ax.twinx()
    #     sax.plot(time, ga_dev_vint_hl, color='tab:blue')
    #     sax.set_ylabel(r'$\langle(\Gamma_m-\Gamma)/\Gamma_m\rangle_{%0.1f}^{%0.1f}$ (%%)' % (vertbnd[0], vertbnd[1]), color='tab:blue')
    #     sax.set_ylim(vmin_ga_dev,vmax_ga_dev)
    #     sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
    #     sax.yaxis.set_minor_locator(MultipleLocator(5))

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.show()
    # tikzplotlib.save('%s.tex' % (plotname))
