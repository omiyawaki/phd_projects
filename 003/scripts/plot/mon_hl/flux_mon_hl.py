import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import *
from misc.load_data import load_flux
# from misc.translate import translate_varname
from misc.means import lat_mean, global_int
from misc import par
from proc.r1 import save_r1
from plot.titles import make_title_sim_time_lat
import os
import pickle
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def flux_mon_hl(sim, **kwargs):

    categ = 'mon_hl'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    latbnd = kwargs.get('latbnd', (80,90))
    latstep = kwargs.get('latstep', 0.25) # latitude step size used for interpolation
    legend = kwargs.get('legend', 0) # draw legend?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    sim_ref = kwargs.get('sim_ref', 'historical')
    timemean_ref = kwargs.get('timemean_ref', 'ymonmean-30')
    yr_span_ref = kwargs.get('yr_span_ref', '186001-200512')

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
        yr_base = 1979

    if latbnd[0] > 0: # NH
        if timemean == 'djfmean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
            vmin = -200
            vmax = 100
            vmin_dev = -50
            vmax_dev = 50
        elif timemean == 'jjamean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
            vmin = -150
            vmax = 50
            vmin_dev = -30
            vmax_dev = 30
    else: # SH
        if timemean == 'djfmean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
            vmin = -200
            vmax = 100
            vmin_dev = -50
            vmax_dev = 50
        elif timemean == 'jjamean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
            vmin = -150
            vmax = 50
            vmin_dev = -30
            vmax_dev = 30

    ##################################
    # LOAD DATA
    ##################################
    
    if isinstance(model, str):
        [flux, grid, datadir, plotdir, modelstr] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [flux, grid, datadir, plotdir, modelstr, flux_mmm] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    if refclim == 'hist-30':
        if isinstance(model, str):
            [flux_ref, _, _, _, _] = load_flux(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        else:
            [flux_ref, _, _, _, _, flux_ref_mmm] = load_flux(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        
        if timemean == 'djfmean':
            for fluxname in flux_ref:
                flux_ref[fluxname] = np.mean(np.roll(flux_ref[fluxname],1,axis=0)[0:3], 0)
        elif timemean == 'jjamean':
            for fluxname in flux_ref:
                flux_ref[fluxname] = np.mean(flux_ref[fluxname][5:8], 0)

        flux_ref_hl = {}

    ############################################
    # AVERAGE FLUXES ONLY AT HIGH LATITUDES
    ############################################
    flux_hl = {}
    flux_dev_hl = {}
    for fluxname in flux:
        flux_hl[fluxname] = lat_mean(flux[fluxname], grid, lat_int, dim=1)
        
        if refclim == 'init':
            flux_dev_hl[fluxname] = flux_hl[fluxname] - flux_hl[fluxname][0]
        else:
            flux_ref_hl[fluxname] = lat_mean(flux_ref[fluxname], grid, lat_int, dim=0)
            flux_dev_hl[fluxname] = flux_hl[fluxname] - flux_ref_hl[fluxname]

    time = yr_base + np.arange(flux_hl['ra'].shape[0]) # create time vector

    ############################################
    # PLOT
    ############################################
    plotname = remove_repdots('%s/flux_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))

    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ra = ax.plot(time, flux_hl['ra'], color='tab:gray', label='$R_a$')
    lp_stg_adv = ax.plot(time, flux_hl['stg_adv'], color='maroon', label='$\partial_t m + \partial_y (vm)$')
    lp_hfls = ax.plot(time, flux_hl['hfls'], color='tab:blue', label='LH')
    lp_hfss = ax.plot(time, flux_hl['hfss'], color='tab:orange', label='SH')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('Energy flux (Wm$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+flux_hl['ra'].shape[0]-1)
    ax.set_ylim(vmin,vmax)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/flux_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ra = ax.plot(time, flux_dev_hl['ra'], color='tab:gray', label='$\Delta R_a$')
    lp_stg_adv = ax.plot(time, flux_dev_hl['stg_adv'], color='maroon', label='$\Delta (\partial_t m + \partial_y (vm))$')
    lp_hfls = ax.plot(time, flux_dev_hl['hfls'], color='tab:blue', label='$\Delta$ LH')
    lp_hfss = ax.plot(time, flux_dev_hl['hfss'], color='tab:orange', label='$\Delta$ SH')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('$\Delta$ Energy flux (Wm$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+flux_dev_hl['ra'].shape[0]-1)
    ax.set_ylim(vmin_dev,vmax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()