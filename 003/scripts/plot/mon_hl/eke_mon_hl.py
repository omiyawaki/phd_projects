import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import *
from misc.load_data import load_flux, load_ke
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

def eke_mon_hl(sim, **kwargs):

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
            vmin = 0.5e6
            vmax = 1.5e6
            vmin_dev = -1e5
            vmax_dev = 1e5
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
    
    if isinstance(model, str) or model is None:
        [ke, grid, datadir, plotdir, modelstr] = load_ke(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [ke, grid, datadir, plotdir, modelstr, ke_mmm] = load_ke(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    if refclim == 'hist-30':
        if isinstance(model, str):
            [ke_ref, _, _, _, _] = load_ke(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        else:
            [ke_ref, _, _, _, _, ke_ref_mmm] = load_ke(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        
        if timemean == 'djfmean':
            for kename in ke_ref:
                ke_ref[kename] = np.mean(np.roll(ke_ref[kename],1,axis=0)[0:3], 0)
        elif timemean == 'jjamean':
            for kename in ke_ref:
                ke_ref[kename] = np.mean(ke_ref[kename][5:8], 0)

        ke_ref_hl = {}

    ############################################
    # AVERAGE FLUXES ONLY AT HIGH LATITUDES
    ############################################
    ke_hl = {}
    ke_dev_hl = {}
    for kename in ke:
        ke_hl[kename] = lat_mean(ke[kename], grid, lat_int, dim=1)
        
        if refclim == 'init':
            ke_dev_hl[kename] = ke_hl[kename] - ke_hl[kename][0]
        else:
            ke_ref_hl[kename] = lat_mean(ke_ref[kename], grid, lat_int, dim=0)
            ke_dev_hl[kename] = ke_hl[kename] - ke_ref_hl[kename]

    time = yr_base + np.arange(ke_hl['tke'].shape[0]) # create time vector

    ############################################
    # PLOT (VEKE, TIME SERIES)
    ############################################

    plotname = remove_repdots('%s/ke_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_tke = ax.plot(time, 1e-6*ke_hl['tke'], color='k', label='Total')
    lp_teke = ax.plot(time, 1e-6*ke_hl['teke'], color='r', label='Transient Eddies')
    lp_seke = ax.plot(time, 1e-6*ke_hl['seke'], color='b', label='Stationary Eddies')
    lp_mke = ax.plot(time, 1e-6*ke_hl['mke'], color='g', label='Mean Meridional Circulation')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('Vertically Integrated Kinetic Energy (MJ m$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+ke_dev_hl['tke'].shape[0]-1)
    # ax.set_ylim(vmin_dev,vmax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    print(plotdir)

    ############################################
    # PLOT (VEKE, DEVIATION FROM INITIAL)
    ############################################

    plotname = remove_repdots('%s/ke_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_tke = ax.plot(time, 1e-6*ke_dev_hl['tke'], color='k', label='Total')
    lp_teke = ax.plot(time, 1e-6*ke_dev_hl['teke'], color='r', label='Transient Eddies')
    lp_seke = ax.plot(time, 1e-6*ke_dev_hl['seke'], color='b', label='Stationary Eddies')
    lp_mke = ax.plot(time, 1e-6*ke_dev_hl['mke'], color='g', label='Mean Meridional Circulation')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('$\Delta$ Vertically Integrated Kinetic Energy (MJ m$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+ke_dev_hl['tke'].shape[0]-1)
    # ax.set_ylim(vmin_dev,vmax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    print(plotdir)

