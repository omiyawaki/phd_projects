import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import *
from misc.load_data import load_flux, load_circ
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

def circ_mon_hl(sim, **kwargs):

    rolling_mean = 80; # smooth data using a rolling mean? (units: yr)

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
    elif sim == 'hist+rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '186001-229912')
        if 'ymonmean' not in timemean:
            yr_base = 1860
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
            vmin_dev = -30
            vmax_dev = 30
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
        [circ, grid, datadir, plotdir, modelstr] = load_circ(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [circ, grid, datadir, plotdir, modelstr, circ_mmm] = load_circ(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    if refclim == 'hist-30':
        if isinstance(model, str):
            [circ_ref, _, _, _, _] = load_circ(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        else:
            [circ_ref, _, _, _, _, circ_ref_mmm] = load_circ(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        
        if timemean == 'djfmean':
            for circname in circ_ref:
                circ_ref[circname] = np.mean(np.roll(circ_ref[circname],1,axis=0)[0:3], 0)
        elif timemean == 'jjamean':
            for circname in circ_ref:
                circ_ref[circname] = np.mean(circ_ref[circname][5:8], 0)

        circ_ref_hl = {}

    ############################################
    # EVALUATE FLUXES AT SPECIFIED LATITUDE BOUND
    ############################################
    # evaluate area of polar cap
    rlat = np.deg2rad(lat_int)
    clat = np.cos(rlat)
    cap_area = 2*np.pi*par.a**2*np.trapz(clat, rlat)

    circ_hl = {}
    circ_dev_hl = {}
    for circname in circ:
        circ_fint = interp1d(grid['lat'], circ[circname], kind='linear', axis=1)
        circ_hl[circname] = circ_fint(latbnd[0]) / cap_area
        
        if refclim == 'init':
            circ_dev_hl[circname] = circ_hl[circname] - circ_hl[circname][0]
        else:
            circ_ref_fint = interp1d(grid['lat'], circ_ref[circname], kind='linear', axis=0)
            circ_ref_hl[circname] = circ_ref_fint(latbnd[0]) / cap_area
            circ_dev_hl[circname] = circ_hl[circname] - circ_ref_hl[circname]

    time = yr_base + np.arange(circ_hl['wap'].shape[0]) # create time vector

    ############################################
    # PLOT (W, EVOLUTION OF ABS MAGNITUDE)
    ############################################
    plotname = remove_repdots('%s/circ_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_tot = ax.plot(time, circ_hl['wap'], color='k', label='$\omega(Pa s$^{-1}$)$')
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
    ax.set_xlim(yr_base,yr_base+circ_dev_hl['wap'].shape[0]-1)
    # ax.set_ylim(vmin_dev,vmax_dev)
    if legend:
        ax.legend()
    fig.set_size_inches(4,3.5)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (W, DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/circ_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_tot = ax.plot(time, circ_dev_hl['wap'], color='k', label='$\Delta\omega(Pa s$^{-1}$)$')
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
    ax.set_xlim(yr_base,yr_base+circ_dev_hl['wap'].shape[0]-1)
    # ax.set_ylim(vmin_dev,vmax_dev)
    if legend:
        ax.legend()
    fig.set_size_inches(4,3.5)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

