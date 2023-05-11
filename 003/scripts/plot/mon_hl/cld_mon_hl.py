import sys
sys.path.append('/cltoject2/tas1/miyawaki/cltojects/003/scripts')
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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib
from plot.preamble import get_predata

def cld_mon_hl(sim, **kwargs):

    categ = 'mon_hl'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    latbnd = kwargs.get('latbnd', (80,90))
    latstep = kwargs.get('latstep', 0.25) # latitude step size used for interpolation
    plotover = kwargs.get('plotover', None) # plot overlay (sic for sea ice, ga_dev for lapse rate deviation, clt for cltecip, decomp for linear decomposition)?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    legend = kwargs.get('legend', 0)
    if plotover == 'ga_dev':
        vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for cltessure, z for height)
        # vertbnd = kwargs.get('vertbnd', (0.7, 0.3)) # sigma bounds of vertical integral
        vertbnd = kwargs.get('vertbnd', (1.0, 0.9)) # sigma bounds of vertical integral

    lat_int = np.arange(latbnd[0], latbnd[1], latstep)

    model, yr_span, yr_base, yr_span_ref, yr_base_show,yr_end_show = get_predata(sim, timemean, kwargs)

    ##########################################
    ## Y AXIS SPECIFICATIONS
    ##########################################
    vmin, vmax = yaxis_r1_mon_hl(sim, timemean=timemean, latbnd=latbnd, plotover=plotover)

    ##################################
    # LOAD DATA
    ##################################

    if isinstance(model, str) or model is None:
        [cld, grid, datadir, plotdir, modelstr] = load_cld(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [cld, grid, datadir, plotdir, modelstr, cld_mmm] = load_cld(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    if refclim == 'hist-30':
        sim_ref = 'historical'
        timemean_ref = 'ymonmean-30'
        if isinstance(model, str):
            [cld_ref, _, _, _, _] = load_cld(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        else:
            [cld_ref, _, _, _, _, cld_ref_mmm] = load_cld(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        
        if timemean == 'djfmean':
            for cldname in cld_ref:
                cld_ref[cldname] = np.mean(np.roll(cld_ref[cldname],1,axis=0)[0:3], 0)
        elif timemean == 'jjamean':
            for cldname in cld_ref:
                cld_ref[cldname] = np.mean(cld_ref[cldname][5:8], 0)

    ############################################
    # take mean in high latitudes
    ############################################

    cld_ref_hl = {}

    clt_hl = lat_mean(cld['clt'], grid, lat_int, dim=1)

    time = yr_base + np.arange(clt_hl.shape[0]) # create time vector

    # multimodel mean statistics if mmm
    if not ( isinstance(model, str) or (model is None) ):
        clt_hl_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in cld_mmm['clt']:
                clt_hl_mmm[i] = lat_mean(cld_mmm['clt'][i], grid, lat_int, dim=1)

    # compute trends
    if 'ymonmean' not in timemean and sim == 'era5':
        A = np.vstack([time, np.ones(len(time))]).T
        m, c = np.linalg.lstsq(A, clt_hl, rcond=None)[0]

    ############################################
    # PLOT CLIVI
    ############################################

    plotname = '%s/clivi_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    ax.axhline(0, time[0], time[-1], color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        ax.fill_between(time, clivi_hl_mmm['clivic25'], clivi_hl_mmm['clivic75'], facecolor='k', alpha=0.2, edgecolor=None)
        ax.fill_between(time, clivic_hl_mmm['clivic25'], clivic_hl_mmm['clivic75'], facecolor='tab:orange', alpha=0.2, edgecolor=None)
        ax.fill_between(time, clivil_hl_mmm['clivic25'], clivil_hl_mmm['clivic75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    ax.plot(time, clivi_hl, color='k')
    if 'ymonmean' not in timemean and sim == 'era5':
        # ax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
        ax.plot(time, m*time + c, '--', color='tab:blue', label='$P_l/P$ trend$ = %.1f$ %% decade$^{-1}$' % (m*10))
    ax.set_xlim(yr_base_show,yr_end_show)
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'Cloud ice (mm d$^{-1}$)')
    ax.set_ylim(vmin['clivi'],vmax['clivi'])
    ax.tick_params(axis='y')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    fig.set_size_inches(5, 4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT CLT
    ############################################

    plotname = '%s/clt_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    ax.axhline(0, time[0], time[-1], color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        ax.fill_between(time, clt_hl_mmm['cltc25'], clt_hl_mmm['cltc75'], facecolor='k', alpha=0.2, edgecolor=None)
        ax.fill_between(time, cltc_hl_mmm['cltc25'], cltc_hl_mmm['cltc75'], facecolor='tab:orange', alpha=0.2, edgecolor=None)
        ax.fill_between(time, cltl_hl_mmm['cltc25'], cltl_hl_mmm['cltc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    ax.plot(time, clt_hl, color='k')
    if 'ymonmean' not in timemean and sim == 'era5':
        # ax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
        ax.plot(time, m*time + c, '--', color='tab:blue', label='$P_l/P$ trend$ = %.1f$ %% decade$^{-1}$' % (m*10))
    ax.set_xlim(yr_base_show,yr_end_show)
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'Cloud area fraction (%)')
    ax.set_ylim(vmin['clt'],vmax['clt'])
    ax.tick_params(axis='y')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    fig.set_size_inches(5, 4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

