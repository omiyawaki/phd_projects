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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def pr_mon_hl(sim, **kwargs):

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
        [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    if refclim == 'hist-30':
        sim_ref = 'historical'
        timemean_ref = 'ymonmean-30'
        yr_span_ref = '186001-200512'
        if isinstance(model, str):
            [hydro_ref, _, _, _, _] = load_hydro(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        else:
            [hydro_ref, _, _, _, _, hydro_ref_mmm] = load_hydro(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        
        if timemean == 'djfmean':
            for hydroname in hydro_ref:
                hydro_ref[hydroname] = np.mean(np.roll(hydro_ref[hydroname],1,axis=0)[0:3], 0)
        elif timemean == 'jjamean':
            for hydroname in hydro_ref:
                hydro_ref[hydroname] = np.mean(hydro_ref[hydroname][5:8], 0)

    ############################################
    # take mean in high latitudes
    ############################################

    hydro_ref_hl = {}

    pr_hl = lat_mean(hydro['pr'], grid, lat_int, dim=1)
    prc_hl = lat_mean(hydro['prc'], grid, lat_int, dim=1)
    prl_hl = lat_mean(hydro['prl'], grid, lat_int, dim=1)
    # prw_hl = lat_mean(hydro['prw'], grid, lat_int, dim=1)
    vhur_hl = lat_mean(hydro['vhur'], grid, lat_int, dim=1)

    time = yr_base + np.arange(pr_hl.shape[0]) # create time vector
    prfrac_hl = lat_mean(100*hydro['prl']/hydro['pr'], grid, lat_int, dim=1)

    # multimodel mean statistics if mmm
    if not ( isinstance(model, str) or (model is None) ):
        pr_hl_mmm = dict()
        prl_hl_mmm = dict()
        prc_hl_mmm = dict()
        # prw_hl_mmm = dict()
        prfrac_hl_mmm = dict()
        vhur_hl_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in hydro_mmm['pr']:
                pr_hl_mmm[i] = lat_mean(hydro_mmm['pr'][i], grid, lat_int, dim=1)
                prl_hl_mmm[i] = lat_mean(hydro_mmm['prl'][i], grid, lat_int, dim=1)
                prc_hl_mmm[i] = lat_mean(hydro_mmm['prc'][i], grid, lat_int, dim=1)
                # prw_hl_mmm[i] = lat_mean(hydro_mmm['prw'][i], grid, lat_int, dim=1)
                prfrac_hl_mmm[i] = lat_mean(100*hydro_mmm['prl'][i]/hydro_mmm['pr'][i], grid, lat_int, dim=1)
                vhur_hl_mmm[i] = lat_mean(hydro_mmm['vhur'][i], grid, lat_int, dim=1)

    # compute trends
    if 'ymonmean' not in timemean and sim == 'era5':
        A = np.vstack([time, np.ones(len(time))]).T
        m, c = np.linalg.lstsq(A, prfrac_hl, rcond=None)[0]

    ############################################
    # PLOT VHUR
    ############################################

    plotname = '%s/vhur_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    ax.axhline(0, time[0], time[-1], color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        ax.fill_between(time, vhur_hl_mmm['prc25'], vhur_hl_mmm['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
    ax.plot(time, vhur_hl, color='k')
    ax.set_xlim(yr_base,yr_base+vhur_hl.shape[0]-1)
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\langle \mathrm{RH} \rangle$ (%)')
    ax.set_ylim(vmin['vhur'],vmax['vhur'])
    ax.tick_params(axis='y', top=True, right=True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    fig.set_size_inches(5, 4)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT PR
    ############################################

    plotname = '%s/pr_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    ax.axhline(0, time[0], time[-1], color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        ax.fill_between(time, pr_hl_mmm['prc25'], pr_hl_mmm['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
        ax.fill_between(time, prc_hl_mmm['prc25'], prc_hl_mmm['prc75'], facecolor='tab:orange', alpha=0.2, edgecolor=None)
        ax.fill_between(time, prl_hl_mmm['prc25'], prl_hl_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    ax.plot(time, pr_hl, color='k', label='$P_{l}+P_{c}$')
    ax.plot(time, prc_hl, color='tab:orange', label='$P_c$')
    ax.plot(time, prl_hl, color='tab:blue', label='$P_l$')
    if 'ymonmean' not in timemean and sim == 'era5':
        # ax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
        ax.plot(time, m*time + c, '--', color='tab:blue', label='$P_l/P$ trend$ = %.1f$ %% decade$^{-1}$' % (m*10))
    ax.set_xlim(yr_base,yr_base+pr_hl.shape[0]-1)
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'Precipitation (mm d$^{-1}$)')
    ax.set_ylim(vmin['pr'],vmax['pr'])
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

    #############################################
    ## PLOT PRW
    #############################################

    #plotname = '%s/prw_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    #fig, ax = plt.subplots()
    #ax.axhline(0, time[0], time[-1], color='k', linewidth=0.5)
    #if not ( isinstance(model, str) or (model is None) ):
    #    ax.fill_between(time, prw_hl_mmm['prc25'], prw_hl_mmm['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
    #ax.plot(time, prw_hl, color='k', label='Precipitable water')
    #if 'ymonmean' not in timemean and sim == 'era5':
    #    # ax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
    #    ax.plot(time, m*time + c, '--', color='tab:blue', label='$P_l/P$ trend$ = %.1f$ %% decade$^{-1}$' % (m*10))
    #ax.set_xlim(yr_base,yr_base+pr_hl.shape[0]-1)
    #ax.set_xlabel('Time (yr)')
    #ax.set_ylabel(r'Precipitable water (kg m$^{-2}$)')
    ## ax.set_ylim(vmin['pr'],vmax['pr'])
    #ax.tick_params(axis='y')
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())

    #fig.set_size_inches(5, 4)
    ## plt.legend()
    #plt.tight_layout()
    #plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    #if viewplt:
    #    plt.show()
    #plt.close()
