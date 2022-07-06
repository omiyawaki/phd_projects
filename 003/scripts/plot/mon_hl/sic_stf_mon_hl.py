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
from plot.preamble import get_predata
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

def sic_stf_mon_hl(sim, **kwargs):

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

    model, yr_span, yr_base = get_predata(sim, timemean, kwargs)

    ##########################################
    ## Y AXIS SPECIFICATIONS
    ##########################################
    vmin={}
    vmax={}
    vmin['stf'] = -10
    vmax['stf'] = 80
    vmin['ra'] = -160
    vmax['ra'] = -120
    vmin['stg_adv'] = -160
    vmax['stg_adv'] = -90

    ##################################
    # LOAD DATA
    ##################################
    if isinstance(model, str) or model is None:
        [flux, grid, datadir, plotdir, modelstr] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [flux, grid, datadir, plotdir, modelstr, flux_mmm] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    if isinstance(model, str) or model is None:
        [seaice, grid, datadir, plotdir, modelstr] = load_seaice(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [seaice, grid, datadir, plotdir, modelstr, seaice_mmm] = load_seaice(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    ############################################
    # AVERAGE R1 ONLY AT HIGH LATITUDES
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

    seaice_hl = {}
    for varname in seaice:
        seaice_hl[varname] = lat_mean(seaice[varname], grid, lat_int, dim=1)

    if not (isinstance(model, str) or model is None):
        seaice_mmm_hl = {}
        for varname in seaice_mmm:
            seaice_mmm_hl[varname] = {}
            for stat in seaice_mmm[varname]:
                seaice_mmm_hl[varname][stat] = lat_mean(seaice_mmm[varname][stat], grid, lat_int, dim=1)

    # pick out vars and sic
    stf = flux_hl['hfls'] + flux_hl['hfss']
    ra = flux_hl['ra']
    stg_adv = flux_hl['stg_adv']
    sic = seaice_hl['sic']

    stf_mmm = {}
    ra_mmm = {}
    stg_adv_mmm = {}
    sic_mmm = {}
    if not (isinstance(model, str) or model is None):
        for stat in flux_mmm_hl['hfls']:
            stf_mmm[stat] = flux_mmm_hl['hfls'][stat] + flux_mmm_hl['hfss'][stat]
            ra_mmm[stat] = flux_mmm_hl['ra'][stat]
            stg_adv_mmm[stat] = flux_mmm_hl['stg_adv'][stat]
        for stat in seaice_mmm_hl['sic']:
            sic_mmm[stat] = seaice_mmm_hl['sic'][stat]

    time = yr_base + np.arange(stf.shape[0]) # create time vector

    ############################################
    # COMPUTE FIRST and LAST 30 YR AVGs (for yaxis settings)
    ############################################
    stf_f30 = np.mean(stf[:30])
    stf_l30 = np.mean(stf[-30:])
    ra_f30 = np.mean(ra[:30])
    ra_l30 = np.mean(ra[-30:])
    stg_adv_f30 = np.mean(stg_adv[:30])
    stg_adv_l30 = np.mean(stg_adv[-30:])
    sic_f30 = np.mean(sic[:30])
    sic_l30 = np.mean(sic[-30:])

    ############################################
    # PLOT STF
    ############################################

    plotname = '%s/sic_stf_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            spr_r1 = ax.fill_between(time, stf_mmm['prc25'], stf_mmm['prc75'], facecolor='tab:orange', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            spr_r1 = ax.fill_between(time, stf_mmm['mmm'] - stf_mmm['std'], stf_mmm['mmm'] + stf_mmm['std'], facecolor='tab:orange', alpha=0.2, edgecolor=None)

    lp_stf = ax.plot(time, stf, '--', color='tab:orange')
    ax.plot(time, stf, ':', color='tab:blue')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('$LH+SH$ (W m$^{-2}$)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+stf.shape[0]-1)
    ax.set_ylim(vmin['stf'],vmax['stf'])

    fig.set_size_inches(4, 3.5)

    ############################################
    # COMPARE WITH SEA ICE CONCENTRATION
    ###########################################
    # first and last 30 years
    m_axis = (sic_l30 - sic_f30)/(stf_l30 - stf_f30)
    vmax_alg = sic_f30 + m_axis*( vmax['stf'] - stf_f30 )
    vmin_alg = sic_l30 + m_axis*( vmin['stf'] - stf_l30 )

    sax = ax.twinx()
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            sax.fill_between(time, sic_mmm['prc25'], sic_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            sax.fill_between(time, sic_mmm['mmm']-sic_mmm['std'], sic_mmm['mmm']+sic_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    sax.plot(time, sic, color='tab:blue')
    sax.set_ylabel('Sea ice fraction (%)', color='tab:blue')
    # sax.set_ylim(vmin['sic'],vmax['sic'])
    sax.set_ylim(vmin_alg,vmax_alg)
    sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
    sax.yaxis.set_minor_locator(MultipleLocator(5))

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT RA
    ############################################

    plotname = '%s/sic_ra_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            spr_r1 = ax.fill_between(time, ra_mmm['prc25'], ra_mmm['prc75'], facecolor='tab:gray', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            spr_r1 = ax.fill_between(time, ra_mmm['mmm'] - ra_mmm['std'], ra_mmm['mmm'] + ra_mmm['std'], facecolor='tab:gray', alpha=0.2, edgecolor=None)

    lp_ra = ax.plot(time, ra, '-', color='tab:gray')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('$R_a$ (W m$^{-2}$)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+ra.shape[0]-1)
    ax.set_ylim(vmin['ra'],vmax['ra'])

    fig.set_size_inches(4, 3.5)

    ############################################
    # COMPARE WITH SEA ICE CONCENTRATION
    ###########################################
    # first and last 30 years
    m_axis = (sic_l30 - sic_f30)/(ra_l30 - ra_f30)
    vmax_alg = sic_f30 + m_axis*( vmax['ra'] - ra_f30 )
    vmin_alg = sic_l30 + m_axis*( vmin['ra'] - ra_l30 )

    sax = ax.twinx()
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            sax.fill_between(time, sic_mmm['prc25'], sic_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            sax.fill_between(time, sic_mmm['mmm']-sic_mmm['std'], sic_mmm['mmm']+sic_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    sax.plot(time, sic, color='tab:blue')
    sax.set_ylabel('Sea ice fraction (%)', color='tab:blue')
    # sax.set_ylim(vmin['sic'],vmax['sic'])
    sax.set_ylim(vmin_alg,vmax_alg)
    sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
    sax.yaxis.set_minor_locator(MultipleLocator(5))

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT STGADV
    ############################################

    plotname = '%s/sic_stgadv_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            spr_r1 = ax.fill_between(time, stg_adv_mmm['prc25'], stg_adv_mmm['prc75'], facecolor='tab:red', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            spr_r1 = ax.fill_between(time, stg_adv_mmm['mmm'] - stg_adv_mmm['std'], stg_adv_mmm['mmm'] + stg_adv_mmm['std'], facecolor='tab:red', alpha=0.2, edgecolor=None)

    lp_stg_adv = ax.plot(time, stg_adv, '-', color='tab:red')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('$\partial_t m + \partial_y(vm)$ (W m$^{-2}$)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+stg_adv.shape[0]-1)
    ax.set_ylim(vmin['stg_adv'],vmax['stg_adv'])

    fig.set_size_inches(4, 3.5)

    ############################################
    # COMPARE WITH SEA ICE CONCENTRATION
    ###########################################
    # first and last 30 years
    m_axis = (sic_l30 - sic_f30)/(stg_adv_l30 - stg_adv_f30)
    vmax_alg = sic_f30 + m_axis*( vmax['stg_adv'] - stg_adv_f30 )
    vmin_alg = sic_l30 + m_axis*( vmin['stg_adv'] - stg_adv_l30 )

    sax = ax.twinx()
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            sax.fill_between(time, sic_mmm['prc25'], sic_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            sax.fill_between(time, sic_mmm['mmm']-sic_mmm['std'], sic_mmm['mmm']+sic_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    sax.plot(time, sic, color='tab:blue')
    sax.set_ylabel('Sea ice fraction (%)', color='tab:blue')
    # sax.set_ylim(vmin['sic'],vmax['sic'])
    sax.set_ylim(vmin_alg,vmax_alg)
    sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
    sax.yaxis.set_minor_locator(MultipleLocator(5))

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

