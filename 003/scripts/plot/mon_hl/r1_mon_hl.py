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

def r1_mon_hl(sim, **kwargs):

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
        if timemean == 'allseas':
            djf = {}; mam = {}; jja = {}; son = {};
            [djf_r1, djf_grid, djf_datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean='djfmean', try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
            [mam_r1, mam_grid, mam_datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean='mammean', try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
            [jja_r1, jja_grid, jja_datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean='jjamean', try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
            [son_r1, son_grid, son_datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean='sonmean', try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
        else:
            [r1, grid, datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
    else:
        if timemean == 'allseas':
            djf = {}; mam = {}; jja = {}; son = {};
            [djf_r1, djf_grid, djf_datadir, plotdir, modelstr, djf_r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean='djfmean', try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
            [mam_r1, mam_grid, mam_datadir, plotdir, modelstr, mam_r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean='mammean', try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
            [jja_r1, jja_grid, jja_datadir, plotdir, modelstr, jja_r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean='jjamean', try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
            [son_r1, son_grid, son_datadir, plotdir, modelstr, son_r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean='sonmean', try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
        else:
            [r1, grid, datadir, plotdir, modelstr, r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 

    if 'sic' in plotover:
        if isinstance(model, str) or model is None:
            [seaice, grid, datadir, plotdir, modelstr] = load_seaice(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        else:
            [seaice, grid, datadir, plotdir, modelstr, seaice_mmm] = load_seaice(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

        if refclim == 'hist-30':
            sim_ref = 'historical'
            timemean_ref = 'ymonmean-30'
            yr_span_ref = '186001-200512'
            if isinstance(model, str):
                [seaice_ref, _, _, _, _] = load_seaice(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
            else:
                [seaice_ref, _, _, _, _, seaice_ref_mmm] = load_seaice(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
            
            if timemean == 'djfmean':
                for seaicename in seaice_ref:
                    seaice_ref[seaicename] = np.mean(np.roll(seaice_ref[seaicename],1,axis=0)[0:3], 0)
            elif timemean == 'jjamean':
                for seaicename in seaice_ref:
                    seaice_ref[seaicename] = np.mean(seaice_ref[seaicename][5:8], 0)

    if 'pr' in plotover or 'cl' in plotover:
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

    if plotover == 'ga_dev':
        if isinstance(model, str) or model is None:
            [ga, grid_ga, datadir, plotdir, modelstr] = load_ga(sim, categ, vertbnd,vertcoord, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        else:
            [ga, grid_ga, datadir, plotdir, modelstr, ga_mmm] = load_ga(sim, categ, vertbnd, vertcoord, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

        if refclim == 'hist-30':
            sim_ref = 'historical'
            timemean_ref = 'ymonmean-30'
            yr_span_ref = '186001-200512'
            if isinstance(model, str):
                [ga_ref, _, _, _, _] = load_ga(sim_ref, categ, vertbnd, vertcoord, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
            else:
                [ga_ref, _, _, _, _, ga_ref_mmm] = load_ga(sim_ref, categ, vertbnd, vertcoord, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
            
            if timemean == 'djfmean':
                for ganame in ga_ref:
                    ga_ref[ganame] = np.mean(np.roll(ga_ref[ganame],1,axis=0)[0:3], 0)
            elif timemean == 'jjamean':
                for ganame in ga_ref:
                    ga_ref[ganame] = np.mean(ga_ref[ganame][5:8], 0)

    ############################################
    # AVERAGE R1 ONLY AT HIGH LATITUDES
    ############################################
    if timemean == 'allseas':
        djf_r1_hl = lat_mean(djf_r1, djf_grid, lat_int, dim=1)
        mam_r1_hl = lat_mean(mam_r1, mam_grid, lat_int, dim=1)
        jja_r1_hl = lat_mean(jja_r1, jja_grid, lat_int, dim=1)
        son_r1_hl = lat_mean(son_r1, son_grid, lat_int, dim=1)

        rolling_mean = 20; # smooth data using a rolling mean? (units: yr)
        djf_r1_hl = uniform_filter(djf_r1_hl, rolling_mean) # apply rolling mean
        mam_r1_hl = uniform_filter(mam_r1_hl, rolling_mean) # apply rolling mean
        jja_r1_hl = uniform_filter(jja_r1_hl, rolling_mean) # apply rolling mean
        son_r1_hl = uniform_filter(son_r1_hl, rolling_mean) # apply rolling mean

        r1_hl = djf_r1_hl

        djf_time = yr_base + np.arange(djf_r1_hl.shape[0]) # create time vector
        mam_time = yr_base + np.arange(mam_r1_hl.shape[0]) # create time vector
        jja_time = yr_base + np.arange(jja_r1_hl.shape[0]) # create time vector
        son_time = yr_base + np.arange(son_r1_hl.shape[0]) # create time vector
    else:
        if not ( isinstance(model, str) or (model is None) ):
            r1_hl_mmm = dict()
            for i in r1_mmm:
                r1_hl_mmm[i] = lat_mean(r1_mmm[i], grid, lat_int, dim=1)

        r1_hl = lat_mean(r1, grid, lat_int, dim=1)

        rolling_mean = 0; # smooth data using a rolling mean? (units: yr)
        r1_filt = uniform_filter(r1, [rolling_mean,0]) # apply rolling mean

        time = yr_base + np.arange(r1_hl.shape[0]) # create time vector

    ############################################
    # COMPUTE REGRESSION
    ############################################
    if 'ymonmean' not in timemean and sim == 'era5':
        A = np.vstack([time, np.ones(len(time))]).T
        m, c = np.linalg.lstsq(A, r1_hl, rcond=None)[0]

    ############################################
    # COMPUTE FIRST and LAST 30 YR AVGs (for yaxis settings)
    ############################################
    r1_f30 = np.mean(r1_hl[:30])
    r1_l30 = np.mean(r1_hl[-30:])

    ############################################
    # PLOT
    ############################################

    plotname = '%s/r1_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    # rae = patches.Rectangle((0,0.9),r1_hl.shape[0],vmax_r1-0.9, alpha=0.5)
    rae = patches.Rectangle((yr_base,0.9),yr_base + r1_hl.shape[0],vmax['r1']-0.9, alpha=0.5)
    ax.add_patch(rae)
    if timemean == 'allseas':
        lp_djf_r1 = ax.plot(djf_time, djf_r1_hl, color='navy', label='DJF')
        lp_mam_r1 = ax.plot(mam_time, mam_r1_hl, color='maroon', label='MAM')
        lp_jja_r1 = ax.plot(jja_time, jja_r1_hl, color='green', label='JJA')
        lp_son_r1 = ax.plot(son_time, son_r1_hl, color='orange', label='SON')
        ax.legend()
    else:
        if not ( isinstance(model, str) or (model is None) ):
            if spread == 'prc':
                spr_r1 = ax.fill_between(time, r1_hl_mmm['prc25'], r1_hl_mmm['prc75'], facecolor='black', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                spr_r1 = ax.fill_between(time, r1_hl_mmm['mmm'] - r1_hl_mmm['std'], r1_hl_mmm['mmm'] + r1_hl_mmm['std'], facecolor='black', alpha=0.2, edgecolor=None)

        lp_r1 = ax.plot(time, r1_hl, color='black')
    if 'ymonmean' not in timemean and sim == 'era5':
        # lp_trend = ax.plot(time, m*time + c, '--k', label='%g decade$^{-1}$' % (m*10))
        lp_trend = ax.plot(time, m*time + c, '--k', label='$R_1$ trend$ = %.1f$ %% decade$^{-1}$' % (m/np.nanmean(r1_hl) *1e3))
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    print(timemean)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('$R_1$ (unitless)')
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.set_xlim(yr_base,yr_base+r1_hl.shape[0]-1)
    ax.set_ylim(vmin['r1'],vmax['r1'])

    fig.set_size_inches(4, 3.5)

    if plotover == 'sic':
        ############################################
        # COMPARE WITH SEA ICE CONCENTRATION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)

        sic_hl_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in seaice_mmm['sic']:
                sic_hl_mmm[i] = lat_mean(seaice_mmm['sic'][i], grid, lat_int, dim=1)

        sic_hl = lat_mean(seaice['sic'], grid, lat_int, dim=1)

        # first and last 30 years
        sic_f30 = np.mean(sic_hl[:30])
        sic_l30 = np.mean(sic_hl[-30:])
        m_axis = (sic_l30 - sic_f30)/(r1_l30 - r1_f30)
        vmax_alg = sic_f30 + m_axis*( vmax['r1'] - r1_f30 )
        vmin_alg = sic_l30 + m_axis*( vmin['r1'] - r1_l30 )

        sax = ax.twinx()
        if timemean == 'jjamean':
            sax.plot(time, 100-sic_hl, color='tab:blue')
            sax.set_ylabel('Ocean fraction (%)', color='tab:blue')
        else:
            if not (isinstance(model, str) or model is None):
                if spread == 'prc':
                    sax.fill_between(time, sic_hl_mmm['prc25'], sic_hl_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
                elif spread == 'std':
                    sax.fill_between(time, sic_hl_mmm['mmm']-sic_hl_mmm['std'], sic_hl_mmm['mmm']+sic_hl_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
            sax.plot(time, sic_hl, color='tab:blue')
            sax.set_ylabel('Sea ice fraction (%)', color='tab:blue')
        # sax.set_ylim(vmin['sic'],vmax['sic'])
        sax.set_ylim(vmin_alg,vmax_alg)
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(MultipleLocator(5))

    elif plotover == 'ga_dev':
        ############################################
        # COMPARE WITH LAPSE RATE DEVIATION
        ###########################################
        # from proc.ga import make_ga_dev_vint

        plotname = '%s.%s.%g.%g' % (plotname, plotover, vertbnd[0], vertbnd[1])

        # ga_dev_vint = make_ga_dev_vint(sim, vertbnd, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

        # ga_dev_vint_hl = lat_mean(ga_dev_vint['ga_dev_vint'], ga_dev_vint['grid'], lat_int, dim=1)

        ga_dev_vint_hl_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in ga_mmm['ga_dev_vint']:
                ga_dev_vint_hl_mmm[i] = lat_mean(ga_mmm['ga_dev_vint'][i], grid, lat_int, dim=1)

        print(ga['ga_dev_vint'].shape)
        ga_dev_vint_hl = lat_mean(ga['ga_dev_vint'], grid_ga, lat_int, dim=1)
        print(ga_dev_vint_hl.shape)

        # first and last 30 years
        ga_dev_vint_f30 = np.mean(ga_dev_vint_hl[:30])
        ga_dev_vint_l30 = np.mean(ga_dev_vint_hl[-30:])
        m_axis = (ga_dev_vint_l30 - ga_dev_vint_f30)/(r1_l30 - r1_f30)
        vmax_alg = ga_dev_vint_f30 + m_axis*( vmax['r1'] - r1_f30 )
        vmin_alg = ga_dev_vint_l30 + m_axis*( vmin['r1'] - r1_l30 )

        sax = ax.twinx()
        if not (isinstance(model, str) or model is None):
            if spread == 'prc':
                sax.fill_between(time, ga_dev_vint_hl_mmm['prc25'], ga_dev_vint_hl_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                sax.fill_between(time, ga_dev_vint_hl_mmm['mmm']-ga_dev_vint_hl_mmm['std'], ga_dev_vint_hl_mmm['mmm']+ga_dev_vint_hl_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        sax.plot(time, ga_dev_vint_hl, color='tab:blue')
        sax.set_ylabel(r'$\langle(\Gamma_m-\Gamma)/\Gamma_m\rangle_{%0.1f}^{%0.1f}$ (%%)' % (vertbnd[0], vertbnd[1]), color='tab:blue')
        # sax.set_ylim(vmin['ga_dev'],vmax['ga_dev'])
        sax.set_ylim(vmin_alg,vmax_alg)
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(MultipleLocator(5))

    elif plotover == 'prfrac':
        ############################################
        # COMPARE WITH FRACTION OF LARGE SCALE PRECIPITATION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)

        hydro_ref_hl = {}

        # take mean in high latitudes
        pr_hl_mmm = dict()
        prl_hl_mmm = dict()
        prc_hl_mmm = dict()
        prfrac_hl_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in hydro_mmm['pr']:
                pr_hl_mmm[i] = lat_mean(hydro_mmm['pr'][i], grid, lat_int, dim=1)
                prl_hl_mmm[i] = lat_mean(hydro_mmm['prl'][i], grid, lat_int, dim=1)
                prc_hl_mmm[i] = lat_mean(hydro_mmm['prc'][i], grid, lat_int, dim=1)
                prfrac_hl_mmm[i] = lat_mean(hydro_mmm['prc'][i]/hydro_mmm['pr'][i], grid, lat_int, dim=1)

        prfrac_hl = lat_mean(hydro['prc']/hydro['pr'], grid, lat_int, dim=1)

        # first and last 30 years
        prfrac_f30 = np.mean(prfrac_hl[:30])
        prfrac_l30 = np.mean(prfrac_hl[-30:])
        m_axis = (prfrac_l30 - prfrac_f30)/(r1_l30 - r1_f30)
        vmin_alg = prfrac_f30 + m_axis*( vmax['r1'] - r1_f30 )
        vmax_alg = prfrac_l30 + m_axis*( vmin['r1'] - r1_l30 )

        # compute trends
        if 'ymonmean' not in timemean and sim == 'era5':
            A = np.vstack([time, np.ones(len(time))]).T
            m, c = np.linalg.lstsq(A, prfrac_hl, rcond=None)[0]

        sax = ax.twinx()
        if not ( isinstance(model, str) or (model is None) ):
            if spread == 'prc':
                sax.fill_between(time, prfrac_hl_mmm['prc25'], prfrac_hl_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                sax.fill_between(time, prfrac_hl_mmm['mmm']-prfrac_hl_mmm['std'], prfrac_hl_mmm['mmm']+prfrac_hl_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        sax.plot(time, prfrac_hl, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P_c/P$ trend$ = %.1f$ %% decade$^{-1}$' % (m*10))
        sax.set_ylabel(r'$P_c/P$ (unitless)', color='tab:blue')
        if model == 'IPSL-CM5A-LR':
            sax.set_ylim(vmin['prfrac'],vmax['prfrac'])
        else:
            sax.set_ylim(vmin_alg,vmax_alg)
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis
        # sax.set_ylim(sax.get_ylim()[::-1]) # invert r1 axis

    elif plotover == 'pr':
        ############################################
        # COMPARE WITH PRECIPITATION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)
        
        hydro_ref_hl = {}

        # take mean in high latitudes
        pr_hl = lat_mean(hydro['pr'], grid, lat_int, dim=1)

        # compute trends
        if 'ymonmean' not in timemean and sim == 'era5':
            A = np.vstack([time, np.ones(len(time))]).T
            m, c = np.linalg.lstsq(A, pr_hl, rcond=None)[0]

        sax = ax.twinx()
        sax.plot(time, pr_hl, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P$ trend$ = %.1f$ %% decade$^{-1}$' % (m/np.nanmean(pr_hl)*1e3))
        sax.set_ylabel(r'$P$ (mm d$^{-1}$)', color='tab:blue')
        sax.set_ylim(vmin['pr'],vmax['pr'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis

    elif plotover == 'prc':
        ############################################
        # COMPARE WITH CONVECTIVE PRECIPITATION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)
        
        hydro_ref_hl = {}

        # take mean in high latitudes
        prc_hl_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in hydro_mmm['prc']:
                prc_hl_mmm[i] = lat_mean(hydro_mmm['prc'][i], grid, lat_int, dim=1)

        prc_hl = lat_mean(hydro['prc'], grid, lat_int, dim=1)

        # compute trends
        if 'ymonmean' not in timemean and sim == 'era5':
            A = np.vstack([time, np.ones(len(time))]).T
            m, c = np.linalg.lstsq(A, prc_hl, rcond=None)[0]

        sax = ax.twinx()
        if not ( isinstance(model, str) or (model is None) ):
            if spread == 'prc':
                spr_r1 = sax.fill_between(time, prc_hl_mmm['prc25'], prc_hl_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                spr_r1 = sax.fill_between(time, prc_hl_mmm['mmm']-prc_hl_mmm['std'], prc_hl_mmm['mmm']+prc_hl_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        sax.plot(time, prc_hl, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P_c$ trend$ = %.1f$ %% decade$^{-1}$' % (m/np.nanmean(prc_hl)*1e3))
        sax.set_ylabel(r'$P_c$ (mm d$^{-1}$)', color='tab:blue')
        sax.set_ylim(vmin['prc'],vmax['prc'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis

    elif plotover == 'prl':
        ############################################
        # COMPARE WITH LARGE-SCALE PRECIPITATION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)

        pr_file = filenames_raw(sim, 'pr', model=model, timemean=timemean, yr_span=yr_span)
        prc_file = filenames_raw(sim, 'prc', model=model, timemean=timemean, yr_span=yr_span)
        if not zonmean:
            prl = pr_file.variables['pr'][:]-prc_file.variables['prc'][:]
        else:
            prl = np.mean(np.squeeze(pr_file.variables['pr'][:]-prc_file.variables['prc'][:]),2)

        prl = 86400*prl # convert kg m**-2 s**-1 to mm d**-1

        prl_hl = lat_mean(prl, grid, lat_int, dim=1)

        sax = ax.twinx()
        sax.plot(time, prl_hl, color='tab:blue')
        sax.set_ylabel(r'$P_l$ (mm d$^{-1}$)', color='tab:blue')
        sax.set_ylim(vmin['prl'],vmax['prl'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())

        ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis

    elif plotover == 'clt':
        ############################################
        # COMPARE WITH TOTAL CLOUD FRACTION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)
        
        hydro_ref_hl = {}

        # take mean in high latitudes
        clt_hl_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in hydro_mmm['clt']:
                clt_hl_mmm[i] = lat_mean(hydro_mmm['clt'][i], grid, lat_int, dim=1)

        clt_hl = lat_mean(hydro['clt'], grid, lat_int, dim=1)

        # compute trends
        if 'ymonmean' not in timemean and sim == 'era5':
            A = np.vstack([time, np.ones(len(time))]).T
            m, c = np.linalg.lstsq(A, clt_hl, rcond=None)[0]

        sax = ax.twinx()
        if not ( isinstance(model, str) or (model is None) ):
            if spread == 'prc':
                spr_r1 = sax.fill_between(time, clt_hl_mmm['prc25'], clt_hl_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                spr_r1 = sax.fill_between(time, clt_hl_mmm['mmm']-clt_hl_mmm['std'], clt_hl_mmm['mmm']+clt_hl_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        sax.plot(time, clt_hl, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P_c$ trend$ = %.1f$ %% decade$^{-1}$' % (m/np.nanmean(clt_hl)*1e3))
        sax.set_ylabel(r'Cloud fraction (%)', color='tab:blue')
        sax.set_ylim(vmin['clt'],vmax['clt'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis

    elif plotover == 'clwvi':
        ############################################
        # COMPARE WITH CLOUD WATER VAPOR (WATER AND ICE)
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)
        
        hydro_ref_hl = {}

        # take mean in high latitudes
        clwvi_hl_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in hydro_mmm['clwvi']:
                clwvi_hl_mmm[i] = lat_mean(hydro_mmm['clwvi'][i], grid, lat_int, dim=1)

        clwvi_hl = lat_mean(hydro['clwvi'], grid, lat_int, dim=1)

        # compute trends
        if 'ymonmean' not in timemean and sim == 'era5':
            A = np.vstack([time, np.ones(len(time))]).T
            m, c = np.linalg.lstsq(A, clwvi_hl, rcond=None)[0]

        sax = ax.twinx()
        if not ( isinstance(model, str) or (model is None) ):
            if spread == 'prc':
                spr_r1 = sax.fill_between(time, clwvi_hl_mmm['prc25'], clwvi_hl_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                spr_r1 = sax.fill_between(time, clwvi_hl_mmm['mmm']- clwvi_hl_mmm['std'], clwvi_hl_mmm['mmm']+clwvi_hl_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        sax.plot(time, clwvi_hl, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P_c$ trend$ = %.1f$ %% decade$^{-1}$' % (m/np.nanmean(clwvi_hl)*1e3))
        sax.set_ylabel(r'Condensed water path (kg m$^{-2}$)', color='tab:blue')
        sax.set_ylim(vmin['clwvi'],vmax['clwvi'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis

    elif plotover == 'clivi':
        ############################################
        # COMPARE WITH CLOUD WATER VAPOR (ICE ONLY)
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)
        
        hydro_ref_hl = {}

        # take mean in high latitudes
        clivi_hl_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in hydro_mmm['clivi']:
                clivi_hl_mmm[i] = lat_mean(hydro_mmm['clivi'][i], grid, lat_int, dim=1)

        clivi_hl = lat_mean(hydro['clivi'], grid, lat_int, dim=1)

        # compute trends
        if 'ymonmean' not in timemean and sim == 'era5':
            A = np.vstack([time, np.ones(len(time))]).T
            m, c = np.linalg.lstsq(A, clivi_hl, rcond=None)[0]

        sax = ax.twinx()
        if not ( isinstance(model, str) or (model is None) ):
            if spread == 'prc':
                spr_r1 = sax.fill_between(time, clivi_hl_mmm['prc25'], clivi_hl_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                spr_r1 = sax.fill_between(time, clivi_hl_mmm['mmm']-clivi_hl_mmm['std'], clivi_hl_mmm['mmm']+clivi_hl_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        sax.plot(time, clivi_hl, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P_c$ trend$ = %.1f$ %% decade$^{-1}$' % (m/np.nanmean(clivi_hl)*1e3))
        sax.set_ylabel(r'Condensed ice path (kg m$^{-2}$)', color='tab:blue')
        sax.set_ylim(vmin['clivi'],vmax['clivi'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis

    elif plotover == 'decomp':
        ############################################
        # DECOMPOSE INTO DYNAMIC AND RADIATIVE COMPONENTS
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)

        if isinstance(model, str) or model is None:
            [r1_dc, grid, datadir, plotdir, modelstr] = load_r1_dc(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span) 
        else:
            [r1_dc, grid, datadir, plotdir, modelstr, r1_dc_mmm] = load_r1_dc(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span) 

        # # location of pickled R1 seasonality data
        # r1_dc_file = remove_repdots('%s/r1_dc.%s.%s.pickle' % (datadir, zonmean, timemean))

        # if not (os.path.isfile(r1_dc_file) and try_load):
        #     save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        # [r1_dc, grid] = pickle.load(open(r1_dc_file, 'rb'))

        ############################################
        # AVERAGE R1 COMPONENTS ONLY AT HIGH LATITUDES
        ############################################
        r1_dc_hl = {}
        for varname in r1_dc:
            r1_dc_hl[varname] = lat_mean(r1_dc[varname], grid, lat_int, dim=1)

        if not (isinstance(model, str) or model is None):
            r1_dc_mmm_hl = {}
            for varname in r1_dc_mmm:
                r1_dc_mmm_hl[varname] = {}
                for stat in r1_dc_mmm[varname]:
                    r1_dc_mmm_hl[varname][stat] = lat_mean(r1_dc_mmm[varname][stat], grid, lat_int, dim=1)

        sax = ax.twinx()
        if not (isinstance(model, str) or model is None):
            if spread == 'prc':
                # sax.fill_between(time, r1_dc_mmm_hl['dr1']['prc25'], r1_dc_mmm_hl['dr1']['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
                sax.fill_between(time, r1_dc_mmm_hl['dyn']['prc25'], r1_dc_mmm_hl['dyn']['prc75'], facecolor='maroon', alpha=0.2, edgecolor=None)
                sax.fill_between(time, r1_dc_mmm_hl['rad']['prc25'], r1_dc_mmm_hl['rad']['prc75'], facecolor='tab:gray', alpha=0.3, edgecolor=None)
                sax.fill_between(time, r1_dc_mmm_hl['res']['prc25'], r1_dc_mmm_hl['res']['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                sax.fill_between(time, r1_dc_mmm_hl['dyn']['mmm']-r1_dc_mmm_hl['dyn']['std'], r1_dc_mmm_hl['dyn']['mmm']+r1_dc_mmm_hl['dyn']['std'], facecolor='maroon', alpha=0.2, edgecolor=None)
                sax.fill_between(time, r1_dc_mmm_hl['rad']['mmm']-r1_dc_mmm_hl['rad']['std'], r1_dc_mmm_hl['rad']['mmm']+r1_dc_mmm_hl['rad']['std'], facecolor='tab:gray', alpha=0.3, edgecolor=None)
                sax.fill_between(time, r1_dc_mmm_hl['res']['mmm']-r1_dc_mmm_hl['res']['std'], r1_dc_mmm_hl['res']['mmm']+r1_dc_mmm_hl['res']['std'], facecolor='k', alpha=0.2, edgecolor=None)
        sax.axhline(0, color='k', linewidth=0.5)
        sax.plot(time, r1_dc_hl['dr1'], color='k', label='$\Delta{R_1}$')
        sax.plot(time, r1_dc_hl['dyn'], color='maroon', label='$\overline{R_1} \dfrac{\Delta (\partial_t m + \partial_y (vm))}{\overline{\partial_t m + \partial_y (vm)}}$')
        sax.plot(time, r1_dc_hl['res'], '-.k', label='Residual')
        sax.plot(time, r1_dc_hl['rad'], color='lightgray', label='$-\overline{R_1} \dfrac{\Delta R_a}{\overline{R_a}}$')
        sax.set_ylabel(r'$\Delta R_1$ (unitless)', color='black')
        sax.set_ylim(vmin['r1']-r1_hl[0],vmax['r1']-r1_hl[0])
        sax.tick_params(axis='y', labelcolor='black', color='black')
        sax.yaxis.set_minor_locator(AutoMinorLocator())

        # add legend
        sax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.2), ncol=2)

        # cut off excess space on the bottom 
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                             box.width, box.height * 0.9])

        # alter figure aspect ratio to accomodate legend
        fig.set_size_inches(4,3.5)

    if legend and sim == 'era5':
        fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()
