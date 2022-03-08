import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import *
from misc.load_data import load_flux, load_dyn
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

def dyn_mon_hl(sim, **kwargs):

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
    spread = kwargs.get('spread', 'std') # zonal mean?
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
        [dyn, grid, datadir, plotdir, modelstr] = load_dyn(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [dyn, grid, datadir, plotdir, modelstr, dyn_mmm] = load_dyn(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    # if refclim == 'hist-30':
    #     if isinstance(model, str):
    #         [dyn_ref, grid_ref, _, _, _] = load_dyn(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
    #     else:
    #         [dyn_ref, grid_ref, _, _, _, dyn_ref_mmm] = load_dyn(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        
    #     if timemean == 'djfmean':
    #         for dynname in dyn_ref:
    #             dyn_ref[dynname] = np.mean(np.roll(dyn_ref[dynname],1,axis=0)[0:3], 0)
    #             if not (isinstance(model, str) or model is None):
    #                 for stat in dyn_ref_mmm[dynname]:
    #                     dyn_ref_mmm[dynname][stat] = np.mean(np.roll(dyn_ref_mmm[dynname][stat],1,axis=0)[0:3], 0)

    #     elif timemean == 'jjamean':
    #         for dynname in dyn_ref:
    #             dyn_ref[dynname] = np.mean(dyn_ref[dynname][5:8], 0)

    #             if not (isinstance(model, str) or model is None):
    #                 for stat in dyn_ref_mmm[dynname]:
    #                     dyn_ref_mmm[dynname][stat] = np.mean(dyn_ref_mmm[dynname][stat][5:8], 0)


    #     dyn_ref_hl = {}
    #     dyn_ref_mmm_hl = {}

    ############################################
    # EVALUATE FLUXES AT SPECIFIED LATITUDE BOUND
    ############################################
    # evaluate area of polar cap
    rlat = np.deg2rad(lat_int)
    clat = np.cos(rlat)
    cap_area = 2*np.pi*par.a**2*np.trapz(clat, rlat)

    dyn_hl = {}
    dyn_dev_hl = {}
    dyn_mmm_hl = {}
    dyn_dev_mmm_hl = {}
    for dynname in dyn:
        dyn_fint = interp1d(grid['lat'], dyn[dynname], kind='linear', axis=1)
        dyn_hl[dynname] = dyn_fint(latbnd[0]) / cap_area

        if not (isinstance(model, str) or model is None):
            dyn_mmm_hl[dynname] = {}
            for stat in dyn_mmm[dynname]:
                dyn_mmm_fint = interp1d(grid['lat'], dyn_mmm[dynname][stat], kind='linear', axis=1)
                dyn_mmm_hl[dynname][stat] = dyn_mmm_fint(latbnd[0]) / cap_area
        
        # if refclim == 'init':
        #     dyn_dev_hl[dynname] = dyn_hl[dynname] - dyn_hl[dynname][0]

        #     if not (isinstance(model, str) or model is None):
        #         for stat in dyn_mmm[dynname]:
        #             dyn_dev_mmm_hl[dynname][stat] = dyn_mmm_hl[dynname][stat] - dyn_mmm_hl[dynname][stat][0]

        # else:
        #     dyn_ref_fint = interp1d(grid['lat'], dyn_ref[dynname], kind='linear', axis=0)
        #     dyn_ref_hl[dynname] = dyn_ref_fint(latbnd[0]) / cap_area
        #     dyn_dev_hl[dynname] = dyn_hl[dynname] - dyn_ref_hl[dynname]

        #     if not (isinstance(model, str) or model is None):
        #         dyn_ref_mmm_hl[dynname] = {}
        #         dyn_dev_mmm_hl[dynname] = {}
        #         for stat in dyn_mmm[dynname]:
        #             dyn_ref_mmm_fint = interp1d(grid_ref['lat'], dyn_ref_mmm[dynname][stat], kind='linear', axis=0)
        #             dyn_ref_mmm_hl[dynname][stat] = dyn_ref_mmm_fint(latbnd[0]) / cap_area
        #             dyn_dev_mmm_hl[dynname][stat] = dyn_mmm_hl[dynname][stat] - dyn_ref_mmm_hl[dynname][stat]

    time = yr_base + np.arange(dyn_hl['aht'].shape[0]) # create time vector

    ############################################
    # PLOT (SE TE DECOMP, ABS MAGNITUDE)
    ############################################
    plotname = remove_repdots('%s/dyn_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(time, dyn_mmm_hl['aht']['prc25'], dyn_mmm_hl['aht']['prc75'], facecolor='maroon', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['vmte']['prc25'], dyn_mmm_hl['vmte']['prc75'], facecolor='r', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['vmse']['prc25'], dyn_mmm_hl['vmse']['prc75'], facecolor='g', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['vmmmc']['prc25'], dyn_mmm_hl['vmmmc']['prc75'], facecolor='b', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(time, dyn_mmm_hl['aht']['mmm']-dyn_mmm_hl['aht']['std'], dyn_mmm_hl['aht']['mmm']+dyn_mmm_hl['aht']['std'], facecolor='maroon', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['vmte']['mmm']-dyn_mmm_hl['vmte']['std'], dyn_mmm_hl['vmte']['mmm']+dyn_mmm_hl['vmte']['std'], facecolor='r', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['vmse']['mmm']-dyn_mmm_hl['vmse']['std'], dyn_mmm_hl['vmse']['mmm']+dyn_mmm_hl['vmse']['std'], facecolor='g', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['vmmmc']['mmm']-dyn_mmm_hl['vmmmc']['std'], dyn_mmm_hl['vmmmc']['mmm']+dyn_mmm_hl['vmmmc']['std'], facecolor='b', alpha=0.2, edgecolor=None)
    lp_tot = ax.plot(time, dyn_hl['aht'], color='maroon', label='$-\Delta \partial_y(vm)$')
    lp_te = ax.plot(time, dyn_hl['vmte'], color='r', label='TE')
    lp_se = ax.plot(time, dyn_hl['vmse'], color='g', label='SE')
    lp_mmc = ax.plot(time, dyn_hl['vmmmc'], color='b', label='MMC')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('Energy flux convergence (Wm$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+dyn_hl['aht'].shape[0]-1)
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
    # PLOT (SE TE DECOMP, DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/dyn_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(time, dyn_mmm_hl['daht']['prc25'], dyn_mmm_hl['daht']['prc75'], facecolor='maroon', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['dvmte']['prc25'], dyn_mmm_hl['dvmte']['prc75'], facecolor='r', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['dvmse']['prc25'], dyn_mmm_hl['dvmse']['prc75'], facecolor='g', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['dvmmmc']['prc25'], dyn_mmm_hl['dvmmmc']['prc75'], facecolor='b', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(time, dyn_mmm_hl['daht']['mmm']-dyn_mmm_hl['daht']['std'], dyn_mmm_hl['daht']['mmm']+dyn_mmm_hl['daht']['std'], facecolor='maroon', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['dvmte']['mmm']-dyn_mmm_hl['dvmte']['std'], dyn_mmm_hl['dvmte']['mmm']+dyn_mmm_hl['dvmte']['std'], facecolor='r', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['dvmse']['mmm']-dyn_mmm_hl['dvmse']['std'], dyn_mmm_hl['dvmse']['mmm']+dyn_mmm_hl['dvmse']['std'], facecolor='g', alpha=0.2, edgecolor=None)
            ax.fill_between(time, dyn_mmm_hl['dvmmmc']['mmm']-dyn_mmm_hl['dvmmmc']['std'], dyn_mmm_hl['dvmmmc']['mmm']+dyn_mmm_hl['dvmmmc']['std'], facecolor='b', alpha=0.2, edgecolor=None)
    lp_tot = ax.plot(time, dyn_hl['daht'], color='maroon', label='$-\Delta \partial_y(vm)$')
    lp_te = ax.plot(time, dyn_hl['dvmte'], color='r', label='TE')
    lp_se = ax.plot(time, dyn_hl['dvmse'], color='g', label='SE')
    lp_mmc = ax.plot(time, dyn_hl['dvmmmc'], color='b', label='MMC')
    ax.plot(time, dyn_hl['daht'], color='maroon')
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
    ax.set_xlim(yr_base,yr_base+dyn_hl['aht'].shape[0]-1)
    if legend:
        ax.legend()
    fig.set_size_inches(4,3.5)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (TOTAL MSE FLUX DIVERGENCE CHANGE ONLY)
    ############################################
    plotname = remove_repdots('%s/dyn_dev_totonly_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_tot = ax.plot(time, dyn_hl['daht'], color='maroon', label='$-\Delta \partial_y(vm)$')
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
    ax.set_xlim(yr_base,yr_base+dyn_hl['aht'].shape[0]-1)
    plt.tight_layout()
    fig.set_size_inches(4,3.5)
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    if sim == 'hist+rcp85':
        #############################################
        ## SMOOTH DATA using RUNNING MEAN
        #############################################
        #dyn_hl_sm={}
        #dyn_mmm_hl_sm={}
        #dyn_dev_mmm_hl_sm={}
        #for varname in dyn_hl:
        #    dyn_hl_sm[varname] = uniform_filter(dyn_hl[varname], rolling_mean)

        #    if not (isinstance(model, str) or model is None):
        #        dyn_mmm_hl_sm[varname] = {}
        #        for stat in dyn_mmm_hl[varname]:
        #            dyn_mmm_hl_sm[varname][stat] = uniform_filter(dyn_mmm_hl[varname][stat], rolling_mean)

        # dyn_dev_hl_sm={}
        # for varname in dyn_dev_hl:
        #     dyn_dev_hl_sm[varname] = uniform_filter(dyn_dev_hl[varname], rolling_mean)

        #     if not (isinstance(model, str) or model is None):
        #         dyn_dev_mmm_hl_sm[varname] = {}
        #         for stat in dyn_dev_mmm_hl[varname]:
        #             dyn_dev_mmm_hl_sm[varname][stat] = uniform_filter(dyn_dev_mmm_hl[varname][stat], rolling_mean)

        ############################################
        # PLOT SMOOTH (SE TE DECOMP, ABS MAGNITUDE)
        ############################################
        plotname = remove_repdots('%s/dyn_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        if not (isinstance(model, str) or model is None):
            if spread == 'prc':
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['aht_sm']['prc25'][:-rolling_mean], dyn_mmm_hl['aht_sm']['prc75'][:-rolling_mean], facecolor='maroon', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['vmte_sm']['prc25'][:-rolling_mean], dyn_mmm_hl['vmte_sm']['prc75'][:-rolling_mean], facecolor='r', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['vmse_sm']['prc25'][:-rolling_mean], dyn_mmm_hl['vmse_sm']['prc75'][:-rolling_mean], facecolor='g', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['vmmmc_sm']['prc25'][:-rolling_mean], dyn_mmm_hl['vmmmc_sm']['prc75'][:-rolling_mean], facecolor='b', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['aht_sm']['mmm'][:-rolling_mean]-dyn_mmm_hl['aht_sm']['std'][:-rolling_mean], dyn_mmm_hl['aht_sm']['mmm'][:-rolling_mean]+dyn_mmm_hl['aht_sm']['std'][:-rolling_mean], facecolor='maroon', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['vmte_sm']['mmm'][:-rolling_mean]-dyn_mmm_hl['vmte_sm']['std'][:-rolling_mean], dyn_mmm_hl['vmte_sm']['mmm'][:-rolling_mean]+dyn_mmm_hl['vmte_sm']['std'][:-rolling_mean], facecolor='r', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['vmse_sm']['mmm'][:-rolling_mean]-dyn_mmm_hl['vmse_sm']['std'][:-rolling_mean], dyn_mmm_hl['vmse_sm']['mmm'][:-rolling_mean]+dyn_mmm_hl['vmse_sm']['std'][:-rolling_mean], facecolor='g', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['vmmmc_sm']['mmm'][:-rolling_mean]-dyn_mmm_hl['vmmmc_sm']['std'][:-rolling_mean], dyn_mmm_hl['vmmmc_sm']['mmm'][:-rolling_mean]+dyn_mmm_hl['vmmmc_sm']['std'][:-rolling_mean], facecolor='b', alpha=0.2, edgecolor=None)
        lp_tot = ax.plot(time[:-rolling_mean], dyn_hl['aht_sm'][:-rolling_mean], color='maroon', label='$-\Delta \partial_y(vm)$')
        lp_te = ax.plot(time[:-rolling_mean], dyn_hl['vmte_sm'][:-rolling_mean], color='r', label='TE')
        lp_se = ax.plot(time[:-rolling_mean], dyn_hl['vmse_sm'][:-rolling_mean], color='g', label='SE')
        lp_mmc = ax.plot(time[:-rolling_mean], dyn_hl['vmmmc_sm'][:-rolling_mean], color='b', label='MMC')
        ax.plot(time[:-rolling_mean], dyn_hl['aht_sm'][:-rolling_mean], color='maroon')
        make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
        ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
        if 'ymonmean' in timemean:
            ax.set_xticks(np.arange(0,12,1))
            ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
        else:
            ax.set_xlabel('Time (yr)')
        ax.set_ylabel('Energy flux convergence (Wm$^{-2}$)')
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(2006,2300)
        if legend:
            ax.legend()
        fig.set_size_inches(4,3.5)
        plt.tight_layout()
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

        ############################################
        # PLOT SMOOTH (SE TE DECOMP, DEVIATION FROM INITIAL)
        ############################################
        plotname = remove_repdots('%s/dyn_dev_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        if not (isinstance(model, str) or model is None):
            if spread == 'prc':
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['daht_sm']['prc25'][:-rolling_mean], dyn_mmm_hl['daht_sm']['prc75'][:-rolling_mean], facecolor='maroon', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['dvmte_sm']['prc25'][:-rolling_mean], dyn_mmm_hl['dvmte_sm']['prc75'][:-rolling_mean], facecolor='r', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['dvmse_sm']['prc25'][:-rolling_mean], dyn_mmm_hl['dvmse_sm']['prc75'][:-rolling_mean], facecolor='g', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['dvmmmc_sm']['prc25'][:-rolling_mean], dyn_mmm_hl['dvmmmc_sm']['prc75'][:-rolling_mean], facecolor='b', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['daht_sm']['mmm'][:-rolling_mean]-dyn_mmm_hl['daht_sm']['std'][:-rolling_mean], dyn_mmm_hl['daht_sm']['mmm'][:-rolling_mean]+dyn_mmm_hl['daht_sm']['std'][:-rolling_mean], facecolor='maroon', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['dvmte_sm']['mmm'][:-rolling_mean]-dyn_mmm_hl['dvmte_sm']['std'][:-rolling_mean], dyn_mmm_hl['dvmte_sm']['mmm'][:-rolling_mean]+dyn_mmm_hl['dvmte_sm']['std'][:-rolling_mean], facecolor='r', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['dvmse_sm']['mmm'][:-rolling_mean]-dyn_mmm_hl['dvmse_sm']['std'][:-rolling_mean], dyn_mmm_hl['dvmse_sm']['mmm'][:-rolling_mean]+dyn_mmm_hl['dvmse_sm']['std'][:-rolling_mean], facecolor='g', alpha=0.2, edgecolor=None)
                ax.fill_between(time[:-rolling_mean], dyn_mmm_hl['dvmmmc_sm']['mmm'][:-rolling_mean]-dyn_mmm_hl['dvmmmc_sm']['std'][:-rolling_mean], dyn_mmm_hl['dvmmmc_sm']['mmm'][:-rolling_mean]+dyn_mmm_hl['dvmmmc_sm']['std'][:-rolling_mean], facecolor='b', alpha=0.2, edgecolor=None)
        lp_tot = ax.plot(time[:-rolling_mean], dyn_hl['daht_sm'][:-rolling_mean], color='maroon', label='$-\Delta \partial_y(vm)$')
        lp_te = ax.plot(time[:-rolling_mean], dyn_hl['dvmte_sm'][:-rolling_mean], color='r', label='TE')
        lp_se = ax.plot(time[:-rolling_mean], dyn_hl['dvmse_sm'][:-rolling_mean], color='g', label='SE')
        lp_mmc = ax.plot(time[:-rolling_mean], dyn_hl['dvmmmc_sm'][:-rolling_mean], color='b', label='MMC')
        ax.plot(time[:-rolling_mean], dyn_hl['daht_sm'][:-rolling_mean], color='maroon')
        make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
        ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
        if 'ymonmean' in timemean:
            ax.set_xticks(np.arange(0,12,1))
            ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
        else:
            ax.set_xlabel('Time (yr)')
        ax.set_ylabel('$\Delta$ Energy flux convergence (Wm$^{-2}$)')
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(2006,2300)
        ax.set_ylim(vmin_dev,vmax_dev)
        if legend:
            ax.legend()
        fig.set_size_inches(4,3.5)
        plt.tight_layout()
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

        ############################################
        # PLOT SMOOTH (STATIONARY TRANSIENT DECOMP, DEVIATION FROM INITIAL)
        ############################################
        plotname = remove_repdots('%s/dyn_dev_st_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        lp_tot = ax.plot(time[:-rolling_mean], dyn_hl['daht_sm'][:-rolling_mean], color='maroon', label='$-\Delta \partial_y(vm)$')
        lp_te = ax.plot(time[:-rolling_mean], dyn_hl['dvmte_sm'][:-rolling_mean], color='r', label='Transient')
        lp_se = ax.plot(time[:-rolling_mean], dyn_hl['dvmse_sm'][:-rolling_mean]+dyn_hl['vmmmc_sm'][:-rolling_mean], color='tab:blue', label='Stationary')
        ax.plot(time[:-rolling_mean], dyn_hl['daht_sm'][:-rolling_mean], color='maroon')
        make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
        ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
        if 'ymonmean' in timemean:
            ax.set_xticks(np.arange(0,12,1))
            ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
        else:
            ax.set_xlabel('Time (yr)')
        ax.set_ylabel('$\Delta$ Energy flux convergence (Wm$^{-2}$)')
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(2006,2300)
        ax.set_ylim(vmin_dev,vmax_dev)
        if legend:
            ax.legend()
        fig.set_size_inches(4,3.5)
        plt.tight_layout()
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

        ############################################
        # PLOT (TOTAL MSE FLUX DIVERGENCE CHANGE ONLY)
        ############################################
        plotname = remove_repdots('%s/dyn_dev_totonly_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        lp_tot = ax.plot(time, dyn_hl['daht_sm'], color='maroon', label='$-\Delta \partial_y(vm)$')
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
        ax.set_xlim(2006,2200)
        ax.set_ylim(vmin_dev,vmax_dev)
        plt.tight_layout()
        fig.set_size_inches(4,3.5)
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

        ############################################
        # PLOT SMOOTH (AHT moist/dry decomp)
        ############################################
        plotname = remove_repdots('%s/aht_qs_dev_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        lp_tot = ax.plot(time, dyn_hl['daht_sm'], color='maroon', label='$-\Delta \partial_y(vm)$')
        lp_q = ax.plot(time, dyn_hl['dqaht_sm'], ':', color='maroon', label='$-\Delta \partial_y(vLq)$')
        lp_s = ax.plot(time, dyn_hl['dsaht_sm'], '--', color='maroon', label='$-\Delta \partial_y(vs)$')
        ax.plot(time, dyn_hl['daht_sm'], color='maroon')
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
        ax.set_xlim(2006,2200)
        ax.set_ylim(vmin_dev,vmax_dev)
        if legend:
            ax.legend()
        fig.set_size_inches(4,3.5)
        plt.tight_layout()
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

        ############################################
        # PLOT SMOOTH (MMC moist/dry decomp)
        ############################################
        plotname = remove_repdots('%s/mmc_qs_dev_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        lp_mmc = ax.plot(time, dyn_hl['dvmmmc_sm'], color='b', label='MMC (MSE)')
        lp_q = ax.plot(time, dyn_hl['dvqmmc_sm'], ':', color='b', label='MMC (latent)')
        lp_s = ax.plot(time, dyn_hl['dvsmmc_sm'], '--', color='b', label='MMC (DSE)')
        ax.plot(time, dyn_hl['dvmmmc_sm'], color='b')
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
        ax.set_xlim(2006,2200)
        ax.set_ylim(vmin_dev,vmax_dev)
        if legend:
            ax.legend()
        fig.set_size_inches(4,3.5)
        plt.tight_layout()
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

        ############################################
        # PLOT SMOOTH (SE moist/dry decomp)
        ############################################
        plotname = remove_repdots('%s/se_qs_dev_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        lp_se = ax.plot(time, dyn_hl['dvmse_sm'], color='g', label='SE (MSE)')
        lp_q = ax.plot(time, dyn_hl['dvqse_sm'], ':', color='g', label='SE (latent)')
        lp_s = ax.plot(time, dyn_hl['dvsse_sm'], '--', color='g', label='SE (DSE)')
        ax.plot(time, dyn_hl['dvmse_sm'], color='g')
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
        ax.set_xlim(2006,2200)
        ax.set_ylim(vmin_dev,vmax_dev)
        if legend:
            ax.legend()
        fig.set_size_inches(4,3.5)
        plt.tight_layout()
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

        ############################################
        # PLOT SMOOTH (STATIONARY MMC+SE moist/dry decomp)
        ############################################
        plotname = remove_repdots('%s/stationary_qs_dev_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        lp_se = ax.plot(time, dyn_hl['dvmmmc_sm']+dyn_hl['dvmse_sm'], color='tab:blue', label='Stationary (MSE)')
        lp_q = ax.plot(time, dyn_hl['dvqmmc_sm']+dyn_hl['dvqse_sm'], ':', color='tab:blue', label='Stationary (latent)')
        lp_s = ax.plot(time, dyn_hl['dvsmmc_sm']+dyn_hl['dvsse_sm'], '--', color='tab:blue', label='Stationary (DSE)')
        ax.plot(time, dyn_hl['dvmmmc_sm']+dyn_hl['dvmse_sm'], color='tab:blue')
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
        ax.set_xlim(2006,2200)
        ax.set_ylim(vmin_dev,vmax_dev)
        if legend:
            ax.legend()
        fig.set_size_inches(4,3.5)
        plt.tight_layout()
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

        ############################################
        # PLOT SMOOTH (TE moist/dry decomp)
        ############################################
        plotname = remove_repdots('%s/te_qs_dev_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        lp_te = ax.plot(time, dyn_hl['dvmte_sm'], color='r', label='TE (MSE)')
        lp_q = ax.plot(time, dyn_hl['dvqte_sm'], ':', color='r', label='TE (latent)')
        lp_s = ax.plot(time, dyn_hl['dvste_sm'], '--', color='r', label='TE (DSE)')
        ax.plot(time, dyn_hl['dvmte_sm'], color='r')
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
        ax.set_xlim(2006,2200)
        ax.set_ylim(vmin_dev,vmax_dev)
        if legend:
            ax.legend()
        fig.set_size_inches(4,3.5)
        plt.tight_layout()
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

        ############################################
        # PLOT SMOOTH (STATIONARY and TRANSIENT OVERLAY moist/dry decomp)
        ############################################
        plotname = remove_repdots('%s/st_qs_dev_sm_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
        fig, ax = plt.subplots()
        ax.axhline(0, color='k', linewidth=0.5)
        # lp_se = ax.plot(time, dyn_hl['dvmmmc_sm']+dyn_hl['dvmse_sm'], color='tab:blue', label='Stationary (MSE)')
        # lp_te = ax.plot(time, dyn_hl['dvmte_sm'], color='r', label='Transient (MSE)')
        lp_qs = ax.plot(time, dyn_hl['dvqmmc_sm']+dyn_hl['dvqse_sm'], ':', color='tab:blue', label='Stationary (latent)')
        lp_qt = ax.plot(time, dyn_hl['dvqte_sm'], ':', color='r', label='Transient (latent)')
        lp_ss = ax.plot(time, dyn_hl['dvsmmc_sm']+dyn_hl['dvsse_sm'], '--', color='tab:blue', label='Stationary (DSE)')
        lp_st = ax.plot(time, dyn_hl['dvste_sm'], '--', color='r', label='Transient (DSE)')
        # ax.plot(time, dyn_hl['dvmmmc_sm']+dyn_hl['dvmse_sm'], color='tab:blue')
        # ax.plot(time, dyn_hl['dvmte_sm'], color='r')
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
        ax.set_xlim(2006,2200)
        ax.set_ylim(vmin_dev,vmax_dev)
        if legend:
            ax.legend()
        fig.set_size_inches(4,3.5)
        plt.tight_layout()
        plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
        if viewplt:
            plt.show()
        plt.close()

