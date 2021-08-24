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
    if plotover == 'ga_dev':
        vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
        vertbnd = kwargs.get('vertbnd', (0.7, 0.3)) # sigma bounds of vertical integral

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
        [r1, grid, datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
    else:
        [r1, grid, datadir, plotdir, modelstr, r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 

    if 'pr' in plotover:
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
    # AVERAGE R1 ONLY AT HIGH LATITUDES
    ############################################
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
    # PLOT
    ############################################

    plotname = '%s/r1_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    # rae = patches.Rectangle((0,0.9),r1_hl.shape[0],vmax_r1-0.9, alpha=0.5)
    rae = patches.Rectangle((yr_base,0.9),yr_base + r1_hl.shape[0],vmax['r1']-0.9, alpha=0.5)
    ax.add_patch(rae)
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

    if plotover == 'sic':
        ############################################
        # COMPARE WITH SEA ICE CONCENTRATION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)

        sic_file = filenames_raw(sim, 'sic', model=model, timemean=timemean, yr_span=yr_span)
        if not zonmean:
            sic = sic_file.variables['sic'][:]
        else:
            sic = np.mean(np.squeeze(sic_file.variables['sic'][:]),2)

        if sim == 'longrun':
            sic = 100*sic

        sic_hl = lat_mean(sic, grid, lat_int, dim=1)

        sax = ax.twinx()
        if timemean == 'jjamean':
            sax.plot(time, 100-sic_hl, color='tab:blue')
            sax.set_ylabel('Ocean fraction (%)', color='tab:blue')
        else:
            sax.plot(time, sic_hl, color='tab:blue')
            sax.set_ylabel('Sea ice fraction (%)', color='tab:blue')
        sax.set_ylim(vmin['sic'],vmax['sic'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(MultipleLocator(5))

    elif plotover == 'ga_dev':
        ############################################
        # COMPARE WITH LAPSE RATE DEVIATION
        ###########################################
        from proc.ga import make_ga_dev_vint

        plotname = '%s.%s.%g.%g' % (plotname, plotover, vertbnd[0], vertbnd[1])

        ga_dev_vint = make_ga_dev_vint(sim, vertbnd, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

        ga_dev_vint_hl = lat_mean(ga_dev_vint['ga_dev_vint'], ga_dev_vint['grid'], lat_int, dim=1)

        sax = ax.twinx()
        sax.plot(time, ga_dev_vint_hl, color='tab:blue')
        sax.set_ylabel(r'$\langle(\Gamma_m-\Gamma)/\Gamma_m\rangle_{%0.1f}^{%0.1f}$ (%%)' % (vertbnd[0], vertbnd[1]), color='tab:blue')
        sax.set_ylim(vmin['ga_dev'],vmax['ga_dev'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(MultipleLocator(5))

    elif plotover == 'prfrac':
        ############################################
        # COMPARE WITH FRACTION OF LARGE SCALE PRECIPITATION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)
        
        hydro_ref_hl = {}

        # take mean in high latitudes
        prfrac_hl = lat_mean(100*hydro['prl']/hydro['pr'], grid, lat_int, dim=1)

        # compute trends
        if 'ymonmean' not in timemean and sim == 'era5':
            A = np.vstack([time, np.ones(len(time))]).T
            m, c = np.linalg.lstsq(A, prfrac_hl, rcond=None)[0]

        sax = ax.twinx()
        sax.plot(time, prfrac_hl, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P_l/P$ trend$ = %.1f$ %% decade$^{-1}$' % (m*10))
        sax.set_ylabel(r'$P_l/P$ (%)', color='tab:blue')
        sax.set_ylim(vmin['prfrac'],vmax['prfrac'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())
        # ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis
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
        prc_hl = lat_mean(hydro['prc'], grid, lat_int, dim=1)

        # compute trends
        if 'ymonmean' not in timemean and sim == 'era5':
            A = np.vstack([time, np.ones(len(time))]).T
            m, c = np.linalg.lstsq(A, prc_hl, rcond=None)[0]

        sax = ax.twinx()
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
            sax.fill_between(time, r1_dc_mmm_hl['dr1']['prc25'], r1_dc_mmm_hl['dr1']['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
            sax.fill_between(time, r1_dc_mmm_hl['dyn']['prc25'], r1_dc_mmm_hl['dyn']['prc75'], facecolor='maroon', alpha=0.2, edgecolor=None)
            sax.fill_between(time, r1_dc_mmm_hl['rad']['prc25'], r1_dc_mmm_hl['rad']['prc75'], facecolor='tab:gray', alpha=0.3, edgecolor=None)
            sax.fill_between(time, r1_dc_mmm_hl['res']['prc25'], r1_dc_mmm_hl['res']['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
        sax.axhline(0, color='k', linewidth=0.5)
        sax.plot(time, r1_dc_hl['dyn'], color='maroon', label='$\overline{R_1} \dfrac{\Delta (\partial_t m + \partial_y (vm))}{\overline{\partial_t m + \partial_y (vm)}}$')
        sax.plot(time, r1_dc_hl['rad'], color='lightgray', label='$-\overline{R_1} \dfrac{\Delta R_a}{\overline{R_a}}$')
        sax.plot(time, r1_dc_hl['res'], '-.k', label='Residual')
        sax.set_ylabel(r'$\Delta R_1$ (unitless)', color='black')
        sax.set_ylim(vmin['r1']-r1_hl[0],vmax['r1']-r1_hl[0])
        sax.tick_params(axis='y', labelcolor='black', color='black')
        sax.yaxis.set_minor_locator(AutoMinorLocator())

    if legend and sim == 'era5':
        fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()