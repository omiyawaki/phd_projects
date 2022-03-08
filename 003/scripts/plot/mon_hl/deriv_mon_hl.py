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

def deriv_mon_hl(sim, **kwargs):

    rolling_mean = 80; # smooth data using a rolling mean? (units: yr)

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
    elif sim == 'hist+rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '186001-229912')
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

    ############################################
    # AVERAGE R1 ONLY AT HIGH LATITUDES
    ############################################
    if timemean == 'allseas':
        djf_r1_hl = lat_mean(djf_r1, djf_grid, lat_int, dim=1)
        mam_r1_hl = lat_mean(mam_r1, mam_grid, lat_int, dim=1)
        jja_r1_hl = lat_mean(jja_r1, jja_grid, lat_int, dim=1)
        son_r1_hl = lat_mean(son_r1, son_grid, lat_int, dim=1)

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

        r1_hl = uniform_filter(r1_hl, rolling_mean) # apply rolling mean

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
    # PLOT (smoothed timeseries)
    ############################################

    plotname = '%s/r1smooth_mon_hl.%g.%g.%s.%g' % (plotdir, latbnd[0], latbnd[1], timemean, rolling_mean)

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
                spr_r1 = ax.fill_between(time[1:-rolling_mean], r1_hl_mmm['prc25'][1:-rolling_mean], r1_hl_mmm['prc75'][1:-rolling_mean], facecolor='black', alpha=0.2, edgecolor=None)
            elif spread == 'std':
                spr_r1 = ax.fill_between(time[1:-rolling_mean], r1_hl_mmm['mmm'][1:-rolling_mean] - r1_hl_mmm['std'][1:-rolling_mean], r1_hl_mmm['mmm'][1:-rolling_mean] + r1_hl_mmm['std'][1:-rolling_mean], facecolor='black', alpha=0.2, edgecolor=None)

        lp_r1 = ax.plot(time[1:-rolling_mean], r1_hl[1:-rolling_mean], color='black')
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

    fig.set_size_inches(4, 3)

    ############################################
    # DECOMPOSE INTO DYNAMIC AND RADIATIVE COMPONENTS
    ###########################################
    plotname = '%s.%s' % (plotname, 'decomp')

    if isinstance(model, str) or model is None:
        [r1_dc, grid, datadir, plotdir, modelstr] = load_r1_dc(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span) 
    else:
        [r1_dc, grid, datadir, plotdir, modelstr, r1_dc_mmm] = load_r1_dc(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span) 

    ############################################
    # AVERAGE R1 COMPONENTS ONLY AT HIGH LATITUDES
    ############################################
    r1_dc_hl = {}
    for varname in r1_dc:
        r1_dc_hl[varname] = lat_mean(r1_dc[varname], grid, lat_int, dim=1)
        r1_dc_hl[varname] = uniform_filter(r1_dc_hl[varname], rolling_mean)

    if not (isinstance(model, str) or model is None):
        r1_dc_mmm_hl = {}
        for varname in r1_dc_mmm:
            r1_dc_mmm_hl[varname] = {}
            for stat in r1_dc_mmm[varname]:
                r1_dc_mmm_hl[varname][stat] = lat_mean(r1_dc_mmm[varname][stat], grid, lat_int, dim=1)
                rc_dc_mmm_hl = uniform_filter(r1_dc_mmm_hl[varname][stat], rolling_mean) # apply rolling mean

    sax = ax.twinx()
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            # sax.fill_between(time, r1_dc_mmm_hl['dr1']['prc25'], r1_dc_mmm_hl['dr1']['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
            sax.fill_between(time[1:-rolling_mean], r1_dc_mmm_hl['dyn']['prc25'][1:-rolling_mean], r1_dc_mmm_hl['dyn']['prc75'][1:-rolling_mean], facecolor='maroon', alpha=0.2, edgecolor=None)
            sax.fill_between(time[1:-rolling_mean], r1_dc_mmm_hl['rad']['prc25'][1:-rolling_mean], r1_dc_mmm_hl['rad']['prc75'][1:-rolling_mean], facecolor='tab:gray', alpha=0.3, edgecolor=None)
            sax.fill_between(time[1:-rolling_mean], r1_dc_mmm_hl['res']['prc25'][1:-rolling_mean], r1_dc_mmm_hl['res']['prc75'][1:-rolling_mean], facecolor='k', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            sax.fill_between(time[1:-rolling_mean], r1_dc_mmm_hl['dyn']['mmm'][1:-rolling_mean]-r1_dc_mmm_hl['dyn']['std'][1:-rolling_mean], r1_dc_mmm_hl['dyn']['mmm'][1:-rolling_mean]+r1_dc_mmm_hl['dyn']['std'][1:-rolling_mean], facecolor='maroon', alpha=0.2, edgecolor=None)
            sax.fill_between(time[1:-rolling_mean], r1_dc_mmm_hl['rad']['mmm'][1:-rolling_mean]-r1_dc_mmm_hl['rad']['std'][1:-rolling_mean], r1_dc_mmm_hl['rad']['mmm'][1:-rolling_mean]+r1_dc_mmm_hl['rad']['std'][1:-rolling_mean], facecolor='tab:gray', alpha=0.3, edgecolor=None)
            sax.fill_between(time[1:-rolling_mean], r1_dc_mmm_hl['res']['mmm'][1:-rolling_mean]-r1_dc_mmm_hl['res']['std'][1:-rolling_mean], r1_dc_mmm_hl['res']['mmm'][1:-rolling_mean]+r1_dc_mmm_hl['res']['std'][1:-rolling_mean], facecolor='k', alpha=0.2, edgecolor=None)
    sax.axhline(0, color='k', linewidth=0.5)
    sax.plot(time[1:-rolling_mean], r1_dc_hl['dr1'][1:-rolling_mean], color='k', label='$\Delta{R_1}$')
    sax.plot(time[1:-rolling_mean], r1_dc_hl['dyn'][1:-rolling_mean], color='maroon', label='$\overline{R_1} \dfrac{\Delta (\partial_t m + \partial_y (vm))}{\overline{\partial_t m + \partial_y (vm)}}$')
    sax.plot(time[1:-rolling_mean], r1_dc_hl['res'][1:-rolling_mean], '-.k', label='Residual')
    sax.plot(time[1:-rolling_mean], r1_dc_hl['rad'][1:-rolling_mean], color='lightgray', label='$-\overline{R_1} \dfrac{\Delta R_a}{\overline{R_a}}$')
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


    ############################################
    # PLOT (derivative)
    ############################################

    plotname = '%s/r1deriv_mon_hl.%g.%g.%s.%g' % (plotdir, latbnd[0], latbnd[1], timemean, rolling_mean)

    if isinstance(model, str) or model is None:
        [dr1, grid, datadir, plotdir, modelstr] = load_dr1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span) 
    else:
        [dr1, grid, datadir, plotdir, modelstr, dr1_mmm] = load_dr1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span) 

    ############################################
    # AVERAGE R1 COMPONENTS ONLY AT HIGH LATITUDES
    ############################################
    dr1_hl = {}
    for varname in dr1:
        dr1_hl[varname] = lat_mean(dr1[varname], grid, lat_int, dim=1)
        dr1_hl[varname] = uniform_filter(dr1_hl[varname], rolling_mean)

    if not (isinstance(model, str) or model is None):
        dr1_mmm_hl = {}
        for varname in dr1_mmm:
            dr1_mmm_hl[varname] = {}
            for stat in dr1_mmm[varname]:
                dr1_mmm_hl[varname][stat] = lat_mean(dr1_mmm[varname][stat], grid, lat_int, dim=1)
                rc_dc_mmm_hl = uniform_filter(dr1_mmm_hl[varname][stat], rolling_mean) # apply rolling mean

    fig, ax = plt.subplots()
    # if not ( isinstance(model, str) or (model is None) ):
    #     if spread == 'prc':
    #         spr_r1 = ax.fill_between(time, r1_hl_mmm['prc25'], r1_hl_mmm['prc75'], facecolor='black', alpha=0.2, edgecolor=None)
    #     elif spread == 'std':
    #         spr_r1 = ax.fill_between(time, r1_hl_mmm['mmm'] - r1_hl_mmm['std'], r1_hl_mmm['mmm'] + r1_hl_mmm['std'], facecolor='black', alpha=0.2, edgecolor=None)

    ax.axhline(0, color='k', linewidth=0.5)
    lp_res, = ax.plot(time[1:-rolling_mean], dr1_hl['res'][:-rolling_mean], '-.k', label='Residual')
    lp_dyn, = ax.plot(time[1:-rolling_mean], dr1_hl['dyn'][:-rolling_mean], color='maroon', label='Advection')
    lp_rad, = ax.plot(time[1:-rolling_mean], dr1_hl['rad'][:-rolling_mean], color='tab:gray', label='Radiative')
    lp_r1, = ax.plot(time[1:-rolling_mean], dr1_hl['dr1'][:-rolling_mean], color='black', label='$dR_1/dt$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    print(timemean)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('$dR_1/dt$ (yr$^{-1}$)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(yr_base,yr_base+r1_hl.shape[0]-1)
    ax.set_xlim(2006,2300)
    ax.set_ylim(vmin['dr1dt'],vmax['dr1dt'])

    # add legend
    ax.legend(handles=[lp_r1, lp_dyn, lp_res, lp_rad], loc='upper center', bbox_to_anchor=(0.5,-0.2), ncol=2)

    # cut off excess space on the bottom 
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.05,
                         box.width, box.height * 0.95])

    # alter figure aspect ratio to accomodate legend
    fig.set_size_inches(4.5,4)

    if legend and sim == 'era5':
        fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # Dr1Dt only
    ############################################
    plotname = '%s/r1deriv_only_mon_hl.%g.%g.%s.%g' % (plotdir, latbnd[0], latbnd[1], timemean, rolling_mean)

    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_r1, = ax.plot(time[1:-rolling_mean][1:], dr1_hl['dr1'][1:-rolling_mean], color='black', label='$dR_1/dt$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    print(timemean)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('$dR_1/dt$ (yr$^{-1}$)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(yr_base,yr_base+r1_hl.shape[0]-1)
    ax.set_xlim(2006,2300)
    ax.set_ylim(vmin['dr1dt'],vmax['dr1dt'])
    fig.set_size_inches(4,3.5)

    if legend and sim == 'era5':
        fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # RELATIVE COMPONENTS
    ############################################
    plotname = '%s/r1deriv_rel_mon_hl.%g.%g.%s.%g' % (plotdir, latbnd[0], latbnd[1], timemean, rolling_mean)

    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axhline(1, color='k', linewidth=0.5)
    lp_res, = ax.plot(time[1:-rolling_mean], dr1_hl['res'][:-rolling_mean] / dr1_hl['dr1'][:-rolling_mean], '-.k', label='Residual')
    lp_dyn, = ax.plot(time[1:-rolling_mean], dr1_hl['dyn'][:-rolling_mean] / dr1_hl['dr1'][:-rolling_mean], color='maroon', label='Advection')
    lp_rad, = ax.plot(time[1:-rolling_mean], dr1_hl['rad'][:-rolling_mean] / dr1_hl['dr1'][:-rolling_mean], color='tab:gray', label='Radiative')
    # lp_r1, = ax.plot(time[1:], dr1_hl['dr1'], color='black', label='$dR_1/dt$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    print(timemean)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Fraction of $dR_1/dt$ (unitless)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(yr_base,yr_base+r1_hl.shape[0]-1)
    ax.set_xlim(2006,2300)
    ax.set_ylim(-0.5,1.5)

    # add legend
    ax.legend(handles=[lp_dyn, lp_rad, lp_res], loc='upper center', bbox_to_anchor=(0.5,-0.2), ncol=3)

    # cut off excess space on the bottom 
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.05,
                         box.width, box.height * 0.95])

    # alter figure aspect ratio to accomodate legend
    fig.set_size_inches(4.5,4)

    if legend and sim == 'era5':
        fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()
