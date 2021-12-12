import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import *
from misc.filenames import *
# from misc.translate import translate_varname
from misc.means import lat_mean, global_int
from misc.load_data import *
from misc import par
from proc.r1 import save_r1
from plot.titles import make_title_sim_time_lat_lon
import os
import pickle
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def r1_mon_ll(sim, **kwargs):

    categ = 'mon_ll'

    zonmean = kwargs.get('zonmean', 0) # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    latbnd = kwargs.get('latbnd', (10,30))
    lonbnd = kwargs.get('lonbnd', (70,90))
    latstep = kwargs.get('latstep', 0.25) # latitude step size used for interpolation
    lonstep = kwargs.get('lonstep', 0.25) # lonitude step size used for interpolonion
    plotover = kwargs.get('plotover', None) # plot overlay (sic for sea ice, ga_dev for lapse rate deviation, pr for precip, decomp for linear decomposition)?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    legend = kwargs.get('legend', 0)
    if plotover == 'ga_dev':
        vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
        # vertbnd = kwargs.get('vertbnd', (0.7, 0.3)) # sigma bounds of vertical integral
        vertbnd = kwargs.get('vertbnd', (1.0, 0.9)) # sigma bounds of vertical integral

    lat_int = np.arange(latbnd[0], latbnd[1], latstep)
    lon_int = np.arange(lonbnd[0], lonbnd[1], lonstep)

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

    ##################################
    # LOAD DATA
    ##################################
    if isinstance(model, str) or model is None:
        [r1, grid, datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
    else:
        [r1, grid, datadir, plotdir, modelstr, r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 

    if 'pr' in plotover or 'cl' in plotover:
        if isinstance(model, str) or model is None:
            [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        else:
            [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    ############################################
    # AVERAGE R1 at selected LAT-LON
    ############################################

    if not ( isinstance(model, str) or (model is None) ):
        r1_ll_mmm = dict()
        for i in r1_mmm:
            r1_lon_mmm[i] = lon_mean(r1_mmm[i], grid, lon_int, dim=2)
            r1_ll_mmm[i] = lat_mean(r1_mmm[i], grid, lat_int, dim=1)

    r1_lon = lon_mean(r1, grid, lon_int, dim=2)
    r1_ll = lat_mean(r1_lon, grid, lat_int, dim=1)

    time = yr_base + np.arange(r1_ll.shape[0]) # create time vector

    ############################################
    # PLOT
    ############################################

    plotname = '%s/r1_mon_ll.%g.%g.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], lonbnd[0], lonbnd[1], timemean)

    fig, ax = plt.subplots()
    # rae = patches.Rectangle((0,0.9),r1_ll.shape[0],vmax_r1-0.9, alpha=0.5)
    # rae = patches.Rectangle((yr_base,0.9),yr_base + r1_ll.shape[0],vmax['r1']-0.9, alpha=0.5)
    # ax.add_patch(rae)
    if not ( isinstance(model, str) or (model is None) ):
        spr_r1 = ax.fill_between(time, r1_ll_mmm['prc25'], r1_ll_mmm['prc75'], facecolor='black', alpha=0.2, edgecolor=None)

    lp_r1 = ax.plot(time, r1_ll, color='black')
    # make_title_sim_time_lat_lon(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1], lon1=lonbnd[0], lon2=lonbnd[1])
    make_title_lat_lon(ax, lat1=latbnd[0], lat2=latbnd[1], lon1=lonbnd[0], lon2=lonbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    print(timemean)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('$R_1$ (unitless)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+r1_ll.shape[0]-1)
    ax.set_ylim(-2,0.5)

    if plotover == 'ga_dev':
        ############################################
        # COMPARE WITH LAPSE RATE DEVIATION
        ###########################################
        from proc.ga import make_ga_dev_vint

        plotname = '%s.%s.%g.%g' % (plotname, plotover, vertbnd[0], vertbnd[1])

        ga_dev_vint = make_ga_dev_vint(sim, vertbnd, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

        ga_dev_vint_ll = lat_mean(ga_dev_vint['ga_dev_vint'], ga_dev_vint['grid'], lat_int, dim=1)

        sax = ax.twinx()
        sax.plot(time, ga_dev_vint_ll, color='tab:blue')
        sax.set_ylabel(r'$\langle(\Gamma_m-\Gamma)/\Gamma_m\rangle_{%0.1f}^{%0.1f}$ (%%)' % (vertbnd[0], vertbnd[1]), color='tab:blue')
        sax.set_ylim(vmin['ga_dev'],vmax['ga_dev'])
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(MultipleLocator(5))

    elif plotover == 'prfrac':
        ############################################
        # COMPARE WITH FRACTION OF LARGE SCALE PRECIPITATION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)

        hydro_ref_ll = {}

        # take mean in high latitudes
        pr_ll_mmm = dict()
        prl_ll_mmm = dict()
        prc_ll_mmm = dict()
        prfrac_ll_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in hydro_mmm['pr']:
                pr_ll_mmm[i] = lat_mean(hydro_mmm['pr'][i], grid, lat_int, dim=1)
                prl_ll_mmm[i] = lat_mean(hydro_mmm['prl'][i], grid, lat_int, dim=1)
                prc_ll_mmm[i] = lat_mean(hydro_mmm['prc'][i], grid, lat_int, dim=1)
                prfrac_ll_mmm[i] = lat_mean(hydro_mmm['prc'][i]/hydro_mmm['pr'][i], grid, lat_int, dim=1)

        prfrac_ll = lat_mean(hydro['prc']/hydro['pr'], grid, lat_int, dim=1)

        sax = ax.twinx()
        if not ( isinstance(model, str) or (model is None) ):
            sax.fill_between(time, prfrac_ll_mmm['prc25'], prfrac_ll_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        sax.plot(time, prfrac_ll, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P_c/P$ trend$ = %.1f$ %% decade$^{-1}$' % (m*10))
        sax.set_ylabel(r'$P_c/P$ (unitless)', color='tab:blue')
        # sax.set_ylim(vmin['prfrac'],vmax['prfrac'])
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
        
        hydro_ref_ll = {}

        # take mean in high latitudes
        pr_lon = lon_mean(hydro['pr'], grid, lon_int, dim=2)
        pr_ll = lat_mean(pr_lon, grid, lat_int, dim=1)

        sax = ax.twinx()
        sax.plot(time, pr_ll, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P$ trend$ = %.1f$ %% decade$^{-1}$' % (m/np.nanmean(pr_ll)*1e3))
        sax.set_ylabel(r'$P$ (mm d$^{-1}$)', color='tab:blue')
        sax.set_ylim(0,12)
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis

    elif plotover == 'prc':
        ############################################
        # COMPARE WITH CONVECTIVE PRECIPITATION
        ###########################################
        plotname = '%s.%s' % (plotname, plotover)
        
        hydro_ref_ll = {}

        # take mean in high latitudes
        prc_ll_mmm = dict()
        if not ( isinstance(model, str) or (model is None) ):
            for i in hydro_mmm['prc']:
                prc_ll_mmm[i] = lat_mean(hydro_mmm['prc'][i], grid, lat_int, dim=1)

        prc_ll = lat_mean(hydro['prc'], grid, lat_int, dim=1)

        # compute trends
        if 'ymonmean' not in timemean and sim == 'era5':
            A = np.vstack([time, np.ones(len(time))]).T
            m, c = np.linalg.lstsq(A, prc_ll, rcond=None)[0]

        sax = ax.twinx()
        if not ( isinstance(model, str) or (model is None) ):
            spr_r1 = sax.fill_between(time, prc_ll_mmm['prc25'], prc_ll_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        sax.plot(time, prc_ll, color='tab:blue')
        if 'ymonmean' not in timemean and sim == 'era5':
            # sax.plot(time, m*time + c, '--', color='tab:blue', label='%g mm d$^{-1}$ decade$^{-1}$' % (m*10))
            sax.plot(time, m*time + c, '--', color='tab:blue', label='$P_c$ trend$ = %.1f$ %% decade$^{-1}$' % (m/np.nanmean(prc_ll)*1e3))
        sax.set_ylabel(r'$P_c$ (mm d$^{-1}$)', color='tab:blue')
        sax.set_ylim(0,15)
        sax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
        sax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis




    # if legend and sim == 'era5':
    #     fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)

    fig.set_size_inches(3, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

def make_title_lat_lon(ax, **kwargs):

    lat1 = kwargs.get('lat1')
    lat2 = kwargs.get('lat2')
    lon1 = kwargs.get('lon1')
    lon2 = kwargs.get('lon2')

    ax.set_title('$\phi=%g^\circ$%s to $%g^\circ$%s\n$\lambda=%g^\circ$%s to $%g^\circ$%s' % (abs(lat1), lat_ns(lat1), abs(lat2), lat_ns(lat2), lon1, lon_ew(lon1), lon2, lon_ew(lon2)))
    
def lat_ns(lat):

    if lat >= 0:
        latstr = 'N'
    else:
        latstr = 'S'

    return latstr

def lon_ew(lon):

    if lon >= 0:
        lonstr = 'E'
    else:
        lonstr = 'W'

    return lonstr
