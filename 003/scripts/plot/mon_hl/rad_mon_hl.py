import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import *
from misc.load_data import load_flux, load_rad, load_hydro
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

def rad_mon_hl(sim, **kwargs):

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
        [rad, grid, datadir, plotdir, modelstr] = load_rad(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [rad, grid, datadir, plotdir, modelstr, rad_mmm] = load_rad(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    if refclim == 'hist-30':
        if isinstance(model, str):
            [rad_ref, _, _, _, _] = load_rad(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
            [hydro_ref, _, _, _, _] = load_hydro(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        else:
            [rad_ref, _, _, _, _, rad_ref_mmm] = load_rad(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
            [hydro_ref, _, _, _, _, hydro_ref_mmm] = load_hydro(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)

        if timemean == 'djfmean':
            for radname in rad_ref:
                rad_ref[radname] = np.mean(np.roll(rad_ref[radname],1,axis=0)[0:3], 0)
            for hydroname in hydro_ref:
                hydro_ref[hydroname] = np.mean(np.roll(hydro_ref[hydroname],1,axis=0)[0:3], 0)
        elif timemean == 'jjamean':
            for radname in rad_ref:
                rad_ref[radname] = np.mean(rad_ref[radname][5:8], 0)
            for hydroname in hydro_ref:
                hydro_ref[hydroname] = np.mean(hydro_ref[hydroname][5:8], 0)

        rad_ref_hl = {}
        hydro_ref_hl = {}

    ############################################
    # AVERAGE FLUXES ONLY AT HIGH LATITUDES
    ############################################
    rad_hl = {}
    rad_dev_hl = {}
    for radname in rad:
        rad_hl[radname] = lat_mean(rad[radname], grid, lat_int, dim=1)

        if refclim == 'init':
            rad_dev_hl[radname] = rad_hl[radname] - rad_hl[radname][0]
        else:
            rad_ref_hl[radname] = lat_mean(rad_ref[radname], grid, lat_int, dim=0)
            rad_dev_hl[radname] = rad_hl[radname] - rad_ref_hl[radname]

    hydro_hl = {}
    hydro_dev_hl = {}
    for hydroname in hydro:
        hydro_hl[hydroname] = lat_mean(hydro[hydroname], grid, lat_int, dim=1)

        if refclim == 'init':
            hydro_dev_hl[hydroname] = hydro_hl[hydroname] - hydro_hl[hydroname][0]
        else:
            hydro_ref_hl[hydroname] = lat_mean(hydro_ref[hydroname], grid, lat_int, dim=0)
            hydro_dev_hl[hydroname] = hydro_hl[hydroname] - hydro_ref_hl[hydroname]

    time = yr_base + np.arange(rad_hl['ra'].shape[0]) # create time vector

    ############################################
    # PLOT (W/ CLEAR SKY RAD, DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/rad_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ra = ax.plot(time, rad_dev_hl['ra'], color='tab:gray', label='$\Delta R_a$')
    lp_ra_cs = ax.plot(time, rad_dev_hl['ra_cs'], '--', color='tab:gray', label='$\Delta R_a$')
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
    ax.set_xlim(yr_base,yr_base+rad_dev_hl['ra'].shape[0]-1)
    ax.set_ylim(vmin_dev,vmax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (FRACTIONAL CLEAR SKY RAD, DEVIATION FROM INITIAL)
    ############################################
    A = np.vstack([hydro_dev_hl['tas'], np.ones(len(hydro_dev_hl['tas']))]).T
    m, c = np.linalg.lstsq(A, 1e2*rad_dev_hl['ra']/rad_hl['ra'], rcond=None)[0]

    t_vec = np.arange(0,35,101)

    plotname = remove_repdots('%s/rad_frac_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ra = ax.plot(hydro_dev_hl['tas'], 1e2*rad_dev_hl['ra']/rad_hl['ra'], '.', color='tab:gray', label='$\Delta R_a / R_a$')
    ax.plot(t_vec, m*t_vec+c, '-k')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('$\Delta T_{2\,m}$ (K)')
    ax.set_ylabel('$\Delta R_a / R_a$ (%)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(yr_base,yr_base+rad_dev_hl['ra'].shape[0]-1)
    # ax.set_ylim(vmin_dev,vmax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (SW/LW DECOMP, DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/radcs_lwsw_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ra = ax.plot(time, rad_dev_hl['ra'], '-', color='tab:gray', label='$\Delta R_a$')
    lp_ra_lw = ax.plot(time, rad_dev_hl['lw'], '-', color='tab:green', label='$\Delta LW$')
    lp_ra_sw = ax.plot(time, rad_dev_hl['sw'], '-', color='yellow', label='$\Delta SW$')
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
    ax.set_xlim(yr_base,yr_base+rad_dev_hl['ra'].shape[0]-1)
    ax.set_ylim(vmin_dev,vmax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (LW CLEAR SKY RAD, DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/rad_lwcs_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ra = ax.plot(time, rad_dev_hl['ra'], '-', color='tab:gray', label='$\Delta R_a$')
    lp_ra_cs_lw = ax.plot(time, rad_dev_hl['lw_cs'], '--', color='tab:red', label='$\Delta \mathrm{LW}_{clear}$')
    lp_ra_cld_lw = ax.plot(time, rad_dev_hl['lw'] - rad_dev_hl['lw_cs'], ':', color='tab:red', label='$\Delta \mathrm{LW}_{cloud}$')
    lp_ra_cs_sw = ax.plot(time, rad_dev_hl['sw_cs'], '--', color='tab:blue', label='$\Delta \mathrm{SW}_{clear}$')
    lp_ra_cld_sw = ax.plot(time, rad_dev_hl['sw'] - rad_dev_hl['sw_cs'], ':', color='tab:blue', label='$\Delta \mathrm{SW}_{cloud}$')
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
    ax.set_xlim(yr_base,yr_base+rad_dev_hl['ra'].shape[0]-1)
    ax.set_ylim(vmin_dev,vmax_dev)
    fig.set_size_inches(4,3)
    plt.tight_layout()
    # save without legend
    plt.savefig(remove_repdots('%s.noleg.pdf' % (plotname)), format='pdf', dpi=300)
    if legend:
        ax.legend()
    plt.tight_layout()
    # save with legend
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    # save with legend outside
    # add legend
    legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.3), ncol=5)
    # cut off excess space on the bottom 
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.05,
        box.width, box.height * 0.95])

    # alter figure aspect ratio to accomodate legend
    fig.set_size_inches(6,3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.legoutside.pdf' % (plotname)), format='pdf', dpi=300)
    # save legend only
    expand = [-5,-5,5,5]
    fi2  = legend.figure
    fi2.canvas.draw()
    bbox  = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fi2.dpi_scale_trans.inverted())
    fi2.savefig(remove_repdots('%s.legonly.pdf' % (plotname)), dpi=300, bbox_inches=bbox)
    if viewplt:
        plt.show()
    plt.close()
