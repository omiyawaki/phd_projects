import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import *
from misc.load_data import load_flux, load_cl_in
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

def cl_in(sim, **kwargs):

    categ = 'cl_in'

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
            vmin_dev = -50
            vmax_dev = 50
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
        [flux, grid_fl, datadir, plotdir, modelstr] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        [cl_in, grid_cl, datadir, plotdir, modelstr] = load_cl_in(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [flux, grid_fl, datadir, plotdir, modelstr, flux_mmm] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        [cl_in, grid_cl, datadir, plotdir, modelstr, cl_in_mmm] = load_cl_in(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    grid = grid_cl

    if refclim == 'hist-30':
        if isinstance(model, str):
            [flux_ref, _, _, _, _] = load_flux(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
            [cl_in_ref, _, _, _, _] = load_cl_in(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        else:
            [flux_ref, _, _, _, _, flux_ref_mmm] = load_flux(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
            [cl_in_ref, _, _, _, _, cl_in_ref_mmm] = load_cl_in(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
        
        if timemean == 'djfmean':
            for fluxname in flux_ref:
                flux_ref[fluxname] = np.mean(np.roll(flux_ref[fluxname],1,axis=0)[0:3], 0)
            for cl_inname in cl_in_ref:
                cl_in_ref[cl_inname] = np.mean(np.roll(cl_in_ref[cl_inname],1,axis=0)[0:3], 0)
        elif timemean == 'jjamean':
            for fluxname in flux_ref:
                flux_ref[fluxname] = np.mean(flux_ref[fluxname][5:8], 0)
            for cl_inname in cl_in_ref:
                cl_in_ref[cl_inname] = np.mean(cl_in_ref[cl_inname][5:8], 0)

        flux_ref_hl = {}
        cl_in_ref_hl = {}

    ############################################
    # AVERAGE cl_inES ONLY AT HIGH LATITUDES
    ############################################
    flux_hl = {}
    flux_dev_hl = {}
    for fluxname in flux:
        flux_hl[fluxname] = lat_mean(flux[fluxname], grid_fl, lat_int, dim=1)
        
        # if refclim == 'init':
        #     flux_dev_hl[fluxname] = flux_hl[fluxname] - flux_hl[fluxname][0]
        # else:
        #     print(grid_fl['lat'].shape)
        #     print(flux_ref[fluxname].shape)
        #     flux_ref_hl[fluxname] = lat_mean(flux_ref[fluxname], grid_fl, lat_int, dim=0)
        #     flux_dev_hl[fluxname] = flux_hl[fluxname] - flux_ref_hl[fluxname]

    cl_in_hl = {}
    cl_in_dev_hl = {}
    for cl_inname in cl_in:
        cl_in_hl[cl_inname] = lat_mean(cl_in[cl_inname], grid_cl, lat_int, dim=2)
        
        # if refclim == 'init':
        #     cl_in_dev_hl[cl_inname] = cl_in_hl[cl_inname] - cl_in_hl[cl_inname][0]
        # else:
        #     cl_in_ref_hl[cl_inname] = lat_mean(cl_in_ref[cl_inname], grid_cl, lat_int, dim=0)
        #     cl_in_dev_hl[cl_inname] = cl_in_hl[cl_inname] - cl_in_ref_hl[cl_inname]


    ############################################
    # PLOT HUS
    ############################################
    plotname = remove_repdots('%s/hus.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))

    fig, ax = plt.subplots()
    lp_hus = ax.plot(cl_in_hl['hus'][0,:], 1e-2*grid['lev'], color='k', label='$R_a$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('q (kg kg$^{-1}$)')
    ax.set_ylabel('p (hPa)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(yr_base,yr_base+flux_hl['ra'].shape[0]-1)
    ax.set_ylim(1e2,1e3)
    ax.invert_yaxis()
    # if legend:
    #     ax.legend()
    fig.set_size_inches(5, 4)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT HUR
    ############################################
    plotname = remove_repdots('%s/hur.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))

    fig, ax = plt.subplots()
    lp_hur = ax.plot(cl_in_hl['hur'][0,:], 1e-2*grid['lev'], color='k', label='$R_a$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('RH (%)')
    ax.set_ylabel('p (hPa)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(yr_base,yr_base+flux_hl['ra'].shape[0]-1)
    ax.set_ylim(1e2,1e3)
    ax.invert_yaxis()
    # if legend:
    #     ax.legend()
    fig.set_size_inches(5, 4)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT TA
    ############################################
    plotname = remove_repdots('%s/ta.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))

    fig, ax = plt.subplots()
    lp_ta = ax.plot(cl_in_hl['ta'][0,:], 1e-2*grid['lev'], color='k', label='$R_a$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('T (K)')
    ax.set_ylabel('p (hPa)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(yr_base,yr_base+flux_hl['ra'].shape[0]-1)
    ax.set_ylim(1e2,1e3)
    ax.invert_yaxis()
    # if legend:
    #     ax.legend()
    fig.set_size_inches(5, 4)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    # save data
    pickle.dump([cl_in_hl, flux_hl, grid_cl], open(remove_repdots('/project2/tas1/miyawaki/projects/003/climlab/input_data/cl_in.%s.%s.pickle' % (zonmean, timemean)), 'wb'))
