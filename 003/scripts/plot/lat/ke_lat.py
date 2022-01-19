import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
import numpy as np
from misc.load_data import *
from misc.filenames import *
from misc.dirnames import *
from misc.translate import *
from plot.titles import make_title_sim_time
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def ke_lat(sim, **kwargs):

    categ = 'lat'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    legend = kwargs.get('legend', 0) # draw legend?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    sim_ref = kwargs.get('sim_ref', 'historical')
    timemean_ref = kwargs.get('timemean_ref', 'ymonmean-30')
    yr_span_ref = kwargs.get('yr_span_ref', '186001-200512')

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

    ##################################
    # LOAD DATA
    ##################################
    if isinstance(model, str) or model is None:
        [r1, grid, datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
    else:
        [r1, grid, datadir, plotdir, modelstr, r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 

    if isinstance(model, str) or model is None:
        [flux, grid, datadir, plotdir, modelstr] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
    else:
        [flux, grid, datadir, plotdir, modelstr, flux_mmm] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 

    if isinstance(model, str) or model is None:
        [ke, grid, datadir, plotdir, modelstr] = load_ke(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [ke, grid, datadir, plotdir, modelstr, ke_mmm] = load_ke(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    ############################################
    # PLOT (TKE ANN, DJF and JJA)
    ############################################
    plotname = remove_repdots('%s/mean_tke.po92.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ann = ax.plot(grid['lat'], np.mean(ke['mtke'], 0), 'k', label='ANN', linewidth=2)
    lp_jja = ax.plot(grid['lat'], np.mean(ke['mtke'][5:8], 0), '--k', label='JJA', linewidth=1)
    lp_jja = ax.plot(grid['lat'], np.mean(np.roll(ke['mtke'],1,axis=0)[0:3], 0), '-k', label='DJF', linewidth=1)
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Total Kinetic Energy (m$^{2}$ s$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    # ax.set_ylim(divin_dev,divax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (TEKE ANN, DJF and JJA)
    ############################################
    plotname = remove_repdots('%s/mean_teke.po92.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ann = ax.plot(grid['lat'], np.mean(ke['mteke'], 0), 'k', label='ANN', linewidth=2)
    lp_jja = ax.plot(grid['lat'], np.mean(ke['mteke'][5:8], 0), '--k', label='JJA', linewidth=1)
    lp_jja = ax.plot(grid['lat'], np.mean(np.roll(ke['mteke'],1,axis=0)[0:3], 0), '-k', label='DJF', linewidth=1)
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Transient Eddy Kinetic Energy (m$^{2}$ s$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    # ax.set_ylim(divin_dev,divax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (SEKE ANN, DJF and JJA)
    ############################################
    plotname = remove_repdots('%s/mean_seke.po92.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ann = ax.plot(grid['lat'], np.mean(ke['mseke'], 0), 'k', label='ANN', linewidth=2)
    lp_jja = ax.plot(grid['lat'], np.mean(ke['mseke'][5:8], 0), '--k', label='JJA', linewidth=1)
    lp_jja = ax.plot(grid['lat'], np.mean(np.roll(ke['mseke'],1,axis=0)[0:3], 0), '-k', label='DJF', linewidth=1)
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Stationary Eddy Kinetic Energy (m$^{2}$ s$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    # ax.set_ylim(divin_dev,divax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (mke ANN, DJF and JJA)
    ############################################
    plotname = remove_repdots('%s/mean_mke.po92.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_ann = ax.plot(grid['lat'], np.mean(ke['mmke'], 0), 'k', label='ANN', linewidth=2)
    lp_jja = ax.plot(grid['lat'], np.mean(ke['mmke'][5:8], 0), '--k', label='JJA', linewidth=1)
    lp_jja = ax.plot(grid['lat'], np.mean(np.roll(ke['mmke'],1,axis=0)[0:3], 0), '-k', label='DJF', linewidth=1)
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Mean Meridional Circulation Kinetic Energy (m$^{2}$ s$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    # ax.set_ylim(divin_dev,divax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (Mean KE)
    ############################################
    plotname = remove_repdots('%s/mean_ke.annmean.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_tke = ax.plot(grid['lat'], np.mean(ke['mtke'], 0), 'k', label='Total')
    lp_teke = ax.plot(grid['lat'], np.mean(ke['mteke'], 0), 'r', label='Transient Eddies')
    lp_seke = ax.plot(grid['lat'], np.mean(ke['mseke'], 0), 'b', label='Stationary Eddies')
    lp_mke = ax.plot(grid['lat'], np.mean(ke['mmke'], 0), 'g', label='Mean Meridional Circulation')
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Vertical Mean Kinetic Energy (m$^{2}$ s$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    # ax.set_ylim(divin_dev,divax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (KE)
    ############################################
    print(plotdir)
    plotname = remove_repdots('%s/ke.annmean.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_tke = ax.plot(grid['lat'], 1e-6*np.mean(ke['tke'], 0), 'k', label='Total')
    lp_teke = ax.plot(grid['lat'], 1e-6*np.mean(ke['teke'], 0), 'r', label='Transient Eddies')
    lp_seke = ax.plot(grid['lat'], 1e-6*np.mean(ke['seke'], 0), 'b', label='Stationary Eddies')
    lp_mke = ax.plot(grid['lat'], 1e-6*np.mean(ke['mke'], 0), 'g', label='Mean Meridional Circulation')
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Vertically Integrated Kinetic Energy (MJ m$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    # ax.set_ylim(divin_dev,divax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

