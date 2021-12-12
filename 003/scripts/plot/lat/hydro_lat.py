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

def hydro_lat(sim, **kwargs):

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
        [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    ############################################
    # PLOT (E/P, PRFRAC)
    ############################################
    plotname = remove_repdots('%s/e.p.annmean.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axhline(1, color='k', linewidth=0.5)
    # if not (isinstance(model, str) or model is None):
    #     ax.fill_between(grid['lat'], np.mean(hydro['prfrac']-hydro_mmm['prfrac']['std'],0), np.mean(hydro['prfrac']+hydro_mmm['prfrac']['std'],0), color='k', alpha=0.1)
    lp_prc = ax.plot(grid['lat'], np.mean(hydro['prfrac'], 0), '-k', label='$P_c/P$')
    # lp_ep = ax.plot(grid['lat'], np.mean(hydro['evspsbl']/hydro['pr'], 0), '--k', label='E/P')
    lp_ep = ax.plot(grid['lat'], np.mean(hydro['evspsbl']/hydro['pr'], 0), '--k', label='$E/P$')
    # make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.set_title('CMIP5 mean, historical')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    # ax.set_ylabel('Energy flux divergence (W m$^{-2}$)')
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
    # PLOT (R2, PRFRAC)
    ############################################
    plotname = remove_repdots('%s/r2.p.annmean.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axhline(1, color='k', linewidth=0.5)
    # if not (isinstance(model, str) or model is None):
    #     ax.fill_between(grid['lat'], np.mean(hydro['prfrac']-hydro_mmm['prfrac']['std'],0), np.mean(hydro['prfrac']+hydro_mmm['prfrac']['std'],0), color='k', alpha=0.1)
    lp_prc = ax.plot(grid['lat'], np.mean(hydro['prfrac'], 0), '-k', label='$P_c/P$')
    # lp_ep = ax.plot(grid['lat'], np.mean(hydro['evspsbl']/hydro['pr'], 0), '--k', label='E/P')
    lp_ep = ax.plot(grid['lat'], np.mean(-(r1-1), 0), '--k', label='$-R_2$')
    # make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.set_title('CMIP5 mean, historical')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    # ax.set_ylabel('Energy flux divergence (W m$^{-2}$)')
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
    # PLOT (SCALED R2, PRFRAC)
    ############################################
    plotname = remove_repdots('%s/r2.scaled.p.annmean.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axhline(1, color='k', linewidth=0.5)
    # if not (isinstance(model, str) or model is None):
    #     ax.fill_between(grid['lat'], np.mean(hydro['prfrac']-hydro_mmm['prfrac']['std'],0), np.mean(hydro['prfrac']+hydro_mmm['prfrac']['std'],0), color='k', alpha=0.1)
    lp_prc = ax.plot(grid['lat'], np.mean(hydro['prfrac'], 0), '-k', label='$P_c/P$')
    # lp_ep = ax.plot(grid['lat'], np.mean(hydro['evspsbl']/hydro['pr'], 0), '--k', label='E/P')
    lp_ep = ax.plot(grid['lat'], 0.6*np.mean(-(r1-1), 0), '--k', label='$-\eta R_2$')
    # make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.set_title('CMIP5 mean, historical')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    # ax.set_ylabel('Energy flux divergence (W m$^{-2}$)')
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
    # PLOT (INFERRED P FROM RCE TERMS ONLY)
    ############################################
    plotname = remove_repdots('%s/rce.p.annmean.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_prc = ax.plot(grid['lat'], np.mean(par.Lv*hydro['prc']/86400, 0), '-k', label='$LP_c$')
    # lp_ep = ax.plot(grid['lat'], np.mean(hydro['evspsbl']/hydro['pr'], 0), '--k', label='E/P')
    lp_ep = ax.plot(grid['lat'], np.mean(-flux['ra']-flux['hfss'], 0), '--k', label='$-R_a-\mathrm{SH}$')
    ax.set_title('CMIP5 mean, historical')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    # ax.set_ylabel('Energy flux divergence (W m$^{-2}$)')
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
    # PLOT (INFERRED P FROM RAE TERMS ONLY)
    ############################################
    plotname = remove_repdots('%s/rae.p.annmean.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_prc = ax.plot(grid['lat'], np.mean(par.Lv*hydro['prc']/86400, 0), '-k', label='$LP_c$')
    # lp_ep = ax.plot(grid['lat'], np.mean(hydro['evspsbl']/hydro['pr'], 0), '--k', label='E/P')
    lp_ep = ax.plot(grid['lat'], np.mean(-flux['ra']+flux['stg_adv_dse'], 0), '--k', label='$-R_a+\partial_y(vs)}$')
    ax.set_title('CMIP5 mean, historical')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    # ax.set_ylabel('Energy flux divergence (W m$^{-2}$)')
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
