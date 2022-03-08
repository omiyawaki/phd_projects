import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
import numpy as np
from misc.load_data import *
from misc.filenames import *
from misc.dirnames import *
from misc.translate import *
from plot.titles import *
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
    seas = kwargs.get('seas', 'djf') # view plot? (plt.show)
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
    if refclim == 'hist-30':
        sim_ref = 'historical'
        timemean_ref = 'ymonmean-30'
        yr_span_ref = '186001-200512'

    if isinstance(model, str) or model is None:
        [r1, grid, datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
        [r1_ref, _, _, _, _] = load_r1(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref, refclim=refclim) 
    else:
        [r1, grid, datadir, plotdir, modelstr, r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
        [r1_ref, _, _, _, _, r1_ref_mmm] = load_r1(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref, refclim=refclim) 

    if isinstance(model, str) or model is None:
        [flux, grid, datadir, plotdir, modelstr] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
        [flux_ref, _, _, _, _] = load_flux(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref, refclim=refclim) 
    else:
        [flux, grid, datadir, plotdir, modelstr, flux_mmm] = load_flux(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
        [flux_ref, _, _, _, _, flux_ref_mmm] = load_flux(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref, refclim=refclim) 

    if isinstance(model, str) or model is None:
        [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        [hydro_ref, _, _, _, _] = load_hydro(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
    else:
        [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        [hydro_ref, _, _, _, _, hydro_mmm] = load_hydro(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)

    if seas == 'djf':
        r1 = np.mean(np.roll(r1, 1, axis=0)[0:3,...], axis=0)
        r1_ref = np.mean(np.roll(r1_ref, 1, axis=0)[0:3,...], axis=0)
        for fluxname in flux:
            flux[fluxname] = np.mean(np.roll(flux[fluxname], 1, axis=0)[0:3,...], axis=0)
            flux_ref[fluxname] = np.mean(np.roll(flux_ref[fluxname], 1, axis=0)[0:3,...], axis=0)
        for hydroname in hydro:
            hydro[hydroname] = np.mean(np.roll(hydro[hydroname], 1, axis=0)[0:3,...], axis=0)
            hydro_ref[hydroname] = np.mean(np.roll(hydro_ref[hydroname], 1, axis=0)[0:3,...], axis=0)
    elif seas == 'mam':
        r1 = np.mean(r1[2:5,...], axis=0)
        r1_ref = np.mean(r1_ref[2:5,...], axis=0)
        for fluxname in flux:
            flux[fluxname] = np.mean(flux[fluxname][2:5,...], axis=0)
            flux_ref[fluxname] = np.mean(flux_ref[fluxname][2:5,...], axis=0)
        for hydroname in hydro:
            hydro[hydroname] = np.mean(hydro[hydroname][2:5,...], axis=0)
            hydro_ref[hydroname] = np.mean(hydro_ref[hydroname][2:5,...], axis=0)
    elif seas == 'jja':
        r1 = np.mean(r1[5:8,...], axis=0)
        r1_ref = np.mean(r1_ref[5:8,...], axis=0)
        for fluxname in flux:
            flux[fluxname] = np.mean(flux[fluxname][5:8,...], axis=0)
            flux_ref[fluxname] = np.mean(flux_ref[fluxname][5:8,...], axis=0)
        for hydroname in hydro:
            hydro[hydroname] = np.mean(hydro[hydroname][5:8,...], axis=0)
            hydro_ref[hydroname] = np.mean(hydro_ref[hydroname][5:8,...], axis=0)
    elif seas == 'son':
        r1 = np.mean(r1[8:11,...], axis=0)
        r1_ref = np.mean(r1_ref[8:11,...], axis=0)
        for fluxname in flux:
            flux[fluxname] = np.mean(flux[fluxname][8:11,...], axis=0)
            flux_ref[fluxname] = np.mean(flux_ref[fluxname][8:11,...], axis=0)
        for hydroname in hydro:
            hydro[hydroname] = np.mean(hydro[hydroname][8:11,...], axis=0)
            hydro_ref[hydroname] = np.mean(hydro_ref[hydroname][8:11,...], axis=0)
    elif seas == 'ann':
        r1 = np.mean(r1, axis=0)
        r1_ref = np.mean(r1_ref, axis=0)
        for fluxname in flux:
            flux[fluxname] = np.mean(flux[fluxname], axis=0)
            flux_ref[fluxname] = np.mean(flux_ref[fluxname], axis=0)
        for hydroname in hydro:
            hydro[hydroname] = np.mean(hydro[hydroname], axis=0)
            hydro_ref[hydroname] = np.mean(hydro_ref[hydroname], axis=0)

    ############################################
    # PLOT (E/P, PRFRAC)
    ############################################
    plotname = remove_repdots('%s/dprfrac.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    ax.plot(grid['lat'], hydro['prfrac']-hydro_ref['prfrac'], '-k', label='$P_c/P$')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
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
    # PLOT (E/P, PRFRAC)
    ############################################
    plotname = remove_repdots('%s/e.p.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axhline(1, color='k', linewidth=0.5)
    lp_prc = ax.plot(grid['lat'], hydro['prfrac'], '-k', label='$P_c/P$')
    lp_ep = ax.plot(grid['lat'], hydro['evspsbl']/hydro['pr'], '--k', label='$E/P$')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
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
    plotname = remove_repdots('%s/r2.p.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axhline(1, color='k', linewidth=0.5)
    lp_prc = ax.plot(grid['lat'], hydro['prfrac'], '-k', label='$P_c/P$')
    lp_ep = ax.plot(grid['lat'], -(r1-1), '--k', label='$-R_2$')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
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
    plotname = remove_repdots('%s/r2.scaled.p.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axhline(1, color='k', linewidth=0.5)
    lp_prc = ax.plot(grid['lat'], hydro['prfrac'], '-k', label='$P_c/P$')
    lp_ep = ax.plot(grid['lat'], 0.6*-(r1-1), '--k', label='$-\eta R_2$')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
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
    plotname = remove_repdots('%s/rce.p.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_prc = ax.plot(grid['lat'], par.Lv*hydro['prc']/86400, '-k', label='$LP_c$')
    lp_ep = ax.plot(grid['lat'], -flux['ra']-flux['hfss'], '--k', label='$-R_a-\mathrm{SH}$')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
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
    plotname = remove_repdots('%s/rae.p.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_prc = ax.plot(grid['lat'], par.Lv*hydro['prc']/86400, '-k', label='$LP_c$')
    lp_ep = ax.plot(grid['lat'], -flux['ra']+flux['stg_adv_dse'], '--k', label='$-R_a+\partial_y(vs)}$')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
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
