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
plt.rcParams["figure.figsize"] = (4, 3)

def diffv_lat(sim, **kwargs):

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
    spread = kwargs.get('spread', 'prc')

    # diffname = 'tdiffv92500'
    # transname = 'aht'
    # translabel = r'$\langle [ \overline{vm} ] \rangle$ (PW)'

    diffname = 'diffv92500'
    transname = 'vmte'
    translabel = r'$\langle [ \overline{v^{\prime} m^{\prime}} ] \rangle$ (PW)'

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
        [diffv, grid, datadir, plotdir, modelstr] = load_diffv(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        [diffv_ref, _, _, _, _] = load_diffv(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)
    else:
        [diffv, grid, datadir, plotdir, modelstr, diffv_mmm] = load_diffv(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        [diffv_ref, _, _, _, _, diffv_mmm] = load_diffv(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model, yr_span=yr_span_ref)

    if seas == 'djf':
        # r1 = np.mean(np.roll(r1, 1, axis=0)[0:3,...], axis=0)
        # r1_ref = np.mean(np.roll(r1_ref, 1, axis=0)[0:3,...], axis=0)
        # for fluxname in flux:
        #     flux[fluxname] = np.mean(np.roll(flux[fluxname], 1, axis=0)[0:3,...], axis=0)
        #     flux_ref[fluxname] = np.mean(np.roll(flux_ref[fluxname], 1, axis=0)[0:3,...], axis=0)
        for diffvname in diffv:
            diffv[diffvname] = np.mean(np.roll(diffv[diffvname], 1, axis=0)[0:3,...], axis=0)
            diffv_ref[diffvname] = np.mean(np.roll(diffv_ref[diffvname], 1, axis=0)[0:3,...], axis=0)
            if not ( isinstance(model, str) or model is None ):
                for statname in diffv_mmm[diffvname]:
                    diffv_mmm[diffvname][statname] = np.mean(np.roll(diffv_mmm[diffvname][statname], 1, axis=0)[0:3,...], axis=0)
    elif seas == 'mam':
        # r1 = np.mean(r1[2:5,...], axis=0)
        # r1_ref = np.mean(r1_ref[2:5,...], axis=0)
        # for fluxname in flux:
        #     flux[fluxname] = np.mean(flux[fluxname][2:5,...], axis=0)
        #     flux_ref[fluxname] = np.mean(flux_ref[fluxname][2:5,...], axis=0)
        for diffvname in diffv:
            diffv[diffvname] = np.mean(diffv[diffvname][2:5,...], axis=0)
            diffv_ref[diffvname] = np.mean(diffv_ref[diffvname][2:5,...], axis=0)
            if not ( isinstance(model, str) or model is None ):
                for statname in diffv_mmm[diffvname]:
                    diffv_mmm[diffvname][statname] = np.mean(diffv_mmm[diffvname][statname][2:5,...], axis=0)
    elif seas == 'jja':
        # r1 = np.mean(r1[5:8,...], axis=0)
        # r1_ref = np.mean(r1_ref[5:8,...], axis=0)
        # for fluxname in flux:
        #     flux[fluxname] = np.mean(flux[fluxname][5:8,...], axis=0)
        #     flux_ref[fluxname] = np.mean(flux_ref[fluxname][5:8,...], axis=0)
        for diffvname in diffv:
            diffv[diffvname] = np.mean(diffv[diffvname][5:8,...], axis=0)
            diffv_ref[diffvname] = np.mean(diffv_ref[diffvname][5:8,...], axis=0)
            if not ( isinstance(model, str) or model is None ):
                for statname in diffv_mmm[diffvname]:
                    diffv_mmm[diffvname][statname] = np.mean(diffv_mmm[diffvname][statname][5:8,...], axis=0)
    elif seas == 'son':
        # r1 = np.mean(r1[8:11,...], axis=0)
        # r1_ref = np.mean(r1_ref[8:11,...], axis=0)
        # for fluxname in flux:
        #     flux[fluxname] = np.mean(flux[fluxname][8:11,...], axis=0)
        #     flux_ref[fluxname] = np.mean(flux_ref[fluxname][8:11,...], axis=0)
        for diffvname in diffv:
            diffv[diffvname] = np.mean(diffv[diffvname][8:11,...], axis=0)
            diffv_ref[diffvname] = np.mean(diffv_ref[diffvname][8:11,...], axis=0)
            if not ( isinstance(model, str) or model is None ):
                for statname in diffv_mmm[diffvname]:
                    diffv_mmm[diffvname][statname] = np.mean(diffv_mmm[diffvname][statname][8:11,...], axis=0)
    elif seas == '':
        # r1 = np.mean(r1, axis=0)
        # r1_ref = np.mean(r1_ref, axis=0)
        # for fluxname in flux:
        #     flux[fluxname] = np.mean(flux[fluxname], axis=0)
        #     flux_ref[fluxname] = np.mean(flux_ref[fluxname], axis=0)
        for diffvname in diffv:
            diffv[diffvname] = np.mean(diffv[diffvname], axis=0)
            diffv_ref[diffvname] = np.mean(diffv_ref[diffvname], axis=0)
            if not ( isinstance(model, str) or model is None ):
                for statname in diffv_mmm[diffvname]:
                    diffv_mmm[diffvname][statname] = np.mean(diffv_mmm[diffvname][statname], axis=0)

    ############################################
    # PLOT diffusivity
    ############################################
    plotname = remove_repdots('%s/%s.%s.%s' % (plotdir, diffname, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(grid['lat'], 1e-9*diffv_mmm[diffname]['prc25'], 1e-9*diffv_mmm[diffname]['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(grid['lat'], 1e-9*diffv_mmm[diffname]['mmm']-1e-9*diffv_mmm[diffname]['std'], 1e-9*diffv_mmm[diffname]['mmm']+1e-9*ga_dev_vint_hl_mmm[diffname]['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    ax.plot(grid['lat'], 1e-9*diffv[diffname], '-k', label='D')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('D ($10^9$ kg s$^{-1}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    ax.set_ylim([0,20])
    ax.set_xticks(np.arange(-90,91,30))
    # ax.set_ylim(divin_dev,divax_dev)
    # if legend:
    #     ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT diffusivity (SH zoom)
    ############################################
    plotname = remove_repdots('%s/%s_sh.%s.%s' % (plotdir, diffname, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(grid['lat'], 1e-9*diffv_mmm[diffname]['prc25'], 1e-9*diffv_mmm[diffname]['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(grid['lat'], 1e-9*diffv_mmm[diffname]['mmm']-1e-9*diffv_mmm[diffname]['std'], 1e-9*diffv_mmm[diffname]['mmm']+1e-9*ga_dev_vint_hl_mmm[diffname]['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    ax.plot(grid['lat'], 1e-9*diffv[diffname], '-k', label='D')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('D ($10^9$ kg s$^{-1}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-60,-20])
    ax.set_xticks(np.arange(-60,-19,5))
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylim([0,20])
    # ax.set_ylim(divin_dev,divax_dev)
    # if legend:
    #     ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT 925 hPa MSE gradient
    ############################################
    plotname = remove_repdots('%s/gmse92500.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(grid['lat'], 1e-3*diffv_mmm['gmse92500']['prc25'], 1e-3*diffv_mmm['gmse92500']['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(grid['lat'], 1e-3*diffv_mmm['gmse92500']['mmm']-1e-3*diffv_mmm['gmse92500']['std'], 1e-3*diffv_mmm['gmse92500']['mmm']+1e-3*ga_dev_vint_hl_mmm['gmse92500']['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    ax.plot(grid['lat'], 1e-3*diffv['gmse92500'], '-k', label='D')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel(r'$\partial_{\phi} m_{925\,\mathrm{hPa}}$ (J g$^{-1}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    ax.set_xticks(np.arange(-90,91,30))
    # ax.set_ylim(divin_dev,divax_dev)
    # if legend:
    #     ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT vmte
    ############################################
    plotname = remove_repdots('%s/%s.%s.%s' % (plotdir, transname, seas, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(grid['lat'], 1e-15*diffv_mmm[transname]['prc25'], 1e-15*diffv_mmm[transname]['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(grid['lat'], 1e-15*diffv_mmm[transname]['mmm']-1e-15*diffv_mmm[transname]['std'], 1e-15*diffv_mmm[transname]['mmm']+1e-15*ga_dev_vint_hl_mmm[transname]['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    ax.plot(grid['lat'], 1e-15*diffv[transname], '-k', label='D')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    # ax.set_ylabel(r'$\langle [ \overline{v^{\prime} m^{\prime}} ] \rangle$ (PW)')
    ax.set_ylabel(translabel)
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    ax.set_xticks(np.arange(-90,91,30))
    # ax.set_ylim(divin_dev,divax_dev)
    # if legend:
    #     ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT 925 hPa MSE 
    ############################################
    plotname = remove_repdots('%s/mse92500.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots()
    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(grid['lat'], 1e-3*diffv_mmm['mse92500']['prc25'], 1e-3*diffv_mmm['mse92500']['prc75'], facecolor='k', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(grid['lat'], 1e-3*diffv_mmm['mse92500']['mmm']-1e-3*diffv_mmm['mse92500']['std'], 1e-3*diffv_mmm['mse92500']['mmm']+1e-3*ga_dev_vint_hl_mmm['mse92500']['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    ax.plot(grid['lat'], 1e-3*diffv['mse92500'], '-k', label='m925')
    make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel(r'$m_{925\,\mathrm{hPa}}$ (J g$^{-1}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    ax.set_xticks(np.arange(-90,91,30))
    # ax.set_ylim(divin_dev,divax_dev)
    # if legend:
    #     ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

