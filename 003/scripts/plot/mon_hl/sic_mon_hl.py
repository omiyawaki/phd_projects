import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import *
from misc.filenames import *
# from misc.translate import translate_varname
from misc.means import lat_mean, global_int
from misc.load_data import *
from misc import par
from plot.preamble import get_predata
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

def sic_mon_hl(sim, **kwargs):

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

    # if sim == 'longrun':
    #     model = kwargs.get('model', 'MPIESM12_abrupt4x')
    #     yr_span = kwargs.get('yr_span', '1000')
    #     yr_base = 0
    # elif sim == 'rcp85':
    #     model = kwargs.get('model', 'MPI-ESM-LR')
    #     yr_span = kwargs.get('yr_span', '200601-230012')
    #     if 'ymonmean' not in timemean:
    #         yr_base = 2006
    #     else:
    #         yr_base = 0
    # elif sim == 'historical':
    #     model = kwargs.get('model', 'MPI-ESM-LR')
    #     yr_span = kwargs.get('yr_span', '186001-200512')
    #     if 'ymonmean' not in timemean:
    #         yr_base = 1860
    #     else:
    #         yr_base = 0
    # elif sim == 'echam':
    #     model = kwargs.get('model', 'rp000140')
    #     yr_span = kwargs.get('yr_span', '0001_0039')
    #     yr_base = 40
    # elif sim == 'era5':
    #     model = None
    #     yr_span = kwargs.get('yr_span', '1979_2019')
    #     if 'ymonmean' not in timemean:
    #         yr_base = int(yr_span[0:4])
    #     else:
    #         yr_base = 0

    model, yr_span, yr_base = get_predata(sim, timemean, kwargs)

    ##########################################
    ## Y AXIS SPECIFICATIONS
    ##########################################
    vmin={}
    vmax={}
    vmin['stf'] = -10
    vmax['stf'] = 80
    vmin['ra'] = -160
    vmax['ra'] = -120
    vmin['stg_adv'] = -150
    vmax['stg_adv'] = -90

    ##################################
    # LOAD DATA
    ##################################
    if isinstance(model, str) or model is None:
        [seaice, grid, datadir, plotdir, modelstr] = load_seaice(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    else:
        [seaice, grid, datadir, plotdir, modelstr, seaice_mmm] = load_seaice(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    ############################################
    # AVERAGE sea ice ONLY AT HIGH LATITUDES
    ############################################
    seaice_hl = {}
    for varname in seaice:
        seaice_hl[varname] = lat_mean(seaice[varname], grid, lat_int, dim=1)

    if not (isinstance(model, str) or model is None):
        seaice_mmm_hl = {}
        for varname in seaice_mmm:
            seaice_mmm_hl[varname] = {}
            for stat in seaice_mmm[varname]:
                seaice_mmm_hl[varname][stat] = lat_mean(seaice_mmm[varname][stat], grid, lat_int, dim=1)

    # pick out vars and sic
    sic = seaice_hl['sic']

    sic_mmm = {}
    if not (isinstance(model, str) or model is None):
        for stat in seaice_mmm_hl['sic']:
            sic_mmm[stat] = seaice_mmm_hl['sic'][stat]

    # pick out sit
    sit = seaice_hl['sit']

    sit_mmm = {}
    if not (isinstance(model, str) or model is None):
        for stat in seaice_mmm_hl['sit']:
            sit_mmm[stat] = seaice_mmm_hl['sit'][stat]

    time = yr_base + np.arange(sic.shape[0]) # create time vector

    ############################################
    # COMPUTE FIRST and LAST 30 YR AVGs (for yaxis settings)
    ############################################
    sic_f30 = np.mean(sic[:30])
    sic_l30 = np.mean(sic[-30:])

    ############################################
    # PLOT SIC
    ############################################

    plotname = '%s/sic_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()

    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(time, sic_mmm['prc25'], sic_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(time, sic_mmm['mmm']-sic_mmm['std'], sic_mmm['mmm']+sic_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.plot(time, sic, color='tab:blue')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Sea ice fraction (%)', color='tab:blue')
    # ax.set_ylim(vmin['sic'],vmax['sic'])
    # ax.set_ylim(vmin_alg,vmax_alg)
    ax.set_xlim(yr_base,yr_base+sic.shape[0]-1)
    ax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    fig.set_size_inches(4, 3.5)

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT SIT
    ############################################

    plotname = '%s/sit_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()

    if not (isinstance(model, str) or model is None):
        if spread == 'prc':
            ax.fill_between(time, sit_mmm['prc25'], sit_mmm['prc75'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            ax.fill_between(time, sit_mmm['mmm']-sit_mmm['std'], sit_mmm['mmm']+sit_mmm['std'], facecolor='tab:blue', alpha=0.2, edgecolor=None)
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.plot(time, sit, color='tab:blue')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Sea ice thickness (m)', color='tab:blue')
    # ax.set_ylim(vmin['sit'],vmax['sit'])
    # ax.set_ylim(vmin_alg,vmax_alg)
    ax.set_xlim(yr_base,yr_base+sit.shape[0]-1)
    ax.tick_params(axis='y', labelcolor='tab:blue', color='tab:blue')
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    fig.set_size_inches(4, 3.5)

    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

