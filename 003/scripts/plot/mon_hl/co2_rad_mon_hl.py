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
from plot.preamble import get_predata
import os
import pickle
import xarray as xr
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def co2_rad_mon_hl(sim, **kwargs):

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
    spread = kwargs.get('spread', 'prc')

    model, yr_span, yr_base = get_predata(sim, timemean, kwargs)

    if sim == 'echam':
        model_ref = refclim
    else:
        model_ref = model

    lat_int = np.arange(latbnd[0], latbnd[1], latstep)

    if latbnd[0] > 0: # NH
        if timemean == 'djfmean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
            vmin = -200
            vmax = 100
            vmin_dev = -30
            vmax_dev = 5
            vmin_dev_olr = -60
            vmax_dev_olr = 5
        elif timemean == 'jjamean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
            vmin = -150
            vmax = 50
            vmin_dev = -30
            vmax_dev = 30
    else: # SH
        if timemean == 'djfmean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
            vmin = -200
            vmax = 100
            vmin_dev = -20
            vmax_dev = 5
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
    else:
        [rad, grid, datadir, plotdir, modelstr, rad_mmm] = load_rad(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    if not refclim == '':
        if isinstance(model, str):
            [rad_ref, _, _, _, _] = load_rad(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model_ref, yr_span=yr_span_ref)
        else:
            [rad_ref, _, _, _, _, rad_ref_mmm] = load_rad(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, try_load=try_load, model=model_ref, yr_span=yr_span_ref)

        if timemean == 'djfmean':
            for radname in rad_ref:
                rad_ref[radname] = np.mean(np.roll(rad_ref[radname],1,axis=0)[0:3], 0)
                if not ( isinstance(model, str) or (model is None) ):
                    for stat in rad_ref_mmm[radname]:
                        rad_ref_mmm[radname][stat] = np.mean(np.roll(rad_ref_mmm[radname][stat],1,axis=0)[0:3], 0)
        elif timemean == 'jjamean':
            for radname in rad_ref:
                rad_ref[radname] = np.mean(rad_ref[radname][5:8], 0)
                if not ( isinstance(model, str) or (model is None) ):
                    for stat in rad_ref_mmm[radname]:
                        rad_ref_mmm[radname][stat] = np.mean(rad_ref_mmm[radname][stat][5:8], 0)

        rad_ref_hl = {}

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

    
    if not ( isinstance(model, str) or (model is None) ):
        rad_hl_mmm = {}
        rad_ref_hl_mmm = {}
        rad_dev_hl_mmm = {}

        for radname in rad:
            rad_hl_mmm[radname] = {}
            rad_ref_hl_mmm[radname] = {}
            rad_dev_hl_mmm[radname] = {}
            
            for stat in rad_mmm[radname]:
                rad_hl_mmm[radname][stat] = lat_mean(rad_mmm[radname][stat], grid, lat_int, dim=1)

                if refclim == 'init':
                    rad_dev_hl_mmm[radname][stat] = rad_hl_mmm[radname][stat] - rad_hl_mmm[radname][stat][0]
                else:
                    rad_ref_hl_mmm[radname][stat] = lat_mean(rad_ref_mmm[radname][stat], grid, lat_int, dim=0)
                    rad_dev_hl_mmm[radname][stat] = rad_hl_mmm[radname][stat] - rad_ref_hl_mmm[radname][stat]

    time = yr_base + np.arange(rad_hl['ra'].shape[0]) # create time vector

    ############################################
    # Load CO2 timeseries
    ############################################
    print(time)
    ghg = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg_rcp85_1765-2500_c100203.nc')
    co2 = 1e-6*ghg.CO2.sel(time=slice(str(time[0]),str(time[-1]))).values
    lco2 = np.log10(co2)

    ############################################
    # Compute axis endpoints
    ############################################
    ra_f30 = np.mean(rad_dev_hl['ra'][:30])
    ra_l30 = np.mean(rad_dev_hl['ra'][-30:])

    # log co2
    lco2_f30 = np.mean(lco2[:30])
    lco2_l30 = np.mean(lco2[-30:])
    m_axis = (lco2_l30 - lco2_f30)/(ra_l30 - ra_f30)
    vmax_alg_lco2 = lco2_f30 + m_axis*( vmax_dev - ra_f30 )
    vmin_alg_lco2 = lco2_l30 + m_axis*( vmin_dev - ra_l30 )

    # co2
    co2_f30 = np.mean(co2[:30])
    co2_l30 = np.mean(co2[-30:])
    m_axis = (co2_l30 - co2_f30)/(ra_l30 - ra_f30)
    vmax_alg_co2 = co2_f30 + m_axis*( vmax_dev - ra_f30 )
    vmin_alg_co2 = co2_l30 + m_axis*( vmin_dev - ra_l30 )

    ############################################
    # Compute pgrayicted ra_dev
    ############################################
    df_tlcl = -1 # W m**-2 K**-1, radiative cooling flux divergence in T coord at LCL
    # dra = df_tlcl * hydro_dev_hl['tas']
    # dra_mod = df_tlcl * hydro_dev_hl['t850']

    ############################################
    # PLOT (W/ LOG CO2 ALL SKY RAD, DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/lco2_ra_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            spr_ra = ax.fill_between(time, rad_dev_hl_mmm['ra']['prc25'], rad_dev_hl_mmm['ra']['prc75'], facecolor='black', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            spr_ra = ax.fill_between(time, rad_dev_hl_mmm['ra']['mmm'] - rad_dev_hl_mmm['ra']['std'], rad_dev_hl_mmm['ra']['mmm'] + rad_dev_hl_mmm['ra']['std'], facecolor='black', alpha=0.2, edgecolor=None)
    lp_ra = ax.plot(time, rad_dev_hl['ra'], color='tab:gray', label='$\Delta R_a$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('$\Delta R_a$ (Wm$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+rad_dev_hl['ra'].shape[0]-1)
    ax.set_ylim(vmin_dev,vmax_dev)

    sax = ax.twinx()
    sax.plot(time, lco2, color='tab:green', label='$\log_{10}(p_{CO_2})$')
    sax.set_ylabel('$\log_{10}(p_{CO_2})$', color='tab:green')
    # sax.set_ylim(vmin['sic'],vmax['sic'])
    sax.set_ylim(vmin_alg_lco2,vmax_alg_lco2)
    sax.tick_params(axis='y', labelcolor='tab:green', color='tab:green')
    sax.yaxis.set_minor_locator(MultipleLocator(5))

    # if legend:
    #     ax.legend()
    fig.set_size_inches(4,3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (W/ CO2, ALL SKY RAD, DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/co2_ra_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            spr_ra = ax.fill_between(time, rad_dev_hl_mmm['ra']['prc25'], rad_dev_hl_mmm['ra']['prc75'], facecolor='black', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            spr_ra = ax.fill_between(time, rad_dev_hl_mmm['ra']['mmm'] - rad_dev_hl_mmm['ra']['std'], rad_dev_hl_mmm['ra']['mmm'] + rad_dev_hl_mmm['ra']['std'], facecolor='black', alpha=0.2, edgecolor=None)
    lp_ra = ax.plot(time, rad_dev_hl['ra'], color='tab:gray', label='$\Delta R_a$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('$\Delta R_a$ (Wm$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+rad_dev_hl['ra'].shape[0]-1)
    ax.set_ylim(vmin_dev,vmax_dev)

    sax = ax.twinx()
    sax.plot(time, co2, color='tab:green', label='$p_{CO_2}$')
    sax.set_ylabel('$p_{CO_2}$', color='tab:green')
    # sax.set_ylim(vmin['sic'],vmax['sic'])
    sax.set_ylim(vmin_alg_co2,vmax_alg_co2)
    sax.tick_params(axis='y', labelcolor='tab:green', color='tab:green')
    sax.yaxis.set_minor_locator(MultipleLocator(5))

    # if legend:
    #     ax.legend()
    fig.set_size_inches(4,3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (W/ LOG CO2 CLEAR SKY RAD, DEVIATION FROM INITIAL)
    ############################################
    ra_cs_f30 = np.mean(rad_dev_hl['ra_cs'][:30])
    ra_cs_l30 = np.mean(rad_dev_hl['ra_cs'][-30:])

    # log co2
    lco2_f30 = np.mean(lco2[:30])
    lco2_l30 = np.mean(lco2[-30:])
    m_axis = (lco2_l30 - lco2_f30)/(ra_cs_l30 - ra_cs_f30)
    vmax_alg_lco2cs = lco2_f30 + m_axis*( vmax_dev - ra_cs_f30 )
    vmin_alg_lco2cs = lco2_l30 + m_axis*( vmin_dev - ra_cs_l30 )

    plotname = remove_repdots('%s/lco2_ra_cs_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            spr_ra = ax.fill_between(time, rad_dev_hl_mmm['ra_cs']['prc25'], rad_dev_hl_mmm['ra_cs']['prc75'], facecolor='tab:gray', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            spr_ra = ax.fill_between(time, rad_dev_hl_mmm['ra_cs']['mmm'] - rad_dev_hl_mmm['ra_cs']['std'], rad_dev_hl_mmm['ra_cs']['mmm'] + rad_dev_hl_mmm['ra_cs']['std'], facecolor='tab:gray', alpha=0.2, edgecolor=None)
    lp_ra = ax.plot(time, rad_dev_hl['ra_cs'], color='tab:gray', label='$\Delta R_{a,\,cs}$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('$\Delta R_{a,\,cs}$ (Wm$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+rad_dev_hl['ra_cs'].shape[0]-1)
    ax.set_ylim(vmin_dev,vmax_dev)

    sax = ax.twinx()
    sax.plot(time, lco2, color='tab:green', label='$\log_{10}(p_{CO_2})$')
    sax.set_ylabel('$\log_{10}(p_{CO_2})$', color='tab:green')
    # sax.set_ylim(vmin['sic'],vmax['sic'])
    sax.set_ylim(vmin_alg_lco2cs,vmax_alg_lco2cs)
    sax.tick_params(axis='y', labelcolor='tab:green', color='tab:green')
    sax.yaxis.set_minor_locator(MultipleLocator(5))

    # if legend:
    #     ax.legend()
    fig.set_size_inches(4,3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (W/ LOG CO2 CLEAR SKY RAD, DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/lco2_ra_cs_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            spr_ra = ax.fill_between(time, rad_dev_hl_mmm['ra_cs']['prc25'], rad_dev_hl_mmm['ra_cs']['prc75'], facecolor='tab:gray', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            spr_ra = ax.fill_between(time, rad_dev_hl_mmm['ra_cs']['mmm'] - rad_dev_hl_mmm['ra_cs']['std'], rad_dev_hl_mmm['ra_cs']['mmm'] + rad_dev_hl_mmm['ra_cs']['std'], facecolor='tab:gray', alpha=0.2, edgecolor=None)
    lp_ra = ax.plot(time, rad_dev_hl['ra_cs'], color='tab:gray', label='$\Delta R_{a,\,cs}$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('$\Delta R_{a,\,cs}$ (Wm$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+rad_dev_hl['ra_cs'].shape[0]-1)
    ax.set_ylim(vmin_dev,vmax_dev)

    sax = ax.twinx()
    sax.plot(time, lco2, color='tab:green', label='$\log_{10}(p_{CO_2})$')
    sax.set_ylabel('$\log_{10}(p_{CO_2})$', color='tab:green')
    # sax.set_ylim(vmin['sic'],vmax['sic'])
    sax.set_ylim(vmin_alg_lco2,vmax_alg_lco2)
    sax.tick_params(axis='y', labelcolor='tab:green', color='tab:green')
    sax.yaxis.set_minor_locator(MultipleLocator(5))

    # if legend:
    #     ax.legend()
    fig.set_size_inches(4,3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (W/ LOG CO2 OLR, DEVIATION FROM INITIAL)
    ############################################
    rad_dev_hl['rlut'] = -rad_dev_hl['rlut']
    rlut_f30 = np.mean(rad_dev_hl['rlut'][:30])
    rlut_l30 = np.mean(rad_dev_hl['rlut'][-30:])

    # log co2
    lco2_f30 = np.mean(lco2[:30])
    lco2_l30 = np.mean(lco2[-30:])
    m_axis = (lco2_l30 - lco2_f30)/(rlut_l30 - rlut_f30)
    vmax_alg_lco2olr = lco2_f30 + m_axis*( vmax_dev - rlut_f30 )
    vmin_alg_lco2olr = lco2_l30 + m_axis*( vmin_dev - rlut_l30 )

    plotname = remove_repdots('%s/lco2_olr_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            spr_ra = ax.fill_between(time, -rad_dev_hl_mmm['rlut']['prc25'], -rad_dev_hl_mmm['rlut']['prc75'], facecolor='tab:gray', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            spr_ra = ax.fill_between(time, -(rad_dev_hl_mmm['rlut']['mmm'] - rad_dev_hl_mmm['rlut']['std']), -(rad_dev_hl_mmm['rlut']['mmm'] + rad_dev_hl_mmm['rlut']['std']), facecolor='tab:gray', alpha=0.2, edgecolor=None)
    lp_ra = ax.plot(time, rad_dev_hl['rlut'], color='tab:gray', label='$\Delta R_{a,\,cs}$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('$\Delta$ OLR (Wm$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+rad_dev_hl['rlut'].shape[0]-1)
    ax.set_ylim(vmin_dev,vmax_dev)

    sax = ax.twinx()
    sax.plot(time, lco2, color='tab:green', label='$\log_{10}(p_{CO_2})$')
    sax.set_ylabel('$\log_{10}(p_{CO_2})$', color='tab:green')
    # sax.set_ylim(vmin['sic'],vmax['sic'])
    sax.set_ylim(vmin_alg_lco2olr,vmax_alg_lco2olr)
    sax.tick_params(axis='y', labelcolor='tab:green', color='tab:green')
    sax.yaxis.set_minor_locator(MultipleLocator(5))

    # if legend:
    #     ax.legend()
    fig.set_size_inches(4,3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (W/ LOG CO2 OLR CS, DEVIATION FROM INITIAL)
    ############################################
    rad_dev_hl['rlutcs'] = -rad_dev_hl['rlutcs']
    rlutcs_f30 = np.mean(rad_dev_hl['rlutcs'][:30])
    rlutcs_l30 = np.mean(rad_dev_hl['rlutcs'][-30:])

    # log co2
    lco2_f30 = np.mean(lco2[:30])
    lco2_l30 = np.mean(lco2[-30:])
    m_axis = (lco2_l30 - lco2_f30)/(rlutcs_l30 - rlutcs_f30)
    vmax_alg_lco2olrcs = lco2_f30 + m_axis*( vmax_dev_olr - rlutcs_f30 )
    vmin_alg_lco2olrcs = lco2_l30 + m_axis*( vmin_dev_olr - rlutcs_l30 )

    plotname = remove_repdots('%s/lco2_olrcs_dev_mon_hl.%g.%g.%s' % (plotdir, latbnd[0], latbnd[1], timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    if not ( isinstance(model, str) or (model is None) ):
        if spread == 'prc':
            spr_ra = ax.fill_between(time, -rad_dev_hl_mmm['rlutcs']['prc25'], -rad_dev_hl_mmm['rlutcs']['prc75'], facecolor='tab:gray', alpha=0.2, edgecolor=None)
        elif spread == 'std':
            spr_ra = ax.fill_between(time, -(rad_dev_hl_mmm['rlutcs']['mmm'] - rad_dev_hl_mmm['rlutcs']['std']), -(rad_dev_hl_mmm['rlutcs']['mmm'] + rad_dev_hl_mmm['rlutcs']['std']), facecolor='tab:gray', alpha=0.2, edgecolor=None)
    lp_ra = ax.plot(time, rad_dev_hl['rlutcs'], color='tab:gray', label='$\Delta R_{a,\,cs}$')
    make_title_sim_time_lat(ax, sim, model=modelstr, timemean=timemean, lat1=latbnd[0], lat2=latbnd[1])
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Time (yr)')
    ax.set_ylabel('$\Delta$ OLR$_{cs}$ (Wm$^{-2}$)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(yr_base,yr_base+rad_dev_hl['rlutcs'].shape[0]-1)
    ax.set_ylim(vmin_dev_olr,vmax_dev_olr)

    sax = ax.twinx()
    sax.plot(time, lco2, color='tab:green', label='$\log_{10}(p_{CO_2})$')
    sax.set_ylabel('$\log_{10}(p_{CO_2})$', color='tab:green')
    # sax.set_ylim(vmin['sic'],vmax['sic'])
    sax.set_ylim(vmin_alg_lco2olrcs,vmax_alg_lco2olrcs)
    sax.tick_params(axis='y', labelcolor='tab:green', color='tab:green')
    sax.yaxis.set_minor_locator(MultipleLocator(5))

    # if legend:
    #     ax.legend()
    fig.set_size_inches(4,3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

