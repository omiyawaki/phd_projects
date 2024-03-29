import os
import sys
import pickle
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from tools import get_modelidx, get_modelstr

lat0=80
lat1=90
models=['rp000134','rp000190f']

for model in models:

    # specify directory for output data and plots
    datadir = '%s/data' % ( os.getcwd() ) # output data directory
    plotdir = '%s/plot/%s' % ( os.getcwd(), model ) # plot directory
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)

    ##################################################
    # LOAD GCM DATA
    ##################################################
    ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/greenhouse_hist+ssp585.nc')
    co2 = 1e-6*ds.CO2.sel(time=slice('1987-01-01', '2200-12-31')).data

    [clim, grid] = pickle.load(open('../climlab/input_data/clima.%s.%g.%g.pickle'%(model,lat0,lat1), 'rb'))
    if model=='rp000134':
        [forc, grid] = pickle.load(open('../climlab/input_data/forcing.%s.%g.%g.pickle'%('rp000188',lat0,lat1), 'rb'))
    elif model=='rp000190f':
        [forc, grid] = pickle.load(open('../climlab/input_data/forcing.%s.%g.%g.pickle'%('rp000191f',lat0,lat1), 'rb'))

    ##################################################
    # LOAD RRTMG DATA
    ##################################################
    # RRTMG climatology
    ref = {}
    [ref['ra'], ref['racs'], ref['diag']] = pickle.load(open('%s/rad.348.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # RRTMG full RCP8.5 response
    full = {}
    [full['ra'], full['racs'], full['diag']] = pickle.load(open('%s/rad.ssp585.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # RRTMG fixed T and q RCP8.5 response
    fixqt = {}
    [fixqt['ra'], fixqt['racs'], fixqt['diag']] = pickle.load(open('%s/rad.ssp585_fixq_fixT.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # RRTMG fixed q and CO2 RCP8.5 response
    fixqco2 = {}
    [fixqco2['ra'], fixqco2['racs'], fixqco2['diag']] = pickle.load(open('%s/rad.ssp585_fixq_fixco2.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # RRTMG fixed q and CO2 (PLANCK) RCP8.5 response
    fixqco2_planck = {}
    [fixqco2_planck['ra'], fixqco2_planck['racs'], fixqco2_planck['diag']] = pickle.load(open('%s/rad.ssp585_fixq_fixco2_planck.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # RRTMG fixed q and CO2 (LAPSE) RCP8.5 response
    fixqco2_lapse = {}
    [fixqco2_lapse['ra'], fixqco2_lapse['racs'], fixqco2_lapse['diag']] = pickle.load(open('%s/rad.ssp585_fixq_fixco2_lapse.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # RRTMG fixed T and CO2 RCP8.5 response
    fixtco2 = {}
    [fixtco2['ra'], fixtco2['racs'], fixtco2['diag']] = pickle.load(open('%s/rad.ssp585_fixT_fixco2.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # RRTMG fixed relative humidity and CO2 RCP8.5 response
    fixrhco2 = {}
    [fixrhco2['ra'], fixrhco2['racs'], fixrhco2['diag']] = pickle.load(open('%s/rad.ssp585_fixrh_fixco2.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # RRTMG fixed T and CO2, devRH RCP8.5 response
    fixtco2devrh = {}
    [fixtco2devrh['ra'], fixtco2devrh['racs'], fixtco2devrh['diag']] = pickle.load(open('%s/rad.ssp585_fixT_fixco2_devrh.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # # RRTMG fixed CO2 RCP8.5 response
    # fixco2 = {}
    # [fixco2['ra'], fixco2['racs'], fixco2['diag']] = pickle.load(open('%s/rad.ssp585_fixco2.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # # RRTMG fixed specific humidity RCP8.5 response
    # fixq = {}
    # [fixq['ra'], fixq['racs'], fixq['diag']] = pickle.load(open('%s/rad.ssp585_fixq.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    # # RRTMG fixed relative humidity RCP8.5 response
    # fixrh = {}
    # [fixrh['ra'], fixrh['racs'], fixrh['diag']] = pickle.load(open('%s/rad.ssp585_fixrh.%s.%g.%g.pickle' % (datadir, model,lat0,lat1), 'rb'))

    ##################################################
    # MISC data
    ##################################################
    plev = 1e-2*grid['lev']
    modelstr = get_modelstr(model)

    ##################################################
    # RCP
    ##################################################
    time = np.arange(1987,2201)

    for varname in clim:
        clim[varname] = np.nanmean(clim[varname][-20:,...], axis=0)

    ref['ra'] = np.nanmean(ref['ra'][-20:])
    ref['racs'] = np.nanmean(ref['racs'][-20:])

    ##################################################
    # FIX q and T
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_fixq_fixT.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, full['racs']-ref['racs'], color='black', label='RRTMG')
    ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='RRTMG (fixq+fixT)')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # FIX q and CO2
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_fixq_fixco2.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, full['racs']-ref['racs'], color='black', label='RRTMG')
    ax.plot(time, fixqco2['racs']-ref['racs'], '-', color='tab:orange', label='RRTMG (fixq+fixCO2)')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # DECOMP T into PLANCK and LAPSE RATE
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_dt_dc_planck_lapse.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    # ax.plot(time, full['racs']-ref['racs'], color='black', label='RRTMG')
    ax.plot(time, fixqco2['racs']-ref['racs'], '-', color='tab:orange', label=r'$\Delta T$')
    ax.plot(time, fixqco2_planck['racs']-ref['racs'], '--', color='tab:orange', label=r'$\Delta T$ (Planck)')
    ax.plot(time, fixqco2_lapse['racs']-ref['racs'], ':', color='tab:orange', label=r'$\Delta T$ (Lapse Rate)')
    ax.plot(time, fixqco2['racs']-fixqco2_planck['racs']-fixqco2_lapse['racs']+ref['racs'], '-.', color='tab:orange', label=r'Residual')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    ax.set_title('RRTMG')
    fig.set_size_inches(4,3.3)
    plt.legend(frameon=False, fontsize=7)
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # FIX T and CO2
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_fixT_fixco2.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, full['racs']-ref['racs'], color='black', label='RRTMG')
    ax.plot(time, fixtco2['racs']-ref['racs'], '-', color='tab:blue', label='RRTMG (fixT+fixCO2)')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # FIX RH and CO2
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_fixrh_fixco2.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, full['racs']-ref['racs'], color='black', label='RRTMG')
    ax.plot(time, fixrhco2['racs']-ref['racs'], '--', color='tab:blue', label='RRTMG (fixRH+fixCO2)')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # DECOMP CO2, T, q (no legend)
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_decomp_noleg.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    lres,=ax.plot(time, full['racs']-fixqt['racs']-fixqco2['racs']-fixtco2['racs']+2*ref['racs'], '-.', color='k', label='Residual')
    lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
    lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='CO$_2$')
    lt,=ax.plot(time, fixqco2['racs']-ref['racs'], '-', color='tab:orange', label='T')
    lq,=ax.plot(time, fixtco2['racs']-ref['racs'], '-', color='tab:blue', label='q')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    if model == 'mmm':
        ax.set_title('RRTMG')
    else:
        ax.set_title('RRTMG(%s)' % (modelstr))
    fig.set_size_inches(4,3.3)
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # DECOMP CO2, T, q (inside legend)
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_decomp.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    lres,=ax.plot(time, full['racs']-fixqt['racs']-fixqco2['racs']-fixtco2['racs']+2*ref['racs'], '-.', color='k', label='Residual')
    lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
    lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='CO$_2$')
    lt,=ax.plot(time, fixqco2['racs']-ref['racs'], '-', color='tab:orange', label='T')
    lq,=ax.plot(time, fixtco2['racs']-ref['racs'], '-', color='tab:blue', label='q')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    if model == 'mmm':
        ax.set_title('RRTMG')
    else:
        ax.set_title('RRTMG(%s)' % (modelstr))
    fig.set_size_inches(4,3.3)
    plt.tight_layout()
    leg=plt.legend(handles=[lfull,lco2,lt,lq,lres], loc='lower left', prop={'size':8})
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_facecolor('none')
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # DECOMP CO2, T, q (w/ legend)
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_decomp_leg.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    lres,=ax.plot(time, full['racs']-fixqt['racs']-fixqco2['racs']-fixtco2['racs']+2*ref['racs'], '-.', color='k', label='Residual')
    lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
    lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='CO$_2$')
    lt,=ax.plot(time, fixqco2['racs']-ref['racs'], '-', color='tab:orange', label='T')
    lq,=ax.plot(time, fixtco2['racs']-ref['racs'], '-', color='tab:blue', label='q')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(6,5)
    plt.tight_layout()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15,
                         box.width, box.height * 0.85])
    ax.legend(handles=[lfull,lco2,lt,lq,lres], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # DECOMP CO2, T, RH (no legend)
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_decomp_rh_noleg.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    lres,=ax.plot(time, full['racs']-fixqt['racs']-fixrhco2['racs']-fixtco2devrh['racs']+2*ref['racs'], '-.', color='k', label='Residual')
    lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
    lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='$\Delta$CO$_2$')
    lcc,=ax.plot(time, fixrhco2['racs']-ref['racs'], '-', color='tab:orange', label=r'$\Delta T$, $\Delta q|_\mathrm{fix\,RH}$')
    ldcc,=ax.plot(time, fixtco2devrh['racs']-ref['racs'], '-', color='tab:blue', label=r'$\Delta RH q^\ast|_\mathrm{clim}$')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    if model == 'mmm':
        ax.set_title('RRTMG')
    else:
        ax.set_title('RRTMG(%s)' % (modelstr))
    fig.set_size_inches(4,3.3)
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # DECOMP CO2, T, RH (inside legend)
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_decomp_rh.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    lres,=ax.plot(time, full['racs']-fixqt['racs']-fixrhco2['racs']-fixtco2devrh['racs']+2*ref['racs'], '-.', color='k', label='Residual')
    lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
    lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='$\Delta$CO$_2$')
    lcc,=ax.plot(time, fixrhco2['racs']-ref['racs'], '-', color='tab:orange', label=r'$\Delta T$, $\Delta q|_\mathrm{fix\,RH}$')
    ldcc,=ax.plot(time, fixtco2devrh['racs']-ref['racs'], '-', color='tab:blue', label=r'$\Delta RH q^\ast|_\mathrm{clim}$')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    if model == 'mmm':
        ax.set_title('RRTMG')
    else:
        ax.set_title('RRTMG(%s)' % (modelstr))
    fig.set_size_inches(4,3.3)
    plt.tight_layout()
    plt.tight_layout()
    leg=plt.legend(handles=[lfull,lco2,lcc,ldcc,lres], prop={'size':8})
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_facecolor('none')
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # DECOMP CO2, T, RH (w/ legend)
    ##################################################
    # compare change in radiative cooling
    plotname='%s/dracs_decomp_rh_leg.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    lres,=ax.plot(time, full['racs']-fixqt['racs']-fixrhco2['racs']-fixtco2devrh['racs']+2*ref['racs'], '-.', color='k', label='Residual')
    lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
    lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='$\Delta$CO$_2$')
    lcc,=ax.plot(time, fixrhco2['racs']-ref['racs'], '-', color='tab:orange', label=r'$\Delta T$, $\Delta q|_\mathrm{fix\,RH}$')
    ldcc,=ax.plot(time, fixtco2devrh['racs']-ref['racs'], '-', color='tab:blue', label=r'$\Delta RH q^\ast|_\mathrm{clim}$')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(5.5,4.5)
    plt.tight_layout()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15,
                         box.width, box.height * 0.85])
    ax.legend(handles=[lfull,lco2,lcc,ldcc,lres], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ###################################################
    ## FIX q and T (compared to CO2)
    ###################################################
    #vmin = -111; vmax=-108;
    #ra_f30 = np.mean(fixqt['racs'][:30])
    #ra_l30 = np.mean(fixqt['racs'][-30:])
    #co2_f30 = np.mean(co2[:30])
    #co2_l30 = np.mean(co2[-30:])
    #maxis = (co2_l30 - co2_f30) / (ra_l30 - ra_f30)
    #vmax_alg = co2_f30 + maxis*(vmax-ra_f30)
    #vmin_alg = co2_l30 + maxis*(vmin-ra_l30)
    ## compare absolute change in radiative cooling
    #plotname='%s/co2_racs_fixq_fixT.pdf' % (plotdir)
    #fig, ax = plt.subplots()
    #ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    #ax.plot(time, fixqt['racs'], '-', color='tab:green', label='RRTMG (fixq+T)')
    #ax.set_xlim([time[0], time[-1]])
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())
    #ax.set_xlabel('Time (yr)')
    #ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    #ax.set_ylim(vmin, vmax)
    #sax = ax.twinx()
    #sax.plot(time, co2, '-', color='tab:purple', label='p_{CO_2}')
    #sax.set_ylim(vmin_alg, vmax_alg)
    #fig.set_size_inches(4,3.3)
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(plotname, format='pdf', dpi=300)
    #plt.close()

    ###################################################
    ## FIX Q
    ###################################################
    ## compare absolute change in radiative cooling
    #plotname='%s/racs_fixq.pdf' % (plotdir)
    #fig, ax = plt.subplots()
    #ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    #ax.plot(time, forc['racs'], color='tab:red', label=modelstr)
    #ax.plot(time, full['racs'], color='black', label='RRTMG')
    #ax.plot(time, fixq['racs'], color='tab:orange', label='RRTMG (fixq)')
    #ax.set_xlim([time[0], time[-1]])
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())
    #ax.set_xlabel('Time (yr)')
    #ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    #fig.set_size_inches(4,3.3)
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(plotname, format='pdf', dpi=300)
    #plt.close()

    ###################################################
    ## FIX RH
    ###################################################
    ## compare absolute change in radiative cooling
    #plotname='%s/racs_fixrh.pdf' % (plotdir)
    #fig, ax = plt.subplots()
    #ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    #ax.plot(time, forc['racs'], color='tab:red', label=modelstr)
    #ax.plot(time, full['racs'], color='black', label='RRTMG')
    #ax.plot(time, fixrh['racs'], '--', color='tab:blue', label='RRTMG (fixrh)')
    #ax.set_xlim([time[0], time[-1]])
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())
    #ax.set_xlabel('Time (yr)')
    #ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    #fig.set_size_inches(4,3.3)
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(plotname, format='pdf', dpi=300)
    #plt.close()

    ###################################################
    ## FIX RH and CO2
    ###################################################
    ## compare absolute change in radiative cooling
    #plotname='%s/racs_fixrh_fixco2.pdf' % (plotdir)
    #fig, ax = plt.subplots()
    #ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    #ax.plot(time, forc['racs'], color='tab:red', label=modelstr)
    #ax.plot(time, full['racs'], color='black', label='RRTMG')
    #ax.plot(time, fixrhco2['racs'], '--', color='tab:blue', label='RRTMG (fixrh+co2)')
    #ax.set_xlim([time[0], time[-1]])
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())
    #ax.set_xlabel('Time (yr)')
    #ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    #fig.set_size_inches(4,3.3)
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(plotname, format='pdf', dpi=300)
    #plt.close()

    ###################################################
    ## FIX CO2
    ###################################################
    ## compare absolute change in radiative cooling
    #plotname='%s/racs_fixco2.pdf' % (plotdir)
    #fig, ax = plt.subplots()
    #ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    #ax.plot(time, forc['racs'], color='tab:red', label=modelstr)
    #ax.plot(time, full['racs'], color='black', label='RRTMG')
    #ax.plot(time, fixco2['racs'], color='tab:blue', label='RRTMG (fixco2)')
    #ax.set_xlim([time[0], time[-1]])
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())
    #ax.set_xlabel('Time (yr)')
    #ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    #fig.set_size_inches(4,3.3)
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(plotname, format='pdf', dpi=300)
    #plt.close()
