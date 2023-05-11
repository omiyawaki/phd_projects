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
    # LOAD ECHAM DATA
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

    print(ref['ra'].shape)
    
    ##################################################
    # MISC data
    ##################################################
    plev = 1e-2*grid['lev']
    modelstr = get_modelstr(model)

    ##################################################
    # HISTORICAL
    ##################################################
    time0 = np.arange(20,41)

    ##################################################
    # TOTAL RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/ra_348.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time0, clim['ra'], color='tab:gray', label=modelstr)
    ax.plot(time0, ref['ra'], color='black', label='RRTMG')
    ax.set_xlim([time0[0], time0[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_a$ (W m$^{-2}$)')
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # CLEAR SKY RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/racs_348.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time0, clim['racs'], color='tab:red', label=modelstr)
    ax.plot(time0, ref['racs'], color='black', label='RRTMG')
    ax.set_xlim([time0[0], time0[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # RCP
    ##################################################
    time = np.arange(1987,2201)

    for varname in clim:
        clim[varname] = np.nanmean(clim[varname][-20:,...], axis=0)

    ref['ra'] = np.nanmean(ref['ra'][-20:])
    ref['racs'] = np.nanmean(ref['racs'][-20:])

    ##################################################
    # ALL COMP DRA
    ##################################################
    # compare relative change in radiative cooling
    plotname='%s/dra_comp.pdf' % (plotdir)

    sw_forc = forc['srad0'] - forc['srads']
    lw_forc = forc['trad0'] - forc['trads']
    lwcs_forc = forc['traf0'] - forc['trafs']
    lwcld_forc = lw_forc - lwcs_forc
    sw_clim = clim['srad0'] - clim['srads']
    lw_clim = clim['trad0'] - clim['trads']
    lwcs_clim = clim['traf0'] - clim['trafs']
    lwcld_clim = lw_clim - lwcs_clim
    dsw = sw_forc - sw_clim
    dlw = lw_forc - lw_clim
    dlwcs = lwcs_forc - lwcs_clim
    dlwcld = lwcld_forc - lwcld_clim

    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)

    ax.plot(time, forc['ra']-clim['ra'], color='tab:gray', label='$\Delta R_a$')
    ax.plot(time, dsw, color='cyan', label='$\Delta SW$')
    ax.plot(time, dlwcs, '-', color='tab:red', label='$\Delta LW_{clear}$')
    ax.plot(time, dlwcld, '-', color='tab:purple', label='$\Delta LW_{cloud}$')
    ax.plot(time, full['racs']-ref['racs'], color='black', label='$\Delta$RRTMG$_{clear}$')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta$ Energy flux (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    leg=plt.legend()
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_facecolor('none')
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    # no leg
    plotname='%s/dra_comp_noleg.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, forc['ra']-clim['ra'], color='tab:gray', label='$\Delta R_a$')
    ax.plot(time, dsw, color='cyan', label='$\Delta SW$')
    ax.plot(time, dlwcs, '-', color='tab:red', label='$\Delta LW_{clear}$')
    ax.plot(time, dlwcld, '-', color='tab:purple', label='$\Delta LW_{cloud}$')
    ax.plot(time, full['racs']-ref['racs'], color='black', label='$\Delta$RRTMG$_{clear}$')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta$ Energy flux (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    # no leg
    plotname='%s/dra_comp_noleg_sw.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, forc['ra']-clim['ra'], color='tab:gray', label='$\Delta R_a$')
    # ax.plot(time, dsw, color='cyan', label='$\Delta SW$')
    # ax.plot(time, dlwcs, '-', color='tab:red', label='$\Delta LW_{clear}$')
    # ax.plot(time, dlwcld, '-', color='tab:purple', label='$\Delta LW_{cloud}$')
    # ax.plot(time, full['racs']-ref['racs'], color='black', label='$\Delta$RRTMG$_{clear}$')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta$ Energy flux (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()


    ##################################################
    # TOTAL RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/ra.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, forc['ra'], color='tab:gray', label=modelstr)
    ax.plot(time, full['ra'], color='black', label='RRTMG')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_a$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    # compare relative (to 348orical) change in radiative cooling
    plotname='%s/dra.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, forc['ra'] - clim['ra'], color='tab:gray', label=modelstr)
    ax.plot(time, full['ra'] - ref['ra'], color='black', label='RRTMG')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_a$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # CLEAR SKY RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/racs.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, forc['racs'], color='tab:red', label=modelstr)
    ax.plot(time, full['racs'], color='black', label='RRTMG')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()


    # compare relative (to 348orical) change in clear sky radiative cooling
    plotname='%s/dracs.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, forc['racs'] - clim['racs'], color='tab:red', label=modelstr)
    # ax.plot(time, full['racs'] - clim['racs'], color='black', label='RRTMG')
    ax.plot(time, full['racs'] - ref['racs'], color='black', label='RRTMG')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

