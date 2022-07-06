import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import *
from misc.filenames import *
from misc.means import *
from misc.load_data import *
import misc.par as par
from proc.circ import save_circ
from plot.titles import *
import os
import pickle
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def circ_lat_lev(sim, **kwargs):

    categ = 'lat_lev'

    linthresh=1e-2

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    refclim = kwargs.get('refclim', '') 
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    seas = kwargs.get('seas', 'djf') # view plot? (plt.show)

    if sim == 'longrun':
        model = kwargs.get('model', 'MPIESM12_abrupt4x')
        yr_span = kwargs.get('yr_span', '1000')
        yr_base = 0
    elif sim == 'historical':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '186001-200512')
        yr_base = 1860
    elif sim == 'rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '200601-230012')
        yr_base = 2006
    elif sim == 'echam':
        model = kwargs.get('model', 'rp000140')
        yr_span = kwargs.get('yr_span', '0001_0039')
        yr_base = 0
    elif sim == 'era5':
        model = None
        yr_span = kwargs.get('yr_span', '1980_2005')
        yr_base = 1980

    vmin = -90
    vmax = 90

    ##################################
    # LOAD DATA
    ##################################
    if refclim == 'hist-30':
        sim_ref = 'historical'
        timemean_ref = 'ymonmean-30'
        yr_span_ref = '186001-200512'

    if isinstance(model, str) or model == None:
        [circ, grid, datadir, plotdir, modelstr] = load_circ(sim, categ, zonmean=zonmean, timemean=timemean, domain=domain, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
        [dyn, _, _, _, _] = load_dyn(sim, categ, zonmean=zonmean, timemean='djfmean', try_load=try_load, model=model, yr_span=yr_span)

        [circ_ref, grid_ref, _, _, _] = load_circ(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, domain=domain, try_load=try_load, model=model, yr_span=yr_span_ref, refclim=refclim) 
        [dyn_ref, _, _, _, _] = load_dyn(sim_ref, categ, zonmean=zonmean, timemean='djfmean', try_load=try_load, model=model, yr_span=yr_span_ref)
    else:
        [circ, grid, datadir, plotdir, modelstr, circ_mmm] = load_circ(sim, categ, zonmean=zonmean, timemean=timemean, domain=domain, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
        [dyn, _, _, _, _, _] = load_dyn(sim, categ, zonmean=zonmean, timemean='djfmean', try_load=try_load, model=model, yr_span=yr_span)

        [circ_ref, grid_ref, _, _, _, _] = load_circ(sim_ref, categ, zonmean=zonmean, timemean=timemean_ref, domain=domain, try_load=try_load, model=model, yr_span=yr_span_ref, refclim=refclim) 
        [dyn_ref, _, _, _, _, _] = load_dyn(sim_ref, categ, zonmean=zonmean, timemean='djfmean', try_load=try_load, model=model, yr_span=yr_span_ref)

    # print(np.reshape(r1, (-1,96,12)).shape)
    # if timemean == '':
    #     circ = np.mean(np.reshape(circ, (-1,12,circ.shape[1])),1)

    if seas == 'djf':
        for circname in circ:
            circ[circname] = np.mean(np.roll(circ[circname], 1, axis=0)[0:3,...], axis=0)
            circ_ref[circname] = np.mean(np.roll(circ_ref[circname], 1, axis=0)[0:3,...], axis=0)
        # for dynname in dyn:
        #     dyn[dynname] = np.mean(np.roll(dyn[dynname], 1, axis=0)[0:3,...], axis=0)
        #     dyn_ref[dynname] = np.mean(np.roll(dyn_ref[dynname], 1, axis=0)[0:3,...], axis=0)
    elif seas == 'mam':
        for circname in circ:
            circ[circname] = np.mean(circ[circname][2:5,...], axis=0)
            circ_ref[circname] = np.mean(circ_ref[circname][2:5,...], axis=0)
        # for dynname in dyn:
        #     dyn[dynname] = np.mean(dyn[dynname][2:5,...], axis=0)
        #     dyn_ref[dynname] = np.mean(dyn_ref[dynname][2:5,...], axis=0)
    elif seas == 'jja':
        for circname in circ:
            circ[circname] = np.mean(circ[circname][5:8,...], axis=0)
            circ_ref[circname] = np.mean(circ_ref[circname][5:8,...], axis=0)
        # for dynname in dyn:
        #     dyn[dynname] = np.mean(dyn[dynname][5:8,...], axis=0)
        #     dyn_ref[dynname] = np.mean(dyn_ref[dynname][5:8,...], axis=0)
    elif seas == 'son':
        for circname in circ:
            circ[circname] = np.mean(circ[circname][8:11,...], axis=0)
            circ_ref[circname] = np.mean(circ_ref[circname][8:11,...], axis=0)
        # for dynname in dyn:
        #     dyn[dynname] = np.mean(dyn[dynname][8:11,...], axis=0)
        #     dyn_ref[dynname] = np.mean(dyn_ref[dynname][8:11,...], axis=0)
    elif seas == 'ann':
        for circname in circ:
            circ[circname] = np.mean(circ[circname], axis=0)
            circ_ref[circname] = np.mean(circ_ref[circname], axis=0)
        # for dynname in dyn:
        #     dyn[dynname] = np.mean(dyn[dynname], axis=0)
        #     dyn_ref[dynname] = np.mean(dyn_ref[dynname], axis=0)

    for dynname in dyn:
        dyn[dynname] = np.mean(dyn[dynname][:-30,...], axis=0)
        dyn_ref[dynname] = np.mean(dyn_ref[dynname][:-30,...], axis=0)

    dcirc = {}
    for circname in circ:
        dcirc[circname] = circ[circname] - circ_ref[circname]
    ddyn = {}
    for dynname in dyn:
        ddyn[dynname] = dyn[dynname] - dyn_ref[dynname]
        
    if 'ymonmean' in timemean:
        yr_base=0
    [mesh_lev, mesh_lat] = np.meshgrid(grid['lev'], grid['lat']) # create mesh

    ###################################
    ## PLOT VERTICAL VELOCITY
    ###################################
    #plotname = remove_repdots('%s/wap_lat_lev.%s.%s' % (plotdir, seas, timemean))
    #fig, ax = plt.subplots(2, sharex=True)
    #vmin = -50
    #vmax = 50
    #vint = 5
    ## csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    #ax[0].axhline(0, color='k', linewidth=0.5)
    #mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    #csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(circ['wap']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    #ax[1].contour(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(circ['wap']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    #make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    #ax[0].set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    #ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    #ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    #ax[1].set_xticks(np.arange(-90,91,30))
    #ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    #ax[1].set_xlabel('Latitude (deg)')
    #ax[1].set_ylim([0, 1000])
    #ax[1].set_ylim(ax[1].get_ylim()[::-1])
    #ax[1].set_yticks(np.arange(0,1001,200))
    #ax[1].set_ylabel('Pressure (hPa)')
    #ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    #fig.set_size_inches(5, 3)
    #plt.tight_layout()
    #cbar = plt.colorbar(csf, ax=ax)
    #cbar.set_label('$\omega$ (hPa d$^{-1}$)')
    #pos = ax[1].get_position()
    #pos0 = ax[0].get_position()
    #ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    #ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    #plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    #if viewplt:
    #    plt.show()
    #plt.close()

    ###################################
    ## PLOT VERTICAL VELOCITY (45-90)
    ###################################
    #plotname = remove_repdots('%s/wap_lat_lev.np.%s.%s' % (plotdir, seas, timemean))
    #fig, ax = plt.subplots(2, sharex=True)
    #vmin = -10
    #vmax = 10
    #vint =1
    #ax[0].axhline(0, color='k', linewidth=0.5)
    #ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    #ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    #mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    #csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(circ['wap']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    #ax[1].contour(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(circ['wap']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    #make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    #ax[0].set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    #ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    #ax[0].set_yscale('symlog',linthresh=linthresh)
    #ax[0].tick_params(axis='y', labelsize=5)
    #ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    #ax[1].set_xlim([40,90])
    #ax[1].set_xticks(np.arange(40,91,10))
    #ax[1].xaxis.set_minor_locator(MultipleLocator(5))
    #ax[1].set_xlabel('Latitude (deg)')
    #ax[1].set_ylim([0, 1000])
    #ax[1].set_ylim(ax[1].get_ylim()[::-1])
    #ax[1].set_yticks(np.arange(0,1001,200))
    #ax[1].set_ylabel('Pressure (hPa)')
    #ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    #fig.set_size_inches(4.5, 3)
    #plt.tight_layout()
    #cbar = plt.colorbar(csf,ax=ax)
    #cbar.set_label('$\omega$ (hPa d$^{-1}$)')
    #pos = ax[1].get_position()
    #pos0 = ax[0].get_position()
    #ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    #ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    #plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    #if viewplt:
    #    plt.show()
    #plt.close()

    ###################################
    ## PLOT CHANGE IN VERTICAL VELOCITY
    ###################################
    #plotname = remove_repdots('%s/dwap_lat_lev.%s.%s' % (plotdir, seas, timemean))
    #fig, ax = plt.subplots(2, sharex=True)
    #vmin_ref = -50
    #vmax_ref = 50
    #vint_ref = 5
    #vmin = -10
    #vmax = 10
    #vint = 1
    ## csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    #ax[0].axhline(0, color='k', linewidth=0.5)
    #mmc = ax[0].plot(grid['lat'], 1e-15*ddyn['vmmmc'], 'r')
    #ax00 = ax[0].twinx()
    #ax00.plot(grid['lat'], 1e-15*dyn_ref['vmmmc'], 'r', linewidth=0.5)
    #csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(dcirc['wap']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    #ax[1].contour(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(dcirc['wap']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    #cs_ref = ax[1].contour(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(circ_ref['wap']), np.arange(vmin_ref, vmax_ref+1, vint_ref), colors='k', vmin=vmin_ref, vmax=vmax_ref, linewidths=0.5)
    #make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    #ax00.set_ylim([-3.3,3.3])
    #ax00.yaxis.set_minor_locator(AutoMinorLocator())
    #ax[0].set_ylim([-1,1])
    #ax[0].set_ylabel('$\Delta [\overline{v}]\,[\overline{m}]$ (PW)')
    #ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ## ax[0].tick_params(axis='y', labelsize=5)
    #ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    #ax[1].set_xticks(np.arange(-90,91,30))
    #ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    #ax[1].set_xlabel('Latitude (deg)')
    #ax[1].set_ylim([0, 1000])
    #ax[1].set_ylim(ax[1].get_ylim()[::-1])
    #ax[1].set_yticks(np.arange(0,1001,200))
    #ax[1].set_ylabel('Pressure (hPa)')
    #ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    #fig.set_size_inches(5, 3)
    #plt.tight_layout()
    #cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    #cbar.set_label('$\Delta \omega$ (hPa d$^{-1}$)')
    #pos = ax[1].get_position()
    #pos0 = ax[0].get_position()
    #ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    #ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    #plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    #if viewplt:
    #    plt.show()
    #plt.close()

    ###################################
    ## PLOT VERTICAL VELOCITY (45-90)
    ###################################
    #plotname = remove_repdots('%s/dwap_lat_lev.np.%s.%s' % (plotdir, seas, timemean))
    #fig, ax = plt.subplots(2, sharex=True)
    #vmin_ref = -10
    #vmax_ref = 10
    #vint_ref = 1
    #vmin = -5
    #vmax = 5
    #vint =0.5
    ## csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    #ax[0].axhline(0, color='k', linewidth=0.5)
    #ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    #ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    #mmc = ax[0].plot(grid['lat'], 1e-15*ddyn['vmmmc'], 'r')
    #ax00 = ax[0].twinx()
    #ax00.plot(grid['lat'], 1e-15*dyn_ref['vmmmc'], 'r', linewidth=0.5)
    #csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(dcirc['wap']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    #ax[1].contour(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(dcirc['wap']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    #cs_ref = ax[1].contour(mesh_lat, 1e-2*mesh_lev, 864*np.transpose(circ_ref['wap']), np.arange(vmin_ref, vmax_ref+1, vint_ref), colors='k', vmin=vmin_ref, vmax=vmax_ref, linewidths=0.5)
    #make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    #ax00.set_ylim([-3.3,3.3])
    #ax00.yaxis.set_minor_locator(AutoMinorLocator())
    #ax00.set_yscale('symlog',linthresh=linthresh)
    #ax00.tick_params(axis='y', labelsize=5)
    #ax[0].set_ylim([-1,1])
    #ax[0].set_ylabel('$\Delta [\overline{v}]\,[\overline{m}]$ (PW)')
    #ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    #ax[0].set_yscale('symlog',linthresh=linthresh)
    #ax[0].tick_params(axis='y', labelsize=5)
    ## ax[0].tick_params(axis='y', labelsize=5)
    #ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    #ax[1].set_xlim([40,90])
    #ax[1].set_xticks(np.arange(40,91,10))
    #ax[1].xaxis.set_minor_locator(MultipleLocator(5))
    #ax[1].set_xlabel('Latitude (deg)')
    #ax[1].set_ylim([0, 1000])
    #ax[1].set_ylim(ax[1].get_ylim()[::-1])
    #ax[1].set_yticks(np.arange(0,1001,200))
    #ax[1].set_ylabel('Pressure (hPa)')
    #ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    #fig.set_size_inches(5, 3)
    #plt.tight_layout()
    #cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    #cbar.set_label('$\Delta \omega$ (hPa d$^{-1}$)')
    #pos = ax[1].get_position()
    #pos0 = ax[0].get_position()
    #ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    #ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    #plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    #if viewplt:
    #    plt.show()
    #plt.close()

    ##################################
    # PLOT MMC Streamfunction
    ##################################
    plotname = remove_repdots('%s/psi_lat_lev.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin = -10
    vmax = 10
    vint = 1
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ['psi']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ['psi']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].set_xticks(np.arange(-90,91,30))
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    ax[1].set_xlim([-90, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax)
    cbar.set_label('$\psi$ ($10^{10}$ kg s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT MMC Streamfunction (40 to 90)
    ##################################
    plotname = remove_repdots('%s/psi_lat_lev.np.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin = -1
    vmax = 1
    vint = 0.1
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ['psi']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ['psi']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].set_yscale('symlog',linthresh=linthresh)
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(40,91,10))
    ax[1].xaxis.set_minor_locator(MultipleLocator(5))
    ax[1].set_xlim([40, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax)
    cbar.set_label('$\psi$ ($10^{10}$ kg s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT MMC Streamfunction (mass corrected)
    ##################################
    plotname = remove_repdots('%s/psic_lat_lev.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin = -10
    vmax = 10
    vint = 1
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ['psic']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ['psic']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].set_xticks(np.arange(-90,91,30))
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    ax[1].set_xlim([-90, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax)
    cbar.set_label('$\psi_d$ ($10^{10}$ kg s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT MMC Streamfunction (40 to 90, mass corrected)
    ##################################
    plotname = remove_repdots('%s/psic_lat_lev.np.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin = -1
    vmax = 1
    vint = 0.1
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ['psic']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ['psic']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].set_yscale('symlog',linthresh=linthresh)
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(40,91,10))
    ax[1].xaxis.set_minor_locator(MultipleLocator(5))
    ax[1].set_xlim([40, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax)
    cbar.set_label('$\psi_d$ ($10^{10}$ kg s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT DELTA MMC Streamfunction
    ##################################
    plotname = remove_repdots('%s/dpsi_lat_lev.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -10
    vmax_ref = 10
    vint_ref = 1
    vmin = -10
    vmax = 10
    vint = 1
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*ddyn['vmmmc'], 'r')
    ax00 = ax[0].twinx()
    ax00.plot(grid['lat'], 1e-15*dyn_ref['vmmmc'], 'r', linewidth=0.5)
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(dcirc['psi']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(dcirc['psi']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psi']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax00.set_ylim([-3.3,3.3])
    ax00.set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].set_ylim([-1,1])
    ax[0].set_ylabel('$\Delta[\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].set_xticks(np.arange(-90,91,30))
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    ax[1].set_xlim([-90, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label('$\Delta\psi$ ($10^{10}$ kg s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT DELTA MMC Streamfunction (40 to 90)
    ##################################
    plotname = remove_repdots('%s/dpsi_lat_lev.np.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -1
    vmax_ref = 1
    vint_ref = 0.1
    vmin = -1
    vmax = 1
    vint = 0.1
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*ddyn['vmmmc'], 'r')
    ax00 = ax[0].twinx()
    ax00.plot(grid['lat'], 1e-15*dyn_ref['vmmmc'], 'r', linewidth=0.5)
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(dcirc['psi']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(dcirc['psi']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psi']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax00.set_ylim([-3.3,3.3])
    ax00.set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    ax00.set_yscale('symlog',linthresh=linthresh)
    ax00.tick_params(axis='y', labelsize=5)
    ax[0].set_ylim([-3.3,3.3])
    ax[0].set_ylabel('$\Delta[\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].set_yscale('symlog',linthresh=linthresh)
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(40,91,10))
    ax[1].xaxis.set_minor_locator(MultipleLocator(5))
    ax[1].set_xlim([40, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label('$\Delta\psi$ ($10^{10}$ kg s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT MSE-streamfunction overlay (-90 to 90)
    ##################################
    plotname = remove_repdots('%s/mse_lat_lev.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -10
    vmax_ref = 10
    vint_ref = 1
    vmin = 250
    vmax = 400
    vint = 5
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-3*np.transpose(circ['mse']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, extend='both')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-3*np.transpose(circ['mse']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psic']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(-90,91,30))
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    ax[1].set_xlim([-90, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label('MSE (J g$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT MSE-streamfunction overlay (40 to 90)
    ##################################
    plotname = remove_repdots('%s/mse_lat_lev.np.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -1
    vmax_ref = 1
    vint_ref = 0.1
    vmin = 250
    vmax = 400
    vint = 5 
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*ddyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-3*np.transpose(circ['mse']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, extend='both')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-3*np.transpose(circ['mse']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psic']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel('$\Delta [\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].set_yscale('symlog',linthresh=linthresh)
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(40,91,10))
    ax[1].xaxis.set_minor_locator(MultipleLocator(5))
    ax[1].set_xlim([40, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label('MSE (J g$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT DELTA MSE-streamfunction overlay (-90 to 90)
    ##################################
    plotname = remove_repdots('%s/dmse_lat_lev.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -10
    vmax_ref = 10
    vint_ref = 1
    vmin = -30
    vmax = 30
    vint = 5
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*ddyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-3*np.transpose(dcirc['mse']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, extend='both', cmap='RdBu_r')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-3*np.transpose(dcirc['mse']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psic']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel('$\Delta [\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(-90,91,30))
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    ax[1].set_xlim([-90, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label('$\Delta$ MSE (J g$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT DELTA MSE-streamfunction overlay (40 to 90)
    ##################################
    plotname = remove_repdots('%s/dmse_lat_lev.np.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -1
    vmax_ref = 1
    vint_ref = 0.1
    vmin = -30
    vmax = 30
    vint = 5 
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-3*np.transpose(dcirc['mse']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, extend='both', cmap='RdBu_r')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-3*np.transpose(dcirc['mse']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psic']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel('$[\overline{v}]\,[\overline{m}]$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].set_yscale('symlog',linthresh=linthresh)
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(40,91,10))
    ax[1].xaxis.set_minor_locator(MultipleLocator(5))
    ax[1].set_xlim([40, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label('$\Delta$ MSE (J g$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT vvmmmc-streamfunction overlay (-90 to 90)
    ##################################
    plotname = remove_repdots('%s/vvmmmc_lat_lev.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -10
    vmax_ref = 10
    vint_ref = 1
    vmin = -5
    vmax = 5
    vint = 0.5
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-4*np.transpose(circ['vvmmmc']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    # ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-4*np.transpose(circ['vvmmmc']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psic']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel(r'$\langle[\overline{v}]\,[\overline{m}]\rangle$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(-90,91,30))
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    ax[1].set_xlim([-90, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label(r'$\widetilde{[\overline{v}]}\,\widetilde{[\overline{m}]}$ ($10^4$ J kg$^{-1}$ m s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT vvmmmc-streamfunction overlay (40 to 90)
    ##################################
    plotname = remove_repdots('%s/vvmmmc_lat_lev.np.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -1
    vmax_ref = 1
    vint_ref = 0.1
    vmin = -1
    vmax = 1
    vint = 0.1
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*dyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-4*np.transpose(circ['vvmmmc']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    # ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-4*np.transpose(circ['vvmmmc']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psic']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel(r'$\langle[\overline{v}]\,[\overline{m}]\rangle$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].set_yscale('symlog',linthresh=linthresh)
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(40,91,10))
    ax[1].xaxis.set_minor_locator(MultipleLocator(5))
    ax[1].set_xlim([40, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label(r'$\widetilde{[\overline{v}]}\,\widetilde{[\overline{m}]}$ ($10^4$ J kg$^{-1}$ m s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT DELTA vvmmmc-streamfunction overlay (-90 to 90)
    ##################################
    plotname = remove_repdots('%s/dvvmmmc_lat_lev.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -10
    vmax_ref = 10
    vint_ref = 1
    vmin = -5
    vmax = 5
    vint = 0.5
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*ddyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-4*np.transpose(dcirc['vvmmmc']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    # ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-4*np.transpose(circ['vvmmmc']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psic']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel(r'$\Delta \langle[\overline{v}]\,[\overline{m}]\rangle$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(-90,91,30))
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    ax[1].set_xlim([-90, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label(r'$\Delta \widetilde{[\overline{v}]}\,\widetilde{[\overline{m}]}$ ($10^4$ J kg$^{-1}$ m s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ##################################
    # PLOT DELTA vvmmmc-streamfunction overlay (40 to 90)
    ##################################
    plotname = remove_repdots('%s/dvvmmmc_lat_lev.np.%s.%s' % (plotdir, seas, timemean))
    fig, ax = plt.subplots(2, sharex=True)
    vmin_ref = -1
    vmax_ref = 1
    vint_ref = 0.1
    vmin = -1
    vmax = 1
    vint = 0.1
    # csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax[0].axhline(0, color='k', linewidth=0.5)
    ax[0].axhline(linthresh, linestyle=':', color='k', linewidth=0.5)
    ax[0].axhline(-linthresh, linestyle=':', color='k', linewidth=0.5)
    mmc = ax[0].plot(grid['lat'], 1e-15*ddyn['vmmmc'], 'r')
    csf = ax[1].contourf(mesh_lat, 1e-2*mesh_lev, 1e-4*np.transpose(dcirc['vvmmmc']), np.arange(vmin, vmax+1, vint), vmin=vmin, vmax=vmax, cmap='RdBu_r', extend='both')
    # ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-4*np.transpose(circ['vvmmmc']), np.arange(10*vmin, 10*vmax+1, vint), linewidths=0.01, colors='white')
    ax[1].contour(mesh_lat, 1e-2*mesh_lev, 1e-10*np.transpose(circ_ref['psic']), np.arange(10*vmin_ref, 10*vmax_ref+1, vint_ref), colors='k', linewidths=0.5)
    make_title_sim_time_seas(ax[0], sim, model=modelstr, timemean=timemean, seasmean=seas)
    ax[0].set_ylabel(r'$\Delta \langle[\overline{v}]\,[\overline{m}]\rangle$ (PW)')
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[0].set_yscale('symlog',linthresh=linthresh)
    ax[0].tick_params(axis='y', labelsize=5)
    ax[1].set_xticks(np.arange(40,91,10))
    ax[1].xaxis.set_minor_locator(MultipleLocator(5))
    ax[1].set_xlim([40, 90])
    ax[1].set_xlabel('Latitude (deg)')
    ax[1].set_yscale('log')
    ax[1].set_ylim([10, 1000])
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    # ax[1].set_yticks(np.arange(10,1001,200))
    ax[1].set_ylabel('Pressure (hPa)')
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', bottom=True, top=True, left=True, right=True)
    fig.set_size_inches(5, 3)
    plt.tight_layout()
    cbar = plt.colorbar(csf, ax=ax, pad=0.15)
    cbar.set_label(r'$\Delta \widetilde{[\overline{v}]}\,\widetilde{[\overline{m}]}$ ($10^4$ J kg$^{-1}$ m s$^{-1}$)')
    pos = ax[1].get_position()
    pos0 = ax[0].get_position()
    ax[0].set_position([pos.x0,pos0.y0+pos0.height/2,pos.width,pos0.height/2])
    ax[1].set_position([pos.x0,pos.y0,pos.width,pos.height+pos0.height/2])
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

