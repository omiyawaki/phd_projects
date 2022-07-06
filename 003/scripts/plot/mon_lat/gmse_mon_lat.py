import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import *
from misc.filenames import *
# from misc.translate import translate_varname
from misc.means import lat_mean, global_int
from misc.load_data import *
from misc import par
from proc.r1 import save_r1
from plot.titles import *
import os
import pickle
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from cartopy.util import add_cyclic_point
# from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                        # LatitudeLocator, LongitudeLocator)

def gmse_mon_lat(sim, **kwargs):

    categ = 'mon_lat'

    zonmean = kwargs.get('zonmean', 0) # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    latbnd = kwargs.get('latbnd', (10,30))
    lonbnd = kwargs.get('lonbnd', (70,90))
    latstep = kwargs.get('latstep', 0.25) # latitude step size used for interpolation
    lonstep = kwargs.get('lonstep', 0.25) # lonitude step size used for interpolonion
    plotover = kwargs.get('plotover', None) # plot overlay (sic for sea ice, ga_dev for lapse rate deviation, pr for precip, decomp for linear decomposition)?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    legend = kwargs.get('legend', 0)
    if plotover == 'ga_dev':
        vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
        # vertbnd = kwargs.get('vertbnd', (0.7, 0.3)) # sigma bounds of vertical integral
        vertbnd = kwargs.get('vertbnd', (1.0, 0.9)) # sigma bounds of vertical integral

    slat = 80

    transname = 'aht'
    translabel = r'$\langle [ \overline{vm} ] \rangle$ (PW)'
    dtranslabel = r'$\Delta \langle [ \overline{vm} ] \rangle$ (PW)'

    # transname = 'vmte'
    # translabel = r'$\langle [ \overline{v^{\prime} m^{\prime}} ] \rangle$ (PW)'
    # dtranslabel = r'$\Delta \langle [ \overline{v^{\prime} m^{\prime}} ] \rangle$ (PW)'

    latlo = 60

    lat_int = np.arange(latbnd[0], latbnd[1], latstep)
    lon_int = np.arange(lonbnd[0], lonbnd[1], lonstep)

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
    elif sim == 'hist+rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '186001-229912')
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
        if 'ymonmean' not in timemean:
            yr_base = int(yr_span[0:4])
        else:
            yr_base = 0

    ##################################
    # LOAD DATA
    ##################################
    prefix = '/project2/tas1/miyawaki/projects/003/data/raw/%s/bcc-csm1-1' % (sim)

    smooth = '_sm'

    file_mse = Dataset('%s/gmse92500%s_Amon_bcc-csm1-1_%s_r1i1p1_186001-229912.zonmean.shsmooth.djfmean.nc' % (prefix, smooth, sim), 'r')
    gmse = np.squeeze(file_mse.variables['gmse%s' % (smooth)][:]) # (mon x lat)
    lato = np.squeeze(file_mse.variables['lat'][:]) # (lat)

    file_vmte = Dataset('%s/%s%s_Amon_bcc-csm1-1_%s_r1i1p1_186001-229912.zonmean.shsmooth.djfmean.nc' % (prefix, transname, smooth, sim), 'r')
    vmte = np.squeeze(file_vmte.variables['%s%s' % (transname, smooth)][:]) # (mon x lat)

    # interpolate to higher latitude grid
    lat = np.arange(-90,91,0.1) # 0.1 deg resolution lat grid
    fgmse = interp1d(lato, gmse, bounds_error=False)
    gmse = fgmse(lat)
    fvmte = interp1d(lato, vmte, bounds_error=False)
    vmte = fvmte(lat)

    ############################################
    # PLOT
    ############################################
    time = yr_base + np.arange(gmse.shape[0]) # create time vector
    plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)
    plotdir = '%s/%s' % (plotdir, transname)
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)

    [mesh_lat, mesh_time] = np.meshgrid(lat, time) # create mesh

    plotname = '%s/gmse92500_time_lat.%s' % (plotdir, timemean)

    fig, ax = plt.subplots()
    csf = plt.contourf(mesh_time, mesh_lat, 1e-3*gmse, np.arange(-30,31,5), cmap='RdBu_r', extend='both')
    plt.contour(mesh_time, mesh_lat, 1e-3*gmse, np.arange(-90,91,5), linewidths=0.5, colors='0.8')
    plt.contour(mesh_time, mesh_lat, 1e-3*gmse, 0, linewidths=1, colors='k')
    plt.title(r'$\partial_{\phi}m_{925\,\mathrm{hPa}}$ [J g$^{-1}$], RCP8.5 DJF')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Latitude (deg)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(2006,2299)
    ax.set_ylim(latlo,90)
    # ax.set_xticks(np.arange(0,361,60))
    ax.set_yticks(np.arange(latlo,91,(90-latlo)/6))
    cbar = plt.colorbar(csf)
    fig.set_size_inches(4, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.close()

    ############################################
    # PLOT DIFFERENCE FROM CLIMATOLOGY
    ############################################
    yrstart = 1975-yr_base
    yrend = 2005-yr_base
    meangmse = np.nanmean(gmse[yrstart:yrend+1,:], axis=0, keepdims=True)
    diffgmse = gmse - meangmse

    plotname = '%s/dgmse92500_time_lat.%s' % (plotdir, timemean)
    fig, ax = plt.subplots()
    csf = plt.contourf(mesh_time, mesh_lat, 1e-3*diffgmse, np.arange(-20,21,2), cmap='RdBu_r', extend='both')
    plt.contour(mesh_time, mesh_lat, 1e-3*gmse, np.arange(-90,91,5), linewidths=0.5, colors='0.8')
    # plt.contour(mesh_time, mesh_lat, 1e-3*gmse, 0, linewidths=1, colors='k')
    plt.title(r'$\Delta \partial_{\phi}m_{925\,\mathrm{hPa}}$ [J g$^{-1}$], RCP8.5 DJF')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Latitude (deg)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(2006,2299)
    ax.set_ylim(latlo,90)
    # ax.set_xticks(np.arange(0,361,60))
    ax.set_yticks(np.arange(latlo,91,(90-latlo)/6))
    cbar = plt.colorbar(csf)
    fig.set_size_inches(4, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.close()

    ############################################
    # COMPUTE DIFFUSIVITY AND COMPONENTS of VMTE CHANGE
    ############################################
    diffv = -vmte / (2*np.pi*gmse)
    meandiffv = np.nanmean(diffv[yrstart:yrend+1,:], axis=0, keepdims=True)
    diffdiffv = diffv - meandiffv

    meanvmte = np.nanmean(vmte[yrstart:yrend+1,:], axis=0, keepdims=True)
    diffvmte = vmte - meanvmte

    # components
    dcgmse = -meandiffv * diffgmse
    dcdiffv = -meangmse * diffdiffv
    dcres = diffvmte - (dcgmse + dcdiffv)

    # Diffusivity
    plotname = '%s/diffv_time_lat.%s' % (plotdir, timemean)
    fig, ax = plt.subplots()
    csf = plt.contourf(mesh_time, mesh_lat, 1e-9*diffv, np.arange(-100,100,10), cmap='RdBu_r', extend='both')
    plt.contour(mesh_time, mesh_lat, 1e-9*diffv, np.arange(-200,200,50), linewidths=0.5, colors='0.8')
    plt.contour(mesh_time, mesh_lat, 1e-9*diffv, 0, linewidths=1, colors='k')
    plt.title(r'$D$ [10$^{9}$ kg s$^{-1}$], RCP8.5 DJF')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Latitude (deg)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(2006,2299)
    ax.set_ylim(latlo,90)
    # ax.set_xticks(np.arange(0,361,60))
    ax.set_yticks(np.arange(latlo,91,(90-latlo)/6))
    cbar = plt.colorbar(csf)
    fig.set_size_inches(4, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.close()

    # TOTAL
    plotname = '%s/d%s_time_lat.%s' % (plotdir, transname, timemean)
    fig, ax = plt.subplots()
    csf = plt.contourf(mesh_time, mesh_lat, 1e-15*diffvmte, np.arange(-0.25,0.25,0.01), cmap='RdBu_r', extend='both')
    plt.contour(mesh_time, mesh_lat, 1e-15*diffvmte, np.arange(-2,2,0.05), linewidths=0.5, colors='0.8')
    plt.contour(mesh_time, mesh_lat, 1e-15*diffvmte, 0, linewidths=1, colors='k')
    plt.title('%s, RCP8.5 DJF' % (translabel) )
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Latitude (deg)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(2006,2299)
    ax.set_ylim(latlo,90)
    # ax.set_xticks(np.arange(0,361,60))
    ax.set_yticks(np.arange(latlo,91,(90-latlo)/6))
    cbar = plt.colorbar(csf)
    fig.set_size_inches(4, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.close()

    # MSE GRADIENT COMPONENT
    plotname = '%s/dcgmse_time_lat.%s' % (plotdir, timemean)
    fig, ax = plt.subplots()
    csf = plt.contourf(mesh_time, mesh_lat, 1e-15*dcgmse, np.arange(-0.25,0.25,0.01), cmap='RdBu_r', extend='both')
    plt.contour(mesh_time, mesh_lat, 1e-15*dcgmse, np.arange(-2,2,0.05), linewidths=0.5, colors='0.8')
    plt.contour(mesh_time, mesh_lat, 1e-15*dcgmse, 0, linewidths=1, colors='k')
    plt.title(r'$-\overline{D} \Delta \partial_{\phi} m_{925\,\mathrm{hPa}}$ [PW], RCP8.5 DJF')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Latitude (deg)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(2006,2299)
    ax.set_ylim(latlo,90)
    # ax.set_xticks(np.arange(0,361,60))
    ax.set_yticks(np.arange(latlo,91,(90-latlo)/6))
    cbar = plt.colorbar(csf)
    fig.set_size_inches(4, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.close()

    # DIFFUSIVITY COMPONENT
    plotname = '%s/dcdiffv_time_lat.%s' % (plotdir, timemean)
    fig, ax = plt.subplots()
    csf = plt.contourf(mesh_time, mesh_lat, 1e-15*dcdiffv, np.arange(-0.25,0.25,0.01), cmap='RdBu_r', extend='both')
    plt.contour(mesh_time, mesh_lat, 1e-15*dcdiffv, np.arange(-2,2,0.05), linewidths=0.5, colors='0.8')
    plt.contour(mesh_time, mesh_lat, 1e-15*dcdiffv, 0, linewidths=1, colors='k')
    plt.title(r'$-\Delta D \overline{ \partial_{\phi} m_{925\,\mathrm{hPa}} }$ [PW], RCP8.5 DJF')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Latitude (deg)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(2006,2299)
    ax.set_ylim(latlo,90)
    # ax.set_xticks(np.arange(0,361,60))
    ax.set_yticks(np.arange(latlo,91,(90-latlo)/6))
    cbar = plt.colorbar(csf)
    fig.set_size_inches(4, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.close()

    # RESIDUAL 
    plotname = '%s/dcres_time_lat.%s' % (plotdir, timemean)
    fig, ax = plt.subplots()
    csf = plt.contourf(mesh_time, mesh_lat, 1e-15*dcres, np.arange(-0.25,0.25,0.01), cmap='RdBu_r', extend='both')
    plt.contour(mesh_time, mesh_lat, 1e-15*dcres, np.arange(-2,2,0.05), linewidths=0.5, colors='0.8')
    plt.contour(mesh_time, mesh_lat, 1e-15*dcres, 0, linewidths=1, colors='k')
    plt.title(r'Residual [PW], RCP8.5 DJF')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Latitude (deg)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(2006,2299)
    ax.set_ylim(latlo,90)
    # ax.set_xticks(np.arange(0,361,60))
    ax.set_yticks(np.arange(latlo,91,(90-latlo)/6))
    cbar = plt.colorbar(csf)
    fig.set_size_inches(4, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.close()

    ############################################
    # EVALUATE DECOMPOSITION AT 80 N
    ############################################
    fvmte = interp1d(lat, diffvmte)
    diffvmtesel = fvmte(slat)
    fdcdiffv = interp1d(lat, dcdiffv)
    dcdiffvsel = fdcdiffv(slat)
    fdcgmse = interp1d(lat, dcgmse)
    dcgmsesel = fdcgmse(slat)
    fdcres = interp1d(lat, dcres)
    dcressel = fdcres(slat)

    # Diffusivity
    plotname = '%s/%s_decomp_time_lat_%g.%s' % (plotdir, transname, slat, timemean)
    fig, ax = plt.subplots()
    plt.axhline(0, color='k', linewidth=0.5)
    plt.plot(time, 1e-15*diffvmtesel, '-k', label=r'$\Delta [\overline{vm}]$')
    plt.plot(time, 1e-15*dcgmsesel, ':k', label=r'-$\overline{D} \Delta \partial_{\phi} m_{925\,\mathrm{hPa}}$')
    plt.plot(time, 1e-15*dcdiffvsel, '--k', label=r'-$\Delta D \overline{ \partial_{\phi} m_{925\,\mathrm{hPa}} }$')
    plt.plot(time, 1e-15*dcressel, '-.k', label=r'Residual')
    plt.title('$\phi=%g^\circ$N, RCP8.5 DJF' % (slat))
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Time (yr)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel(dtranslabel)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(2006,2299)
    ax.legend(loc='lower left', prop={'size':8})
    fig.set_size_inches(4, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.close()

