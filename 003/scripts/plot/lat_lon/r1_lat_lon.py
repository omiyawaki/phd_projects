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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                        LatitudeLocator, LongitudeLocator)

def r1_lat_lon(sim, **kwargs):

    categ = 'lat_lon'

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

    lat_int = np.arange(latbnd[0], latbnd[1], latstep)
    lon_int = np.arange(lonbnd[0], lonbnd[1], lonstep)

    vmin=-1.7
    vmax=1.7

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
        if 'ymonmean' not in timemean:
            yr_base = int(yr_span[0:4])
        else:
            yr_base = 0

    ##################################
    # LOAD DATA
    ##################################
    if isinstance(model, str) or model is None:
        [r1, grid, datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
    else:
        [r1, grid, datadir, plotdir, modelstr, r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 

    # ANN
    r1_ann = np.nanmean(r1, axis=0)
    r1_djf = np.nanmean(np.roll(r1,1,axis=0)[0:3,...], axis=0)
    r1_mam = np.nanmean(r1[2:5,...], axis=0)
    r1_jja = np.nanmean(r1[5:8,...], axis=0)
    r1_son = np.nanmean(r1[8:11,...], axis=0)

    # if 'pr' in plotover or 'cl' in plotover:
    #     if isinstance(model, str) or model is None:
    #         [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
    #     else:
    #         [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

    ##################################
    # CATEGORIZE REGIONS INTO ANNUAL RCE, SEASONAL RCE
    ##################################

    categ = np.empty([r1.shape[1], r1.shape[2]])

    categ[np.logical_and(np.min(r1, axis=0)<0.1, np.max(r1, axis=0)>0.9)] = 5 # label seasonal RCE/RCAE/RAE as 5

    categ[np.logical_and(np.min(r1, axis=0)<0.1, np.max(r1, axis=0)<0.9)] = 3 # label seasonal RCE/RCAE as 3
    categ[np.logical_and(np.min(r1, axis=0)>0.1, np.max(r1, axis=0)>0.9)] = 4 # label seasonal RAE/RCAE as 4

    categ[np.max(r1, axis=0)<0.1] = 0 # label yearround RCE as 0
    categ[np.min(r1, axis=0)>0.9] = 1 # label yearround RAE as 1
    categ[np.logical_and(np.min(r1, axis=0)>0.1, np.max(r1, axis=0)<0.9)] = 2 # label yearround RCAE as 2

    n_cats = 6

    ############################################
    # PLOT R1 ANNMEAN
    ############################################
    [mesh_lat, mesh_lon] = np.meshgrid(grid['lat'], grid['lon']) # create mesh
    plotname = '%s/r1_lat_lon.ann.%s' % (plotdir, timemean)
    fig = plt.figure()
    axc = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    # axc.set_extent([-180, 180, -30, 30], ccrs.PlateCarree())
    csf = plt.contourf(mesh_lon, mesh_lat, np.transpose(r1_ann),np.arange(vmin,vmax,0.1), vmin=vmin, vmax=vmax, cmap='RdBu', extend='both', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_ann), [0.1], linewidths=0.5, colors='sandybrown', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_ann), [0.9], linewidths=0.5, colors='royalblue', transform=ccrs.PlateCarree())
    make_title_sim_time_seas(axc, sim, model=modelstr, timemean=timemean, seasmean='')
    cbar = plt.colorbar(csf)
    cbar.set_label('$R_1$')
    axc.coastlines(color='gray', linewidth=0.5)
    gl = axc.gridlines(draw_labels=True, linewidth=0.5, alpha=0.2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = LongitudeLocator()
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    fig.set_size_inches(6, 3)
    # plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT R1 DJF
    ############################################
    [mesh_lat, mesh_lon] = np.meshgrid(grid['lat'], grid['lon']) # create mesh
    plotname = '%s/r1_lat_lon.djf.%s' % (plotdir, timemean)
    fig = plt.figure()
    axc = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    # axc.set_extent([-180, 180, -30, 30], ccrs.PlateCarree())
    csf = plt.contourf(mesh_lon, mesh_lat, np.transpose(r1_djf),np.arange(vmin,vmax,0.1), vmin=vmin, vmax=vmax, cmap='RdBu', extend='both', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_djf), [0.1], linewidths=0.5, colors='sandybrown', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_djf), [0.9], linewidths=0.5, colors='royalblue', transform=ccrs.PlateCarree())
    make_title_sim_time_seas(axc, sim, model=modelstr, timemean=timemean, seasmean='djf')
    cbar = plt.colorbar(csf)
    cbar.set_label('$R_1$')
    axc.coastlines(color='gray', linewidth=0.5)
    gl = axc.gridlines(draw_labels=True, linewidth=0.5, alpha=0.2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = LongitudeLocator()
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    fig.set_size_inches(6, 3)
    # plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT R1 mam
    ############################################
    [mesh_lat, mesh_lon] = np.meshgrid(grid['lat'], grid['lon']) # create mesh
    plotname = '%s/r1_lat_lon.mam.%s' % (plotdir, timemean)
    fig = plt.figure()
    axc = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    # axc.set_extent([-180, 180, -30, 30], ccrs.PlateCarree())
    csf = plt.contourf(mesh_lon, mesh_lat, np.transpose(r1_mam),np.arange(vmin,vmax,0.1), vmin=vmin, vmax=vmax, cmap='RdBu', extend='both', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_mam), [0.1], linewidths=0.5, colors='sandybrown', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_mam), [0.9], linewidths=0.5, colors='royalblue', transform=ccrs.PlateCarree())
    make_title_sim_time_seas(axc, sim, model=modelstr, timemean=timemean, seasmean='mam')
    cbar = plt.colorbar(csf)
    cbar.set_label('$R_1$')
    axc.coastlines(color='gray', linewidth=0.5)
    gl = axc.gridlines(draw_labels=True, linewidth=0.5, alpha=0.2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = LongitudeLocator()
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    fig.set_size_inches(6, 3)
    # plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT R1 jja
    ############################################
    [mesh_lat, mesh_lon] = np.meshgrid(grid['lat'], grid['lon']) # create mesh
    plotname = '%s/r1_lat_lon.jja.%s' % (plotdir, timemean)
    fig = plt.figure()
    axc = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    # axc.set_extent([-180, 180, -30, 30], ccrs.PlateCarree())
    csf = plt.contourf(mesh_lon, mesh_lat, np.transpose(r1_jja),np.arange(vmin,vmax,0.1), vmin=vmin, vmax=vmax, cmap='RdBu', extend='both', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_jja), [0.1], linewidths=0.5, colors='sandybrown', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_jja), [0.9], linewidths=0.5, colors='royalblue', transform=ccrs.PlateCarree())
    make_title_sim_time_seas(axc, sim, model=modelstr, timemean=timemean, seasmean='jja')
    cbar = plt.colorbar(csf)
    cbar.set_label('$R_1$')
    axc.coastlines(color='gray', linewidth=0.5)
    gl = axc.gridlines(draw_labels=True, linewidth=0.5, alpha=0.2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = LongitudeLocator()
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    fig.set_size_inches(6, 3)
    # plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT R1 son
    ############################################
    [mesh_lat, mesh_lon] = np.meshgrid(grid['lat'], grid['lon']) # create mesh
    plotname = '%s/r1_lat_lon.son.%s' % (plotdir, timemean)
    fig = plt.figure()
    axc = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    # axc.set_extent([-180, 180, -30, 30], ccrs.PlateCarree())
    csf = plt.contourf(mesh_lon, mesh_lat, np.transpose(r1_son),np.arange(vmin,vmax,0.1), vmin=vmin, vmax=vmax, cmap='RdBu', extend='both', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_son), [0.1], linewidths=0.5, colors='sandybrown', transform=ccrs.PlateCarree())
    plt.contour(mesh_lon, mesh_lat, np.transpose(r1_son), [0.9], linewidths=0.5, colors='royalblue', transform=ccrs.PlateCarree())
    make_title_sim_time_seas(axc, sim, model=modelstr, timemean=timemean, seasmean='son')
    cbar = plt.colorbar(csf)
    cbar.set_label('$R_1$')
    axc.coastlines(color='gray', linewidth=0.5)
    gl = axc.gridlines(draw_labels=True, linewidth=0.5, alpha=0.2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = LongitudeLocator()
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    fig.set_size_inches(6, 3)
    # plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()


    ############################################
    # PLOT
    ############################################
    [mesh_lat, mesh_lon] = np.meshgrid(grid['lat'], grid['lon']) # create mesh

    plotname = '%s/categ_lat_lon.%s' % (plotdir, timemean)

    fig, ax = plt.subplots()
    axc = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    # axc.set_extent([-180, 180, -30, 30], ccrs.PlateCarree())

    csf = plt.contourf(mesh_lon, mesh_lat, np.transpose(categ), levels=np.arange(0,n_cats+1,1)-0.5, cmap='Pastel1', transform=ccrs.PlateCarree())
    make_title_sim_time(axc, sim, model=modelstr, timemean=timemean)
    # ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    # ax.set_xlabel('Longitude (deg)')
    # ax.xaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_ylabel('Latitude (deg)')
    # ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(0,360)
    # ax.set_ylim(-30,30)
    # ax.set_yticks(np.arange(-90,91,30))
    cbar = plt.colorbar(csf)
    tick_locs = np.arange(n_cats+1)#*(n_cats-1)/n_cats
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(['Yearround\nRCE', 'Yearround\nRAE', 'Yearround\nRCAE', 'Seasonal\nRCE/RCAE', 'Seasonal\nRAE/RCAE', 'Seasonal\nRCE/RCAE/RAE'])

    axc.coastlines(color='gray', linewidth=0.5)
    gl = axc.gridlines(draw_labels=True, linewidth=0.5, alpha=0.2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = LongitudeLocator()
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}


    fig.set_size_inches(6, 3)
    # plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # PLOT (with monsoon boxes)
    ############################################
    [mesh_lat, mesh_lon] = np.meshgrid(grid['lat'], grid['lon']) # create mesh

    plotname = '%s/categ_lat_lon_box.%s' % (plotdir, timemean)

    fig, ax = plt.subplots()
    axc = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    # axc.set_extent([-180, 180, -30, 30], ccrs.PlateCarree())

    csf = plt.contourf(mesh_lon, mesh_lat, np.transpose(categ), levels=np.arange(0,n_cats+1,1)-0.5, cmap='Pastel1', transform=ccrs.PlateCarree())
    make_title_sim_time(axc, sim, model=modelstr, timemean=timemean)
    # ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    # ax.set_xlabel('Longitude (deg)')
    # ax.xaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_ylabel('Latitude (deg)')
    # ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(0,360)
    # ax.set_ylim(-30,30)
    # ax.set_yticks(np.arange(-90,91,30))
    cbar = plt.colorbar(csf)
    tick_locs = np.arange(n_cats+1)#*(n_cats-1)/n_cats
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(['Yearround\nRCE', 'Yearround\nRAE', 'Yearround\nRCAE', 'Seasonal\nRCE/RCAE', 'Seasonal\nRAE/RCAE', 'Seasonal\nRCE/RCAE/RAE'])

    safr=patches.Rectangle((10,-20),40,15,fill=False,transform=ccrs.PlateCarree())
    axc.add_patch(safr)

    aus=patches.Rectangle((120,-20),30,15,fill=False,transform=ccrs.PlateCarree())
    axc.add_patch(aus)

    sam=patches.Rectangle((285,-20),30,15,fill=False,transform=ccrs.PlateCarree())
    axc.add_patch(sam)

    sam=patches.Rectangle((0,5),30,10,fill=False,transform=ccrs.PlateCarree())
    axc.add_patch(sam)

    wafr=patches.Rectangle((0,5),30,10,fill=False,transform=ccrs.PlateCarree())
    axc.add_patch(wafr)

    nam=patches.Rectangle((250,5),25,15,fill=False,transform=ccrs.PlateCarree())
    axc.add_patch(nam)

    sas=patches.Rectangle((70,10),30,20,fill=False,transform=ccrs.PlateCarree())
    axc.add_patch(sas)

    axc.coastlines(color='gray', linewidth=0.5)
    gl = axc.gridlines(draw_labels=True, linewidth=0.5, alpha=0.2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = LongitudeLocator()
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}


    fig.set_size_inches(6, 3)
    # plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)

    if viewplt:
        plt.show()
    plt.close()
