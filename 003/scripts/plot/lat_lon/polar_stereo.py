import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import *
from misc.filenames import *
# from misc.translate import translate_varname
from misc.means import lat_mean, global_int
from misc.load_data import *
from misc import par
from proc.r1 import save_r1
from plot.titles import make_title_sim_time_lat
import os
import pickle
import numpy as np
from cb_info import cb_info
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

def polar_stereo(sim, **kwargs):

    categ = 'polar_stereo'

    annmean = kwargs.get('annmean', 0)
    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    latbnd = kwargs.get('latbnd', (80,90))
    latstep = kwargs.get('latstep', 0.25) # latitude step size used for interpolation
    plotvar = kwargs.get('plotvar', None) # plot overlay (sic for sea ice, ga_dev for lapse rate deviation, pr for precip, decomp for linear decomposition)?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    legend = kwargs.get('legend', 0)

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

    # load convection precipitation
    if plotvar == 'prc':
        if isinstance(model, str) or model is None:
            [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        else:
            [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

        # save prc
        binvar = hydro['prc']

        # if multimodel mean save statistics
        if not ( isinstance(model, str) or model is None ):
            binvar_mmm = hydro_mmm['prc']

    # convective precipitation fraction
    if plotvar == 'prfrac':
        if isinstance(model, str) or model is None:
            [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        else:
            [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

        # save prfrac
        binvar = hydro['prfrac']

        # if multimodel mean save statistics
        if not ( isinstance(model, str) or model is None ):
            binvar_mmm = hydro_mmm['prfrac']

    # sea ice
    if plotvar == 'sic':
        if isinstance(model, str) or model is None:
            [seaice, grid, datadir, plotdir, modelstr] = load_seaice(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        else:
            [seaice, grid, datadir, plotdir, modelstr, seaice_mmm] = load_seaice(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

        # save sea ice
        binvar = seaice['sic']

        # if multimodel mean save statistics
        if not ( isinstance(model, str) or model is None ):
            binvar_mmm = seaice_mmm['sic']

    # # take annual mean?
    # if annmean:
    #     r1 = np.mean(r1, axis=0, keepdims=True)
    #     binvar = np.mean(binvar, axis=0, keepdims=True)

    # # limit data to select region?
    # lat_sel = np.where( np.logical_and(grid['lat'] >= latbnd[0], grid['lat'] <= latbnd[1]) )
    # r1 = r1[:,lat_sel[0]]
    # binvar = binvar[:,lat_sel[0]]

    # # flatten data
    # r1_flat = np.ndarray.flatten(r1)
    # binvar_flat = np.ndarray.flatten(binvar)

    # import colorbar info
    [cb_levs, cb_label] = cb_info(plotvar)

    ############################################
    # POLAR STEREOGRAPHIC PLOT
    ############################################

    if annmean:
        plotname = '%s/%s.polste.annmean.%s' % (plotdir, plotvar, timemean)
    else:
        plotname = '%s/%s.polste.%s' % (plotdir, plotvar, timemean)

    fig, ax = plt.subplots()
    # axc = plt.axes(projection=ccrs.Stereographic())
    axc = plt.axes(projection=ccrs.NorthPolarStereo())
    axc.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())

    cyclic_data, cyclic_lon = add_cyclic_point(binvar[0,...], axis=1, coord=grid['lon'])
    csf = plt.contourf(cyclic_lon, grid['lat'], cyclic_data, cb_levs, extend='both', transform=ccrs.PlateCarree())
    # plt.contourf(grid['lon'], grid['lat'], binvar[0,...], cb_levs, transform=ccrs.PlateCarree())
    # axc.add_feature(cfeature.LAND)
    # axc.add_feature(cfeature.OCEAN)
    cbar = plt.colorbar(csf)
    cbar.set_label(cb_label)
    axc.coastlines(color='gray', linewidth=0.5)
    gl = axc.gridlines(draw_labels=False, xlocs=[], ylocs=[40,50,60,70,80], linewidth=0.5, alpha=0.2)
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    

    if viewplt:
        plt.show()
    plt.close()

