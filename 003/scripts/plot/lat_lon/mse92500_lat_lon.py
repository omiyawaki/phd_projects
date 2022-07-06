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

def mse92500_lat_lon(sim, **kwargs):

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
    prefix = '/project2/tas1/miyawaki/projects/003/data/raw/historical/bcc-csm1-1'

    file_mse = Dataset('%s/mse92500_Amon_bcc-csm1-1_historical_r1i1p1_186001-200512.ymonmean-30.djfmean.nc' % (prefix), 'r')                                       
    mse = np.squeeze(file_mse.variables['mse'][:]) # (lat x lon)
    lat = np.squeeze(file_mse.variables['lat'][:]) # (lat)
    lon = np.squeeze(file_mse.variables['lon'][:]) # (lon)

    ############################################
    # PLOT
    ############################################
    plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)
    [mesh_lat, mesh_lon] = np.meshgrid(lat, lon) # create mesh

    plotname = '%s/mse92500_lat_lon.%s' % (plotdir, timemean)

    fig, ax = plt.subplots()
    csf = plt.contourf(mesh_lon, mesh_lat, 1e-3*np.transpose(mse))
    plt.title(r'$m_{925\,\mathrm{hPa}}$ [J g$^{-1}$], 1979-2005 DJF Climatology')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Longitude (deg)')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('Latitude (deg)')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(0,360)
    ax.set_ylim(-90,90)
    ax.set_xticks(np.arange(0,361,60))
    ax.set_yticks(np.arange(-90,91,30))
    cbar = plt.colorbar(csf)
    fig.set_size_inches(4, 3)
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.close()

