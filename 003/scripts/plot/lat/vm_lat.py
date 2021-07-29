import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
import numpy as np
from misc.filenames import *
from misc.dirnames import *
from misc.translate import *
from plot.titles import make_title_sim_time
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def vm_lat(sim, **kwargs):

    categ = 'lat'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    legend = kwargs.get('legend', 0) # draw legend?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    sim_ref = kwargs.get('sim_ref', 'historical')
    timemean_ref = kwargs.get('timemean_ref', 'ymonmean-30')
    yr_span_ref = kwargs.get('yr_span_ref', '186001-200512')

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
    
    vm = {}
    grid = {}
    varnames = ['aht', 'vmmmc', 'vmse', 'vmte']
    for varname in varnames:
        file = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)
        if varname == 'aht':
            varname_str = 'vE'
        else:
            varname_str = varname
        vm[varname] = np.nanmean(file.variables[varname_str][:],axis=0)

    grid['lat'] = file.variables[latetrans_grid(sim, 'lat')][:]
        
    plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

    ############################################
    # PLOT (CLIMATOLOGY)
    ############################################
    plotname = remove_repdots('%s/vm.%s' % (plotdir, timemean))
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linewidth=0.5)
    lp_vm = ax.plot(grid['lat'], 1e-15*vm['aht'], color='k', label='Total')
    lp_vmmmc = ax.plot(grid['lat'], 1e-15*vm['vmmmc'], color='g', label='MOC')
    lp_vmse = ax.plot(grid['lat'], 1e-15*vm['vmse'], color='b', label='Stationary Eddy')
    lp_vmte = ax.plot(grid['lat'], 1e-15*vm['vmte'], color='r', label='Transients')
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Energy transport (PW)')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim([-90,90])
    # ax.set_ylim(vmin_dev,vmax_dev)
    if legend:
        ax.legend()
    plt.tight_layout()
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()

    # if refclim is not None: