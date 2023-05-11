import sys
import os
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.translate import *
from misc.dirnames import *
from misc.filenames import *
from misc.means import *
from misc import par
from proc.r1 import save_r1
from proc.rad import save_rad
from proc.ke import save_ke
from proc.dyn import save_dyn
from proc.diffv import save_diffv
from proc.circ import save_circ
from proc.hydro import save_hydro
from proc.cld import save_cld
from proc.cl_in import save_cl_in
from proc.co2 import save_co2
from proc.seaice import save_seaice
from proc.ga import make_ga_dev_vint
import pickle

def load_var(sim, varname, **kwargs):

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', 'zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    yr_span = kwargs.get('yr_span') # considered span of years
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

    # directory to save pickled data
    datadir = get_datadir(sim, model=model, yr_span=yr_span)

    # location of pickled data if available
    file = remove_repdots('%s/%s.%s.%s.pickle' % (datadir, varname, zonmean, timemean))

    if (os.path.isfile(file) and try_load):
        alldata = pickle.load(open(file, 'rb'))
    else:
        file_raw = {}
        grid = {}
        alldata = {}

        file_raw[varname] = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)

        grid['lat'] = file_raw[varname].variables[latetrans_grid(sim, 'lat')][:]
        grid['lon'] = file_raw[varname].variables[latetrans_grid(sim, 'lon')][:]
        try: # try to load vertical grid  data
            grid['lev'] = file_raw[varname].variables[latetrans_grid(sim, 'lev')][:]
            if sim == 'era5':
                grid['lev'] = 100 * grid['lev'] # convert hPa to Pa
        except:
            print('%s does not contain vertical level data.' % (varname))

        vardata = file_raw[varname].variables[varname][:]

        if sim == 'era5' and varname in ['z', 'zs']:
            vardata = 1/par.g * vardata # convert geopotential (m**2 s**-2) to geopotential height (m)

        alldata[translate_varname(varname)] = vardata
        alldata['grid'] = grid

        pickle.dump(alldata, open(remove_repdots('%s/%s.%s.%s.pickle' % (datadir, varname, zonmean, timemean)), 'wb'), protocol=4)

    return alldata

#
# R1
#

def load_r1(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        r1_file = remove_repdots('%s/r1.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(r1_file) and try_load):
            save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span, refclim=refclim)

        [r1, grid] = pickle.load(open(r1_file, 'rb'))

        return r1, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        r1list = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            r1_file = remove_repdots('%s/r1.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(r1_file) and try_load):
                save_r1(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [r1list[i], gridlist[i]] = pickle.load(open(r1_file, 'rb'))
        r1_mmm = mmm_mon_lat(r1list, gridlist, grid)
        r1 = r1_mmm['mmm']
        
        return r1, grid, datadir, plotdir, modelstr, r1_mmm

#
# R1a
#

def load_r1a(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        r1a_file = remove_repdots('%s/r1a.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(r1a_file) and try_load):
            save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span, refclim=refclim)

        [r1a, grid] = pickle.load(open(r1a_file, 'rb'))

        return r1a, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        r1alist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            r1a_file = remove_repdots('%s/r1a.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(r1a_file) and try_load):
                save_r1(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [r1alist[i], gridlist[i]] = pickle.load(open(r1a_file, 'rb'))
        r1a_mmm = mmm_mon_lat(r1alist, gridlist, grid)
        r1a = r1a_mmm['mmm']
        
        return r1a, grid, datadir, plotdir, modelstr, r1a_mmm

#
# R1 DECOMP
#

def load_r1_dc(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str): # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        r1_dc_file = remove_repdots('%s/r1_dc.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(r1_dc_file) and try_load):
            save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [r1_dc, grid] = pickle.load(open(r1_dc_file, 'rb'))

        return r1_dc, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        r1_dclist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            r1_dc_file = remove_repdots('%s/r1_dc.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(r1_dc_file) and try_load):
                save_r1(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [r1_dclist[i], gridlist[i]] = pickle.load(open(r1_dc_file, 'rb'))
            
        r1_dc_mmm = {}
        r1_dc = {}
        for varname in r1_dclist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = r1_dclist[i][varname]
            r1_dc_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            r1_dc[varname] = r1_dc_mmm[varname]['mmm']
        
        return r1_dc, grid, datadir, plotdir, modelstr, r1_dc_mmm
    
#
# R1 DERIVATIVE
#

def load_dr1(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str): # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        dr1_file = remove_repdots('%s/dr1.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(dr1_file) and try_load):
            save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [dr1, grid] = pickle.load(open(dr1_file, 'rb'))

        return dr1, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        dr1list = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            dr1_file = remove_repdots('%s/dr1.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(dr1_file) and try_load):
                save_r1(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [dr1list[i], gridlist[i]] = pickle.load(open(dr1_file, 'rb'))
            
        dr1_mmm = {}
        dr1 = {}
        for varname in dr1list[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = dr1list[i][varname]
            dr1_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            dr1[varname] = dr1_mmm[varname]['mmm']
        
        return dr1, grid, datadir, plotdir, modelstr, dr1_mmm
    
#
# ENERGY FLUXES
#

def load_flux(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str): # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        flux_file = remove_repdots('%s/flux.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(flux_file) and try_load):
            save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [flux, grid] = pickle.load(open(flux_file, 'rb'))

        return flux, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        print(sim)
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        else:
            modelstr = ''
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        fluxlist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            flux_file = remove_repdots('%s/flux.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(flux_file) and try_load):
                save_r1(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [fluxlist[i], gridlist[i]] = pickle.load(open(flux_file, 'rb'))
        
        flux_mmm = {}
        flux = {}
        for varname in fluxlist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = fluxlist[i][varname]
            flux_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            flux[varname] = flux_mmm[varname]['mmm']
        
        return flux, grid, datadir, plotdir, modelstr, flux_mmm

#
# CL IN
#

def load_cl_in(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str): # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        cl_in_file = remove_repdots('%s/cl_in.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(cl_in_file) and try_load):
            save_cl_in(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [cl_in, grid] = pickle.load(open(cl_in_file, 'rb'))

        return cl_in, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid0['lev'] = par.p50
        grid = kwargs.get('gridstd', grid0)

        cl_inlist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            cl_in_file = remove_repdots('%s/cl_in.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(cl_in_file) and try_load):
                save_cl_in(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [cl_inlist[i], gridlist[i]] = pickle.load(open(cl_in_file, 'rb'))
        
        cl_in_mmm = {}
        cl_in = {}
        for varname in cl_inlist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = cl_inlist[i][varname]
            cl_in_mmm[varname] = mmm_mon_lat_lev(varnamelist, gridlist, grid)
            cl_in[varname] = cl_in_mmm[varname]['mmm']
        
        return cl_in, grid, datadir, plotdir, modelstr, cl_in_mmm

#
# RADIATION
#

def load_rad(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        rad_file = remove_repdots('%s/rad.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(rad_file) and try_load):
            save_rad(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [rad, grid] = pickle.load(open(rad_file, 'rb'))

        return rad, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        else:
            modelstr = ''
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        radlist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            rad_file = remove_repdots('%s/rad.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(rad_file) and try_load):
                save_rad(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [radlist[i], gridlist[i]] = pickle.load(open(rad_file, 'rb'))
        
        rad_mmm = {}
        rad = {}
        for varname in radlist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = radlist[i][varname]
            rad_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            rad[varname] = rad_mmm[varname]['mmm']
        
        return rad, grid, datadir, plotdir, modelstr, rad_mmm
    
#
# KINETIC ENERGY
#

def load_ke(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        ke_file = remove_repdots('%s/ke.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(ke_file) and try_load):
            save_ke(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [ke, grid] = pickle.load(open(ke_file, 'rb'))

        return ke, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        kelist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            ke_file = remove_repdots('%s/ke.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(ke_file) and try_load):
                save_ke(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [kelist[i], gridlist[i]] = pickle.load(open(ke_file, 'rb'))
        
        ke_mmm = {}
        ke = {}
        for varname in kelist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = kelist[i][varname]
            ke_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            ke[varname] = ke_mmm[varname]['mmm']
        
        return ke, grid, datadir, plotdir, modelstr, ke_mmm
    
#
# HYDRO
#

def load_hydro(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        hydro_file = remove_repdots('%s/hydro.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(hydro_file) and try_load):
            save_hydro(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [hydro, grid] = pickle.load(open(hydro_file, 'rb'))

        return hydro, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        else:
            modelstr = ''
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        hydrolist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            hydro_file = remove_repdots('%s/hydro.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(hydro_file) and try_load):
                save_hydro(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [hydrolist[i], gridlist[i]] = pickle.load(open(hydro_file, 'rb'))
        
        hydro_mmm = {}
        hydro = {}
        for varname in hydrolist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = hydrolist[i][varname]
            hydro_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            hydro[varname] = hydro_mmm[varname]['mmm']
        
        return hydro, grid, datadir, plotdir, modelstr, hydro_mmm
    
#
# CLOUD
#

def load_cld(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        cld_file = remove_repdots('%s/cld.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(cld_file) and try_load):
            save_cld(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [cld, grid] = pickle.load(open(cld_file, 'rb'))

        return cld, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        else:
            modelstr = ''
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        cldlist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            cld_file = remove_repdots('%s/cld.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(cld_file) and try_load):
                save_cld(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [cldlist[i], gridlist[i]] = pickle.load(open(cld_file, 'rb'))
        
        cld_mmm = {}
        cld = {}
        for varname in cldlist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = cldlist[i][varname]
            cld_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            cld[varname] = cld_mmm[varname]['mmm']
        
        return cld, grid, datadir, plotdir, modelstr, cld_mmm
    
#
# CO2
#

def load_co2(sim, categ, **kwargs):

    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled co2 data
        co2_file = remove_repdots('%s/co2.%s.pickle' % (datadir, timemean))

        if not (os.path.isfile(co2_file) and try_load):
            save_co2(sim, model=model, timemean=timemean, yr_span=yr_span)

        co2 = pickle.load(open(co2_file, 'rb'))

        return co2, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        co2list = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            co2_file = remove_repdots('%s/co2.%s.%s.pickle' % (datadir, timemean))

            if not (os.path.isfile(co2_file) and try_load):
                save_co2(sim, model=modelname, timemean=timemean, yr_span=yr_span)

            co2list[i] = pickle.load(open(co2_file, 'rb'))
        
        co2_mmm = {}
        co2 = {}
        for varname in co2list[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = co2list[i][varname]
            co2_mmm[varname] = mmm_mon(varnamelist)
            co2[varname] = co2_mmm[varname]['mmm']
        
        return co2, datadir, plotdir, modelstr, co2_mmm
    
#
# SEA ICE
#

def load_seaice(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        seaice_file = remove_repdots('%s/seaice.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(seaice_file) and try_load):
            save_seaice(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [seaice, grid] = pickle.load(open(seaice_file, 'rb'))

        return seaice, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        else:
            modelstr = ''
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lon'] = par.lon360
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        seaicelist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            seaice_file = remove_repdots('%s/seaice.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(seaice_file) and try_load):
                save_seaice(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [seaicelist[i], gridlist[i]] = pickle.load(open(seaice_file, 'rb'))
        
        seaice_mmm = {}
        seaice = {}
        for varname in seaicelist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                filledvar = seaicelist[i][varname].filled(0)
                varnamelist[i] = filledvar
            seaice_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            seaice[varname] = seaice_mmm[varname]['mmm']
        
        return seaice, grid, datadir, plotdir, modelstr, seaice_mmm

#
# LARGE SCALE CIRCULATION (e.g., wap, psi)
#

def load_circ(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        circ_file = remove_repdots('%s/circ.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(circ_file) and try_load):
            save_circ(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [circ, grid] = pickle.load(open(circ_file, 'rb'))

        return circ, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lat'] = par.lat180
        grid0['lev'] = par.p50
        grid = kwargs.get('gridstd', grid0)

        circlist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            circ_file = remove_repdots('%s/circ.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(circ_file) and try_load):
                save_circ(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [circlist[i], gridlist[i]] = pickle.load(open(circ_file, 'rb'))
        
        circ_mmm = {}
        circ = {}
        for varname in circlist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = circlist[i][varname]
            circ_mmm[varname] = mmm_mon_lat_lev(varnamelist, gridlist, grid)
            circ[varname] = circ_mmm[varname]['mmm']
        
        return circ, grid, datadir, plotdir, modelstr, circ_mmm
    
#
# DYNAMICS
#

def load_dyn(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        dyn_file = remove_repdots('%s/dyn.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(dyn_file) and try_load):
            save_dyn(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [dyn, grid] = pickle.load(open(dyn_file, 'rb'))

        return dyn, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        dynlist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            dyn_file = remove_repdots('%s/dyn.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(dyn_file) and try_load):
                save_dyn(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [dynlist[i], gridlist[i]] = pickle.load(open(dyn_file, 'rb'))
        
        dyn_mmm = {}
        dyn = {}
        for varname in dynlist[0]:
            print(varname)
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                print(modellist[i])
                print(dynlist[i][varname].shape)
                if sim == 'historical' and not modellist[i] == 'HadGEM2-ES':
                    trunc = dynlist[i][varname][:-1,...]
                    print(trunc.shape)
                else:
                    trunc = dynlist[i][varname]
                varnamelist[i] = trunc
                # varnamelist[i] = dynlist[i][varname]
            dyn_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            dyn[varname] = dyn_mmm[varname]['mmm']
        
        return dyn, grid, datadir, plotdir, modelstr, dyn_mmm

#
# DIFFUSIVITY
#

def load_diffv(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        diffv_file = remove_repdots('%s/diffv.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(diffv_file) and try_load):
            save_diffv(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [diffv, grid] = pickle.load(open(diffv_file, 'rb'))

        return diffv, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lat'] = par.lat180
        grid = kwargs.get('gridstd', grid0)

        diffvlist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            diffv_file = remove_repdots('%s/diffv.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(diffv_file) and try_load):
                save_diffv(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [diffvlist[i], gridlist[i]] = pickle.load(open(diffv_file, 'rb'))
        
        diffv_mmm = {}
        diffv = {}
        for varname in diffvlist[0]:
            print(varname)
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                print(modellist[i])
                print(diffvlist[i][varname].shape)
                # if sim == 'historical' and not modellist[i] == 'HadGEM2-ES':
                #     trunc = diffvlist[i][varname][:-1,...]
                #     print(trunc.shape)
                # else:
                trunc = diffvlist[i][varname]
                varnamelist[i] = trunc
                # varnamelist[i] = diffvlist[i][varname]
            diffv_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            diffv[varname] = diffv_mmm[varname]['mmm']
        
        return diffv, grid, datadir, plotdir, modelstr, diffv_mmm

#
# LAPSE RATE
#

def load_ga(sim, categ, vertbnd, vertcoord, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        ga_file = remove_repdots('%s/ga_dev_vint.%g.%g.%s.%s.%s.pickle' % (datadir, vertbnd[0], vertbnd[1], vertcoord, zonmean, timemean))

        print(ga_file)

        # if not (os.path.isfile(ga_file) and try_load):
        #     make_ga_dev_vint(sim, vertbnd, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)
        make_ga_dev_vint(sim, vertbnd, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

        [ga, grid] = pickle.load(open(ga_file, 'rb'))

        return ga, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        else:
            modelstr = ''
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lat'] = par.lat180
        grid0['lon'] = par.lon360
        grid = kwargs.get('gridstd', grid0)

        galist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            ga_file = remove_repdots('%s/ga_dev_vint.%g.%g.%s.%s.%s.pickle' % (datadir, vertbnd[0], vertbnd[1], vertcoord, zonmean, timemean))

            if not (os.path.isfile(ga_file) and try_load):
                make_ga_dev_vint(sim, vertbnd, model=modelname, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

            [galist[i], gridlist[i]] = pickle.load(open(ga_file, 'rb'))
        
        # ga_mmm = {}
        # ga = {}
        # varnamelist = [ [] for _ in range(len(modellist)) ]
        # # gridlist = [ [] for _ in range(len(modellist)) ]
        # for i in range(len(modellist)):
        #     print(i)
        #     print(modellist[i])
        #     varnamelist[i] = galist[i]['ga_dev_vint']
        #     # varnamelist[i] = np.nanmean(galist[i]['ga_dev_vint'],2)
        #     # gridlist[i] = galist[i]['grid']
        # ga_mmm['ga_dev_vint'] = mmm_mon_lat(varnamelist, gridlist, grid)
        # ga['ga_dev_vint'] = ga_mmm['ga_dev_vint']['mmm']
        # ga['grid'] = grid

        ga_mmm = {}
        ga = {}
        for varname in galist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = galist[i]['ga_dev_vint']
            ga_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
            ga[varname] = ga_mmm[varname]['mmm']

#         hydro_mmm = {}
#         hydro = {}
#         for varname in hydrolist[0]:
#             varnamelist = [ [] for _ in range(len(modellist)) ]
#             for i in range(len(modellist)):
#                 varnamelist[i] = hydrolist[i][varname]
#             hydro_mmm[varname] = mmm_mon_lat(varnamelist, gridlist, grid)
#             hydro[varname] = hydro_mmm[varname]['mmm']

        return ga, grid, datadir, plotdir, modelstr, ga_mmm
    
#
# TEMPERATURE (e.g., tas, ta)
#

def load_dta(sim, categ, **kwargs):

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    if isinstance(model, str) or model is None: # single model
        modelstr = model

        # load data and plot directories
        datadir = get_datadir(sim, model=model, yr_span=yr_span)
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        # location of pickled R1 data
        dta_file = remove_repdots('%s/dta.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(dta_file) and try_load):
            save_dta(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [dta, grid] = pickle.load(open(dta_file, 'rb'))

        return dta, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        if sim in ['rcp85', 'hist+rcp85']:
            modelstr = 'CMIP5 mean'
        elif sim in ['ssp585', 'hist+ssp585']:
            modelstr = 'CMIP6 mean'
        plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

        grid0 = {}
        grid0['lat'] = par.lat180
        grid0['lev'] = par.p50
        grid = kwargs.get('gridstd', grid0)

        dtalist = [ [] for _ in range(len(modellist)) ]
        gridlist = [ [] for _ in range(len(modellist)) ]
        for i in range(len(modellist)):
            modelname = modellist[i]
            # load data and plot directories
            datadir = get_datadir(sim, model=modelname, yr_span=yr_span)

            # location of pickled R1 data
            dta_file = remove_repdots('%s/dta.%s.%s.pickle' % (datadir, zonmean, timemean))

            if not (os.path.isfile(dta_file) and try_load):
                save_dta(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

            [dtalist[i], gridlist[i]] = pickle.load(open(dta_file, 'rb'))
        
        dta_mmm = {}
        dta = {}
        for varname in dtalist[0]:
            varnamelist = [ [] for _ in range(len(modellist)) ]
            for i in range(len(modellist)):
                varnamelist[i] = dtalist[i][varname]
            dta_mmm[varname] = mmm_mon_lat_lev(varnamelist, gridlist, grid)
            dta[varname] = dta_mmm[varname]['mmm']
        
        return dta, grid, datadir, plotdir, modelstr, dta_mmm
    
