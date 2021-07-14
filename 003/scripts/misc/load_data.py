import sys
import os
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.translate import *
from misc.dirnames import *
from misc.filenames import *
from misc.means import *
from misc import par
from proc.r1 import save_r1
import pickle

def load_var(sim, varname, **kwargs):

    model = kwargs.get('model') # name of model
    zonmean = kwargs.get('zonmean', 'zonmean') # do zonal mean? (bool)
    timemean = kwargs.get('timemean', '') # do annual mean? (bool)
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

def load_r1(sim, categ, **kwargs):

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
        r1_file = remove_repdots('%s/r1.%s.%s.pickle' % (datadir, zonmean, timemean))

        if not (os.path.isfile(r1_file) and try_load):
            save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [r1, grid] = pickle.load(open(r1_file, 'rb'))

        return r1, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        modelstr = 'CMIP5 mean'
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
            save_r1_dc(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [r1_dc, grid] = pickle.load(open(r1_dc_file, 'rb'))

        return r1_dc, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        modelstr = 'CMIP5 mean'
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
                save_r1_dc(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

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
            save_flux(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

        [flux, grid] = pickle.load(open(flux_file, 'rb'))

        return flux, grid, datadir, plotdir, modelstr

    else: # multi-model mean
        modellist = model
        model = 'mmm'
        modelstr = 'CMIP5 mean'
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
                save_flux(sim, model=modelname, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

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
