import sys
from tqdm import tqdm
from isvertvar import isvertvar
import numpy as np
import pickle
from scipy.interpolate import interp1d
from netCDF4 import Dataset

prefix = '/project2/tas1/miyawaki/projects/003/data/raw/historical'
suffix = '.zonmean.djfmean.nc'

ntime = 147
grid = {}
grid['lat'] = np.linspace(80,90,101)
grid['lev'] = np.logspace(3,5,101)

models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MPI-ESM-LR','GISS-E2-H','GISS-E2-R']
# models = []
varnames = ['rsdt', 'rsut', 'rsutcs', 'rlut', 'rlutcs', 'rsds', 'rsdscs', 'rsus', 'rsuscs', 'rlus','rlds', 'rldscs', 'adv', 'ftoa', 'ftoacs', 'fsfc','tas', 'ts', 'ra', 'racs','hfls','hfss', 'tendv', 'hur', 'hus', 'ta']
clim = 'historical'
yr_span = '186001-200512'

mmm = {}
modeldata = {}
for varname in varnames:
    if isvertvar(varname):
        # mmm[varname] = np.empty([len(grid['lev'])])
        modeldata[varname] = np.empty([len(models), ntime, len(grid['lev'])])
    else:
        modeldata[varname] = np.empty([len(models), ntime])

    for m in tqdm(range(len(models))):
        model = models[m]

        f = Dataset('%s/%s/%s_Amon_%s_%s_r1i1p1_%s%s' % (prefix, model, varname, model, clim, yr_span, suffix) )
        raw = np.squeeze(f.variables[varname][:])
        if raw.shape[0] == 146:
            raw = np.append(raw, raw[-1,...][None,:], axis=0)
        lat = np.squeeze(f.variables['lat'][:])
        if isvertvar(varname):
            lev = np.squeeze(f.variables['plev'][:])
            lataxis = 2
            clat = np.cos(np.radians(grid['lat']))[None, None, :]
        else:
            lataxis = 1
            clat = np.cos(np.radians(grid['lat']))[None, :]

        # area average high latitudes
        fint = interp1d(lat, raw, axis=lataxis, bounds_error=False)
        rawi = fint(grid['lat'])
        idxnan = np.invert(np.isnan(rawi)).astype(int)
        rawa = np.nansum(idxnan*clat*rawi, axis=lataxis) / np.nansum(idxnan*clat, axis=lataxis)

        # interpolate to standard vertical grid if relevant
        if isvertvar(varname):
            if model in ['CCSM4','GISS-E2-H','GISS-E2-R']:
                for itime in tqdm(range(ntime)):
                    rawa0 = rawa[itime,:]
                    nanfilt = ~np.isnan(rawa0)
                    lev0 = lev[nanfilt]
                    rawa0 = rawa0[nanfilt]
                    fint = interp1d(lev0, rawa0, kind='linear', bounds_error=False, fill_value='extrapolate')
                    modeldata[varname][m,itime,:] = fint(grid['lev'])
            else:
                fint = interp1d(lev, rawa, kind='linear', bounds_error=False)
                modeldata[varname][m,:,:] = fint(grid['lev'])
        else:
            modeldata[varname][m,:] = rawa

    # multimodel mean
    mmm[varname] = np.nanmean(modeldata[varname], axis=0)

    print(varname)
    print(mmm[varname])
    print(modeldata[varname].shape)

pickle.dump([mmm, grid, modeldata], open('./clima.pickle', 'wb'))
