import sys
from tqdm import tqdm
from isvertvar import isvertvar
import numpy as np
import pickle
from scipy.interpolate import interp1d
# from netCDF4 import Dataset
import xarray as xr

prefix = '/project2/tas1/miyawaki/projects/003/data/raw/historical'
suffix = '.zonmean.djfmean.nc'

ntime = 156
grid = {}
grid['lat'] = np.linspace(-30,30,101)
grid['lev'] = np.logspace(3,5,101)

models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5', 'CESM2-WACCM', 'IPSL-CM6A-LR', 'MRI-ESM2-0','MIROC-ES2L','GISS-E2-1-G','GISS-E2-1-H','UKESM1-0-LL']
# models = ['GISS-E2-1-H']

varnames = ['rsdt', 'rsut', 'rsutcs', 'rlut', 'rlutcs', 'rsds', 'rsdscs', 'rsus', 'rsuscs', 'rlus','rlds', 'rldscs', 'adv', 'ftoa', 'ftoacs', 'fsfc','tas', 'ts', 'ra', 'racs','hfls','hfss', 'tendv', 'hur', 'hus', 'ta']
# varnames = ['hus']
clim = 'historical'
yr_span = '186001-201412'

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

        # print('%s/%s/%s_Amon_%s_%s_*_%s%s' % (prefix, model, varname, model, clim, yr_span, suffix))
        f = xr.open_mfdataset('%s/%s/%s_Amon_%s_%s_*_%s%s' % (prefix, model, varname, model, clim, yr_span, suffix) )
        raw = np.squeeze(f[varname].data)
        if raw.shape[0] == 146:
            raw = np.append(raw, raw[-1,...][None,:], axis=0)
        lat = np.squeeze(f['lat'].data)
        if isvertvar(varname):
            lev = np.squeeze(f['plev'].data)
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
            if model in ['CanESM5']:
                fint = interp1d(lev, rawa, kind='linear', bounds_error=False)
                modeldata[varname][m,:,:] = fint(grid['lev'])
            else:
                for itime in tqdm(range(ntime)):
                    rawa0 = rawa[itime,:]
                    nanfilt = ~np.isnan(rawa0)
                    lev0 = lev[nanfilt]
                    rawa0 = rawa0[nanfilt]
                    fint = interp1d(lev0, rawa0, kind='linear', bounds_error=False, fill_value='extrapolate')
                    modeldata[varname][m,itime,:] = fint(grid['lev'])
        else:
            modeldata[varname][m,:] = rawa

    # multimodel mean
    mmm[varname] = np.nanmean(modeldata[varname], axis=0)

    # print(varname)
    # print(mmm[varname])

pickle.dump([mmm, grid, modeldata], open('./clima.cmip6.pickle', 'wb'))
