import sys
from tqdm import tqdm
from isvertvar import isvertvar
import numpy as np
import pickle
from scipy.interpolate import interp1d
from netCDF4 import Dataset

prefix = '/project2/tas1/miyawaki/projects/003/data/raw/historical'
suffix = '.zonmean.ymonmean-30.djfmean.nc'

grid = {}
grid['lat'] = np.linspace(80,90,101)
grid['lev'] = np.logspace(3,5,101)
print(grid['lev'])

models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MPI-ESM-LR']
varnames = ['rsdt', 'rsut', 'rsutcs', 'rlut', 'rlutcs', 'rsds', 'rsdscs', 'rsus', 'rsuscs', 'rlus','rlds', 'rldscs', 'adv', 'ftoa', 'ftoacs', 'fsfc','tas', 'ts', 'ra', 'racs','hfls','hfss', 'tendv', 'hur', 'hus', 'ta']
clim = 'historical'
yr_span = '186001-200512'

mmm = {}
modeldata = {}
for varname in varnames:
    if isvertvar(varname):
        # mmm[varname] = np.empty([len(grid['lev'])])
        modeldata[varname] = np.empty([len(models), len(grid['lev'])])
    else:
        modeldata[varname] = np.empty([len(models)])

    for m in tqdm(range(len(models))):
        model = models[m]

        f = Dataset('%s/%s/%s_Amon_%s_%s_r1i1p1_%s%s' % (prefix, model, varname, model, clim, yr_span, suffix) )
        raw = np.squeeze(f.variables[varname][:])
        lat = np.squeeze(f.variables['lat'][:])
        if isvertvar(varname):
            lev = np.squeeze(f.variables['plev'][:])
            lataxis = 1
            clat = np.cos(np.radians(grid['lat']))[None, :]
        else:
            lataxis = 0
            clat = np.cos(np.radians(grid['lat']))

        # area average high latitudes
        fint = interp1d(lat, raw, axis=lataxis, bounds_error=False)
        rawi = fint(grid['lat'])
        rawa = np.nansum(clat*rawi, axis=lataxis) / np.nansum(clat, axis=lataxis)

        # interpolate to standard vertical grid if relevant
        if isvertvar(varname):
            fint = interp1d(lev, rawa, kind='linear', bounds_error=False)
            modeldata[varname][m,:] = fint(grid['lev'])
        else:
            modeldata[varname][m] = rawa

    # multimodel mean
    mmm[varname] = np.nanmean(modeldata[varname], axis=0)

    print(varname)
    print(mmm[varname])

pickle.dump([mmm, grid], open('./indata.pickle', 'wb'))
