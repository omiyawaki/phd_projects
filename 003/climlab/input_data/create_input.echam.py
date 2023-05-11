import sys
from tqdm import tqdm
from isvertvar import isvertvar
import numpy as np
import pickle
from scipy.interpolate import interp1d
from netCDF4 import Dataset

prefix = '/project2/tas1/miyawaki/projects/003/data/raw/echam'
suffix = '.zonmean.ymonmean-20.djfmean.nc'

lat0=70
lat1=90

grid = {}
grid['lat'] = np.linspace(lat0,lat1,101)
grid['lev'] = np.logspace(3,5,101)

models = ['rp000134','rp000190f']
varnames = ['srad0', 'trad0', 'sraf0', 'traf0', 'srads','trads','srafs','trafs', 'adv', 'ftoa', 'ftoacs', 'fsfc','temp2', 'tsurf', 'ra', 'racs','ahfl','ahfs', 'tendv', 'rhumidity', 'q', 't']
yr_span = '0020_0039'

for m in tqdm(range(len(models))):
    model = models[m]

    modeldata = {}
    for varname in varnames:
        if isvertvar(varname):
            modeldata[varname] = np.empty(len(grid['lev']))

            f = Dataset('%s/%s/%s_%s_%s%s' % (prefix, model, varname, model, yr_span, suffix) )
            raw = np.squeeze(f.variables[varname][:])
            lat = np.squeeze(f.variables['lat'][:])
            if isvertvar(varname):
                try:
                    lev = np.squeeze(f.variables['plev'][:])
                except:
                    lev = np.squeeze(f.variables['lev'][:])
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
                modeldata[varname] = fint(grid['lev'])
            else:
                modeldata[varname] = rawa

    pickle.dump([modeldata, grid], open('./indata.%s.%g.%g.pickle'%(model,lat0,lat1), 'wb'))
