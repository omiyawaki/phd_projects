import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/ebm/python')
import pickle
import numpy as np
from netCDF4 import Dataset

ens = 'r1i1p1'
seas = 'djf'
mean = '.ymonmean-30.zonmean'
clim = 'historical'
date = '186001-200512'
model = 'MPI-ESM-LR'
vlist = ['rsdt', 'rsut', 'rsds', 'rsus', 'rlds', 'rlus', 'hfls', 'hfss', 'tend']

inpath = '/project2/tas1/miyawaki/projects/003/data/raw/historical/MPI-ESM-LR'
outpath = '/project2/tas1/miyawaki/projects/003/ebm/python/sims/%s/in' % (clim)

f = {}
for var in vlist:
    fname = '%s_Amon_%s_%s_%s_%s%s.nc' % (var, model, clim, ens, date, mean)
    file = Dataset('%s/%s' % (inpath, fname), 'r')
    f[var] = file[var][:]

    if seas == 'djf':
        f[var] = np.squeeze( np.mean( np.roll(f[var], 1, axis=0)[0:3,...], axis=0) )
    elif seas == 'ann':
        f[var] = np.squeeze( np.mean( f[var], axis=0) )

# save input data
pickle.dump(f, open('%s/in_data.%s.pickle' % (outpath, seas), 'wb'))
