import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (4, 3) 
from netCDF4 import Dataset

ens = 'r1i1p1'
seas = 'ann'
mean = '.ymonmean-30.zonmean'
clim = 'historical'
date = '186001-200512'
model = 'MPI-ESM-LR'
vlist = ['rlut', 'tas']

inpath = '/project2/tas1/miyawaki/projects/003/data/raw/historical/MPI-ESM-LR'

f = {}
for var in vlist:
    fname = '%s_Amon_%s_%s_%s_%s%s.nc' % (var, model, clim, ens, date, mean)
    file = Dataset('%s/%s' % (inpath, fname), 'r')
    f[var] = file[var][:]

    if seas == 'djf':
        f[var] = np.squeeze( np.mean( np.roll(f[var], 1, axis=0)[0:3,...], axis=0) )
    elif seas == 'ann':
        f[var] = np.squeeze( np.mean( f[var], axis=0) )

# regress OLR and Ts
A = np.vstack([f['tas'], np.ones(len(f['tas']))]).T
a, b = np.linalg.lstsq(A, f['rlut'], rcond=None)[0]

print(a)
print(b)

tv = np.linspace(230,300,101)

# temperature                                                                       
plt.figure()                                                                        
plt.plot(f['tas'][:],f['rlut'][:],'k.',label=None)
reg = plt.plot(tv,a*tv+b,'k--', label=r'$OLR=%.2fT%g$' % (a,b))
plt.xlabel(r'$T_s$ (K)')
plt.ylabel(r'OLR (W m$^{-2}$)') 
plt.tight_layout()                                                                  
plt.legend()
plt.savefig('./ab.pdf', format='pdf', dpi=300)

sys.exit()
