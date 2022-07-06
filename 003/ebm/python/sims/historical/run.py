import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/ebm/python')
import pickle
from src.ebm_orig import ebm_orig
import numpy as np
from netCDF4 import Dataset

ens = 'r1i1p1'
seas = 'djf'
mean = '.ymonmean-30.zonmean'
clim = 'historical'
date = '186001-200512'
model = 'MPI-ESM-LR'

inpath = '/project2/tas1/miyawaki/projects/003/ebm/python/sims/%s/in' % (clim)
outpath = '/project2/tas1/miyawaki/projects/003/ebm/python/sims/%s/out' % (clim)

f = pickle.load(open('%s/in_data.%s.pickle' % (inpath, seas), 'rb'))

sw = f['rsdt'] - f['rsut'] # positive heats atmos
sfc = f['rsds'] - f['rsus'] + f['rlds'] - f['rlus'] - f['hfls'] - f['hfss'] # positive cools atmos
stg = f['tend'] # positive cools atmos (because we move it to the RHS)

r = {}
[r['To'], r['lat'], r['qo'], r['mo'], r['trn'], r['lath'], r['nl']] = ebm_orig(sw, sfc=sfc, stg=stg, vb=True)

# save data
pickle.dump(r, open('%s/raw_data.%s.pickle' % (outpath, seas), 'wb'))
