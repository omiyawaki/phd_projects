import numpy as np
from netCDF4 import Dataset

path_t1 = '/glade/scratch/miyawaki/archive/arctic_test.FSCAM.padv/atm/hist/arctic_test.FSCAM.padv.cam.h0.0001-01.nc'
path_t2 = '/glade/scratch/miyawaki/archive/arctic_test.FSCAM.padv/atm/hist/arctic_test.FSCAM.padv.cam.h0.0002-01.nc'

file_t1 = Dataset(path_t1, 'r')
file_t2 = Dataset(path_t2, 'r')

ta1 = file_t1.variables['T'][:] # (lev)
ta2 = file_t2.variables['T'][:] # (lev)

print(ta1)
print(ta2)
