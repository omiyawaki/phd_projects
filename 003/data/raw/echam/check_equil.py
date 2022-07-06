import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_ftoa = sys.argv[1]
plotname = sys.argv[2]

# open files
file_ftoa = Dataset(path_ftoa, 'r')

# read data
ftoa = np.squeeze(file_ftoa.variables['ftoa'][:]) # (mon x lat x lon)

fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(ftoa, '-k')
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$F_{TOA}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig(plotname, format='pdf', dpi=300)
