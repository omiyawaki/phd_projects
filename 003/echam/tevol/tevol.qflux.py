import os,sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

mean='.yearmean'
varn='tsurf'
plotdir='./plot/%s' % ('32xcomp')
plotn='%s/%s%s.qflux.pdf' % (plotdir, varn, mean)

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

runqpi='rp000200'
trangeqpi='0040_0057'
datadirqpi='/project2/tas1/miyawaki/projects/003/data/raw/echam/%s' % runqpi
filenqpi='%s/%s_%s_%s%s.nc' % (datadirqpi, varn, runqpi, trangeqpi, mean)
dsqpi=xr.open_dataset(filenqpi)
# take global mean
weightsqpi=np.cos(np.deg2rad(dsqpi.lat))
tsqpi=dsqpi.tsurf.weighted(weightsqpi)
tsqpi=tsqpi.mean(('lon','lat'))
timeqpi=np.arange(1,len(dsqpi.tsurf.data)+1)

runq32='rp000196'
trangeq32='0040_0065'
datadirq32='/project2/tas1/miyawaki/projects/003/data/raw/echam/%s' % runq32
filenq32='%s/%s_%s_%s%s.nc' % (datadirq32, varn, runq32, trangeq32, mean)
dsq32=xr.open_dataset(filenq32)
# take global mean
weightsq32=np.cos(np.deg2rad(dsq32.lat))
tsq32=dsq32.tsurf.weighted(weightsq32)
tsq32=tsq32.mean(('lon','lat'))
timeq32=np.arange(1,len(dsq32.tsurf.data)+1)

# MPI
dsMPI=xr.open_dataset('/project2/tas1/miyawaki/projects/003/data/raw/longrun/MPIESM12_abrupt32x/ts_mon_MPIESM12_abrupt32x_103.yearmean.nc')
# take global mean
weightsMPI=np.cos(np.deg2rad(dsMPI.lat))
tsMPI=dsMPI.ts.weighted(weightsMPI)
tsMPI=tsMPI.mean(('lon','lat'))
timeMPI=np.arange(1,len(dsMPI.ts.data)+1)

plt.figure(figsize=(4,3))
plt.plot(timeqpi, tsqpi,label='$Q_{pi}$ m AQUA')
plt.plot(timeq32, tsq32,label='$Q_{32}$ m AQUA')
plt.plot(timeMPI, tsMPI,label='MPI-ESM')
plt.xlabel('Time (yr)')
plt.ylabel(r'Global mean $T_s$ (K)')
plt.title('abrupt32xCO2')
plt.legend()
plt.tight_layout()
plt.savefig(plotn, format='pdf', dpi=300)
