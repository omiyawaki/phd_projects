import os,sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

mean='.yearmean'
varn='r1'
plotdir='./plot/%s' % ('32xcomp')
plotn='%s/%s%s.qflux.pdf' % (plotdir, varn, mean)
lamin=-75
lamax=-55

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

runqpi='rp000200'
trangeqpi='0040_0057'
datadirqpi='/project2/tas1/miyawaki/projects/003/data/raw/echam/%s' % runqpi
filenqpi='%s/%s_%s_%s%s.nc' % (datadirqpi, varn, runqpi, trangeqpi, mean)
dsqpi=xr.open_dataset(filenqpi)
# take SO mean
dsqpi=dsqpi.sel(lat=slice(lamax,lamin))
weighr1qpi=np.cos(np.deg2rad(dsqpi.lat))
r1qpi=dsqpi.r1.weighted(weighr1qpi)
r1qpi=r1qpi.mean(('lat'))
timeqpi=np.arange(1,len(r1qpi.data)+1)

runq32='rp000196'
trangeq32='0040_0065'
datadirq32='/project2/tas1/miyawaki/projects/003/data/raw/echam/%s' % runq32
filenq32='%s/%s_%s_%s%s.nc' % (datadirq32, varn, runq32, trangeq32, mean)
dsq32=xr.open_dataset(filenq32)
# take global mean
dsq32=dsq32.sel(lat=slice(lamax,lamin))
weighr1q32=np.cos(np.deg2rad(dsq32.lat))
r1q32=dsq32.r1.weighted(weighr1q32)
r1q32=r1q32.mean(('lat'))
timeq32=np.arange(1,len(r1q32.data)+1)

# MPI
dsMPI=xr.open_dataset('/project2/tas1/miyawaki/projects/003/data/raw/longrun/MPIESM12_abrupt32x/r1_mon_MPIESM12_abrupt32x_103.yearmean.nc')
# take global mean
dsMPI=dsMPI.sel(lat=slice(lamin,lamax))
weighr1MPI=np.cos(np.deg2rad(dsMPI.lat))
r1MPI=dsMPI.r1.weighted(weighr1MPI)
r1MPI=r1MPI.mean(('lat'))
timeMPI=np.arange(1,len(dsMPI.r1.data)+1)

plt.figure(figsize=(4,3))
plt.axhline(0.9,linewidth=0.5)
plt.plot(timeqpi, r1qpi,label='$Q_{pi}$ AQUA')
plt.plot(timeq32, r1q32,label='$Q_{32}$ AQUA')
plt.plot(timeMPI, r1MPI,label='MPI-ESM')
plt.xlabel('Time (yr)')
plt.ylabel(r'Southern Ocean $R_1$ (unitless)')
plt.title('abrupt32xCO2')
plt.legend()
plt.tight_layout()
plt.savefig(plotn, format='pdf', dpi=300)
