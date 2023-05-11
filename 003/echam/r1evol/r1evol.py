import os,sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

mean='.yearmean'
varn='r1'
plotdir='./plot/%s' % ('32xcomp')
plotn='%s/%s%s.pdf' % (plotdir, varn, mean)
lamin=-75
lamax=-55

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

run40='rp000193'
trange40='0040_0051'
datadir40='/project2/tas1/miyawaki/projects/003/data/raw/echam/%s' % run40
filen40='%s/%s_%s_%s%s.nc' % (datadir40, varn, run40, trange40, mean)
ds40=xr.open_dataset(filen40)
# take SO mean
ds40=ds40.sel(lat=slice(lamax,lamin))
weighr140=np.cos(np.deg2rad(ds40.lat))
r140=ds40.r1.weighted(weighr140)
r140=r140.mean(('lat'))
time40=np.arange(1,len(r140.data)+1)

run80='rp000196'
trange80='0040_0065'
datadir80='/project2/tas1/miyawaki/projects/003/data/raw/echam/%s' % run80
filen80='%s/%s_%s_%s%s.nc' % (datadir80, varn, run80, trange80, mean)
ds80=xr.open_dataset(filen80)
# take global mean
ds80=ds80.sel(lat=slice(lamax,lamin))
weighr180=np.cos(np.deg2rad(ds80.lat))
r180=ds80.r1.weighted(weighr180)
r180=r180.mean(('lat'))
time80=np.arange(1,len(r180.data)+1)

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
plt.plot(time40, r140,label='40 m AQUA')
plt.plot(time80, r180,label='80 m AQUA')
plt.plot(timeMPI, r1MPI,label='MPI-ESM')
plt.xlabel('Time (yr)')
plt.ylabel(r'Southern Ocean $R_1$ (unitless)')
plt.title('abrupt32xCO2')
plt.legend()
plt.tight_layout()
plt.savefig(plotn, format='pdf', dpi=300)
