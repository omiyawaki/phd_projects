import os,sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

mean='.yearmean'
varn='tsurf'
plotdir='./plot/%s' % ('32xcomp')
plotn='%s/%s%s.pdf' % (plotdir, varn, mean)

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

run40='rp000193'
trange40='0040_0051'
datadir40='/project2/tas1/miyawaki/projects/003/data/raw/echam/%s' % run40
filen40='%s/%s_%s_%s%s.nc' % (datadir40, varn, run40, trange40, mean)
ds40=xr.open_dataset(filen40)
# take global mean
weights40=np.cos(np.deg2rad(ds40.lat))
ts40=ds40.tsurf.weighted(weights40)
ts40=ts40.mean(('lon','lat'))
time40=np.arange(1,len(ds40.tsurf.data)+1)

run80='rp000196'
trange80='0040_0065'
datadir80='/project2/tas1/miyawaki/projects/003/data/raw/echam/%s' % run80
filen80='%s/%s_%s_%s%s.nc' % (datadir80, varn, run80, trange80, mean)
ds80=xr.open_dataset(filen80)
# take global mean
weights80=np.cos(np.deg2rad(ds80.lat))
ts80=ds80.tsurf.weighted(weights80)
ts80=ts80.mean(('lon','lat'))
time80=np.arange(1,len(ds80.tsurf.data)+1)

# MPI
dsMPI=xr.open_dataset('/project2/tas1/miyawaki/projects/003/data/raw/longrun/MPIESM12_abrupt32x/ts_mon_MPIESM12_abrupt32x_103.yearmean.nc')
# take global mean
weightsMPI=np.cos(np.deg2rad(dsMPI.lat))
tsMPI=dsMPI.ts.weighted(weightsMPI)
tsMPI=tsMPI.mean(('lon','lat'))
timeMPI=np.arange(1,len(dsMPI.ts.data)+1)

plt.figure(figsize=(4,3))
plt.plot(time40, ts40,label='40 m AQUA')
plt.plot(time80, ts80,label='80 m AQUA')
plt.plot(timeMPI, tsMPI,label='MPI-ESM')
plt.xlabel('Time (yr)')
plt.ylabel(r'Global mean $T_s$ (K)')
plt.title('abrupt32xCO2')
plt.legend()
plt.tight_layout()
plt.savefig(plotn, format='pdf', dpi=300)
