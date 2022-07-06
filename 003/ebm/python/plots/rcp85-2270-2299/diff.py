import sys
import pickle
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (4, 3)

seas = 'djf'

# control 
clim = 'historical'
inpath = '/project2/tas1/miyawaki/projects/003/ebm/python/sims/%s/in' % (clim)
outpath = '/project2/tas1/miyawaki/projects/003/ebm/python/sims/%s/out' % (clim)
hist = pickle.load(open('%s/raw_data.%s.pickle' % (outpath, seas), 'rb'))

# future
clim = 'rcp85-2270-2299'
inpath = '/project2/tas1/miyawaki/projects/003/ebm/python/sims/%s/in' % (clim)
outpath = '/project2/tas1/miyawaki/projects/003/ebm/python/sims/%s/out' % (clim)
fut = pickle.load(open('%s/raw_data.%s.pickle' % (outpath, seas), 'rb'))

# temperature
plt.figure()
plt.plot(hist['lat'],fut['To']-hist['To'],'k-')
plt.xlim([-90,90])
plt.xlabel('latitude (deg)')
plt.xticks(np.arange(-90,91,30))
plt.ylabel(r'$\Delta T$ (K)')
plt.tight_layout()
plt.savefig('./diff_t.%s.pdf' % (seas), format='pdf', dpi=300)

# mse
plt.figure()
plt.plot(hist['lat'],1e-3*(fut['mo']-hist['mo']),'k-')
plt.xlim([-90,90])
plt.xlabel('latitude (deg)')
plt.xticks(np.arange(-90,91,30))
plt.ylabel(r'$\Delta m$ (J g$^{-1}$)')
plt.tight_layout()
plt.savefig('./diff_m.%s.pdf' % (seas), format='pdf', dpi=300)

# mse transport
plt.figure()
plt.axhline(0, linewidth=0.5, color='k')
plt.plot(hist['lath'],1e-15*(fut['trn']-hist['trn']),'k-')
plt.xlim([-90,90])
plt.xlabel('latitude (deg)')
plt.xticks(np.arange(-90,91,30))
plt.ylabel(r'$\Delta F_a$ (PW)')
plt.tight_layout()
plt.savefig('./diff_fa.%s.pdf' % (seas), format='pdf', dpi=300)
