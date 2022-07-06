import sys
import pickle
from tqdm import tqdm
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (4, 3)

forcings = np.arange(0,11,1)

plotpath = sys.argv[1]

c = 0
# load data
for f in tqdm(range(len(forcings))):
    forcing = forcings[f]
    djf=pickle.load(open('%s/frad_%g/djfmean.pickle' % (plotpath, forcing), 'rb'))

    if c == 0:
        x = djf['x']
        cos = (1-x**2)**(1/2)
        Lf = djf['Lf']

        # select poleward of 80N
        # extra = np.where(x < np.sin(np.deg2rad(80)))
        a = np.where(x >= np.sin(np.deg2rad(80)))
        xa = np.take(x, a[0])
        ca = np.take(cos, a[0])

        # initialize arrays
        Ea = np.empty([len(forcings)])
        Ta = np.empty([len(forcings)])
        DTa = np.empty([len(forcings)])
        Sa = np.empty([len(forcings)])

        c = c+1

    # compute Arctic mean (poleward of 80N)
    # select second to last season (last season only has Dec)
    Esel = np.take(djf['E120'][:,-2], a[0], axis=0)
    Tsel = np.take(djf['T120'][:,-2], a[0], axis=0)
    DTsel = np.take(djf['DT120'][:,-2], a[0], axis=0)

    ##compute seasonal ice fraction
    sic = np.zeros(len(xa))
    sic[Esel<0] = 1

    Ea[f] = np.sum(ca * Esel, axis=0) / np.sum(ca)
    Ta[f] = np.sum(ca * Tsel, axis=0) / np.sum(ca)
    DTa[f] = np.sum(ca * DTsel, axis=0) / np.sum(ca)
    Sa[f] = np.sum(ca * sic, axis=0) / np.sum(ca)

# ------------------------------------------------------------------------
# DJF plots ------------------------
# ------------------------------------------------------------------------

# temperature
plt.figure()
plt.plot(forcings,Ta,'k-')
plt.xlabel('F (W m$^{-2}$)')
plt.ylabel(r'T ($^\circ$C)')
plt.tight_layout()
plt.savefig('%s/agg/arctic_t_djf.pdf' % (plotpath), format='pdf', dpi=300)

# enthalpy
plt.figure()
plt.plot(forcings,Ea,'k-')
plt.xlabel('F (W m$^{-2}$)')
plt.ylabel(r'E (J m$^{-2}$)')
plt.tight_layout()
plt.savefig('%s/agg/arctic_e_djf.pdf' % (plotpath), format='pdf', dpi=300)

# ice thickness
plt.figure()
Ha = -Ea/Lf*(Ea<0)
plt.plot(forcings,Ha,'k-')
plt.xlabel('F (W m$^{-2}$)')
plt.ylabel('h (m)')
plt.tight_layout()
plt.savefig('%s/agg/arctic_h_djf.pdf' % (plotpath), format='pdf', dpi=300)

# ice fraction
plt.figure()
plt.plot(forcings,Sa,'k-')
plt.xlabel('F (W m$^{-2}$)')
plt.ylabel('Sea ice fraction (unitless)')
plt.tight_layout()
plt.savefig('%s/agg/arctic_s_djf.pdf' % (plotpath), format='pdf', dpi=300)

# diffused heat flux convergence
plt.figure()
plt.plot(forcings,DTa,'k-')
plt.xlabel('F (W m$^{-2}$)')
plt.ylabel(r'$D\nabla^2 T$ (W m$^{-2}$)')
plt.tight_layout()
plt.savefig('%s/agg/arctic_dt_djf.pdf' % (plotpath), format='pdf', dpi=300)

