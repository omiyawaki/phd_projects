import sys
import pickle
from tqdm import tqdm
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (4, 3)

forcings = np.arange(-10,11,1)

plotpath = sys.argv[1]

winter = 91 #time of coldest <T>
summer = 31 #time of warmest <T>
c = 0
# load data
for f in tqdm(range(len(forcings))):
    forcing = forcings[f]
    raw=pickle.load(open('%s/frad_%g/raw_data.pickle' % (plotpath, forcing), 'rb'))

    if c == 0:
        x = raw['x']
        Lf = raw['Lf']

        Efin = np.empty([len(x), len(forcings)])
        Tfin = np.empty([len(x), len(forcings)])

        Ewin = np.empty([len(x), len(forcings)])
        Esum = np.empty([len(x), len(forcings)])

        c = c+1
    #-------------------------------------------------------------------------
    #output only converged, final year
    Efin[:,f] = np.mean(raw['E120'][:,-121:-1], axis=1)
    Tfin[:,f] = np.mean(raw['T120'][:,-121:-1], axis=1)

    Ewin[:,f] = raw['E120'][:,-(120-winter)]
    Esum[:,f] = raw['E120'][:,-(120-summer)]

# ------------------------------------------------------------------------
#WE15, Figure 3: sensitivity to forcing ------------------------
# ------------------------------------------------------------------------
#compute seasonal ice edge
xann = np.zeros(len(forcings))
xwin = np.zeros(len(forcings))
xsum = np.zeros(len(forcings))
#if isempty(find(E<0,1))==0:
for j in range(0,len(forcings)):
    Eann = Efin[:,j]
    if any(Eann<0):
        ice = np.where(Eann<0)[0]
        xann[j] = x[ice[0]];
    else:
        xann[j] = max(x);

    E = Ewin[:,j]
    if any(E<0):
        ice = np.where(E<0)[0]
        xwin[j] = x[ice[0]];
    else:
        xwin[j] = max(x);

    E = Esum[:,j]
    if any(E<0):
        ice = np.where(E<0)[0]
        xsum[j] = x[ice[0]];
    else:
        xsum[j] = max(x);

plt.figure()
#plot temperature (Fig 3a)
clevsT = np.arange(-20.001,35.,5.)
plt.contourf(forcings,x,Tfin,clevsT,cmap='turbo')
plt.colorbar()
#plot ice edge on T
plt.plot(forcings,xwin,'k--')
plt.plot(forcings,xsum,'k--')
#plot T=0 contour (the region between ice edge and T=0 contour is the
plt.xlabel('F (W m$^{-2}$)')
plt.ylabel('$\sin(\phi)$')
plt.title(r'T ($^\circ$C)')
plt.tight_layout()
plt.savefig('%s/agg/f3a.pdf' % (plotpath), format='pdf', dpi=300)

plt.figure()
#plot enthalpy (Fig 3b)
clevsE = np.arange(-150.001,350,50)
plt.contourf(forcings,x,Efin,clevsE,cmap='turbo',extend='both')
plt.colorbar()
#plot ice edge on E
plt.plot(forcings,xwin,'k--')
plt.plot(forcings,xsum,'k--')
plt.xlabel('F (W m$^{-2}$)')
plt.ylabel('$\sin(\phi)$')
plt.title(r'E (J m$^{-2}$)')
plt.tight_layout()
plt.savefig('%s/agg/f3b.pdf' % (plotpath), format='pdf', dpi=300)

plt.figure()
#plot the ice thickness (Fig 2c)
clevsh = [1e-5,0.4,0.8,2,3,4,5,6]
hfin = -Efin/Lf*(Efin<0)
plt.contourf(forcings,x,hfin,clevsh,extend='max')
plt.colorbar()
#plot ice edge on h
plt.plot(forcings,xwin,'k--')
plt.plot(forcings,xsum,'k--')
plt.xlabel('F (W m$^{-2}$)')
plt.ylabel('$\sin(\phi)$')
plt.title('h (m)')
plt.tight_layout()
plt.savefig('%s/agg/f3c.pdf' % (plotpath), format='pdf', dpi=300)
