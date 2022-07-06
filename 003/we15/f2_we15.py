import sys
import pickle
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (4, 3)

f_in = float(sys.argv[1])
plotpath = sys.argv[2]

# load data
raw=pickle.load(open('%s/raw_data.pickle' % (plotpath), 'rb'))

#-------------------------------------------------------------------------
#output only converged, final year
x = raw['x']
Lf = raw['Lf']
tfin = np.linspace(0,1,120)
Efin = raw['E120'][:,-121:-1]
Tfin = raw['T120'][:,-121:-1]

# ------------------------------------------------------------------------
#WE15, Figure 2: Default Steady State Climatology ------------------------
# ------------------------------------------------------------------------
winter = 31 #time of coldest <T>
summer = 91 #time of warmest <T>
#compute seasonal ice edge
xi = np.zeros(120)
#if isempty(find(E<0,1))==0:
for j in range(0,len(tfin)):
    E = Efin[:,j]
    if any(E<0):
        ice = np.where(E<0)[0]
        xi[j] = x[ice[0]];
    else:
        xi[j] = max(x);

plt.figure()
#plot enthalpy (Fig 2a)
clevsE = np.append(np.arange(-40,20,20),np.arange(50,350,50))
plt.contourf(tfin,x,Efin,clevsE)
plt.colorbar()
#plot ice edge on E
plt.plot(tfin,xi,'k')
plt.xlabel('t (final year)')
plt.ylabel('$\sin(\phi)$')
plt.title(r'E (J m$^{-2}$)')
plt.tight_layout()
plt.savefig('%s/f2a.pdf' % (plotpath), format='pdf', dpi=300)

plt.figure()
#plot temperature (Fig 2b)
clevsT = np.arange(-30.001,35.,5.)
plt.contourf(tfin,x,Tfin,clevsT)
plt.colorbar()
#plot ice edge on T
plt.plot(tfin,xi,'k')
#plot T=0 contour (the region between ice edge and T=0 contour is the
#region of summer ice surface melt)
plt.contour(tfin,x,Tfin,[-0.001],colors='r',linestyles='-')
plt.xlabel('t (final year)')
plt.ylabel('$\sin(\phi)$')
plt.title(r'T ($^\circ$C)')
plt.tight_layout()
plt.savefig('%s/f2b.pdf' % (plotpath), format='pdf', dpi=300)

plt.figure()
#plot the ice thickness (Fig 2c)
clevsh = np.arange(0.00001,4.5,.5)
hfin = -Efin/Lf*(Efin<0)
plt.contourf(tfin,x,hfin,clevsh)
plt.colorbar()
#plot ice edge on h
# plt.contour(tfin,x,hfin,[0],colors='k')
plt.axvline(tfin[winter], color='k')
plt.axvline(tfin[summer], linestyle='--', color='k')
plt.xlabel('t (final year)')
plt.ylabel('$\sin(\phi)$')
plt.title('h (m)')
plt.tight_layout()
plt.savefig('%s/f2c.pdf' % (plotpath), format='pdf', dpi=300)

plt.figure()
#plot temperature profiles (Fig 2d)
Summer, = plt.plot(x,Tfin[:,summer],'k--',label='summer')
Winter, = plt.plot(x,Tfin[:,winter],'k',label='winter')
plt.axhline(0, color='k')
plt.xlabel('$\sin(\phi)$')
plt.ylabel(r'T ($^\circ$C)')
plt.legend(handles = [Summer,Winter],loc=0)
plt.tight_layout()
plt.savefig('%s/f2d.pdf' % (plotpath), format='pdf', dpi=300)

plt.figure()
#plot ice thickness profiles (Fig 2e)
plt.plot(x,hfin[:,summer],'k--')
plt.plot(x,hfin[:,winter],'k')
plt.axhline(0, color='k')
plt.xlim([0.7,1])
plt.xlabel('$\sin(\phi)$')
plt.ylabel('h (m)')
plt.legend(handles = [Summer,Winter],loc=0)
plt.tight_layout()
plt.savefig('%s/f2e.pdf' % (plotpath), format='pdf', dpi=300)

plt.figure()
#plot seasonal thickness cycle at pole (Fig 2f)
plt.plot(tfin,hfin[-1,:],'k')
plt.xlabel('t (final year)')
plt.ylabel(r'h$_p$ (m)')
plt.ylim([2, 1.1*max(hfin[-1,:])])
plt.tight_layout()
plt.savefig('%s/f2f.pdf' % (plotpath), format='pdf', dpi=300)

plt.figure()
#plot ice edge seasonal cycle (Fig 2g)
xideg = np.degrees(np.arcsin(xi));
plt.plot(tfin,xideg,'k-')
plt.ylim([0,90])
plt.xlabel('t (final year)')
plt.ylabel(r'$\theta_i$ (deg)');
plt.tight_layout()
plt.savefig('%s/f2g.pdf' % (plotpath), format='pdf', dpi=300)
