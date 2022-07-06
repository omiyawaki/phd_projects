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
yearmean=pickle.load(open('%s/yearmean.pickle' % (plotpath), 'rb'))

#-------------------------------------------------------------------------
x = yearmean['x']
rcos = (1-x**2)**(1/2)

F = f_in
Fb = raw['Fb']
R = yearmean['R120']
DT = yearmean['DT120']
N = R + DT + Fb + F

# compute global mean net flux of the atmosphere + surface system
Ng = np.sum(rcos[:,None] * N, axis=0) / np.sum(rcos)

# ------------------------------------------------------------------------
# Plot trajectory toward equilibrium ------------------------
# ------------------------------------------------------------------------

plt.figure()
plt.plot(Ng, 'k')
plt.axhline(0, linewidth=0.5, color='k')
plt.xlabel('t (year)')
plt.ylabel('$\partial E_g / \partial t$ (W m$^{-2}$)')
plt.tight_layout()
plt.savefig('%s/equil.pdf' % (plotpath), format='pdf', dpi=300)
