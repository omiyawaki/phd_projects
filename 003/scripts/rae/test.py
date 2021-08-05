import sys
from model import model
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

plev = np.linspace(1e5,1e3,50)
ps = 1000e2
tau0 = 3
Fs = 30
Fa = 150
b = 1
beta = 0
n = 2

tau0s = np.linspace(1.5,10,50)

invstr = np.empty_like(tau0s)
OLR = np.empty_like(tau0s)
for i in tqdm(range(len(tau0s))):
    _, _, invstr[i], OLR[i] = model(tau0s[i], b, beta, Fs, Fa, plev, ps, n)
    
print(OLR)
sys.exit()

fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(tau0s, invstr)
ax.set_xlabel('$\\tau_0$ (unitless)')
ax.set_ylabel('$T(\\tau_0) - T_s$')
plt.show()
plt.close()

# fig, ax = plt.subplots()
# ax.plot(T, plev*1e-2)
# ax.plot(Ts, ps*1e-2, '.')
# ax.set_ylim(ax.get_ylim()[::-1]) # invert r1 axis
# plt.show()
# plt.close()