# Till Wagner's python version:
#
# Reference: "How Model Complexity Influences Sea Ice Stability",
# T.J.W. Wagner & I. Eisenman, J Clim 28,10 (2015)
#
# This script numerically solves the system described by eqns (2), (8) and
# (9) in the article above (henceforth WE15).
# For computational convenience a diffusive 'ghost layer' is invoked, as
# described in WE15, Appendix A. This allows us to solve the system using
# an Implicit Euler scheme on the ghost layer (which effectively solves the
# diffusion equation) and a Forward Euler scheme for the evolution of the
# surface enthalpy. For further detailed documentation to go with this
# script, see WE15_NumericIntegration.pdf.
#
# The present script uses the default parameter values of WE15 (Table 1 and
# Section 2d) with no climate forcing (F=0) and the final part of the code
# produces the plot corresponding essentially to Figure 2 in WE15.
#
# The default configuration here runs a simulation for 30 years at 1000
# timesteps/year and a spatial resolution of 100 gridboxes, equally spaced
# between equator and pole.
#
# The computational time of this code is comparable to that of the
# corresponding Matlab code (without having performed proper performance checks)
#
# Till Wagner (tjwagner@ucsd.edu) & Ian Eisenman (eisenman@ucsd.edu), 
# created Oct 2015, minor bug fix Jan 2022 (in eq.A1, S[:,i] -> S[:,i+1]).
#
#--------------------------------------------------------------------------
import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

f_in = float(sys.argv[1])
plotpath = sys.argv[2]

print(f_in)
print(plotpath)

D = 0.6 # diffusivity for heat transport (W m^-2 K^-1)
A = 193 # OLR when T = 0 (W m^-2)
B = 2.1 # OLR temperature dependence (W m^-2 K^-1)
cw = 9.8 # ocean mixed layer heat capacity (W yr m^-2 K^-1)
S0 = 420 # insolation at equator (W m^-2)
S1 = 338; # insolation seasonal dependence (W m^-2)
S2 = 240 # insolation spatial dependence (W m^-2)
a0 = 0.7 # ice-free co-albedo at equator
a2 = 0.1 # ice=free co-albedo spatial dependence
ai = 0.4 # co-albedo where there is sea ice
Fb = 4; # heat flux from ocean below (W m^-2)
k = 2; # sea ice thermal conductivity (W m^-2 K^-1)
Lf = 9.5; # sea ice latent heat of fusion (W yr m^-3)
F = f_in; # forcing (W m**-2)
cg = 0.01*cw; # ghost layer heat capacity(W yr m^-2 K^-1)
tau = 1e-5; # ghost layer coupling timescale (yr)
##The default run in WE15, Fig 2 uses the time-stepping parameters: -------
n=400; # of evenly spaced latitudinal gridboxes (equator to pole)
# nt=1000; # of timesteps per year (approx lower limit of stability)
dur=200; # of years for the whole run
##For a quicker computation, use the parameters: --------------------------
# n = 100;
nt = 1200; #[OM: changed 1000 -> 1200 for cleaner monthly means]
# dur= 30;
# dur= 200;
dt = 1/nt;
#Spatial Grid -------------------------------------------------------------
dx = 1.0/n #grid box width
x = np.arange(dx/2,1+dx/2,dx) #native grid
xb = np.arange(dx,1,dx)
##Diffusion Operator (WE15, Appendix A) -----------------------------------
lam = D/dx**2*(1-xb**2)
L1=np.append(0, -lam); L2=np.append(-lam, 0); L3=-L1-L2
diffop = - np.diag(L3) - np.diag(L2[:n-1],1) - np.diag(L1[1:n],-1);
##Definitions for implicit scheme on Tg
cg_tau = cg/tau;
dt_tau = dt/tau;
dc = dt_tau*cg_tau;
kappa = (1+dt_tau)*np.identity(n)-dt*diffop/cg;
##Seasonal forcing (WE15 eq.3)
ty = np.arange(dt/2,1+dt/2,dt)
S = (np.tile(S0-S2*x**2,[nt,1])-np.tile(S1*np.cos(2*np.pi*ty),[n,1]).T*np.tile(x,[nt,1]));
S = np.vstack((S,S[0,:]));

##Further definitions
M = B+cg_tau;
aw = a0-a2*x**2 # open water albedo
kLf = k*Lf;

# collect all output into this dictionary
dt120 = dt*10
td = np.arange(dt120/2,dur+dt120/2,dt120) # diagnostic time
raw={}
raw['Lf'] = Lf; raw['Fb'] = Fb; # constants/parameters
raw['x'] = x; raw['td'] = td; # independent variables
raw['varnames'] = ['E120', 'T120', 'R120', 'DT120'] # list of dependent variable names
#Set up output arrays, saving 120 timesteps/year
raw['E120'] = np.zeros([n,dur*120]); raw['T120'] = np.zeros([n,dur*120])
# [OM: initialize additional diagnostics]
raw['R120'] = np.zeros([n,dur*120]); raw['DT120'] = np.zeros([n,dur*120]);
p = -1; m = -1

#Initial conditions ------------------------------------------------------
T = 7.5+20*(1-2*x**2);
Tg = T; E = cw*T;
# [OM: solve initial alpha for R120]
alpha = aw*(E>0) + ai*(E<0)

#Integration (see WE15_NumericIntegration.pdf)----------------------------
#Loop over Years ---------------------------------------------------------
for years in tqdm(range(0,dur)):
    #Loop within One Year-------------------------------------------------
    for i in range(0,int(nt)):
        m = m+1

        # DIAGNOSTICS
        # store output data 120 timesteps per year
        if (p+1)*10 == m:
            p = p+1
            raw['E120'][:,p] = E
            raw['T120'][:,p] = T

            # [OM: ADDITIONAL DIAGNOSTICS]
            # [OM: compute D \nabla^2 T to diagnose diffusive heat flux convergence]
            raw['DT120'][:,p] = np.dot(diffop,Tg)
            # [OM: TOA net Radiation, R, diagnostic]
            raw['R120'][:,p] = alpha*S[i,:] - A - B*T 

        # FORWARD STEPPING
        #forcing
        alpha = aw*(E>0) + ai*(E<0) #WE15, eq.4
        C = alpha*S[i,:]+cg_tau*Tg-A
        #surface temperature
        T0 = C/(M-kLf/E) #WE15, eq.A3
        T = E/cw*(E>=0)+T0*(E<0)*(T0<0); #WE15, eq.9
        #Forward Euler on E
        E = E+dt*(C-M*T+Fb+F); #WE15, eq.A2 [OM: added +F]
        #Implicit Euler on Tg
        Tg = np.linalg.solve(kappa-np.diag(dc/(M-kLf/E)*(T0<0)*(E<0)),
                Tg+(dt_tau*(E/cw*(E>=0)+(ai*S[i+1,:]-A)/(M-kLf/E)*(T0<0)*(E<0)))) #WE15, eq.A1

# save data
pickle.dump(raw, open('%s/raw_data.pickle' % (plotpath), 'wb'))
