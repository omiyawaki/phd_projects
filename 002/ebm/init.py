from ebm2d import Ebm2D
import climlab
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

def initialize(N, mld):
    # create directories if they doesn't exist
    plotpath = './figures/d%g' % (mld)
    if not os.path.isdir(plotpath):
        os.makedirs(plotpath)
    datapath = './data/d%g' % (mld)
    if not os.path.isdir(datapath):
        os.makedirs(datapath)

    ebm_init = Ebm2D(d_ml = mld, conv_scheme=6.5, sol_state='seasonal')
    D = 0.04 # thermal diffusivity in W/m**2/degC
    ebm_init.add_diff(D)
    Cd = 3e-3
    ebm_init.add_stf(Cd)

    # steady state test
    t_end = N # integrate for 5 years
    dt = 1/12 # one month steps
    Nt = int(t_end/dt)
    t = np.linspace(0, t_end, Nt)

    imbalance = np.empty([Nt, 1])

    plt.figure()
    for i in np.arange(0,Nt):
        ebm_init.integrate(dt) # integrate for 1 year
        imbalance[i] = ebm_init.comp_toa_imbalance()

        plt.cla()
        plt.plot(t[:i], imbalance[:i])
        plt.xlabel('t (yr)')
        plt.ylabel('$F_{TOA}$ (Wm$^{-2}$)')
        plt.savefig('%s/imbalance.png' % (plotpath))
        plt.close()

    t_yr = np.arange(0.5,0.5+t_end)
    imbalance_yr = np.empty([t_end, 1])
    for yr in np.arange(0,t_end):
        idx = np.arange(0,12) + 12*yr
        imbalance_yr[yr] = imbalance[idx].sum()/12

    plt.figure()
    plt.plot(t, imbalance)
    plt.plot(t_yr, imbalance_yr)
    plt.xlabel('t (yr)')
    plt.ylabel('$F_{TOA}$ (Wm$^{-2}$)')
    plt.savefig('%s/imbalance.png' % (plotpath))
    plt.close()

    outfile = open('%s/ebm_init.pickle' % (datapath), 'wb')
    pickle.dump(ebm_init, outfile)
    outfile.close()

# initialize EBM for N years and save the file
N = 5
mld = 5 # mixed layer depth
initialize(N, mld)
