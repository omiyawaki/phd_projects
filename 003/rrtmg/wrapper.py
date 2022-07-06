import os
import sys
import pickle
import climlab
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from rrtmg import rrtmg
from tqdm import tqdm
from tools import calc_esat, e2q, get_modelidx

# specify directory for output data and plots
datadir = '%s/data' % ( os.getcwd() ) # output data directory
if not os.path.exists(datadir):
    os.mkdir(datadir)

lat=83.3
day=0
alb=0.6
cld=True
# clims=['hist', 'rcp'] # hist, rcp
# clims=['hist', 'rcp', 'rcp_fixq_fixT', 'rcp_fixq_fixco2', 'rcp_fixT_fixco2', 'rcp_fixrh_fixco2'] # hist, rcp
clims=['rcp_fixT_fixco2_devrh'] # hist, rcp
models=['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MPI-ESM-LR']
# models=['CNRM-CM5', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MPI-ESM-LR']
# models=['mmm']

for model in models:
    for clim in clims:

        if 'rcp' in clim:
            ntime = 295
            [mmm, grid, indiv] = pickle.load(open('../climlab/input_data/forcing.pickle', 'rb'))
            [mmm0, _, indiv0] = pickle.load(open('../climlab/input_data/clima.pickle', 'rb'))
            if model == 'mmm':
                gcm_forc = mmm
                gcm_clim = mmm0
            else:
                gcm_forc = {}
                gcm_clim = {}
                imod = get_modelidx(model)
                for varname in indiv:
                    gcm_forc[varname] = indiv[varname][imod,...]
                    gcm_clim[varname] = indiv0[varname][imod,...]

            T = gcm_forc['ta']
            Ts = gcm_forc['ts']
            q = gcm_forc['hus']
            ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/ghg_rcp85_1765-2500_c100203.nc')
            co2 = 1e-6*ds.CO2.sel(time=slice('2006-01-01', '2300-12-31')).data

            if 'fixq' in clim:
                # q[1:,:] = q[0,:]
                q[:,:] = np.nanmean(gcm_clim['hus'][-30:,:],axis=0)[None,:]

            if 'fixT' in clim:
                # T[1:,:] = T[0,:]
                # Ts[1:] = Ts[0]
                T[:,:] = np.nanmean(gcm_clim['ta'][-30:,:],axis=0)[None,:]
                Ts[:] = np.nanmean(gcm_clim['ts'][-30:])

            if 'fixco2' in clim:
                # co2[1:] = co2[0]
                co2[:] = 1e-6*np.nanmean(ds.CO2.sel(time=slice('1976-01-01', '2005-12-31')).data)

            if 'fixrh' in clim:
                # rh = 1e-2*gcm_forc['hur'][0,:]
                rh = 1e-2*np.nanmean(gcm_clim['hur'][-30:,:],axis=0)
                esat = calc_esat(T)
                e = rh[None,:]*esat
                p = grid['lev'][None,:]
                q = e2q(p,e)

            if 'devrh' in clim:
                rh = 1e-2*gcm_forc['hur']
                esat = calc_esat(T)
                e = rh*esat
                p = grid['lev'][None,:]
                q = e2q(p,e)

        elif clim == 'hist':
            ntime = 147
            [mmm, grid, indiv] = pickle.load(open('../climlab/input_data/clima.pickle', 'rb'))

            if model == 'mmm':
                gcm_clim = mmm
            else:
                gcm_clim = {}
                imod = get_modelidx(model)
                for varname in indiv:
                    gcm_clim[varname] = indiv[varname][imod,...]

            T = gcm_clim['ta']
            Ts = gcm_clim['ts'];
            q = gcm_clim['hus']
            ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/ghg_rcp85_1765-2500_c100203.nc')
            co2 = 1e-6*ds.CO2.sel(time=slice('1860-01-01', '2006-12-31')).data

            # [gcm_clim, grid] = pickle.load(open('../climlab/input_data/indata.pickle', 'rb'))
            # T = gcm_clim['ta']; T=T[None,:]
            # Ts = np.array([gcm_clim['ts']]);
            # q = gcm_clim['hus']; q=q[None,:]
            # co2 = np.array([348e-6]);

        plev = 1e-2*grid['lev']
        # compute rrtmg radiation
        ra = np.empty(ntime)
        racs = np.empty(ntime)
        raddiag = []
        for i in tqdm(range(ntime)):
            rad = rrtmg(co2[i], T[i,:], Ts[i], q[i,:], plev, lat=lat, day=day, alb=alb, cld=cld)
            raddiag.append(rad)
            ra[i] = rad.ASR - rad.SW_sfc + rad.LW_sfc - rad.OLR
            racs[i] = rad.ASRclr - rad.SW_sfc_clr + rad.LW_sfc_clr - rad.OLRclr
            print(ra[i])

        print(ra[-1]-ra[0])

        # save data
        pickle.dump([ra, racs, raddiag], open('%s/rad.%s.%s.pickle' % (datadir, clim, model), 'wb'))

        # # #  Check for energy balance
        # print('SCM RTOA = %g W m**-2' % (rad.ASR - rad.OLR) )
        # print('GCM RTOA = %g W m**-2\n' % gcm_clim['ftoacs'] )

        # print('SCM ASR = %g W m**-2' % (rad.ASR) )
        # print('GCM ASR = %g W m**-2' % (gcm_clim['rsdt']-gcm_clim['rsutcs']) )

        # print('SCM OLR = %g W m**-2' % (rad.OLR) )
        # print('GCM OLR = %g W m**-2' % (gcm_clim['rlutcs']) )

        # # total radiative cooling
        # Ra = rad.ASR - rad.SW_sfc + rad.LW_sfc - rad.OLR
        # Rsfc = rad.SW_sfc - rad.LW_sfc
        # print('SCM Ra = %g W m**-2' % (Ra))
        # print('GCM Ra = %g W m**-2\n' % gcm_clim['racs'])
        # print('Rsfc = %g W m**-2' % (Rsfc))

