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
from tools import calc_esat, e2q, get_model6idx

# specify directory for output data and plots
datadir = '%s/data' % ( os.getcwd() ) # output data directory
if not os.path.exists(datadir):
    os.mkdir(datadir)

lat=83.3
day=0
alb=0.6
cld=True
# clims=[] # hist, rcp
clims=['hist','hist+ssp585','hist+ssp585_fixq_fixT', 'hist+ssp585_fixq_fixco2', 'hist+ssp585_fixT_fixco2', 'hist+ssp585_fixrh_fixco2', 'hist+ssp585_fixq_fixco2_planck', 'hist+ssp585_fixq_fixco2_lapse','hist+ssp585_fixT_fixco2_devrh'] # hist, hist+ssp585
# clims=['hist','hist+ssp585'] # hist, rcp
# models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5', 'CESM2-WACCM', 'IPSL-CM6A-LR', 'MRI-ESM2-0','MIROC-ES2L','GISS-E2-1-G','GISS-E2-1-H','UKESM1-0-LL']
models = ['mmm']

for model in models:
    for clim in clims:

        if 'hist+ssp585' in clim:
            ntime = 441
            [mmm, grid, indiv] = pickle.load(open('../climlab/input_data/forcing.cmip6.pickle', 'rb'))
            [mmm0, _, indiv0] = pickle.load(open('../climlab/input_data/clima.cmip6.pickle', 'rb'))
            if model == 'mmm':
                gcm_forc = mmm
                gcm_clim = mmm0
            else:
                gcm_forc = {}
                gcm_clim = {}
                imod = get_model6idx(model)
                for varname in indiv:
                    gcm_forc[varname] = indiv[varname][imod,...]
                    gcm_clim[varname] = indiv0[varname][imod,...]

            T = gcm_forc['ta']
            Ts = gcm_forc['ts']
            q = gcm_forc['hus']
            ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/greenhouse_hist+ssp585.nc')
            co2 = 1e-6*ds.CO2.sel(time=slice('1859-12-31', '2300-12-31')).data

            if 'fixq' in clim:
                # q[1:,:] = q[0,:]
                # q[:,:] = np.nanmean(gcm_clim['hus'][-30:,:],axis=0)[None,:]
                q[:,:] = np.nanmean(gcm_clim['hus'][-39:-9,:],axis=0)[None,:]

            if 'fixT' in clim:
                # T[1:,:] = T[0,:]
                # Ts[1:] = Ts[0]
                # T[:,:] = np.nanmean(gcm_clim['ta'][-30:,:],axis=0)[None,:]
                T[:,:] = np.nanmean(gcm_clim['ta'][-39:-9,:],axis=0)[None,:]
                # Ts[:] = np.nanmean(gcm_clim['ts'][-30:])
                Ts[:] = np.nanmean(gcm_clim['ts'][-39:-9])

            if 'fixco2' in clim:
                # co2[1:] = co2[0]
                co2[:] = 1e-6*np.nanmean(ds.CO2.sel(time=slice('1976-01-01', '2005-12-31')).data)

            if 'fixrh' in clim:
                # rh = 1e-2*gcm_forc['hur'][0,:]
                # rh = 1e-2*np.nanmean(gcm_clim['hur'][-30:,:],axis=0)
                rh = 1e-2*np.nanmean(gcm_clim['hur'][-39:-9,:],axis=0)
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

            if 'planck' in clim:
                # tbar = np.nanmean(gcm_clim['ta'][-30:,:],axis=0)
                tbar = np.nanmean(gcm_clim['ta'][-39:-9,:],axis=0)
                dt = T - tbar
                dp = grid['lev'][1:]-grid['lev'][:-1]
                dp = np.append(dp, dp[-1])
                vmean_dt = np.nansum( dp*dt, axis=1, keepdims=True) / np.nansum(dp)
                T = tbar + vmean_dt

                # tsbar = np.nanmean(gcm_clim['ts'][-30:])
                tsbar = np.nanmean(gcm_clim['ts'][-39:-9])
                dts = Ts - tsbar
                Ts = tsbar + np.squeeze(vmean_dt)

            if 'lapse' in clim:
                # tbar = np.nanmean(gcm_clim['ta'][-30:,:],axis=0)
                tbar = np.nanmean(gcm_clim['ta'][-39:-9,:],axis=0)
                dt = T - tbar
                dp = grid['lev'][1:]-grid['lev'][:-1]
                dp = np.append(dp, dp[-1])
                vmean_dt = np.nansum( dp*dt , axis=1, keepdims=True) / np.nansum(dp)
                T = tbar + dt-vmean_dt

                # tsbar = np.nanmean(gcm_clim['ts'][-30:])
                tsbar = np.nanmean(gcm_clim['ts'][-39:-9])
                dts = Ts - tsbar
                Ts = tsbar + dts - np.squeeze(vmean_dt)

        elif clim == 'hist':
            ntime = 156
            [mmm, grid, indiv] = pickle.load(open('../climlab/input_data/clima.cmip6.pickle', 'rb'))

            if model == 'mmm':
                gcm_clim = mmm
            else:
                gcm_clim = {}
                imod = get_model6idx(model)
                for varname in indiv:
                    gcm_clim[varname] = indiv[varname][imod,...]

            T = gcm_clim['ta']
            Ts = gcm_clim['ts'];
            q = gcm_clim['hus']
            ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/greenhouse_hist+ssp585.nc')
            co2 = 1e-6*ds.CO2.sel(time=slice('1859-12-30', '2015-12-31')).data

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
            # print(co2[i])
            # print(T[i,:])
            # print(Ts[i])
            # print(q[i,:])
            # print(plev)
            # sys.exit()
            rad = rrtmg(co2[i], T[i,:], Ts[i], q[i,:], plev, lat=lat, day=day, alb=alb, cld=cld)
            raddiag.append(rad)
            ra[i] = rad.ASR - rad.SW_sfc + rad.LW_sfc - rad.OLR
            racs[i] = rad.ASRclr - rad.SW_sfc_clr + rad.LW_sfc_clr - rad.OLRclr
            print(ra[i])

        print(ra[-1]-ra[0])

        # save data
        pickle.dump([ra, racs, raddiag], open('%s/rad.%s.%s.cmip6.pickle' % (datadir, clim, model), 'wb'))

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

