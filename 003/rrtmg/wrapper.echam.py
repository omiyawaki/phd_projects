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

lat0=80
lat1=90
# lat=76.7 # cos-weighted average of latitude 70-90 deg
lat=83.3 # " " 80-90 deg
day=0
alb=0.6
cld=True
# runs=['348','rcp'] # hist, rcp
# runs=['348', 'rcp', 'rcp_fixq_fixT', 'rcp_fixq_fixco2', 'rcp_fixT_fixco2', 'rcp_fixrh_fixco2'] # hist, rcp
# runs=['rcp_fixq_fixT', 'rcp_fixq_fixco2', 'rcp_fixT_fixco2', 'rcp_fixrh_fixco2','rcp_fixq_fixco2_planck','rcp_fixq_fixco2_lapse','rcp_fixT_fixco2_devrh'] # hist, rcp
# runs=['rcp_fixT_fixco2_devrh','rcp_fixrh_fixco2']
# runs=['rcp_fixq_fixco2_planck'] # hist, rcp
runs=['348','ssp585','ssp585_fixq_fixT', 'ssp585_fixq_fixco2', 'ssp585_fixT_fixco2', 'ssp585_fixrh_fixco2', 'ssp585_fixq_fixco2_planck', 'ssp585_fixq_fixco2_lapse','ssp585_fixT_fixco2_devrh'] # hist, ssp585
# runs=['ssp585','ssp585_fixq_fixT', 'ssp585_fixq_fixco2', 'ssp585_fixT_fixco2', 'ssp585_fixrh_fixco2', 'ssp585_fixq_fixco2_planck', 'ssp585_fixq_fixco2_lapse','ssp585_fixT_fixco2_devrh'] # hist, ssp585
# runs=['348'] # hist, ssp585
models=['rp000134','rp000190f']
# models=['rp000190f']

for model in models:
    for run in runs:

        if 'ssp' in run:
            ntime = 214
            [clim, _] = pickle.load(open('../climlab/input_data/clima.%s.%g.%g.pickle'%(model,lat0,lat1), 'rb'))
            if model=='rp000134':
                [forc, grid] = pickle.load(open('../climlab/input_data/forcing.%s.%g.%g.pickle'%('rp000188',lat0,lat1), 'rb'))
            elif model=='rp000190f':
                [forc, grid] = pickle.load(open('../climlab/input_data/forcing.%s.%g.%g.pickle'%('rp000191f',lat0,lat1), 'rb'))

            T = forc['t']
            Ts = forc['tsurf']
            q = forc['q']
            ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/greenhouse_hist+ssp585.nc')
            # co2 = 1e-6*ds.CO2.sel(time=slice('1987-01-01', '2088-12-31')).data
            co2 = 1e-6*ds.CO2.sel(time=slice('1987-01-01', '2201-12-31')).data
            # co2 = 1e-6*ds.CO2.sel(time=slice('1987-01-01', '2288-12-31')).data

            if 'fixq' in run:
                # q[1:,:] = q[0,:]
                q[:,:] = np.nanmean(clim['q'][-20:,:],axis=0)[None,:]

            if 'fixT' in run:
                # T[1:,:] = T[0,:]
                # Ts[1:] = Ts[0]
                T[:,:] = np.nanmean(clim['t'][-20:,:],axis=0)[None,:]
                Ts[:] = np.nanmean(clim['tsurf'][-20:])

            if 'fixco2' in run:
                # co2[1:] = co2[0]
                co2[:] = 1e-6*348

            if 'fixrh' in run:
                # rh = 1e-2*forc['hur'][0,:]
                rh = np.nanmean(clim['rhumidity'][-20:,:],axis=0)
                esat = calc_esat(T)
                e = rh[None,:]*esat
                p = grid['lev'][None,:]
                q = e2q(p,e)

            if 'devrh' in run:
                rh = forc['rhumidity']
                esat = calc_esat(T)
                e = rh*esat
                p = grid['lev'][None,:]
                q = e2q(p,e)

            if 'planck' in run:
                tbar = np.nanmean(clim['t'][-20:,:],axis=0)
                dt = T - tbar
                dp = grid['lev'][1:]-grid['lev'][:-1]
                dp = np.append(dp, dp[-1])
                vmean_dt = np.nansum( dp*dt, axis=1, keepdims=True) / np.nansum(dp)
                T = tbar + vmean_dt

                tsbar = np.nanmean(clim['tsurf'][-20:])
                dts = Ts - tsbar
                Ts = tsbar + np.squeeze(vmean_dt)

            if 'lapse' in run:
                tbar = np.nanmean(clim['t'][-20:,:],axis=0)
                dt = T - tbar
                dp = grid['lev'][1:]-grid['lev'][:-1]
                dp = np.append(dp, dp[-1])
                vmean_dt = np.nansum( dp*dt , axis=1, keepdims=True) / np.nansum(dp)
                T = tbar + dt-vmean_dt

                tsbar = np.nanmean(clim['tsurf'][-20:])
                dts = Ts - tsbar
                Ts = tsbar + dts - np.squeeze(vmean_dt)

        elif run == '348':
            ntime = 21
            [clim, grid] = pickle.load(open('../climlab/input_data/clima.%s.%g.%g.pickle'%(model,lat0,lat1), 'rb'))

            T = clim['t']
            Ts = clim['tsurf'];
            q = clim['q']
            ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/greenhouse_hist+ssp585.nc')
            # co2 = 1e-6*ds.CO2.sel(time=slice('1860-01-01', '2006-12-31')).data
            co2=1e-6*348*np.ones(ntime)

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
        pickle.dump([ra, racs, raddiag], open('%s/rad.%s.%s.%g.%g.pickle' % (datadir, run, model,lat0,lat1), 'wb'))

        # # #  Check for energy balance
        # print('SCM RTOA = %g W m**-2' % (rad.ASR - rad.OLR) )
        # print('GCM RTOA = %g W m**-2\n' % clim['ftoacs'] )

        # print('SCM ASR = %g W m**-2' % (rad.ASR) )
        # print('GCM ASR = %g W m**-2' % (clim['rsdt']-clim['rsutcs']) )

        # print('SCM OLR = %g W m**-2' % (rad.OLR) )
        # print('GCM OLR = %g W m**-2' % (clim['rlutcs']) )

        # # total radiative cooling
        # Ra = rad.ASR - rad.SW_sfc + rad.LW_sfc - rad.OLR
        # Rsfc = rad.SW_sfc - rad.LW_sfc
        # print('SCM Ra = %g W m**-2' % (Ra))
        # print('GCM Ra = %g W m**-2\n' % clim['racs'])
        # print('Rsfc = %g W m**-2' % (Rsfc))

