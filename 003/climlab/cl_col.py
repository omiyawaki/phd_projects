import os
import sys
import pickle
import climlab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# specify directory for output data and plots
datadir = '%s/data' % ( os.getcwd() ) # output data directory
if not os.path.exists(datadir):
    os.mkdir(datadir)

lat=83.3
day=0
# co2_lev=280e-6 # volume mixing ratio
co2_lev=348e-6 # volume mixing ratio
# co2_lev=4*280e-6 # volume mixing ratio
n_yr = 1 # number of years to run

par={}
par['n']=2
par['b']=0.5
par['g']=9.81 # m s**-2
par['cp']=1004 # J kg**-1 K**-1
par['ps']=1000 # hPa
par['d']=10 # mixed layer depth [m]
par['rhow']=1e3 # density of water [kg/m^3]
par['cw']=4.182e3 # specific heat capacity of water [J/kg/K]

# load GCM data
# [cl_hl, flux_hl, grid] = pickle.load(open('./input_data/cl_in.zonmean.ymonmean-30.pickle', 'rb'))
# take DJF mean
# gcm_clim = {} # gcm climatology data
# gcm_clim['ta_djf'] = np.nanmean(np.roll(cl_hl['ta'], 1, axis=0)[0:3], 0)
# gcm_clim['hur_djf'] = np.nanmean(np.roll(cl_hl['hur'], 1, axis=0)[0:3], 0)
# gcm_clim['hus_djf'] = np.nanmean(np.roll(cl_hl['hus'], 1, axis=0)[0:3], 0)
# gcm_clim['adv_djf'] = np.nanmean(np.roll(flux_hl['stg_adv'], 1, axis=0)[0:3], 0)
# Fadv = -gcm_clim['adv_djf'] # W m **-2
# gcm_clim['plev'] = 1e-2*grid['lev']

[gcm_clim, grid] = pickle.load(open('./input_data/indata.pickle', 'rb'))
Fadv = -gcm_clim['adv']
gcm_clim['plev'] = 1e-2*grid['lev']

# Fadv=130 # W m **-2


alb = 0.6
Q = climlab.solar.insolation.daily_insolation( lat, day )

#  State variables (Air and surface temperature)
state = climlab.column_state(lev=gcm_clim['plev'], water_depth=par['d'])

# #  Initialize a nearly dry column (small background stratospheric humidity)
q = np.ones_like(state.Tatm) * 5.E-6
# q = gcm_clim['hus_djf']
#  Add specific_humidity to the state dictionary
state['q'] = q

#  Parent model process
rcm = climlab.TimeDependentProcess(state=state)

# prescribe vertical structure of advective heat flux (expressed in heating tendency, units K/s)
Fadv_lev = par['n']*par['b']*par['g']/par['cp'] * Fadv / (1e2*par['ps']) * (rcm.lev / par['ps'])**(par['n']*par['b']-1)

# atmospheric storage tendency (units K/s)
Fatm = gcm_clim['tendv']/par['cp']

# surface storage/heat flux conv tendency (units K/s)
Fsfc = gcm_clim['fsfc']/(par['rhow']*par['cw']*par['d'])

# #  Fixed relative humidity
# h2o = climlab.radiation.ManabeWaterVapor(state=state)
h2o = climlab.radiation.water_vapor.FixedRelativeHumidity(state=state)
h2o.RH_profile = 1e-2*gcm_clim['hur']

#  RADIATION
abs_vmr = climlab.radiation.radiation.default_absorbers(rcm.Tatm)
abs_vmr['CO2']=co2_lev
#  Daily insolation as a function of latitude and time of year
# rad = climlab.radiation.RRTMG(state=state, specific_humidity=h2o.q, albedo=alb, absorber_vmr=abs_vmr)
rad = climlab.radiation.RRTMG(state=state, specific_humidity=h2o.q, albedo=alb, insolation=Q, absorber_vmr=abs_vmr)
# rad = climlab.radiation.RRTMG(state=state, specific_humidity=gcm_clim['hus'], albedo=alb, insolation=Q, absorber_vmr=abs_vmr)

#  Convective adjustment
conv = climlab.convection.ConvectiveAdjustment(state=state, adj_lapse_rate='MALR')

#  LH and SH
shf = climlab.surface.SensibleHeatFlux(state=state, Cd=3e-3)
lhf = climlab.surface.LatentHeatFlux(state=state, Cd=3e-5)

# ADVECTIVE HEATING
adv = climlab.process.ExternalForcing(state=state)
adv.forcing_tendencies['Tatm'] = Fadv_lev + Fatm
adv.forcing_tendencies['Ts'] = Fsfc

#  Couple everything together
rcm.add_subprocess('Radiation', rad)
rcm.add_subprocess('WaterVapor', h2o)
rcm.add_subprocess('Convection', conv)
rcm.add_subprocess('SHF', shf)
rcm.add_subprocess('LHF', lhf)
rcm.add_subprocess('Advection', adv)

#  Run the model
rcm.integrate_years(n_yr)
print(rcm.Tatm[-1])
print('SCM Ts = %g K' % rcm.Ts)
print('GCM Ts = %g K\n' % gcm_clim['tas'])
# for i in range(10):
#     rcm.step_forward()
#     print(rcm.Tatm[-1])
#     print(rcm.Ts)
# #  Check for energy balance
print('SCM RTOA = %g W m**-2' % (rcm.ASR - rcm.OLR) )
print('GCM RTOA = %g W m**-2\n' % gcm_clim['ftoa'] )

# total radiative cooling
Ra = rcm.ASR - rcm.SW_sfc + rcm.LW_sfc - rcm.OLR
Rsfc = rcm.SW_sfc - rcm.LW_sfc
print('SCM Ra = %g W m**-2' % (Ra))
print('GCM Ra = %g W m**-2\n' % gcm_clim['ra'])
print('SCM SH = %g W m**-2' % (rcm.SHF))
print('GCM SH = %g W m**-2\n' % gcm_clim['hfss'])
print('SCM LH = %g W m**-2' % (rcm.LHF))
print('GCM LH = %g W m**-2\n' % gcm_clim['hfls'])
print('Rsfc = %g W m**-2' % (Rsfc))
print('Fa = %g W m**-2' % (Fadv))

# save data
pickle.dump([rcm, Ra, Rsfc, Fadv_lev, gcm_clim], open('%s/rcm.%g.%g.pickle' % (datadir, 1e6*co2_lev, n_yr), 'wb'))
