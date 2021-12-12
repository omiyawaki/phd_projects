import os
import sys
import pickle
import climlab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

lat=85
day=0
# co2_lev=280e-6 # volume mixing ratio
co2_lev=355e-6 # volume mixing ratio
# co2_lev=4*280e-6 # volume mixing ratio
try_load=0

par={}
par['n']=2
par['b']=1
par['g']=9.81 # m s**-2
par['cp']=1004 # J kg**-1 K**-1
par['ps']=1000 # hPa

# load GCM data
[cl_hl, flux_hl, grid] = pickle.load(open('./input_data/cl_in.zonmean.ymonmean-30.pickle', 'rb'))

# take DJF mean
ta_djf = np.nanmean(np.roll(cl_hl['ta'], 1, axis=0)[0:3], 0)
hur_djf = np.nanmean(np.roll(cl_hl['hur'], 1, axis=0)[0:3], 0)
hus_djf = np.nanmean(np.roll(cl_hl['hus'], 1, axis=0)[0:3], 0)
adv_djf = np.nanmean(np.roll(flux_hl['stg_adv'], 1, axis=0)[0:3], 0)

# Fadv=130 # W m **-2
Fadv = -adv_djf # W m **-2

plev_gcm = 1e-2*grid['lev']

alb = 0.5

def run_model(bool_sh, co2):
    Q = climlab.solar.insolation.daily_insolation( lat, day )

    #  State variables (Air and surface temperature)
    state = climlab.column_state(lev=plev_gcm, water_depth=1)

    # #  Initialize a nearly dry column (small background stratospheric humidity)
    # q = np.ones_like(state.Tatm) * 5.E-6
    #  Add specific_humidity to the state dictionary
    # state['q'] = q

    #  Parent model process
    rcm = climlab.TimeDependentProcess(state=state)

    # #  Fixed relative humidity
    # h2o = climlab.radiation.ManabeWaterVapor(state=state)
    h2o = climlab.radiation.water_vapor.FixedRelativeHumidity(state=state)
    h2o.RH_profile = 1e-2*hur_djf

    #  RADIATION
    abs_vmr = climlab.radiation.radiation.default_absorbers(rcm.Tatm)
    abs_vmr['CO2']=co2
    #  Daily insolation as a function of latitude and time of year
    rad = climlab.radiation.RRTMG(state=state, specific_humidity=h2o.q, albedo=alb, insolation=Q, absorber_vmr=abs_vmr, ozone_file='apeozone_cam3_5_54.nc')
    # rad = climlab.radiation.RRTMG(state=state, specific_humidity=hus_djf, albedo=alb, insolation=Q, absorber_vmr=abs_vmr, ozone_file='apeozone_cam3_5_54.nc')

    #  Convective adjustment
    conv = climlab.convection.ConvectiveAdjustment(state=state, adj_lapse_rate='MALR')

    #  LH and SH
    if bool_sh:
        shf = climlab.surface.SensibleHeatFlux(state=state, Cd=0.5e-3)
    # lhf = climlab.surface.LatentHeatFlux(state=state, Cd=0.5e-3)

    # ADVECTIVE HEATING
    adv = climlab.process.ExternalForcing(state=state)
    adv.forcing_tendencies['Tatm'] = par['n']*par['b']*par['g']/par['cp'] * Fadv / (1e2*par['ps']) * (rcm.lev / par['ps'])**(par['n']*par['b']-1)

    #  Couple everything together
    rcm.add_subprocess('Radiation', rad)
    rcm.add_subprocess('WaterVapor', h2o)
    rcm.add_subprocess('Convection', conv)
    if bool_sh:
        rcm.add_subprocess('SHF', shf)
    # rcm.add_subprocess('LHF', lhf)
    rcm.add_subprocess('Advection', adv)

    #  Run the model
    rcm.integrate_years(1)
    #  Check for energy balance
    print(rcm.ASR - rcm.OLR)

    # total radiative cooling
    Ra = rcm.ASR - rcm.SW_sfc + rcm.LW_sfc - rcm.OLR
    Rsfc = rcm.SW_sfc - rcm.LW_sfc
    print('Ra = %g W m**-2' % (Ra))
    if bool_sh:
        print('SH = %g W m**-2' % (rcm.SHF))
    # print('LH = %g W m**-2' % (rcm.LHF))
    print('Rsfc = %g W m**-2' % (Rsfc))

    if bool_sh:
        return Ra, rcm.SHF, rcm.Ts, rcm.Tatm, rcm, state
    else:
        return Ra, rcm.Ts, rcm.Tatm, rcm, state


# run with and without sh
co2_list = 1e-6*np.linspace(355,2000,11)

ra = np.empty([len(co2_list),1])
sh = np.empty([len(co2_list),1])
ts = np.empty([len(co2_list),1])
tatm = np.empty([len(co2_list),len(plev_gcm)])

# location of pickled data
co2file = 'co2sweep.pickle'

if not (os.path.isfile(co2file) and try_load):
    counter = 0
    for co2 in co2_list:
        # [ra[counter], sh[counter], ts[counter], tatm[counter,:],_, _] = run_model(1, co2)
        [ra[counter], ts[counter], tatm[counter,:],_, _] = run_model(0, co2)
        counter = counter + 1

    # pickle.dump([co2_list, ra, sh, ts, tatm], open(co2file, 'wb'))
    pickle.dump([co2_list, ra, ts, tatm], open(co2file, 'wb'))
else:
    # [co2_list, ra, sh, ts, tatm] = pickle.load(open(co2file, 'rb'))
    [co2_list, ra, ts, tatm] = pickle.load(open(co2file, 'rb'))

# PLOT CO2 vs Ra
fig, ax = plt.subplots()
ax.plot(1e6*co2_list, ra, '.k')
fig.set_size_inches(5,4)
# ax.set_ylim([1e2,1e3])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('$p_{CO_2}$ (ppmv)')
ax.set_ylabel('$R_a$ (W m$^{-2}$)')
plt.savefig('co2_ra.pdf', format='pdf', dpi=300)
plt.close()

# # PLOT CO2 vs Delta Ra and SH
# fig, ax = plt.subplots()
# ax.axhline(0, color='k', linewidth=0.5)
# ax.plot(1e6*co2_list, sh-sh[0], '.', color='orange', label='SH')
# ax.plot(1e6*co2_list, ra-ra[0], '.k', label='$R_a$')
# fig.set_size_inches(5,4)
# # ax.set_ylim([1e2,1e3])
# ax.xaxis.set_minor_locator(AutoMinorLocator())
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.set_xlabel('$p_{CO_2}$ (ppmv)')
# ax.set_ylabel('$\Delta$ Energy flux (W m$^{-2}$)')
# ax.legend()
# plt.savefig('co2_dflux.pdf', format='pdf', dpi=300)
# plt.close()

# PLOT CO2 vs Ts
fig, ax = plt.subplots()
ax.plot(1e6*co2_list, tatm[:,-1], '.k')
fig.set_size_inches(5,4)
# ax.set_ylim([1e2,1e3])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('$p_{CO_2}$ (ppmv)')
ax.set_ylabel('$T$ (K)')
plt.savefig('co2_ts.pdf', format='pdf', dpi=300)
plt.close()

# PLOT temp responses at various co2
fig, ax = plt.subplots()
ax.axvline(0, color='k', linewidth=0.5)
for co2 in range(len(co2_list[1:])):
    ax.plot(tatm[co2+1,:]-tatm[0,:], plev_gcm, label='%.0f ppmv' % (1e6*co2_list[co2+1]))
fig.set_size_inches(5,4)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('$p_{CO_2}$ (ppmv)')
ax.set_ylabel('$\Delta T$ (K)')
ax.set_ylim([1e2,1e3])
ax.invert_yaxis()
ax.legend()
plt.savefig('tatm_response.pdf', format='pdf', dpi=300)
plt.close()
