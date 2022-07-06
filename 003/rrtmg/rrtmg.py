import os
import sys
import pickle
import climlab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def rrtmg(pco2,T,Ts,q,plev,**kwargs):

    lat=kwargs.get('lat', 83.3)
    day=kwargs.get('day', 0)
    alb=kwargs.get('alb', 0.6)
    cld=kwargs.get('cld', False)

    par={}
    par['g']=9.81 # m s**-2
    par['cp']=1004 # J kg**-1 K**-1
    par['ps']=1000 # hPa
    par['d']=1 # mixed layer depth [m]

    Q = climlab.solar.insolation.daily_insolation( lat, day )

    #  State variables (Air and surface temperature)
    state = climlab.column_state(lev=plev, water_depth=par['d'])
    state['Tatm'][:] = T
    state['Ts'][:] = Ts

    # #  Initialize a nearly dry column (small background stratospheric humidity)
    # q = np.ones_like(state.Tatm) * 5.E-6
    #  Add specific_humidity to the state dictionary
    # state['q'] = q
    # state['q'][:] = gcm_clim['hus']

    #  Parent model process
    rcm = climlab.TimeDependentProcess(state=state)

    #  RADIATION
    abs_vmr = climlab.radiation.radiation.default_absorbers(state['Tatm'])
    abs_vmr['CO2'] = pco2
    # abs_vmr['CO2'] = 700e-6

    # abs_vmr = { 'CO2':co2_lev,
    #         'CH4':1662.2e-9, # historical 1987
    #         'N2O':306.2625e-9, # historical 1987
    #         'O2':0.21, # historical 1987
    #         'CFC11':491.2003e-12, # historical 1987
    #         'CFC12':408.65e-12, # historical 1987
    #         'CFC22':0.,
    #         'CCL4':0.,
    #         'O3':0.}

    #  Daily insolation as a function of latitude and time of year
    # rad = climlab.radiation.RRTMG(state=state, specific_humidity=h2o.q, albedo=alb, absorber_vmr=abs_vmr)
    rad = climlab.radiation.RRTMG(state=state, icld=cld, specific_humidity=q, albedo=alb, insolation=Q, absorber_vmr=abs_vmr)

    rad.compute_diagnostics()

    # # #  Check for energy balance
    # print('SCM RTOA = %g W m**-2' % (rad.ASR - rad.OLR) )
    # print('GCM RTOA = %g W m**-2\n' % gcm_clim['ftoacs'] )

    # print('SCM ASR = %g W m**-2' % (rad.ASR) )
    # print('GCM ASR = %g W m**-2' % (gcm_clim['rsdt']-gcm_clim['rsutcs']) )

    # print('SCM OLR = %g W m**-2' % (rad.OLR) )
    # print('GCM OLR = %g W m**-2' % (gcm_clim['rlutcs']) )

    # print(q)
    # print(T)
    # print(Ts)

    # # total radiative cooling
    # Ra = rad.ASR - rad.SW_sfc + rad.LW_sfc - rad.OLR
    # Rsfc = rad.SW_sfc - rad.LW_sfc
    # print('Ra = %g W m**-2' % (Ra))
    # print('GCM Ra = %g W m**-2\n' % gcm_clim['racs'])
    # print('Rsfc = %g W m**-2' % (Rsfc))

    return rad
