import numpy as np
from scipy import integrate as sciint
from scipy import interpolate
import matplotlib.pyplot as plt
import climlab
from climlab import constants as cn

class Ebm2D:

    def __init__(self, tstep=None, tsteps=None, N_lat=None, N_lev=None, sol_state=None, eccent=None, precess=None, obliq=None, d_ml=None, albedo=None, RH_surf=None, conv_scheme=None, transport=None):
        # DEFAULT VALUES
        # TIMESTEP
        if tstep is None: tstep = 24*3600 # 1 day time step for radiation
        if tsteps is None: tsteps = 3600 # 1 hourly time step for everything else
        # CONFIGURE DOMAIN
        if N_lat is None: N_lat = 40
        if N_lev is None: N_lev = 30
        if albedo is None: albedo = 0.07
        # CONFIGURE INSOLATION
        if sol_state is None: sol_state = 'seasonal'
        if eccent is None: eccent = 0.
        if precess is None: precess = 0.
        if obliq is None: obliq = 23.446
        # PRESCRIBE BOUNDARY CONDITIONS
        if d_ml is None: d_ml = 5. # mixed layer depth used to calculate surface heat capacity [m]
        if RH_surf is None: RH_surf = 0.8 # surface relative humidity [unitless]
        # CONFIGURE CONVECTION SCHEME
        if conv_scheme is None: conv_scheme = 6.5 # hard adjustment convection scheme [MALR, DALR, or set custom lapse rate]
        # CONFIGURE HEAT TRANSPORT SCHEME
        if transport is None: transport = 'dry_diff' # diffuse internal energy

        self.tstep = tstep
        self.tsteps = tsteps
        self.N_lat = N_lat
        self.N_lev = N_lev
        self.albedo = albedo
        self.sol_state = sol_state
        self.orb = {'ecc':eccent, 'long_peri':precess, 'obliquity':obliq}
        self.d_ml = d_ml
        self.RH_surf = RH_surf
        self.conv_scheme = conv_scheme
        self.transport = transport

        # initialize state with mixed layer depth d_ml
        self.full_state = climlab.column_state(num_lev=self.N_lev, num_lat=self.N_lat, water_depth=self.d_ml)
        self.temp_state = {'Tatm':self.full_state.Tatm, 'Ts':self.full_state.Ts}

        self.h2o = climlab.radiation.ManabeWaterVapor(state=self.full_state, timestep=self.tsteps)
        if self.conv_scheme == 'emanuel':
            # self.full_state['q'] = np.ones_like(self.full_state.Tatm) * 5.E-6
            self.full_state['q'] = self.h2o.q
        # construct the model
        self.model = climlab.TimeDependentProcess(state=self.full_state, timestep=self.tsteps) #num_steps_per_year=12)

        # insolation scheme
        if self.sol_state == 'seasonal':
            self.sol = climlab.radiation.DailyInsolation(domains=self.model.Ts.domain, orb=self.orb, timestep=self.tstep)
        elif self.sol_state == 'annual':
            self.sol = climlab.radiation.AnnualMeanInsolation(domains=self.model.Ts.domain, timestep=self.tstep)

        # these processes depend on the chosen convection scheme
        if self.conv_scheme == 'emanuel':
            self.conv = climlab.convection.EmanuelConvection(state=self.full_state, timestep=self.tsteps)
            self.rad = climlab.radiation.RRTMG(state=self.temp_state, insolation=self.sol.insolation, specific_humidity=self.full_state.q, coszen=self.sol.coszen, timestep=self.tstep, albedo=self.albedo, ozone_file=None)
        else:
            self.conv = climlab.convection.ConvectiveAdjustment(state=self.temp_state, adj_lapse_rate=self.conv_scheme, timestep=self.tsteps)
            self.rad = climlab.radiation.RRTMG(state=self.temp_state, insolation=self.sol.insolation, specific_humidity=self.h2o.q, coszen=self.sol.coszen, timestep=self.tstep, albedo=self.albedo, ozone_file=None)

        self.model.add_subprocess('radiation', self.rad)
        self.model.add_subprocess('insolation', self.sol)
        self.model.add_subprocess('convection', self.conv)
        if self.conv_scheme != 'emanuel':
            self.model.add_subprocess('water', self.h2o)

        [self.mesh_lat, self.mesh_lev] = np.meshgrid(self.model.lat, self.model.lev, indexing='ij')

    def integrate(self, N_yr):
        self.model.integrate_years(N_yr)

    def add_diff(self, D, U=None):
        if U is None: U=0

        if self.transport == 'dry_diff':
            # meridional diffusivity in m**2/s
            K = D / self.model.Tatm.domain.heat_capacity[0] * cn.a**2
            self.diff = climlab.dynamics.MeridionalDiffusion(state={'Tatm': self.model.state['Tatm']}, K=K, **self.model.param)
        elif self.transport == 'adv_diff':
            K = D / self.model.Tatm.domain.heat_capacity[0] * cn.a**2
            self.diff = climlab.dynamics.MeridionalAdvectionDiffusion(state={'Tatm': self.model.state['Tatm']}, K=K, U=U, **self.model.param)
        elif self.transport == 'moist_diff':
            self.diff = climlab.dynamics.MeridionalMoistDiffusion(state={'Tatm': self.model.state['Tatm'], 'Ts': self.model.state['Ts']}, **self.model.param)
        self.model.add_subprocess('diffusion', self.diff)

    def add_stf(self, Cd):
        #  Add surface heat fluxes
        self.shf = climlab.surface.SensibleHeatFlux(state=self.temp_state, Cd=Cd, timestep=self.tsteps)
        self.lhf = climlab.surface.LatentHeatFlux(state=self.full_state, Cd=Cd, timestep=self.tsteps)
        if self.conv_scheme != 'emanuel':
            # set the water vapor input field for LHF
            self.lhf.q = self.h2o.q
        self.model.add_subprocess('sh', self.shf)
        self.model.add_subprocess('lh', self.lhf)

    def add_co2(self, nxco2):
        self.model.subprocess.radiation.input['absorber_vmr']['CO2'] = nxco2 * self.model.subprocess.radiation.input['absorber_vmr']['CO2']

    def comp_ra(self):
        self.ra = self.model.ASR - self.model.SW_sfc - (self.model.OLR - self.model.LW_sfc)
        return self.ra

    def comp_r1(self):
        self.r1 = self.infer_divfm()/self.comp_ra()
        return self.r1

    def infer_divfm(self):
        self.divfm = self.comp_ra() + self.model.LHF + self.model.SHF
        return self.divfm

    def infer_heat_transport(self):
        '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
        self.Rtoa = np.squeeze(self.model.timeave['ASR'] - self.model.timeave['OLR'])
        lat_rad = np.deg2rad( self.model.lat )
        return ( 1E-15 * 2 * np.math.pi * cn.a**2 *
                sciint.cumtrapz( np.cos(lat_rad)*self.Rtoa,
                x=lat_rad, initial=0. ) )

    def comp_toa_imbalance(self):
        # calculate globally averaged TOA energy flux imbalance
        trans = self.infer_heat_transport()
        self.toa_net = 1e15*trans[-1]/(4*np.math.pi*cn.a**2)
        return self.toa_net

    def plot_T(self, filename):
        plt.figure()
        plt.contourf(self.mesh_lat, self.mesh_lev, self.model.Tatm, cmap=plt.cm.coolwarm)
        plt.xlabel('Latitude (deg)')
        plt.ylabel('Pressure (hPa)')
        plt.colorbar(label='Temperature (K)')
        plt.gca().invert_yaxis()
        plt.savefig(filename)
        plt.close()

    def plot_Ts(self, filename):
        plt.figure()
        plt.plot(self.model.lat, self.model.Ts, 'k')
        plt.xlabel('Latitude (deg)')
        plt.ylabel('$T_s$ (K)')
        plt.savefig(filename)
        plt.close()

    def plot_T_sel(self, filename):
        ps = self.model.lev_bounds[-1]
        pfull = np.append(self.model.lev, ps)
        tfull = np.append(self.model.Tatm, self.model.Ts, axis=1)
        t_lat = interpolate.interp1d(self.model.lat, tfull, axis=0)
        plt.figure()
        plt.plot(t_lat(0), pfull, 'C1')
        plt.plot(t_lat(45), pfull, 'C7')
        plt.plot(t_lat(85), pfull, 'C0')
        plt.xlabel('T (K)')
        plt.ylabel('p (hPa)')
        plt.gca().invert_yaxis()
        plt.savefig(filename)
        plt.close()

    def plot_energy_fluxes(self, filename):
        plt.figure()
        plt.plot(self.model.lat, self.comp_ra(), 'C7', label='$R_a$')
        plt.plot(self.model.lat, self.infer_divfm(), 'C3', label=r'$\nabla \cdot F_m$')
        plt.plot(self.model.lat, self.model.LHF, 'C0', label='LH')
        plt.plot(self.model.lat, self.model.SHF, 'C1', label='SH')
        plt.xlabel('Latitude (deg)')
        plt.ylabel('Energy flux (Wm$^{-2}$)')
        plt.legend()
        plt.savefig(filename)
        plt.close()

    def plot_r1(self, filename):
        plt.figure()
        plt.plot(self.model.lat, self.comp_r1(), 'k')
        plt.xlabel('Latitude (deg)')
        plt.ylabel('$R_1$ (unitless)')
        plt.savefig(filename)
        plt.close()

    def plot_stf(self, filename):
        plt.figure()
        plt.plot(self.model.lat, self.model.timeave['LHF'], label='LH')
        plt.plot(self.model.lat, self.model.timeave['SHF'], label='SH')
        plt.xlabel('Latitude (deg)')
        plt.ylabel('Energy flux (Wm$^{-2}$)')
        plt.legend()
        plt.savefig(filename)
        plt.close()

    def plot_trans(self, filename):
        plt.figure()
        plt.plot(self.model.lat, self.infer_heat_transport())
        plt.xlabel('Latitude (deg)')
        plt.ylabel('Energy transport (PW)')
        plt.savefig(filename)
        plt.close()
