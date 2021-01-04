from ebm2d import Ebm2D
import climlab
from climlab import constants as cn
import numpy as np
import scipy.interpolate as interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
from calc_ma import calc_ma
import os
import pickle

mpl.rcParams['font.size'] = 10
mpl.rcParams['mathtext.fontset'] = 'dejavusans'

class SeasonalEbm:

    def __init__(self, plotpath, datapath, tstep=None, tsteps=None, N_lat=None, N_lev=None, mld=None, albedo=None, conv_scheme=None, transport=None, sol_state=None, D=None, U=None, Cd=None, Ninit=None, forcing=None):
        # default values if none is specified
        if tstep is None: tstep = 24*3600
        if tsteps is None: tsteps = 3600
        if N_lat is None: N_lat = 40
        if N_lev is None: N_lev = 30
        if mld is None: mld = 5
        if albedo is None: albedo = 0.07
        if conv_scheme is None: conv_scheme = 6.5
        if transport is None: transport = 'dry_diff'
        if sol_state is None: sol_state = 'seasonal'
        if D is None: D = 0.04
        if U is None: U = 0
        if Cd is None: Cd = 3e-3
        if Ninit is None: Ninit = 5

        self.plotpath = plotpath
        self.datapath = datapath
        self.tstep = tstep
        self.tsteps = tsteps
        self.N_lat = N_lat
        self.N_lev = N_lev
        self.mld = mld
        self.albedo = albedo
        self.conv_scheme = conv_scheme
        self.transport = transport
        self.sol_state = sol_state
        self.D = D
        self.U = U
        self.Cd = Cd
        self.Ninit = Ninit
        self.forcing = forcing

        self.ebm = Ebm2D(tstep = self.tstep, tsteps=self.tsteps, N_lat=self.N_lat, N_lev=self.N_lev, d_ml = self.mld, albedo=self.albedo, conv_scheme=self.conv_scheme, transport=self.transport ,sol_state=self.sol_state)
        self.ebm.add_diff(self.D, self.U)
        self.ebm.add_stf(self.Cd)

    def run_nyears(self):
        # steady state test
        t_end = self.Ninit # integrate for Ninit years
        dt = 1/12 # one month steps
        Nt = int(t_end/dt)
        t = np.linspace(0, t_end, Nt)

        if self.forcing == '4xco2':
            imbalance = np.empty([2*Nt, 1])
        else:
            imbalance = np.empty([Nt, 1])

        plt.figure()
        for i in np.arange(0,Nt):
            self.ebm.integrate(dt)
            imbalance[i] = self.ebm.comp_toa_imbalance()

            plt.cla()
            plt.hlines(0,0,t[i], 'k')
            plt.plot(t[:i], imbalance[:i])
            plt.xlabel('t (yr)')
            plt.ylabel('$F_{TOA}$ (Wm$^{-2}$)')
            plt.savefig('%s/imbalance.png' % (self.plotpath))
            plt.close()

        if self.forcing == '4xco2':
            self.ebm.add_co2(4)
            t_end = 2*self.Ninit
            Nt = int(t_end/dt)
            t = np.linspace(0, t_end, Nt)
            for i in np.arange(0,Nt)+Nt:
                self.ebm.integrate(dt)
                imbalance[i] = self.ebm.comp_toa_imbalance()

                plt.cla()
                plt.hlines(0,0,t[i], 'k')
                plt.plot(t[:i], imbalance[:i])
                plt.xlabel('t (yr)')
                plt.ylabel('$F_{TOA}$ (Wm$^{-2}$)')
                plt.savefig('%s/imbalance.png' % (self.plotpath))
                plt.close()

        t_yr = np.arange(0.5,0.5+t_end)
        imbalance_yr = np.empty([t_end, 1])
        for yr in np.arange(0,t_end):
            idx = np.arange(0,12) + 12*yr
            imbalance_yr[yr] = imbalance[idx].sum()/12

        plt.figure()
        plt.hlines(0,0,t[i], 'k')
        plt.plot(t, imbalance)
        plt.plot(t_yr, imbalance_yr)
        plt.xlabel('t (yr)')
        plt.ylabel('$F_{TOA}$ (Wm$^{-2}$)')
        plt.savefig('%s/imbalance.png' % (self.plotpath))
        plt.close()


    def run_year(self):
        t_end = 1 # integrate for 1 year to obtain seasonality
        dt = 1/12 # one month steps
        Nt = int(t_end/dt)
        t = np.linspace(0, t_end, Nt)

        self.r1 = np.empty([len(self.ebm.model.lat), Nt])
        self.ra = np.empty([len(self.ebm.model.lat), Nt])
        self.divfm = np.empty([len(self.ebm.model.lat), Nt])
        self.lh = np.empty([len(self.ebm.model.lat), Nt])
        self.sh = np.empty([len(self.ebm.model.lat), Nt])

        self.Ts = np.empty([len(self.ebm.model.lat), Nt])
        self.Tatm = np.empty([len(self.ebm.model.lat), len(self.ebm.model.lev), Nt])

        [self.mesh_mon, self.mesh_lat] = np.meshgrid(np.arange(1,13), self.ebm.model.lat, indexing='ij')

        for i in np.arange(0,Nt):
            self.ebm.integrate(dt)
            self.r1[:,[i]] = self.ebm.comp_r1()
            self.ra[:,[i]] = self.ebm.comp_ra()
            self.divfm[:,[i]] = self.ebm.infer_divfm()
            self.lh[:,[i]] = self.ebm.model.LHF
            self.sh[:,[i]] = self.ebm.model.SHF

            self.Ts[:,[i]] = self.ebm.model.Ts
            self.Tatm[:,:,[i]] = np.expand_dims(self.ebm.model.Tatm, axis=2)

    def plot_energy_fluxes(self, filename):
        lh = np.nanmean(self.lh, axis=1)
        sh = np.nanmean(self.sh, axis=1)
        ra = np.nanmean(self.ra, axis=1)
        divfm = np.nanmean(self.divfm, axis=1)

        fig, ax=plt.subplots(1, figsize=(10/3, 10/3), dpi=300)
        plt.hlines(0,-90,90, 'k')
        ax.plot(self.ebm.model.lat, ra, 'C7', label='$R_a$')
        ax.plot(self.ebm.model.lat, divfm, 'C3', label=r'$\nabla\cdot F_m$')
        ax.plot(self.ebm.model.lat, lh, 'C0', label='LH')
        ax.plot(self.ebm.model.lat, sh, 'C1', label='SH')
        ax.set_xlabel('Latitude (deg)')
        # ax.set_ylim([-150, 200])
        ax.set_ylabel('Energy flux (Wm$^{-2}$)')
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        ax.set_xlim([-90,90])
        ax.set_xticks(np.arange(-90,91,30))
        ax.tick_params(which='both', top='on', right='on')
        plt.title('climlab 2D EBM, %g m' % (self.mld), fontsize=10)
        plt.tight_layout()
        plt.savefig('%s/%s.png' % (self.plotpath, filename))

        # fig.set_size_inches(21/3, 7/3)
        # plt.tight_layout(rect=[0,0,0.6,1])
        # plt.legend(bbox_to_anchor=(1.2, 0.5), loc='center left', ncol=2)
        # plt.savefig('%s/%s_leg.png' % (self.plotpath, filename))
        # plt.close()

    def plot_r1_ann(self, eps, filename):
        r1_mon = self.divfm/self.ra
        r1 = np.nanmean(self.r1, axis=1)

        fig, ax=plt.subplots(1, figsize=(10/3, 7/3), dpi=300)
        ax.add_patch(patches.Rectangle([-90,-0.6], 180, eps+0.6, facecolor='C1', alpha=0.5))
        ax.add_patch(patches.Rectangle([-90,1-eps], 180, 1.2, facecolor='C0', alpha=0.5))
        ax.plot(self.ebm.model.lat, r1, 'k', label='$R_1$')
        ax.set_xlabel('Latitude (deg)')
        ax.set_ylabel('Energy flux (Wm$^{-2}$)')
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        ax.set_xlim([-90,90])
        ax.set_xticks(np.arange(-90,91,30))
        ax.set_ylim([-0.6,1])
        ax.tick_params(which='both', top='on', right='on')
        plt.title('climlab 2D EBM, %g m' % (self.mld), fontsize=10)
        plt.tight_layout()
        plt.savefig('%s/%s.png' % (self.plotpath, filename))

    def plot_mon_lat(self, eps, varname, clabel, filename):
        fig,ax = plt.subplots(1, figsize=(13/3, 7/3), dpi=300)
        if varname == 'r1':
            levels=np.arange(-2,2,0.2)
            norm=colors.Normalize(vmin=-1, vmax=1)
            # norm=colors.BoundaryNorm(np.arange(-1,1,0.2), 10)
        h = ax.contourf(self.mesh_mon, self.mesh_lat, np.transpose(getattr(self, varname)), levels=levels, norm=norm, cmap=plt.cm.bwr_r, extend='both')
        ax.contour(self.mesh_mon, self.mesh_lat, np.transpose(getattr(self, varname)), levels=[eps], colors='C1')
        ax.contour(self.mesh_mon, self.mesh_lat, np.transpose(getattr(self, varname)), levels=[1-eps], colors='C0')
        ax.set_ylim([-90,90])
        ax.set_yticks(np.arange(-90,91,30))
        ax.set_ylabel('Latitude (deg)')
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        ax.set_xlim([1,12])
        ax.set_xticks(np.arange(1,13))
        ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
        fig.colorbar(h, label=clabel)
        plt.tight_layout()
        plt.savefig('%s/%s.png' % (self.plotpath, filename))
        plt.close()

    def plot_mse(self, loc, filename):
        if loc == 'nhmid':
            lat_lo = 40; lat_up = 50
            lat = np.arange(lat_lo,lat_up,0.25)
            lat = np.expand_dims(lat, axis=-1)
            clat = np.cos(lat*np.pi/180)

        lhlat = interpolate.interp1d(self.ebm.model.lat, self.lh, axis=0)
        lh = np.sum(np.squeeze(lhlat(lat))*clat, axis=0)/clat.sum(axis=0)
        shlat = interpolate.interp1d(self.ebm.model.lat, self.sh, axis=0)
        sh = np.sum(np.squeeze(shlat(lat))*clat, axis=0)/clat.sum(axis=0)
        divfmlat = interpolate.interp1d(self.ebm.model.lat, self.divfm, axis=0)
        divfm = np.sum(np.squeeze(divfmlat(lat))*clat, axis=0)/clat.sum(axis=0)
        ralat = interpolate.interp1d(self.ebm.model.lat, self.ra, axis=0)
        ra = np.sum(np.squeeze(ralat(lat))*clat, axis=0)/clat.sum(axis=0)

        fig, ax=plt.subplots(1, figsize=(13/3, 7/3), dpi=300)
        plt.hlines(0,1,12, 'k')
        ax.plot(np.arange(1,13), ra, 'C7', label='$R_a$')
        ax.plot(np.arange(1,13), divfm, 'C3', label=r'$\nabla\cdot F_m$')
        ax.plot(np.arange(1,13), lh, 'C0', label='LH')
        ax.plot(np.arange(1,13), sh, 'C1', label='SH')
        ax.set_ylim([-150, 100])
        ax.set_ylabel('Energy flux (Wm$^{-2}$)')
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        ax.set_xlim([1,12])
        ax.set_xticks(np.arange(1,13))
        ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
        ax.tick_params(which='both', top='on', right='on')
        plt.title('climlab 2D EBM, %g m, $\phi=%g^\circ$ to $%g^\circ$' % (self.mld, lat_lo, lat_up), fontsize=10)
        plt.tight_layout()
        plt.savefig('%s/%s.png' % (self.plotpath, filename))

        fig.set_size_inches(21/3, 7/3)
        plt.tight_layout(rect=[0,0,0.6,1])
        plt.legend(bbox_to_anchor=(1.2, 0.5), loc='center left', ncol=2)
        plt.savefig('%s/%s_leg.png' % (self.plotpath, filename))
        plt.close()

    def plot_dr1(self, eps, loc, filename):
        if loc == 'nhmid':
            lat_lo = 40; lat_up = 50
            lat = np.arange(lat_lo,lat_up,0.25)
            lat = np.expand_dims(lat, axis=-1)
            clat = np.cos(lat*np.pi/180)

        r1lat = interpolate.interp1d(self.ebm.model.lat, self.r1, axis=0)
        r1 = np.sum(np.squeeze(r1lat(lat))*clat, axis=0)/clat.sum(axis=0)

        divfmlat = interpolate.interp1d(self.ebm.model.lat, self.divfm, axis=0)
        divfm = np.sum(np.squeeze(divfmlat(lat))*clat, axis=0)/clat.sum(axis=0)
        ralat = interpolate.interp1d(self.ebm.model.lat, self.ra, axis=0)
        ra = np.sum(np.squeeze(ralat(lat))*clat, axis=0)/clat.sum(axis=0)

        r1_ann = np.nanmean(r1)
        ra_ann = np.nanmean(ra)
        divfm_ann = np.nanmean(divfm)

        dr1 = r1-r1_ann
        dra = ra-ra_ann
        ddivfm = divfm-divfm_ann

        comp1 = ddivfm/ra_ann
        comp2 = -divfm_ann/ra_ann**2 * dra
        res = dr1 - comp1 - comp2

        ylim_lo = -0.2;
        ylim_up = 0.5;

        fig,ax = plt.subplots(1, figsize=(13/3,7/3), dpi=300)
        plt.hlines(r1_ann, 1,12, 'k')
        ax.add_patch(patches.Rectangle([1,-0.6], 11, eps+0.6, facecolor='C1', alpha=0.5))
        ax.plot(np.arange(1,13), r1, 'k', label='$R_1$')
        ax.set_ylabel('$R_1$ (unitless)')
        ax.set_xlim([1,12])
        ax.set_ylim([min(np.append(r1, ylim_lo)), max(np.append(r1, ylim_up))])
        ax.set_xticks(np.arange(1,13))
        ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        ax.tick_params(which='both', top='on')
        ax2 = ax.twinx()
        ax2.plot(np.arange(1,13), dr1, 'k', label='$\Delta R_1$')
        ax2.plot(np.arange(1,13), comp1, 'C3', label=r'$\frac{\Delta (\nabla \cdot F_m)}{\overline{R_a}}$')
        ax2.plot(np.arange(1,13), res, '-.k', label='Residual')
        ax2.plot(np.arange(1,13), comp2, 'C7', label=r'$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a}^2}\Delta R_a$')
        ax2.set_ylabel('$\Delta R_1$ (unitless)')
        ax2.set_xlim([1,12])
        ax2.set_ylim([min(np.append(r1, ylim_lo))-r1_ann, max(np.append(r1, ylim_up))-r1_ann])
        ax2.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        plt.title('climlab 2D EBM, %g m, $\phi=%g^\circ$ to $%g^\circ$' % (self.mld, lat_lo, lat_up), fontsize=10)
        plt.tight_layout()
        plt.savefig('%s/%s.png' % (self.plotpath, filename))

        fig.set_size_inches(21/3, 7/3)
        plt.tight_layout(rect=[0,0,0.6,1])
        plt.legend(bbox_to_anchor=(1.3, 0.5), loc='center left', ncol=2)
        plt.savefig('%s/%s_leg.png' % (self.plotpath, filename))
        plt.close()

    def plot_temp_ann(self, filename):
        ps = self.ebm.model.lev_bounds[-1]
        pfull = np.append(self.ebm.model.lev, ps)
        tfull = np.append(self.ebm.model.Tatm, self.ebm.model.Ts, axis=1)
        t_lat = interpolate.interp1d(self.ebm.model.lat, tfull, axis=0)

        fig,ax = plt.subplots(1, figsize=(10/7,10/7), dpi=300)
        ax.plot(t_lat(0), pfull, 'C1')
        ax.plot(t_lat(45), pfull, 'C7')
        ax.plot(t_lat(85), pfull, 'C0')
        ax.set_xlabel('T (K)')
        ax.set_ylabel('p (hPa)')
        plt.savefig('%s/%s.png' % (self.plotpath, filename))
        plt.close()

    def plot_temp_sel(self, filename):
        pinit = 950 # initial parcel level for moist adiabat (hPa)
        ps = self.ebm.model.lev_bounds[-1]
        pfull = np.append(self.ebm.model.lev, ps)
        si = pfull/ps

        tfull = np.append(self.Tatm, np.expand_dims(self.Ts, axis=1), axis=1)
        # tjan = np.squeeze(tfull[:,:,1]); tjun = np.squeeze(tfull[:,:,6]);
        t_lat = interpolate.interp1d(self.ebm.model.lat, tfull, axis=0)

        monlist = [1,6] # months to plot
        monlabel = ['January','June'] # months to plot
        counter = 0
        for mon in monlist:
            t_lev=interpolate.interp1d(pfull, t_lat(45)[:,mon], axis=0)
            ma=calc_ma(pinit*100, t_lev(pinit), 100, pfull*100, frz=0); # initiate saturated moist adiabat from level pinit
            fig,ax = plt.subplots(1, figsize=(10/3,10/3), dpi=300)
            # ax.plot(t_lat(0)[:,mon], si, 'C1')
            # ax.plot(ma, si, ':C1')
            if mon == 1:
                ax.plot(t_lat(45)[:,mon], si, 'C7')
                ax.plot(ma, si, ':C7')
            elif mon == 6:
                ax.plot(t_lat(45)[:,mon], si, 'C1')
                ax.plot(ma, si, ':C1')
                # ax.plot(t_lat(85)[:,mon], si, 'C0')
            ax.set_xlabel('T (K)')
            ax.set_ylabel('$\sigma$ (unitless)')
            ax.set_xlim([210,290])
            ax.set_ylim([1,0.2])
            plt.title('climlab 2D EBM, %g m, %s' % (self.mld, monlabel[counter]))
            plt.tight_layout()
            plt.savefig('%s/%s_%s.png' % (self.plotpath, filename, monlabel[counter]))
            plt.close()
            counter =+ 1
