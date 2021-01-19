clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));
addpath(genpath('./matlab'));

gcm_info
echam_info

figure_params

%% set parameters
% lat grid type
if 1
par.echam_clims = {'echr0001'}; % par.echam.all_mld; % par.echam.sel; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
% par.echam_clims = par.echam.noice_mld; % {'echr0001'}; % par.echam.sel; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
par.era5.yr_span = '1979_2005';
par.era5c.yr_span = '1979_2005';
par.jra55.yr_span = '1979_2005';
par.merra2.yr_span = '1980_2005';
par.ep_swp = 0.1; %[0.25 0.3 0.35]; % threshold value for determining RCE
par.ga_swp = 0.9; % threshold for determining RAE
par.si_eval = [0.8 0.85 0.9]; % sigma level for calculating inversion strength
par.ma_init = 0.95; % 'surf' for initializing with 2 m data, otherwise enter starting sigma level for moist adiabat
par.pa_eval = 500e2; % pressure level for calculating ma_diff
par.si_bl_swp = [0.85 0.9 0.95]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_up = 0.4; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.ta_thresh = 6.5; % criteria for temperature difference in K of plotting contour for closeness to moist adiabat
par.ga_thresh = 10; % criteria for lapse rate difference in % of plotting contour for closeness to moist adiabat
par.ga_bl_thresh = 90; % criteria for lapse rate difference in % of plotting contour for closeness to moist adiabat
par.inv_thresh = -4; % criteria for temperature difference in K of plotting contour for inversion strength
par.albedo_thresh = 0.6; % criteria for surface albedo of plotting contour for inversion strength
par.sn_thresh = 1e-1; % criteria for surface albedo of plotting contour for inversion strength
par.alb_thresh = 0.8; % criteria for surface albedo of plotting contour for inversion strength
par.albcs_thresh = 0.8; % criteria for surface albedo of plotting contour for inversion strength
par.r1_bins = [-0.55:0.1:1.35]; % bins for sorting temperature profiles according to r1 values
% par.era.fw = {'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt', 'div79'};
% par.era.fw = {'div79', 'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt'};
par.era.fw = {'mse', 'mse_ac', 'mse_sc', 'mse_ac_ra', 'mse_sc_ra', 'dse'};
par.jra55.fw = {'mse', 'dse'};
par.merra2.fw = {'mse', 'dse'};
par.gcm.fw = {'mse', 'dse'};
par.echam.fw = {'mse', 'dse'};
par.pa = linspace(1000,10,100)*1e2; % high resolution vertical grid to interpolate to
par.z = [0:500:40e3]';
% set how to close energy budget
% if == teten, use TETEN data from Donohoe to close energy budget (option for ERA-Interim)
% if == stf, use SH and LH data from ERA-Interim to close energy budget (option for ERA-Interim)
% if == era5, use ERA5 radiative cooling and surface turbulent fluxes to close energy budget (only option for ERA5)
par.closure = 'era5';
par.do_surf = 0; % whether or not to calculate temperature profile in pressure grids including ts and ps data
par.cpd = 1005.7; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
end

%% call functions
% plot_rad_lat(par)
% plot_rad_lon_lat(par)
% plot_tediv_lat(par)

% type = 'jra55';
% par.lat_interp = 'native';
% choose_plots(type, par);
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    type='echam';
    par.lat_interp = 'native';
    disp(par.echam.clim)
    choose_plots(type, par);
end
for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
    % type = 'gcm';
    % disp(par.model)
    % choose_plots(type, par);
end

% % sweep through various boundary layer heights
% for i = 1:length(par.si_bl_swp); par.si_bl = par.si_bl_swp(i);
%     type = 'merra2';
%     % par.lat_interp = 'native';
%     % choose_plots_si_bl(type, par)
%     for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
%         type='echam';
%         % par.lat_interp = 'native';
%         % disp(par.echam.clim)
%         % choose_plots_si_bl(type, par);
%     end
%     for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
%         type = 'gcm';
%         disp(par.model)
%         choose_plots_si_bl(type, par)
%     end
% end

% % sweep through various threshold values
% for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = par.ga_swp(i);
%     % type = 'era5c'; par.lat_interp = 'native';
%     % choose_plots_ep(type, par)
%     for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
%         type='echam'; par.lat_interp = 'native';
%         disp(par.echam.clim)
%         choose_plots_ep(type, par);
%     end
%     for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
%         % type = 'gcm';
%         % disp(par.model)
%         % choose_plots_ep(type, par)
%     end
% end

function choose_plots(type, par)
    % plot_temp_zon_select(type, par) % plot temperature profiles at specific latitudes
    % plot_temp_binned_r1(type, par) % plot temperature profiles at specific latitudes
    % plot_dmse_midlatitude_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dmse_polar_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dlh_polar_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dlh_polar_line_so(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    plot_srfc_polar_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_siced(type, par) % sea ice depth
    % plot_friac(type, par) % sea ice fraction
    % plot_ahfres(type, par) % ice melt
    % plot_ahfliac(type, par) % LH over ice
    % plot_ahfllac(type, par) % LH over land
    % plot_ahflwac(type, par) % LH over water
    % plot_alb(type, par) % surface albedo
    % plot_sftlf(type, par) % land fraction

    % plot_dmse_toasfc_midlatitude_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dmse_polar_line_asym(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dmse_polar_line_topocomp(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dmse_so_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_tas(type, par) % 2m temperature
    % plot_tas_midlatitude_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_sol_midlatitude_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_ts(type, par) % skin temperature
    % plot_dmse_toasfc_midlatitude_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dra_polar_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_ma_diff(type, par) % plot difference of temperature profile from moist adiabat
    % plot_ga_diff(type, par) % plot difference of temperature profile from moist adiabat
    % plot_dtdz_binned_r1(type, par) % plot temperature profiles at specific latitudes
    % plot_thetaeq_zon_select(type, par) % plot eq pot temp profiles at specific latitudes
    % plot_dtdz_zon_select(type, par) % plot temperature profiles at specific latitudes
    % plot_ga_malr_si_diff(type, par) % plot difference of temperature profile from moist adiabat in sigma
    % plot_temp_zon(type, par) % plot temperature profiles at specific latitudes

end % select which functions to run at a time
function choose_plots_si_bl(type, par)
    % plot_ga_malr_diff_mon_lat(type, par) % plot ga diff
    % plot_ga_malr_diff_lon_lat(type, par) % plot ga diff
end
function choose_plots_ep(type, par)
    % plot_energy_lat(type, par); % plot all energy fluxes vs latitude a la Fig. 6.1 in Hartmann (2016)
    % plot_energy_lat_comp(type, par); % plot all energy fluxes vs latitude a la Fig. 6.1 in Hartmann (2016)
    % plot_r1z_lat(type, par); % compare r1 line plot with ERA5
    % plot_flux(type, par) % plot various energy fluxes in mon x lat and lon x lat space
    % plot_flux_comp(type, par) % plot various energy fluxes in mon x lat and lon x lat space
    % plot_temp(type, par) % plot temperature profiles
    % plot_temp_ann(type, par) % plot temperature profiles
    plot_dr1_midlatitude_line(type, par) % plot decomposition of R1 in mon at specific latitudes
    % plot_dr1_polar_line(type, par) % plot decomposition of R1 in mon at specific latitudes

    % plot_dr1_so_line(type, par) % plot decomposition of R1 in mon at specific latitudes
    % plot_dr1_polar_line_repl(type, par) % plot decomposition of R1 in mon at specific latitudes
    % plot_dr1_polar_line_topocomp(type, par) % plot decomposition of R1 in mon at specific latitudes
end % select which ep-functions to run at a time
