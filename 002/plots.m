clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));
addpath(genpath('./matlab'));

gcm_info
echam_info
hahn_info

figure_params

%% set parameters
% lat grid type
if 1
% par.echam_clims = {'rp000135', 'rp000134'}; % par.echam.all_mld; % par.echam.sel; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
% par.echam_clims = {'rp000134', 'rp000135'}; % par.echam.all_mld; % par.echam.sel; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
par.echam_clims = {'rp000144'}; % par.echam.all_mld; % par.echam.sel; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
% par.echam_clims = {"rp000046",... % 50 m
%                        "rp000149",... % 45 m
%                        "rp000135",... % 40 m
%                        "rp000147",... % 35 m
%                        "rp000131",... % 30 m
%                        "rp000145",... % 25 m
%                        "rp000133",... % 20 m
%                        "rp000141",... % 15 m
%                        "rp000034",... % 10 m
%                        "rp000086",... % 5 m
%                        "rp000172"}; % 3 m

% par.echam_clims = {            "rp000126",... % 50 m
%                                  "rp000148",... % 45 m
%                                  "rp000134",... % 40 m
%                                  "rp000146",... % 35 m
%                                  "rp000130",... % 30 m
%                                  "rp000144",... % 25 m
%                                  "rp000132",... % 20 m
%                                  "rp000140",... % 15 m
%                                  "rp000124"};   % 10 m
% par.echam_clims = par.echam.noice_mld; % {'echr0001'}; % par.echam.sel; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
par.hahn_clims = {'Flat1850', 'Control1850'}; % Control1850, Flat1850, Control2xCO2, Flat2xCO2
par.rea.yr_span = '1980_2005';
par.erai.yr_span = '1980_2005';
par.era5.yr_span = '1980_2005';
par.era5c.yr_span = '1980_2005';
par.jra55.yr_span = '1980_2005';
par.merra2c.yr_span = '1980_2005';
par.levtype = 'pl'; % analyze model level (ml) or pressure level (pl) data?
par.ep_swp = 0.1; %[0.25 0.3 0.35]; % threshold value for determining RCE
par.ga_swp = 0.9; % threshold for determining RAE
par.si_eval = [0.8 0.85 0.9]; % sigma level for calculating inversion strength
par.ma = 1; % plot moist adiabat?
par.frz = 0; % consider effects of freezing on moist adiabat, computation of saturation vapor pressure, etc
par.ma_init = 0.95; % 'surf' for initializing with 2 m data, otherwise enter starting sigma level for moist adiabat
par.pa_eval = 500e2; % pressure level for calculating ma_diff
% par.si_bl_swp = [0.85 0.9 0.95]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_bl_swp = [0.7 0.9]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_up_list = [0.3]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
% par.si_up_list = [0.1]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.ta_thresh = 6.5; % criteria for temperature difference in K of plotting contour for closeness to moist adiabat
par.ga_thresh = 15; % criteria for lapse rate difference in % of plotting contour for closeness to moist adiabat
par.ga_bl_thresh = 100; % criteria for lapse rate difference in % of plotting contour for closeness to moist adiabat
par.inv_thresh = -4; % criteria for temperature difference in K of plotting contour for inversion strength
par.albedo_thresh = 0.6; % criteria for surface albedo of plotting contour for inversion strength
par.sn_thresh = 1e-1; % criteria for surface albedo of plotting contour for inversion strength
par.alb_thresh = 0.8; % criteria for surface albedo of plotting contour for inversion strength
par.albcs_thresh = 0.8; % criteria for surface albedo of plotting contour for inversion strength
par.r1_bins = [-0.55:0.1:1.5]; % bins for sorting temperature profiles according to r1 values
par.r1_bins_45 = [-0.55:0.05:1.5]; % bins for sorting temperature profiles according to r1 values
par.r1_bins_hl = [0.8-0.025/2:0.025:1.5+0.025/2]; % bins for sorting temperature profiles according to r1 values
% par.era.fw = {'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt', 'div79'};
% par.era.fw = {'div79', 'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt'};
% par.era.fw = {'ceresrad'};
par.land_list = {'lo'};
% par.rea.fw = {'mse', 'mse_old'};
par.rea.fw = {'mse', 'mse_old'};
par.era.fw = {'mse', 'mse_old'};
par.era5c.fw = {'mse', 'mse_old'};
% par.era5c.fw = {'mse', 'mse_old', 'mse_lat'};
par.jra55.fw = {'mse', 'mse_old'};
par.merra2c.fw = {'mse', 'mse_old'};
par.gcm.fw = {'mse', 'mse_old'};
par.echam.fw = {'mse', 'mse_old'};
par.hahn.fw = {'mse_old'};
par.pa = linspace(1000,10,100)*1e2; % high resolution vertical grid to interpolate to
par.z = [0:500:40e3]';
par.make_tikz = 0; % save figures as latex tikz files?
% set how to close energy budget
% if == teten, use TETEN data from Donohoe to close energy budget (option for ERA-Interim)
% if == stf, use SH and LH data from ERA-Interim to close energy budget (option for ERA-Interim)
% if == era5, use ERA5 radiative cooling and surface turbulent fluxes to close energy budget (only option for ERA5)
par.closure = 'era5';
par.do_surf = 0; % whether or not to calculate temperature profile in pressure grids including ts and ps data
par.cpd = 1005.7; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.Ls = 2.834e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
end

%% call functions
% plot_rad_lat(par)
% plot_rad_lon_lat(par)
% plot_tediv_lat(par)

type = 'era5c'; par.lat_interp = 'native';
% type = 'rea'; par.lat_interp = '1.00';
% choose_plots(type, par);
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    % type='echam';
    % par.lat_interp = 'native';
    % disp(par.echam.clim)
    % choose_plots(type, par);
end
for k=1:length(par.hahn_clims); par.hahn.clim=par.hahn_clims{k};
    type='hahn';
    par.lat_interp = 'native';
    disp(par.hahn.clim)
    choose_plots(type, par);
end
for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
     % type = 'gcm';
     % disp(par.model)
     % choose_plots(type, par);
end

% sweep through various boundary layer heights
for i = 1:length(par.si_bl_swp); par.si_bl = par.si_bl_swp(i);
    for i = 1:length(par.si_up_list); par.si_up = par.si_up_list(i);
        % type = 'era5c'; par.lat_interp = 'native';
        % type = 'rea'; par.lat_interp = '1.00';
        % choose_plots_si_bl(type, par)
        for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
            % type='echam';
            % par.lat_interp = 'native';
            % disp(par.echam.clim)
            % choose_plots_si_bl(type, par);
        end
        for k=1:length(par.hahn_clims); par.hahn.clim=par.hahn_clims{k};
            % type='hahn';
            % par.lat_interp = 'native';
            % disp(par.hahn.clim)
            % choose_plots_si_bl(type, par);
        end
        for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
            % type = 'gcm';
            % disp(par.model)
            % choose_plots_si_bl(type, par)
        end
    end
end

% sweep through various threshold values
for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = par.ga_swp(i);
    % type = 'era5c'; par.lat_interp = 'native';
    % type = 'rea'; par.lat_interp = '1.00';
    % choose_plots_ep(type, par)
    for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
        % type='echam'; par.lat_interp = 'native';
        % disp(par.echam.clim)
        % choose_plots_ep(type, par);
    end
    for k=1:length(par.hahn_clims); par.hahn.clim=par.hahn_clims{k};
        % type='hahn'; par.lat_interp = 'native';
        % disp(par.hahn.clim)
        % choose_plots_ep(type, par);
    end
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        % type = 'gcm';
        % disp(par.model)
        % choose_plots_ep(type, par)
    end
end

function choose_plots(type, par)
    % plot_ga_zon_select(type, par) % plot temperature profiles at specific latitudes
    % plot_ga_frac_binned_r1(type, par) % plot lapse rate profiles binned by R1
    % plot_dtempsi_binned_r1(type, par) % plot temperature responses binned by R1

    % plot_dmse_midlatitude_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dmse_polar_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dmse_midlatitude_line_comp(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dmse_polar_line_comp(type, par) % plot decomposition of R1 in mon x lat and lon x lat space

    % plot_dlh_polar_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dlh_polar_line_basic(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dlh_polar_line_decomp(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dlh_polar_line_asymdecomp(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dlh_polar_line_so(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_srfc_polar_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_srfc_polar_line_asym(type, par) % plot decomposition of R1 in mon x lat and lon x lat space

    % plot_alb(type, par) % surface albedo
    % plot_alb_comp(type, par) % surface albedo
    % plot_alb_mld(type, par) % surface albedo
    plot_alb_hemi(type, par) % surface albedo
    % plot_sice(type, par) % sea ice depth
    % plot_sice_rea(type, par) % sea ice depth
    % plot_siced(type, par) % sea ice depth

end % select which functions to run at a time
function choose_plots_si_bl(type, par)
    % plot_ga_malr_diff_mon_lat(type, par) % plot ga diff
    % plot_ga_malr_diff_lon_lat(type, par) % plot ga diff

    % plot_r1_ga_lat_line(type, par) % plot ga diff

    % plot_r1_ga_midlatitude_line(type, par) % plot R1 and lapse rate deviation together at specific latitudes
    plot_r1_ga_polar_line(type, par) % plot R1 and lapse rate deviation together at specific latitudes
end
function choose_plots_ep(type, par)
    % plot_energy_lat(type, par); % plot all energy fluxes vs latitude a la Fig. 6.1 in Hartmann (2016)
    % plot_r1z_lat(type, par); % compare r1 line plot with ERA5
    % plot_flux(type, par) % plot various energy fluxes in mon x lat and lon x lat space
    % plot_flux_comp(type, par) % plot various energy fluxes in mon x lat and lon x lat space

    % plot_ga_frac_ann(type, par) % plot temperature profiles

    % plot_ga_frac_midlatitude(type, par)
    % plot_ga_frac_polar(type, par)

    % plot_dr1_midlatitude_line(type, par) % plot decomposition of R1 in mon at specific latitudes
    % plot_dr1_polar_line(type, par) % plot decomposition of R1 in mon at specific latitudes

    % plot_dlh_polar_asymline_basic(type, par) % plot decomposition of R1 in mon x lat and lon x lat space

    % plot_dr1_comp_midlatitude_line(type, par) % plot decomposition of R1 in mon at specific latitudes
    % plot_dr1_comp_polar_line(type, par) % plot decomposition of R1 in mon at specific latitudes

    % plot_dr1_so_line(type, par) % plot decomposition of R1 in mon at specific latitudes
    % plot_dr1_polar_line_repl(type, par) % plot decomposition of R1 in mon at specific latitudes
    % plot_dr1_polar_line_topocomp(type, par) % plot decomposition of R1 in mon at specific latitudes

    plot_hahn_comp(type, par) % plot decomposition of R1 in mon at specific latitudes
end % select which ep-functions to run at a time
