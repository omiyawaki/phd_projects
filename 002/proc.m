clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));
addpath(genpath('./matlab'));

gcm_info
echam_info
hahn_info

%% set parameters
if 1
% par.erai.yr_span = '2000_2012'; % spanning years for ERA-Interim
par.erai.yr_span = '1980_2005'; % spanning years for ERA-Interim
par.era5.yr_span = '1980_2005'; % spanning years for ERA5
par.jra55.yr_span = '1980_2005'; % spanning years for JRA55
par.era5c.yr_span = par.era5.yr_span;
par.merra2c.yr_span = '1980_2005'; % spanning years for MERRA2
% par.echam_clims = {'rp000134', 'rp000135', 'rp000141'}; % par.echam.all_mld; % par.echam.sel; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
% par.echam_clims = {'rp000149'};
par.echam_clims = {"rp000046",... % 50 m
                       "rp000149",... % 45 m
                       "rp000135",... % 40 m
                       "rp000147",... % 35 m
                       "rp000131",... % 30 m
                       "rp000145",... % 25 m
                       "rp000133",... % 20 m
                       "rp000141",... % 15 m
                       "rp000034",... % 10 m
                       "rp000086",... % 5 m
                       "rp000172"}; % 3 m

% par.echam_clims = {            "rp000126",... % 50 m
%                                  "rp000148",... % 45 m
%                                  "rp000134",... % 40 m
%                                  "rp000146",... % 35 m
%                                  "rp000130",... % 30 m
%                                  "rp000144",... % 25 m
%                                  "rp000132",... % 20 m
%                                  "rp000140",... % 15 m
%                                  "rp000124"};   % 10 m
% par.echam_clims = par.echam.noice_mld; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
par.hahn_clims = {'Control1850'}; % Control1850, Flat1850, Control2xCO2, Flat2xCO2
% par.hahn_clims = {'Flat1850', 'Control1850'}; % Control1850, Flat1850, Control2xCO2, Flat2xCO2
par.lat_interp = 'native'; % which latitudinal grid to interpolate to: native (no interpolation), don (donohoe, coarse), era (native ERA-Interim, fine), or std (custom, very fine)
par.levtype = 'pl'; % analyze model level (ml) or pressure level (pl) data?
par.lat_std = transpose(-90:0.25:90); % define standard latitude grid for 'std' interpolation
par.ep_swp = 0.1; %[0.25 0.3 0.35]; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ga_swp = 0.9; % optional threshold for RAE. If undefined, the default value is 1-par.ep
par.ep_cp = 0.5; % additional flag for RCE definition using convective precipitation. RCE is defined as where lsp/cp < ep_cp
par.ma_type = 'std'; % choose the type of moist adiabat: reversible, pseudo, or std
par.ma_init = 0.95; % initial starting level for moist adiabat ('surf' = start with dry adiabat until LCL or enter sigma level for starting saturated adiabat)
par.frz = 0; % consider latent heat of fusion in moist adiabat?
par.pa_span = [1000 10]*100; % pressure range for calculating moist adiabat
par.pa = linspace(1000,10,100)*1e2; % high resolution vertical grid to interpolate to
par.si = 1e-5*par.pa; % high resolution vertical grid to interpolate to
par.dpa = -10; % pressure increment for integrating dry adiabat section of moist adiabat (matters for how accurately the LCL is computed)
par.z_span = [0 25]*10^3; % height range for calculating moist adiabat
par.dz = 10; % pressure increment for integrating dry adiabat section of moist adiabat (matters for how accurately the LCL is computed)
par.do_surf = 0; % whether or not to calculate temperature profile in pressure grids including ts and ps data
par.si_eval = [0.8 0.9]; % sigma level for evaluating inversion strength (T(si_eval) - T(surface))
par.si_bl_swp = [0.7 0.8 0.9]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_up_list = [0.3]; % sigma level for upper boundary of vertical average for close to moist adiabatic
% par.si_up_list = [0.1 0.2 0.3 0.4]; % sigma level for upper boundary of vertical average for close to moist adiabatic
% par.era.fw = {'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt', 'div79'};
% par.era.fw = {'div79', 'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt'};
par.land_list = {'lo'};
par.era.fw = {'mse', 'mse_old'};
% par.era.fw = {'mse_old'};
par.jra55.fw = {'mse', 'mse_old'};
par.merra2c.fw = {'mse', 'mse_old'};
par.gcm.fw = {'mse', 'mse_old'};
par.echam.fw = {'mse', 'mse_old'};
par.hahn.fw = {'mse_old'};
par.cpd = 1005.7; par.cpv = 1870; par.cpl = 4186; par.cpi = 2108; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
end

%% call functions
% comp_flux(par)
% ceres_flux(par)
% choose_disp(par)

% par.rea_models = {'era5c'};
par.rea_models = {'era5c', 'merra2c', 'jra55'};
for k=1:length(par.rea_models); type = par.rea_models{k};
    % choose_proc(type, par)
end
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    % type='echam';
    % disp(par.echam.clim)
    % choose_proc(type, par);
end
for k=1:length(par.hahn_clims); par.hahn.clim=par.hahn_clims{k};
    % type='hahn';
    % disp(par.hahn.clim)
    % choose_proc(type, par);
end
for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    disp(par.model)
    choose_proc(type, par)
end

for i=1:length(par.si_bl_swp); par.si_bl = par.si_bl_swp(i);
    for i=1:length(par.si_up_list); par.si_up = par.si_up_list(i);
        % par.rea_models = {'era5c'};
        par.rea_models = {'era5c', 'jra55', 'merra2c'};
        for k=1:length(par.rea_models); type = par.rea_models{k};
            % choose_proc_si_bl(type, par)
        end
        for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
            % type='echam';
            % disp(par.echam.clim)
            % choose_proc_si_bl(type, par);
        end
        for k=1:length(par.hahn_clims); par.hahn.clim=par.hahn_clims{k};
            % type='hahn';
            % disp(par.hahn.clim)
            % choose_proc_si_bl(type, par);
        end
        for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
            % type = 'gcm';
            % disp(par.model)
            % choose_proc_si_bl(type, par)
        end
    end
end

for i=1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = par.ga_swp(i);
    % type = 'jra55';
    % choose_proc_ep(type, par)
    for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
        % type='echam';
        % disp(par.echam.clim)
        % choose_proc_ep(type, par);
    end
    for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
        % type = 'gcm';
        % disp(par.model)
        % choose_proc_ep(type, par)
    end
end

function choose_proc(type, par)
    % proc_flux(type, par) % calculate energy fluxes in the vertically-integrated MSE budget using ERA-Interim data
    % proc_temp_mon_lat(type, par) % calculate mon x lat temperature profiles
    % proc_temp_pl_mon_lat(type, par) % calculate mon x lat temperature profiles
    make_tai(type, par) % calculate interpolated temperature profile in p coordinates with 2 m temp insert in lon x lat x mon
    % proc_tai_mon_lat(type, par) % calculate mon x lat temperature profiles
    % make_masi(type, par) % calculate moist adiabats at every lon x lat x mon
    % proc_ma_mon_lat(type, par) % calculate mon x lat moist adiabats
    make_dtdzsi(type, par) % calculate model lapse rate and interpolate to sigma coordinates
    % make_dtdzsi_alt(type, par) % calculate model lapse rate and interpolate to sigma coordinates
    make_malrsi(type, par) % calculate moist adiabatic lapse rate of model temperature sigma coordinates
    % save_mask(type, par) % save land and ocean masks once (faster than creating mask every time I need it)

    proc_ga_malr_mon_lat(type, par) % calculate mon x lat MALR profiles
    proc_ga_frac_mon_lat(type, par) % calculate mon x lat lapse rate deviation from a MALR profiles
    % proc_dtempsi_mon_lat(type, par) % calculate mon x lat temperature response profiles
    % proc_gad_frac_mon_lat(type, par) % calculate mon x lat lapse rate deviation from a DALR profiles

    % proc_thetaeq_mon_lat(type, par) % calculate mon x lat equivalent potential temperature profiles
    % make_dthedzsi(type, par) % calculate eq pot temp lapse rate and interpolate to sigma coordinates
    % proc_dthedz_mon_lat(type, par) % calculate mon x lat eq pot temp lapse rate profiles

    % make_dthedpasi(type, par) % calculate eq pot temp lapse rate and interpolate to sigma coordinates
    % proc_dthedpa_mon_lat(type, par) % calculate mon x lat eq pot temp lapse rate profiles

    % proc_mse_mon_lat(type, par) % calculate mon x lat equivalent potential temperature profiles
    % make_dmsedzsi(type, par) % calculate MSE lapse rate and interpolate to sigma coordinates
    % proc_dmsedz_mon_lat(type, par) % calculate mon x lat eq pot temp lapse rate profiles
    
    % proc_dmse_midlatitude_line(type, par);
    % proc_dmse_polar_line(type, par);
    % proc_dr1_midlatitude_line(type, par);
    % proc_dr1_polar_line(type, par);

    % proc_ga_frac_midlatitude(type, par);
    % proc_ga_frac_polar(type, par);

    % proc_temp_mon_lat_interp(type, par) % calculate mon x lat temperature profiles
    % proc_temp_mon_lat_interp_mean(type, par) % calculate mon x lat temperature profiles
    % proc_dtdz_mon_lat(type, par) % calculate mon x lat lapse rate
    % proc_malr_mon_lat(type, par) % calculate mon x lat moist adiabats
    % proc_net_flux(type, par) % calculate net energy fluxes at TOA and surface
end % select which functions to run at a time
function choose_proc_si_bl(type, par)
    % proc_ga_malr_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    % proc_ga_dalr_bl_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    proc_ga_malr_bl_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    % proc_dthedpa_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    % proc_ga_trop_malr_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference

    % proc_ga_malr_diff_midlatitude_line(type, par);
    proc_ga_malr_bl_diff_polar_line(type, par);
end
function choose_proc_ep(type, par)
    % proc_rcae(type, par) % calculate RCE and RAE regimes
    % proc_rcae_alt(type, par) % calculate RCE and RAE regimes (all divergence allowed for RCE)
    % proc_ta_si(type, par) % calculate RCE and RAE temperature profiles
    % proc_ma_si(type, par) % calculate moist adiabats corresponding to RCE profiles
end % select which ep-functions to run at a time
function choose_disp(par)
    % disp_global_rad(par)
    disp_global_stf(par)
end % display globally-averaged radiation fluxes
