clc; clear variables; close all;

addpath(genpath('/project2/tas1/miyawaki/matlab'));
addpath(genpath('./matlab'));

gcm_info
echam_info
raw_varnames

%% set parameters
% par.erai.yr_span = '2000_2012'; % spanning years for ERA-Interim
% par.erai.yr_span = '1979_2018'; % spanning years for ERA-Interim
par.erai.yr_span = '1980_2005'; % spanning years for ERA-Interim
par.era5.yr_span = '1980_2005'; % spanning years for ERA5
% par.era5c.yr_span_list = {'1979_2005','2000_2018'}; % spanning years for ERA5
par.era5c.yr_span = par.era5.yr_span; % spanning years for ERA5
par.merra2c.yr_span = '1980_2005'; % spanning years for MERRA2
par.jra55.yr_span = '1980_2005'; % spanning years for JRA-55
% par.gcm.yr_span = '198001-200512'; % number of years that I am considering in the GCM climatology
% par.gcm.yr_span = 30; % number of years that I am considering in the GCM climatology
% par.echam_clims = par.echam.noice_mld; %{'echr0001'}; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), echr0001 (AMIP), echr0023 (AMIP no elevation), or rp000*** (various mixed layer depth and with/without sea ice)
% par.echam_clims = {'rp000126'}; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), echr0001 (AMIP), echr0023 (AMIP no elevation), or rp000*** (various mixed layer depth and with/without sea ice)

% par.echam_clims = {            "rp000126",... % 50 m
%                                  "rp000148",... % 45 m
%                                  "rp000134",... % 40 m
%                                  "rp000146",... % 35 m
%                                  "rp000130",... % 30 m
%                                  "rp000144",... % 25 m
%                                  "rp000132",... % 20 m
%                                  "rp000140",... % 15 m
%                                  "rp000124"};   % 10 m

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

% par.hahn_clims = {'Control1850'}; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), echr0001 (AMIP), echr0023 (AMIP no elevation), or rp000*** (various mixed layer depth and with/without sea ice)
par.hahn_clims = {'Flat1850', 'Control1850'}; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), echr0001 (AMIP), echr0023 (AMIP no elevation), or rp000*** (various mixed layer depth and with/without sea ice)
par.ceres.yr_span = '200003-201802'; % spanning years for CERES data
par.echam.exceptions = {'rp000092', 'rp000172'}; % ECHAM data that are not in the data11 archive yet
% standard p coordinate for interpolation
par.pa = 1e2*linspace(1000,10,100);
% low res grid
par.pa_lo = 1e2*[1000 925 850 775 700 600 500 400 300 250 200 150 100 70 50 30 20 10 7 5 3 2 1 0.5 0.2 0.1];
% standard z coordinate for interpolation
par.z = [0:500:40e3]';
par.z_hires = linspace(0,par.z(end),1001); % high resolution grid for computing tropopause
par.si = linspace(1,1e-2,1e2);

% consider freezing for moist processes?
par.frz = 0;

% useful constants
par.cpd = 1005.7; par.Rd = 287; par.Rv = 461; par.L = 2.501e6; par.g = 9.81; par.a = 6357e3; par.eps = par.Rd/par.Rv;

% call functions
par.rea_models = {'era5c'};
% par.rea_models = {'era5c', 'merra2c', 'jra55'};
for k=1:length(par.rea_models); type=par.rea_models{k};
    % run_func(type, par);
end
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    type='echam';
    disp(par.echam.clim)
    run_func(type, par);
end
for k=1:length(par.hahn_clims); par.hahn.clim=par.hahn_clims{k};
    % type='hahn';
    % disp(par.hahn.clim)
    % run_func(type, par);
end
for k=1:length(par.gcm_models); par.model=par.gcm_models{k};
    % type='gcm';
    % disp(par.model)
    % run_func(type, par);
end

function run_func(type, par)
    % read_grid(type, 'ymonmean', par)
    % read_srfc(type, 'ymonmean', par) % other surface variables, e.g. 2-m temperature, surface pressure
    % read_sfcWind(type, 'ymonmean', par) % other surface variables, e.g. 2-m temperature, surface pressure
    % read_hus_ml0(type, 'ymonmean', par) % other surface variables, e.g. 2-m temperature, surface pressure
    % read_tend(type, 'ymonmean', par) % mse tendency
    % make_tempsi(type, par) % convert temp from plev to sigma
    % make_dtempsi(type, par) % take temperature difference (e.g., rcp85 - historical)
    % make_zgsi(type, 'ymonmean', par) % convert zg from plev to sigma
    % make_psi(type, par) % compute plev in si coords
    % read_lfrac(type, par) % land fraction (%)
    % read_sice(type, par) % sea ice cover (1)
    % read_siced(type, par) % sea ice cover (1)
    read_albedo(type, par) % surface albedo (1)
    
    % read_rad(type, 'mon', par) % radiation fluxes
    % read_hydro(type, 'mon', par) % hydrological variables, e.g. precip, evap
    % read_stf(type, 'mon', par) % surface turbulent fluxes
    % read_srfc(type, 'mon', par) % other surface variables, e.g. 2-m temperature, surface pressure

    % read_radcs(type, 'ymonmean', par) % clear sky radiation fluxes
    % read_orog(type, par) % orography (m)
    % for varname = {'siced', 'friac', 'ahfres', 'ahfliac', 'ahfllac', 'ahflwac'}; % varnames: siced, friac, ahfres, ahfliac, ahfllac, ahflwac
    %     read_echam(varname{1}, type, par);
    % end

end

