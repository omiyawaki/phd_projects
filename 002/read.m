clc; clear variables; close all;

addpath(genpath('/project2/tas1/miyawaki/matlab'));
addpath(genpath('./matlab'));

gcm_info
echam_info

%% set parameters
if 1
% par.erai.yr_span = '2000_2012'; % spanning years for ERA-Interim
% par.erai.yr_span = '1979_2018'; % spanning years for ERA-Interim
par.erai.yr_span = '1979_2005'; % spanning years for ERA-Interim
par.era5.yr_span = '1979_2005'; % spanning years for ERA5
% par.era5c.yr_span_list = {'1979_2005','2000_2018'}; % spanning years for ERA5
par.era5c.yr_span = par.era5.yr_span; % spanning years for ERA5
par.merra2.yr_span = '1980_2005'; % spanning years for MERRA2
par.jra55.yr_span = '1979_2005'; % spanning years for JRA-55
par.gcm.yr_span = 30; % number of years that I am considering in the GCM climatology
% par.echam_clims = par.echam.noice_mld; %{'echr0001'}; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), echr0001 (AMIP), echr0023 (AMIP no elevation), or rp000*** (various mixed layer depth and with/without sea ice)
par.echam_clims = {'rp000086'}; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), echr0001 (AMIP), echr0023 (AMIP no elevation), or rp000*** (various mixed layer depth and with/without sea ice)
par.hahn_clims = {'Control1850'}; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), echr0001 (AMIP), echr0023 (AMIP no elevation), or rp000*** (various mixed layer depth and with/without sea ice)
par.ceres.yr_span = '200003-201802'; % spanning years for CERES data
par.era.vars.rad = {'ssr', 'str', 'tsr', 'ttr'}; % radiation variables to read
par.era.vars.radcs = {'ssrc', 'strc', 'tsrc', 'ttrc'}; % radiation variables to read
par.era.vars.hydro = {'cp', 'lsp', 'e'}; % radiation variables to read
par.era.vars.div = {'p85.162', 'p84.162', 'p83.162'}; % radiation variables to read
par.era.vars.div_txt = {'divg', 'divq', 'divt'}; % radiation variables to read
par.era.vars.stf = {'sshf', 'slhf'}; % surface turbulent flux variables to read
par.era.vars.vert = {'t'}; % 3d variables to read (t = temp)
par.era.vars.srfc = {'sp', 't2m', 'd2m', 'zs'}; % surface variables to read (sp = surface pressure, t2m = 2m temp, d2m = 2m dew point temp)
par.era.vars.tend = {'p62.162'}; % 3d variables to read (t = temp)
par.era.vars.tend_txt = {'tend'}; % 3d variables to read (t = temp)
par.hahn.vars.rad = {'FLNT', 'FLNS', 'FSNT', 'FSNS'}; % radiation variables to read
par.hahn.vars.hydro = {'PRECC', 'PRECL', 'PRECSC', 'PRECSL'}; % hydrology variables
par.hahn.vars.stf = {'SHFLX', 'LHFLX'};
par.hahn.vars.vert = {'T'};
par.hahn.vars.srfc = {'PS', 'TREFHT', 'TS', 'zs'};
par.merra2.vars.rad = {'SWTNT', 'SWGNT', 'LWTUP', 'LWGNT'}; % radiation variables to read
par.merra2.vars.hydro = {'PRECTOT', 'PRECCON', 'EVAP'}; % hydrology variables
par.merra2.vars.stf = {'HFLUX', 'EFLUX'};
par.merra2.vars.vert = {'T'};
par.merra2.vars.srfc = {'PS', 'T2M', 'QV2M', 'zs'};
% par.jra55.vars.rad = {'dswrf', 'uswrf', 'dlwrf', 'ulwrf'}; % radiation variables to read
par.jra55.vars.rad = {'rsus', 'rsds', 'rlus', 'rlds', 'rsdt', 'rsut', 'rlut'}; % radiation variables to read
par.jra55.vars.hydro = {'pr', 'prc', 'evspsbl'}; % hydrology variables
% par.jra55.vars.stf = {'lhtfl', 'shtfl'};
par.jra55.vars.stf = {'hfss', 'hfls'}; % surface turbulent flux variables to read
par.jra55.vars.vert = {'ta'};
par.jra55.vars.srfc = {'ps', 'tas', 'hurs', 'zs'};
par.gcm.vars.rad = {'rsus', 'rsds', 'rlus', 'rlds', 'rsdt', 'rsut', 'rlut'}; % radiation variables to read
par.gcm.vars.radcs = {'rsuscs', 'rsdscs', 'rldscs', 'rsutcs', 'rlutcs'}; % radiation variables to read
par.gcm.vars.hydro = {'pr', 'evspsbl'}; % radiation variables to read
par.gcm.vars.stf = {'hfss', 'hfls'}; % surface turbulent flux variables to read
par.gcm.vars.vert = {'ta'}; % 3d variables to read (removed va)
par.gcm.vars.srfc = {'ps', 'ts', 'tas', 'hurs', 'zs'}; % surface variables to read
par.echam.vars.rad = {'srads', 'trads', 'srad0', 'trad0', 'tradsu', 'sradsu'}; % radiation variables to read
par.echam.vars.radcs = {'srafs', 'trafs', 'sraf0', 'traf0'}; % radiation variables to read
par.echam.vars.hydro = {'aprc', 'aprl', 'evap'}; % radiation variables to read
par.echam.vars.stf = {'ahfl', 'ahfs'}; % surface turbulent flux variables to read
par.echam.vars.vert = {'t', 'v'}; % 3d variables to read
par.echam.vars.srfc = {'aps', 'tsurf', 'temp2', 'dew2'}; % surface variables to read
par.ceres.vars.rad = {'sfc_net_sw_all_mon', 'sfc_net_lw_all_mon', 'toa_sw_all_mon', 'solar_mon', 'toa_lw_all_mon'}; % radiation variables to read
par.ceres.vars.rad_txt = {'ssr', 'str', 'tsur', 'tsdr', 'ttr'}; % radiation variables to read
% standard p coordinate for interpolation
par.pa = 1e2*linspace(1000,10,100);
% low res grid
par.pa_lo = 1e2*[1000 925 850 775 700 600 500 400 300 250 200 150 100 70 50 30 20 10 7 5 3 2 1 0.5 0.2 0.1];
% standard z coordinate for interpolation
par.z = [0:500:40e3]';
par.z_hires = linspace(0,par.z(end),1001); % high resolution grid for computing tropopause
par.si = linspace(1,1e-2,1e2);
% useful constants
par.cpd = 1005.7; par.Rd = 287; par.Rv = 461; par.L = 2.501e6; par.g = 9.81; par.a = 6357e3; par.eps = par.Rd/par.Rv;
end

% call functions
type='erai';
run_func(type, par);
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    %type='echam';
    %disp(par.echam.clim)
    %run_func(type, par);
end
for k=1:length(par.hahn_clims); par.hahn.clim=par.hahn_clims{k};
    %type='hahn';
    %disp(par.hahn.clim)
    %run_func(type, par);
end
for k=1:length(par.gcm_models); par.model=par.gcm_models{k};
    % type='gcm';
    % disp(par.model)
    % run_func(type, par);
end

function run_func(type, par)
    %read_grid(type, par) % grid, i.e. lon, lat, plev
    %read_rad(type, 'ymonmean', par) % radiation fluxes
    %read_hydro(type, 'ymonmean', par) % hydrological variables, e.g. precip, evap
    %read_stf(type, 'ymonmean', par) % surface turbulent fluxes
    %read_srfc(type, 'ymonmean', par) % other surface variables, e.g. 2-m temperature, surface pressure
    %read_tend(type, par) % mse tendency
    read_tempml(type, par); % read model level data and convert to standard sigma coord.
    %make_tempsi(type, par) % convert temp from plev to sigma
    %make_zgsi(type, par) % convert zg from plev to sigma
    %make_psi(type, par) % compute plev in si coords
    %read_lfrac(type, par) % land fraction (%)
    
    % read_rad(type, 'mon', par) % radiation fluxes
    % read_hydro(type, 'mon', par) % hydrological variables, e.g. precip, evap
    % read_stf(type, 'mon', par) % surface turbulent fluxes
    % read_srfc(type, 'mon', par) % other surface variables, e.g. 2-m temperature, surface pressure

    % read_radcs(type, 'ymonmean', par) % clear sky radiation fluxes
    % make_tempz(type, par) % convert temp from plev to z
    % make_pz(type, par) % compute plev in z coords
    % read_orog(type, par) % orography (m)
    % for varname = {'siced', 'friac', 'ahfres', 'ahfliac', 'ahfllac', 'ahflwac'}; % varnames: siced, friac, ahfres, ahfliac, ahfllac, ahflwac
    %     read_echam(varname{1}, type, par);
    % end
    % read_alb(type, par) % surface albedo (1)

    % make_thetaeqsi(type, par) % convert temp from plev to sigma
    % read_div(type, par) % divergence terms to calculate MSE flux divergence
    % read_dondiv79(type, par) % compute mse flux divergence from Donohoe mse transport, only for ERA-I data
    % read_dondiv00(type, par) % compute mse flux divergence from Donohoe mse transport, only for ERA-I data

    % make_tempsi_from_ml(type, par) % convert temp from ml to sigma
end

