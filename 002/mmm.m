% This script imports all CMIP data and creates the multi-model mean
clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));
addpath(genpath('./matlab'));

mmm_info
raw_varnames

%% set parameters
if 1
par.rea.yr_span = '1980_2005'; % spanning years for ERA-Interim
par.erai.yr_span = '1980_2005'; % spanning years for ERA-Interim
par.era5.yr_span = '1980_2005'; % spanning years for ERA5
par.jra55.yr_span = '1980_2005'; % spanning years for JRA55
par.gcm.yr_span = '198001-200512'; % spanning years for JRA55
par.era5c.yr_span = par.era5.yr_span;
par.merra2c.yr_span = '1980_2005'; % spanning years for MERRA2
par.outname = 'mmm'; % use mmm for standard multimodel mean (all models available) and mmm_subset where subset=name of climate where you are taking the mmm of the subset of models in common with the subset climate
par.lat_interp = '1.00'; % latitudinal grid spacing to interpolate to (deg)
par.lat = -90:str2num(par.lat_interp):90; % define standard latitude grid for 'std' interpolation
par.lon_interp = '1.00'; % longitudinal grid spacing to interpolone to (deg)
par.lon = 0:str2num(par.lon_interp):360; % define standard lonitude grid for 'std' interpolonion
par.ep_swp = 0.1; %[0.25 0.3 0.35]; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ga_swp = 0.9; % optional threshold for RAE. If undefined, the default value is 1-par.ep
par.ma_init = 0.95; % initial starting level for moist adiabat ('surf' = start with dry adiabat until LCL or enter sigma level for starting saturated adiabat)
par.pa = linspace(1000,10,100)*1e2; % high resolution vertical grid to interpolate to
par.si = 1e-5*par.pa; % high resolution vertical grid to interpolate to
% par.si_eval = [0.8 0.85 0.9]; % sigma level for evaluating inversion strength (T(si_eval) - T(surface))
% par.si_bl_swp = [0.85 0.9 0.95]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_eval = [0.9]; % sigma level for evaluating inversion strength (T(si_eval) - T(surface))
par.si_bl_swp = [0.7 0.8 0.9]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_up_list = [0.3]; % sigma level for upper boundary of vertical average for close to moist adiabatic
% par.si_up_list = [0.1]; % sigma level for upper boundary of vertical average for close to moist adiabatic
par.z = [0:500:40e3]';
par.land_list = {'lo'};
par.frz = 0;
par.rea.fw = {'mse', 'mse_old'};
par.gcm.fw = {'mse', 'mse_old'};
par.cpd = 1005.7; par.cpv = 1870; par.cpl = 4186; par.cpi = 2108; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
end

type = 'gcm';
par.model_list = par.gcm_models;
par.clim = par.gcm.clim;

% type = 'rea';
% par.model_list = {'era5c', 'jra55', 'merra2c'};
% par.clim = '1980_2005';

make_grid(type, par);
% choose_mmm_raw_2d(type, par);
choose_mmm_lat_mon_lev(type, par); % make mmm of mon x lat x lev data
% choose_mmm_lat_mon(type, par); % make mmm of mon x lat data
% choose_mmm_lon_lat(type, par); % make mmm of lon x lat data
% choose_mmm_lat(type, par); % make mmm of lat data
% choose_mmm_mon(type, par); % make mmm of lat data
for i=1:length(par.si_bl_swp); par.si_bl = par.si_bl_swp(i);
    for i=1:length(par.si_up_list); par.si_up = par.si_up_list(i);
        % choose_mmm_lat_mon_bl(type, par);
        % choose_mmm_lat_bl(type, par);
    end
end

function make_grid(type, par)
    filename = 'grid.mat';
    grid.dim2.lon = par.lon';
    grid.dim3.lon = par.lon';
    grid.dim2.lat = par.lat';
    grid.dim3.lat = par.lat';
    grid.dim3.plev = par.pa;
    grid.dim3.z = par.z;
    grid.dim3.si = 1e-5*par.pa;
    newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/%s', type, par.outname, par.clim, par.(type).yr_span);
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    save(sprintf('%s/%s', newdir, filename), 'grid');
    if strcmp(type, 'rea')
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid');
    end
end

%% Choose functions to run

function choose_mmm_raw_2d(type, par)
    % mmm_srfc(type, par);
    mmm_sfcWind(type, par);
    % mmm_e(type, par);
end

function choose_mmm_lat_mon_lev(type, par)
    % mmm_ta_mon_lat(type, par);
    % mmm_tai_mon_lat(type, par);
    % mmm_ta_pl_mon_lat(type, par);

    % mmm_dtempsi_mon_lat(type, par);

    % mmm_ma_mon_lat(type, par);
    mmm_ga_diff_mon_lat(type, par);
    mmm_ga_frac_mon_lat(type, par);
    mmm_gad_frac_mon_lat(type, par);

    % mmm_ga_frac_midlatitude(type, par);
    % mmm_ga_frac_polar(type, par);
end

function choose_mmm_lat_mon(type, par)
    mmm_flux_z(type, par);
    mmm_vh_mon(type, par);
end

function choose_mmm_lat_mon_bl(type, par)
    mmm_ga_dalr_bl_diff_si_mon_lat(type, par);
    % mmm_ga_malr_bl_diff_si_mon_lat(type, par);
    % mmm_ga_malr_diff_si_mon_lat(type, par);

    % mmm_ga_malr_diff_midlatitude_line(type, par);
    % mmm_ga_malr_bl_diff_polar_line(type, par);
end

function choose_mmm_lon_lat(type, par)
    mmm_flux_t(type, par);
end

function choose_mmm_lat(type, par)
    % mmm_flux_zt(type, par);
    % mmm_vh(type, par);
end

function choose_mmm_lat_bl(type, par)
    mmm_ga_dalr_bl_diff_si_lat(type, par);
    % mmm_ga_malr_bl_diff_si_lat(type, par);
    % mmm_ga_malr_diff_si_lat(type, par);
end

function choose_mmm_mon(type, par)
    mmm_dmse_midlatitude_line(type, par);
    mmm_dmse_polar_line(type, par);
    mmm_dr1_midlatitude_line(type, par);
    mmm_dr1_polar_line(type, par);
    mmm_dr2_midlatitude_line(type, par);
    mmm_dr2_polar_line(type, par);
end

