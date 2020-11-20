clc; clear variables; close all;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

gcm_info
echam_info

%% set parameters
if 1
% par.erai.yr_span = '2000_2012'; % spanning years for ERA-Interim
% par.erai.yr_span = '1979_2018'; % spanning years for ERA-Interim
par.erai.yr_span = '2000_2018'; % spanning years for ERA-Interim
par.era5.yr_span = '2000_2018'; % spanning years for ERA5
par.merra2.yr_span = '2000_2018'; % spanning years for MERRA2
par.gcm.yr_span = 30; % number of years that I am considering in the GCM climatology
par.gcm.clim = 'piControl'; % choose from piControl or abrupt4xCO2 or historical
par.echam_clims = par.echam.ice_mld; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
par.ceres.yr_span = '200003-201802'; % spanning years for CERES data
par.era.vars.rad = {'ssr', 'str', 'tsr', 'ttr'}; % radiation variables to read
par.era.vars.hydro = {'cp', 'lsp', 'e'}; % radiation variables to read
par.era.vars.div = {'p85.162', 'p84.162', 'p83.162'}; % radiation variables to read
par.era.vars.div_txt = {'divg', 'divq', 'divt'}; % radiation variables to read
par.era.vars.stf = {'sshf', 'slhf'}; % surface turbulent flux variables to read
par.era.vars.vert = {'t'}; % 3d variables to read (t = temp)
par.era.vars.srfc = {'sp', 't2m', 'd2m', 'zs'}; % surface variables to read (sp = surface pressure, t2m = 2m temp, d2m = 2m dew point temp)
par.era.vars.tend = {'p62.162'}; % 3d variables to read (t = temp)
par.era.vars.tend_txt = {'tend'}; % 3d variables to read (t = temp)
par.merra2.vars.rad = {'SWTNT', 'SWGNT', 'LWTUP', 'LWGNT'}; % radiation variables to read
par.merra2.vars.hydro = {'PRECTOT', 'PRECCON', 'EVAP'}; % hydrology variables
par.merra2.vars.stf = {'HFLUX', 'EFLUX'};
par.merra2.vars.vert = {'T'};
par.merra2.vars.srfc = {'PS', 'T2M', 'QV2M', 'zs'};
par.gcm.vars.rad = {'rsus', 'rsds', 'rlus', 'rlds', 'rsdt', 'rsut', 'rlut'}; % radiation variables to read
par.gcm.vars.radcs = {'rsuscs', 'rsdscs', 'rldscs', 'rsutcs', 'rlutcs'}; % radiation variables to read
par.gcm.vars.hydro = {'prc', 'pr', 'evspsbl'}; % radiation variables to read
par.gcm.vars.stf = {'hfss', 'hfls'}; % surface turbulent flux variables to read
par.gcm.vars.vert = {'ta'}; % 3d variables to read (removed va)
par.gcm.vars.srfc = {'ps', 'ts', 'tas', 'hurs', 'zs'}; % surface variables to read
par.echam.vars.rad = {'srads', 'trads', 'srad0', 'trad0'}; % radiation variables to read
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

%% call functions
type='era5';
% run_func(type, par);
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    type='echam';
    % disp(par.echam.clim)
    % run_func(type, par);
end
for k=1:length(par.gcm_models); par.model=par.gcm_models{k};
    type='gcm';
    disp(par.model)
    run_func(type, par);
end

function run_func(type, par)
    % read_grid(type, par) % grid, i.e. lon, lat, plev
    % read_rad(type, par) % radiation fluxes
    % read_hydro(type, par) % hydrological variables, e.g. precip, evap
    % read_stf(type, par) % surface turbulent fluxes
    % read_srfc(type, par) % other surface variables, e.g. 2-m temperature, surface pressure
    % read_lfrac(type, par) % land fraction (%)
    % read_orog(type, par) % orography (m)
    % read_siced(type, par) % sea ice thickness (m)
    % make_tempsi(type, par) % convert temp from plev to sigma
    make_thetaeqsi(type, par) % convert temp from plev to sigma
    % make_psi(type, par) % compute plev in si coords
    % make_zgsi(type, par) % convert zg from plev to sigma

    % read_radcs(type, par) % clear sky radiation fluxes
    % make_tempz(type, par) % convert temp from plev to z
    % make_pz(type, par) % compute plev in z coords

    % read_div(type, par) % divergence terms to calculate MSE flux divergence
    % read_tend(type, par) % mse tendency, only for ERA data
    % read_dondiv79(type, par) % compute mse flux divergence from Donohoe mse transport, only for ERA-I data
    % read_dondiv00(type, par) % compute mse flux divergence from Donohoe mse transport, only for ERA-I data

    % make_tempsi_from_ml(type, par) % convert temp from ml to sigma
end

%% define functions
function read_grid(type, par)
    filename = 'grid.mat';
    % read data net SW and LW radiation data downloaded from Era5
    % first read lon and lat vectors since this is different from the Donohoe grid
    if any(strcmp(type, {'era5', 'erai'}))
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'longitude'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'latitude'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 10^2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'level')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type), 'grid')
    elseif strcmp(type, 'merra2')
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lon'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lat'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 1e2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lev')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid')
    elseif strcmp(type, 'gcm')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.nc', par.model, 'tas', par.model, par.gcm.clim));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.nc', par.model, 'ta', par.model, par.gcm.clim));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=ncread(fullpath.dim2, 'lon');
        grid.dim3.lon=ncread(fullpath.dim3, 'lon');
        grid.dim2.lat=ncread(fullpath.dim2, 'lat');
        grid.dim3.lat=ncread(fullpath.dim3, 'lat');
        grid.dim3.plev=ncread(fullpath.dim3, 'plev');
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'echam')
        if contains(par.echam.clim, 'rp000')
            file.dim2=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            file.dim3=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT_rjg_%s_*.ymonmean.nc', par.echam.clim));
            file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.plev=double(ncread(fullpath.dim3, 'lev'));
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'echam_ml')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_*.ymonmean.nc'));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.a=double(ncread(fullpath.dim3, 'hyam'));
        grid.dim3.b=double(ncread(fullpath.dim3, 'hybm'));
        grid.dim3.z = par.z;
        grid.dim3.plev = par.pa;
        % grid.dim3.si = 1e-5*([grid.dim3.a+grid.dim3.b*1e5; 1e5]);
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_ml');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'echam_pl')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_*.ymonmean.nc'));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.plev=double(ncread(fullpath.dim3, 'lev'));
        grid.dim3.z = par.z;
        % grid.dim3.si = 1e-5*grid.dim3.plev;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_pl');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'ceres')
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s.ymonmean.nc', type, par.(type).yr_span), 'lon'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s.ymonmean.nc', type, par.(type).yr_span), 'lat'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type), 'grid')
    end

    % save grid
end
function read_rad(type, par)
    if strcmp(type, 'erai') | strcmp(type, 'era5')
        rad_vars=par.era.vars.rad;
        for i=1:length(rad_vars)
            % dimensions are (lon x lat x time)
            % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
            rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), rad_vars{i}));
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            rad.(rad_vars{i}) = rad.(rad_vars{i})/86400;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/rad.mat', type), 'rad', 'rad_vars');
        if strcmp(type, 'era5')
            for i=1:length(rad_vars)
                rad.(rad_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.erai.yr_span), rad_vars{i});
                rad.(rad_vars{i}) = rad.(rad_vars{i})/86400;
            end
            save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/rad_2000_2012.mat', type), 'rad', 'rad_vars');
        end
    elseif strcmp(type, 'merra2')
        rad_vars=par.merra2.vars.rad;
        for i=1:length(rad_vars)
            % dimensions are (lon x lat x time)
            rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), rad_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/rad.mat', type), 'rad', 'rad_vars');
    elseif strcmp(type, 'gcm')
        rad_vars=par.gcm.vars.rad;
        for i=1:length(par.gcm.vars.rad); var = par.gcm.vars.rad{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            rad.(var)=ncread(fullpath, var);
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='rad.mat';
        save(sprintf('%s/%s', newdir, filename), 'rad', 'rad_vars');
    elseif strcmp(type, 'echam')
        rad_vars=par.echam.vars.rad;
        for i=1:length(par.echam.vars.rad); var = par.echam.vars.rad{i};
            if contains(par.echam.clim, 'rp000')
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            else
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT_rjg_%s_*.ymonmean.nc', par.echam.clim));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            rad.(var)=double(ncread(fullpath, var));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s/rad.mat', par.echam.clim), 'rad', 'rad_vars');
    elseif strcmp(type, 'ceres')
        rad_vars = par.ceres.vars.rad;
        rad_vars_txt = par.ceres.vars.rad_txt;
        for i=1:length(rad_vars)
            rad.(rad_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s.ymonmean.nc', type, par.(type).yr_span), rad_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/rad.mat', type), 'rad');

        for i=1:length(rad_vars)
            rad.(rad_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_200101-200912.ymonmean.nc', type), rad_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/rad_2001_2009.mat', type), 'rad');
    end
end
function read_radcs(type, par)
    if strcmp(type, 'erai') | strcmp(type, 'era5')
        radcs_vars=par.era.vars.radcs;
        for i=1:length(radcs_vars)
            % dimensions are (lon x lat x time)
            % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
            radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_radcs_%s.ymonmean.nc', type, type, par.(type).yr_span), radcs_vars{i}));
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            radcs.(radcs_vars{i}) = radcs.(radcs_vars{i})/86400;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/radcs.mat', type), 'radcs', 'radcs_vars');
        if strcmp(type, 'era5')
            for i=1:length(radcs_vars)
                radcs.(radcs_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_radcs_%s.ymonmean.nc', type, type, par.erai.yr_span), radcs_vars{i});
                radcs.(radcs_vars{i}) = radcs.(radcs_vars{i})/86400;
            end
            save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/radcs_2000_2012.mat', type), 'radcs', 'radcs_vars');
        end
    elseif strcmp(type, 'gcm')
        radcs_vars=par.gcm.vars.radcs;
        for i=1:length(par.gcm.vars.radcs); var = par.gcm.vars.radcs{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            radcs.(var)=ncread(fullpath, var);
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='radcs.mat';
        save(sprintf('%s/%s', newdir, filename), 'radcs', 'radcs_vars');
    elseif strcmp(type, 'echam')
        radcs_vars=par.echam.vars.radcs;
        for i=1:length(par.echam.vars.radcs); var = par.echam.vars.radcs{i};
            if contains(par.echam.clim, 'rp000')
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            else
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT_rjg_%s_*.ymonmean.nc', par.echam.clim));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            radcs.(var)=double(ncread(fullpath, var));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='radcs.mat';
        save(sprintf('%s/%s', newdir, filename), 'radcs', 'radcs_vars');
    elseif strcmp(type, 'ceres')
        radcs_vars = par.ceres.vars.radcs;
        radcs_vars_txt = par.ceres.vars.radcs_txt;
        for i=1:length(radcs_vars)
            radcs.(radcs_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s.ymonmean.nc', type, par.(type).yr_span), radcs_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/radcs.mat', type), 'radcs');

        for i=1:length(radcs_vars)
            radcs.(radcs_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_200101-200912.ymonmean.nc', type), radcs_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/radcs_2001_2009.mat', type), 'radcs');
    end
end
function read_hydro(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        hydro_vars=par.era.vars.hydro;
        for i=1:length(hydro_vars)
            % dimensions are (lon x lat x time)
            hydro.(hydro_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hydro/%s_hydro_%s.ymonmean.nc', type, type, par.(type).yr_span), hydro_vars{i}));
            % the data is originally reported as m (depth) hydror day, so
            % divide by 86400 s and multiply by 1000 kg/m^3 to get the
            % conventional kg/m^2/s mass flux over the full day
            hydro.(hydro_vars{i}) = hydro.(hydro_vars{i})/86400*1e3;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/hydro.mat', type), 'hydro', 'hydro_vars');

    elseif strcmp(type, 'merra2')
        hydro_vars=par.merra2.vars.hydro;
        for i=1:length(hydro_vars)
            % dimensions are (lon x lat x time)
            hydro.(hydro_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hydro/%s_hydro_%s.ymonmean.nc', type, type, par.(type).yr_span), hydro_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/hydro.mat', type), 'hydro', 'hydro_vars');

    elseif strcmp(type, 'gcm')
        hydro_vars=par.gcm.vars.hydro;
        for i=1:length(par.gcm.vars.hydro); var = par.gcm.vars.hydro{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            hydro.(var)=ncread(fullpath, var);
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='hydro.mat';
        save(sprintf('%s/%s', newdir, filename), 'hydro', 'hydro_vars');
    elseif strcmp(type, 'echam')
        hydro_vars=par.echam.vars.hydro;
        for i=1:length(par.echam.vars.hydro); var = par.echam.vars.hydro{i};
            if contains(par.echam.clim, 'rp000')
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            else
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT_rjg_%s_*.ymonmean.nc', par.echam.clim));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            hydro.(var)=double(ncread(fullpath, var));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='hydro.mat';
        save(sprintf('%s/%s', newdir, filename), 'hydro', 'hydro_vars');
    end
end
function read_div(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        div_vars=par.era.vars.div;
        div_vars_txt=par.era.vars.div_txt;
        for i=1:length(div_vars)
            % dimensions are (lon x lat x time)
            div.(div_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/div/%s_div_%s.ymonmean.nc', type, type, par.(type).yr_span), div_vars{i}));
            div.(div_vars_txt{i}) = div.(div_vars_txt{i});
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/div.mat', type), 'div', 'div_vars_txt');

    else
        error('Divergence data are only available for ERA.');
    end
end
function read_stf(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        stf_vars=par.era.vars.stf;
        for i=1:length(stf_vars)
            % dimensions are (lon x lat x time)
            stf.(stf_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/stf/%s_stf_%s.ymonmean.nc', type, type, par.(type).yr_span), stf_vars{i}));
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            stf.(stf_vars{i}) = stf.(stf_vars{i})/86400;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/stf.mat', type), 'stf', 'stf_vars');

        if strcmp(type, 'era5')
            for i=1:length(stf_vars)
                stf.(stf_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/stf/%s_stf_%s.ymonmean.nc', type, type, par.erai.yr_span), stf_vars{i});
                stf.(stf_vars{i}) = stf.(stf_vars{i})/86400;
            end
            save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/stf_2000_2012.mat', type), 'stf', 'stf_vars');
        end
    elseif strcmp(type, 'merra2')
        stf_vars=par.merra2.vars.stf;
        for i=1:length(stf_vars)
            % dimensions are (lon x lat x time)
            stf.(stf_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/stf/%s_stf_%s.ymonmean.nc', type, type, par.(type).yr_span), stf_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/stf.mat', type), 'stf', 'stf_vars');
    elseif strcmp(type, 'gcm')
        stf_vars=par.gcm.vars.stf;
        for i=1:length(par.gcm.vars.stf); var = par.gcm.vars.stf{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            stf.(var)=ncread(fullpath, var);
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='stf.mat';
        save(sprintf('%s/%s', newdir, filename), 'stf', 'stf_vars');
    elseif strcmp(type, 'echam')
        stf_vars=par.echam.vars.stf;
        for i=1:length(par.echam.vars.stf); var = par.echam.vars.stf{i};
            if contains(par.echam.clim, 'rp000')
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            else
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT_rjg_%s_*.ymonmean.nc', par.echam.clim));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            stf.(var)=double(ncread(fullpath, var));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='stf.mat';
        save(sprintf('%s/%s', newdir, filename), 'stf', 'stf_vars');
    end
end
function read_srfc(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        srfc_vars=par.era.vars.srfc;
        for i=1:length(srfc_vars); var = srfc_vars{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/srfc/%s_srfc_%s.ymonmean.nc', type, type, par.(type).yr_span));
            fullpath=sprintf('%s/%s', file.folder, file.name);

            if strcmp(var, 'zs'); % use orography as surface geopotential if data exists
                if strcmp(type, 'erai')
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/orog/interim_%s.nc', type, 'orog'));
                elseif strcmp(type, 'era5')
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/orog/%s_%s_%s.ymonmean.nc', type, type, 'orog', par.(type).yr_span));
                end
                fullpath=sprintf('%s/%s', file.folder, file.name);
                if exist(fullpath, 'file')
                    srfc.(var)=double(ncread(fullpath, 'z')); srfc.(var) = srfc.(var)/par.g; % divide by g to get height
                    if strcmp(type, 'erai')
                        srfc.(var)=repmat(srfc.(var), [1 1 12]);
                    end
                else % create surface geopotential height using surface pressure data
                    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
                    load(sprintf('%s/grid.mat', prefix)); % read grid data
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
                    fullpath=sprintf('%s/%s', file.folder, file.name);
                    zg = ncread(fullpath, 'z');
                    zg = permute(zg, [3 1 2 4]);
                    pb=CmdLineProgressBar("Calculating zs..."); % track progress of this loop
                    for lo = 1:length(grid.dim2.lon)
                        pb.print(lo, length(grid.dim2.lon));
                        for la = 1:length(grid.dim2.lat)
                            for mo = 1:12
                                srfc.zs(lo,la,mo) = interp1(grid.dim3.plev, zg(:,lo,la,mo), srfc.sp(lo,la,mo), 'linear', 'extrap');
                            end
                        end
                    end
                end
            else
                srfc.(var) = double(squeeze(ncread(fullpath, srfc_vars{i})));
            end

        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/srfc.mat', type), 'srfc', 'srfc_vars');

    elseif strcmp(type, 'merra2')
        srfc_vars=par.merra2.vars.srfc;
        for i=1:length(srfc_vars); var = srfc_vars{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/srfc/%s_srfc_%s.ymonmean.nc', type, type, par.(type).yr_span));
            fullpath=sprintf('%s/%s', file.folder, file.name);

            if strcmp(var, 'zs'); % use orography as surface geopotential if data exists
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/orog/MERRA2_*.nc4', type));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                if exist(fullpath, 'file')
                    srfc.(var) = double(ncread(fullpath, 'PHIS')); srfc.(var)=srfc.(var)/par.g; % convert geopotential to height
                    srfc.(var) = repmat(srfc.(var), [1 1 12]);
                else % create surface geopotential height using surface pressure data
                    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
                    load(sprintf('%s/grid.mat', prefix)); % read grid data
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
                    fullpath=sprintf('%s/%s', file.folder, file.name);
                    zg = ncread(fullpath, 'H');
                    zg = permute(zg, [3 1 2 4]);
                    pb=CmdLineProgressBar("Calculating zs..."); % track progress of this loop
                    for lo = 1:length(grid.dim2.lon)
                        pb.print(lo, length(grid.dim2.lon));
                        for la = 1:length(grid.dim2.lat)
                            for mo = 1:12
                                srfc.zs(lo,la,mo) = interp1(grid.dim3.plev, zg(:,lo,la,mo), srfc.PS(lo,la,mo), 'linear', 'extrap');
                            end
                        end
                    end
                end
            else
                srfc.(var) = double(squeeze(ncread(fullpath, srfc_vars{i})));
            end

        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/srfc.mat', type), 'srfc', 'srfc_vars');

    elseif strcmp(type, 'gcm')
        srfc_vars=par.gcm.vars.srfc;
        for i=1:length(par.gcm.vars.srfc); var = par.gcm.vars.srfc{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            if ~exist(fullpath)
                if strcmp(var, 'hurs')
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.nc', par.model, 'hur', par.model, par.gcm.clim));
                    fullpath=sprintf('%s/%s', file.folder, file.name);
                    hur=ncread(fullpath, 'hur');
                    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s/grid.mat', par.model, par.gcm.clim));
                    if ~isequal(grid.dim2.lon, grid.dim3.lon); hur=interp1(grid.dim3.lon, hur, grid.dim2.lon); end; % interpolate to 2D lon if different from 3D
                    hur=permute(hur, [2 1 3 4]); % bring lat to first dim
                    if ~isequal(grid.dim2.lat, grid.dim3.lat); hur=interp1(grid.dim3.lat, hur, grid.dim2.lat); end;
                    hur=permute(hur, [3 2 1 4]); % bring plev to first dim
                    pb=CmdLineProgressBar("Calculating hurs..."); % track progress of this loop
                    for id_lon=1:length(grid.dim2.lon)
                        pb.print(id_lon, length(grid.dim2.lon));
                        for id_lat=1:length(grid.dim2.lat)
                            for id_time=1:size(srfc.ps, 3)
                                srfc.hurs(id_lon, id_lat, id_time)=interp1(grid.dim3.plev, hur(:,id_lon,id_lat,id_time), srfc.ps(id_lon, id_lat, id_time), 'linear', 'extrap');
                            end
                        end
                    end
                elseif strcmp(var, 'zs')
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s_raw/%s/%s_*.nc', par.gcm.clim, par.model, 'orog'));
                    fullpath=sprintf('%s/%s', file.folder, file.name);
                    if exist(fullpath, 'file')
                        srfc.(var)=ncread(fullpath, 'orog');
                        srfc.(var)=repmat(srfc.(var), [1 1 12]);
                    else
                        % create surface geopotential height using surface pressure data
                        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
                        load(sprintf('%s/grid.mat', prefix)); % read grid data
                        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, 'zg', par.model, par.gcm.clim));
                        fullpath=sprintf('%s/%s', file.folder, file.name);
                        zg = ncread(fullpath, 'zg');
                        zg = permute(zg, [3 1 2 4]);
                        pb=CmdLineProgressBar("Calculating zs..."); % track progress of this loop
                        for lo = 1:length(grid.dim2.lon)
                            pb.print(lo, length(grid.dim2.lon));
                            for la = 1:length(grid.dim2.lat)
                                for mo = 1:12
                                    srfc.zs(lo,la,mo) = interp1(grid.dim3.plev, zg(:,lo,la,mo), srfc.ps(lo,la,mo), 'linear', 'extrap');
                                end
                            end
                        end
                    end
                else
                    error(sprintf('The file for variable %s does not exist. Check in the raw data folder to see if you forgot to download the file.'))
                end
            else
                srfc.(var)=ncread(fullpath, var);

            end
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='srfc.mat';
        save(sprintf('%s/%s', newdir, filename), 'srfc', 'srfc_vars');
    elseif strcmp(type, 'echam')
        srfc_vars=par.echam.vars.srfc;
        for i=1:length(par.echam.vars.srfc); var = par.echam.vars.srfc{i};
            if contains(par.echam.clim, 'rp000')
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            else
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT_rjg_%s_*.ymonmean.nc', par.echam.clim));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            srfc.(var)=double(ncread(fullpath, var));

            if strcmp(var, 'aps'); % create surface geopotential height using surface pressure data
                prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
                load(sprintf('%s/grid.mat', prefix)); % read grid data
                if contains(par.echam.clim, 'rp000')
                    file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
                else
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
                end
                fullpath=sprintf('%s/%s', file.folder, file.name);
                zg = double(ncread(fullpath, 'geopoth'));
                zg = permute(zg, [3 1 2 4]);
                pb=CmdLineProgressBar("Calculating zs..."); % track progress of this loop
                for lo = 1:length(grid.dim2.lon)
                    pb.print(lo, length(grid.dim2.lon));
                    for la = 1:length(grid.dim2.lat)
                        for mo = 1:12
                            srfc.zs(lo,la,mo) = interp1(grid.dim3.plev, zg(:,lo,la,mo), srfc.aps(lo,la,mo), 'linear', 'extrap');
                        end
                    end
                end
            end

        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='srfc.mat';
        save(sprintf('%s/%s', newdir, filename), 'srfc', 'srfc_vars');
    elseif strcmp(type, 'echam_ml')
        srfc_vars=par.echam.vars.srfc;
        for i=1:length(par.echam.vars.srfc); var = par.echam.vars.srfc{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_*.ymonmean.nc'));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            srfc.(var)=double(ncread(fullpath, var));

            if strcmp(var, 'aps'); % create surface geopotential height using surface pressure data
                prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
                load(sprintf('%s/grid.mat', prefix)); % read grid data
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                zg = double(ncread(fullpath, 'geopoth'));
                srfc.zs(:,:,:) = squeeze(zg(:,:,1,:));
            end

        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_ml');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='srfc.mat';
        save(sprintf('%s/%s', newdir, filename), 'srfc', 'srfc_vars');
    elseif strcmp(type, 'echam_pl')
        srfc_vars=par.echam.vars.srfc;
        for i=1:length(par.echam.vars.srfc); var = par.echam.vars.srfc{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_*.ymonmean.nc'));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            srfc.(var)=double(ncread(fullpath, var));

            if strcmp(var, 'aps'); % create surface geopotential height using surface pressure data
                prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
                load(sprintf('%s/grid.mat', prefix)); % read grid data
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                zg = double(ncread(fullpath, 'geopoth'));
                srfc.zs(:,:,:) = squeeze(zg(:,:,1,:));
            end

        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_pl');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='srfc.mat';
        save(sprintf('%s/%s', newdir, filename), 'srfc', 'srfc_vars');
    end
end
function read_lfrac(type, par)
% land fraction
    if strcmp(type, 'erai')
        sftlf = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/lmask/interim_lmask.nc', type), 'lsm'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/sftlf.mat', type), 'sftlf');
    elseif strcmp(type, 'gcm')
        if any(strcmp(par.gcm.clim, {'piControl', 'abrupt4xCO2'}))
            file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sftlf_*.nc', par.model));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s_raw/%s/sftlf_*.nc', par.gcm.clim, par.model));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sftlf=ncread(fullpath, 'sftlf');
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sftlf.mat';
        save(sprintf('%s/%s', newdir, filename), 'sftlf');
    elseif strcmp(type, 'echam')
        file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sftlf_*.nc', 'MPI-ESM-LR'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sftlf=ncread(fullpath, 'sftlf');
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sftlf.mat';
        save(sprintf('%s/%s', newdir, filename), 'sftlf');
    end
end
function read_tend(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        tend_vars=par.era.vars.tend;
        tend_vars_txt=par.era.vars.tend_txt;
        for i=1:length(tend_vars)
            % dimensions are (lon x lat x time)
            tend.(tend_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/tend/%s_tend_%s.deltat.ymonmean.nc', type, type, par.(type).yr_span), tend_vars{i}));
            % the data is originally reported as J m^-2, so
            % divide by 6 hr (because data is 6 hourly) to
            % convert to W m^-2.
            tend.(tend_vars_txt{i}) = tend.(tend_vars_txt{i})/(6*3600);
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/tend.mat', type), 'tend', 'tend_vars_txt');

    else
        error('MSE tendency data are read only for ERA data.');
    end
end
function read_dondiv79(type, par)
% from 1979-10 to 2018-09
    if strcmp(type, 'erai')

        rootdir = "/project2/tas1/miyawaki/projects/002/data/raw/don/ERA_MHT"; % root directory of Donohoe MSE transport data
        means = load(sprintf('%s/means/1979_10means.mat', rootdir));
        lat = means.lat;
        latr = deg2rad(lat);
        dlat = latr(2)-latr(1);
        clat = cos(latr); clat(1)=nan; clat(end)=nan;

        startdate = datetime(1979,10,1);
        enddate = datetime(2018,9,1);
        dates = datetime(startdate:calmonths(1):enddate, 'format', 'yyyy_M');
        nfiles = length(dates);
        pb=CmdLineProgressBar("Reading Donohoe heat transport data..."); % track progress of this loop
        div_orig = nan([floor(nfiles/13)+1 length(lat) 12]);
        for ifile = 1:nfiles
            pb.print(ifile, nfiles);
            trans_orig = load(sprintf('%s/heat_transport/%s_heattrans.mat', rootdir, dates(ifile)));

            imonth = month(dates(ifile));
            iyear = year(dates(ifile)) - year(startdate) + 1;
            fmse = trans_orig.MME + trans_orig.SE + trans_orig.TE;
            % div_arg = fmse'.*clat;
            % div_orig(year,:,month) = 1./(2*pi*par.a^2*clat.^2).*gradient(div_arg, dlat);
            div_arg = fmse';
            div_orig(iyear,:,imonth) = 1./(2*pi*par.a^2*clat).*gradient(div_arg, dlat);
            % figure; clf; hold all;
            % plot(lat, fmse, 'k');
            % print(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/test_folder/fmse_test', type), '-dpng', '-r300')
            % figure; clf; hold all;
            % plot(lat, div_orig(year,:,month), 'k');
            % print(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/test_folder/test', type), '-dpng', '-r300')
            % clear fmse div_arg trans_orig
            % return
        end

        % div_orig = circshift(div_orig, -3, 3); % shift months so that January is the first entry (note that Donohoe ERA-I begins on Oct 1979)

        dondiv = squeeze(nanmean(div_orig, 1)); % take climatology
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/dondiv79.mat', type), 'dondiv', 'lat');

    else
        error('Donohoe MSE transport data are available only for ERA-I data.');
    end
end
function read_dondiv00(type, par)
% from 2000-03 to 2018-02
    if strcmp(type, 'erai')

        rootdir = "/project2/tas1/miyawaki/projects/002/data/raw/don/ERA_MHT"; % root directory of Donohoe MSE transport data
        means = load(sprintf('%s/means/1979_10means.mat', rootdir));
        lat = means.lat;
        latr = deg2rad(lat);
        dlat = latr(2)-latr(1);
        clat = cos(latr); clat(1)=nan; clat(end)=nan;

        startdate = datetime(2000,3,1);
        enddate = datetime(2018,2,1);
        dates = datetime(startdate:calmonths(1):enddate, 'format', 'yyyy_M');
        nfiles = length(dates);
        pb=CmdLineProgressBar("Reading Donohoe heat transport data..."); % track progress of this loop
        div_orig = nan([floor(nfiles/13)+1 length(lat) 12]);
        for ifile = 1:nfiles
            pb.print(ifile, nfiles);
            trans_orig = load(sprintf('%s/heat_transport/%s_heattrans.mat', rootdir, dates(ifile)));

            imonth = month(dates(ifile));
            iyear = year(dates(ifile)) - year(startdate) + 1;
            fmse = trans_orig.MME + trans_orig.SE + trans_orig.TE;
            % div_arg = fmse'.*clat;
            % div_orig(year,:,month) = 1./(2*pi*par.a^2*clat.^2).*gradient(div_arg, dlat);
            div_arg = fmse';
            div_orig(iyear,:,imonth) = 1./(2*pi*par.a^2*clat).*gradient(div_arg, dlat);
            % figure; clf; hold all;
            % plot(lat, fmse, 'k');
            % print(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/test_folder/fmse_test', type), '-dpng', '-r300')
            % figure; clf; hold all;
            % plot(lat, div_orig(iyear,:,imonth), 'k');
            % print(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/test_folder/test', type), '-dpng', '-r300')
            % clear fmse div_arg trans_orig
            % return
        end

        % div_orig = circshift(div_orig, -3, 3); % shift months so that January is the first entry (note that Donohoe ERA-I begins on Oct 1979)

        dondiv = squeeze(nanmean(div_orig, 1)); % take climatology
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/dondiv00.mat', type), 'dondiv', 'lat');

    else
        error('Donohoe MSE transport data are available only for ERA-I data.');
    end
end
function read_orog(type, par)
% orography
    if any(strcmp(type, {'era5', 'erai'}))
        % dimensions are (lon x lat x time)
        % albedo = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/albedo/%s_albedo_%s.ymonmean.nc', type, type, par.(type).yr_span), 'fal'));
        % save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/albedo.mat', type), 'albedo');
    elseif strcmp(type, 'gcm')
        for i=1:length(par.gcm.vars.rad); var = par.gcm.vars.rad{i};
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
            file=dir(sprintf('/project2/tas1/CMIP5_%s/%s/orog_*.nc', par.gcm.clim, par.model));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            orog=ncread(fullpath, 'orog');
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='orog.mat';
        save(sprintf('%s/%s', newdir, filename), 'orog');
    elseif strcmp(type, 'echam')
        file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/orog_*.nc', 'MPI-ESM-LR'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        orog=ncread(fullpath, 'orog');
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='orog.mat';
        save(sprintf('%s/%s', newdir, filename), 'orog');
    end
end
function read_siced(type, par)
% read sea ice depth
    if any(strcmp(type, {'era5', 'erai'})) % only for MPI-ESM-LR
        % dimensions are (lon x lat x time)
        sn = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sn/%s_sn_%s.ymonmean.nc', type, type, par.(type).yr_span), 'fal'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/sn.mat', type), 'sn');
    elseif strcmp(type, 'gcm') & strcmp(par.model, 'MPI-ESM-LR') % only for MPI-ESM-LR
        for i=1:length(par.gcm.vars.rad); var = par.gcm.vars.rad{i};
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
            file=dir(sprintf('/project2/tas1/echam6/piControl/BOT_piControl_r1i1p1-P_1900_1949.nc'));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            sn=ncread(fullpath, 'siced');
            lat=ncread(fullpath, 'lat');
            lon=ncread(fullpath, 'lon');
        end
        load(sprintf('%s/grid.mat', prefix)); % read grid data
        sn = interp1(lon, sn, grid.dim2.lon); % interpolate to MPI-ESM-LR longitude grid
        sn = permute(sn, [2 1 3]); % bring lat to 1st
        sn = interp1(lat, sn, grid.dim2.lat); % interpolate to MPI-ESM-LR latitude grid
        sn = permute(sn, [2 1 3]); % bring lat back to 2nd
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sn.mat';
        save(sprintf('%s/%s', newdir, filename), 'sn');
    elseif strcmp(type, 'echam')
        file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        siced=double(ncread(fullpath, 'siced'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='siced.mat';
        save(sprintf('%s/%s', newdir, filename), 'siced');
    end
end

function make_tempz(type, par)
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_ml')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data

    temp = permute(temp, [3 1 2 4]);
    zg = permute(zg, [3 1 2 4]);

    pb=CmdLineProgressBar("Calculating tempz..."); % track progress of this loop
    for lo = 1:length(grid.dim3.lon)
        pb.print(lo, length(grid.dim3.lon));
        for la = 1:length(grid.dim3.lat)
            for mo = 1:12
                % only keep nonnan data and interpolate
                notnan = find(~isnan(zg(:,lo,la,mo)));

                tempz(:,lo,la,mo) = interp1(zg(notnan,lo,la,mo), temp(notnan,lo,la,mo), grid.dim3.z);
            end
        end
    end

    tempz = permute(tempz, [2 3 1 4]); % reorder to lon x lat x z x mon

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif strcmp(type, 'echam_ml'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_ml');
    elseif strcmp(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_pl'); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='tempz.mat';
    save(sprintf('%s/%s', newdir, filename), 'tempz', '-v7.3');
end
function make_tempsi(type, par)
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 'T');
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif contains(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    elseif contains(type, 'echam_pl')
        ps_vert = repmat(srfc.aps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    sm = nan(size(ta_orig));
    sm(pa < 0.9961*ps_vert) = 1;
    ta_sm = ta_orig.*sm; % filter ta with surface mask

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    ta_sm = permute(ta_sm, [3 1 2 4]);

    pb = CmdLineProgressBar("Sorting and interpolating temperature to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);

                % add surface data
                if strcmp(type, 'era5') | strcmp(type, 'erai')
                    tmp(end) = srfc.t2m(lo,la,mo); % surface is the last element in era plev
                elseif strcmp(type, 'merra2')
                    tmp(1) = srfc.T2M(lo,la,mo);
                elseif strcmp(type, 'gcm')
                    tmp(1) = srfc.tas(lo,la,mo);
                elseif contains(type, 'echam')
                    tmp(1) = srfc.temp2(lo,la,mo);
                end

                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(tmp)));

                % ta_si.lin(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'linear', nan);
                % ta_si.cub(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'pchip', nan);
                ta_si.spl(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);
                % ta_si.mak(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'makima', nan);

                clear tmp

            end
        end
    end

    % ta_si.lin = permute(ta_si.lin, [2 3 1 4]); % reorder to lon x lat x si x mon
    % ta_si.cub = permute(ta_si.cub, [2 3 1 4]); % reorder to lon x lat x si x mon
    ta_si.spl = permute(ta_si.spl, [2 3 1 4]); % reorder to lon x lat x si x mon
    % ta_si.mak = permute(ta_si.mak, [2 3 1 4]); % reorder to lon x lat x si x mon

    % ta_si = permute(ta_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif contains(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ta_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'ta_si', '-v7.3');
end
function make_tempsi_from_ml(type, par) % model level to sigma temperature
    if strcmp(type, 'echam_ml')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
        var = 'aps';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ps_orig = double(ncread(fullpath, var));
    else
        error('This code only works for data output in the model vertical grid.')
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    % compute sigma from a and b
    ps_vert = repmat(ps_orig, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
    ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    a = permute(repmat(grid.dim3.a, [1 size(ps_orig)]), [2 3 1 4]);
    b = permute(repmat(grid.dim3.b, [1 size(ps_orig)]), [2 3 1 4]);
    si = a./ps_vert + b;

    % interpolate to standard sigma levels
    si = permute(si, [3 1 2 4]); % bring si front
    ta_orig = permute(ta_orig, [3 1 2 4]);
    pb = CmdLineProgressBar("Interpolating temperature to new standard grid...");
    for lo=1:size(si,2)
        pb.print(lo, size(si,2));
        for la=1:size(si,3)
            for mo=1:size(si,4)
                tmp_ta = ta_orig(:,lo,la,mo);
                tmp_ta(end+1) = squeeze(srfc.temp2(lo,la,mo));

                tmp_si = si(:,lo,la,mo);
                tmp_si(end+1) = 1;

                ta_si.lin(:,lo,la,mo) = interp1(tmp_si, tmp_ta, grid.dim3.si, 'linear', nan);
                ta_si.cub(:,lo,la,mo) = interp1(tmp_si, tmp_ta, grid.dim3.si, 'pchip', nan);
                ta_si.spl(:,lo,la,mo) = interp1(tmp_si, tmp_ta, grid.dim3.si, 'spline', nan);
                ta_si.mak(:,lo,la,mo) = interp1(tmp_si, tmp_ta, grid.dim3.si, 'makima', nan);

                clear tmp_ta tmp_si

            end
        end
    end

    ta_si.lin = permute(ta_si.lin, [2 3 1 4]); % reorder to lon x lat x si x mon
    ta_si.cub = permute(ta_si.cub, [2 3 1 4]); % reorder to lon x lat x si x mon
    ta_si.spl = permute(ta_si.spl, [2 3 1 4]); % reorder to lon x lat x si x mon
    ta_si.mak = permute(ta_si.mak, [2 3 1 4]); % reorder to lon x lat x si x mon

    % ta_si = permute(ta_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    if strcmp(type, 'echam_ml'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_ml'); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ta_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'ta_si', '-v7.3');
end
function make_zgsi(type, par)
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = ncread(fullpath, 'z');
    elseif contains(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = ncread(fullpath, 'H');
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 'geopoth';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = double(ncread(fullpath, var));
    elseif contains(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(zg_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(zg_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(zg_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif contains(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(zg_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    sm = nan(size(zg_orig));
    sm(pa < ps_vert) = 1;
    zg_sm = zg_orig.*sm; % filter zg with surface mask

    % % add tsurf dazg and interpolate to higher resolution vertical grid
    % [pa_plus zg_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    % pa_plus(:,:,1:end-1,:) = pa; % populate with szgndard pressure grid
    % zg_plus(:,:,1:end-1,:) = zg_sm; % populate with standard geopotential
    % pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    % zg_plus(:,:,end,:) = srfc.zs(:,:,:); % add surface height data
    % pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    % zg_plus = permute(zg_plus, [3 1 2 4]); % bring plev dimension to front
    % [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    % zgi_sm = nan(length(par.si), size(pa, 1), size(pa, 2), size(pa, 4));
    % pb = CmdLineProgressBar("Sorting and interpolating zg to new standard grid...");
    % for lo=1:size(pa_plus,2)
    %     pb.print(lo, size(pa_plus,2));
    %     for la=1:size(pa_plus,3)
    %         for mo=1:size(pa_plus,4)
    %             zg_plus(:,lo,la,mo) = zg_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort zg (has to be in loop because sort_index works for vector calls only)
    %             zgsi(:,lo,la,mo) = interp1(pa_plus(:,lo,la,mo)/ps_vert(lo,la,1,mo), zg_plus(:,lo,la,mo), grid.dim3.si);
    %         end
    %     end
    % end
    % clear pa_plus zg_plus; % clear unneeded variables

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    zg_sm = permute(zg_sm, [3 1 2 4]);

    pb = CmdLineProgressBar("Sorting and interpolating temperature to new szgndard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);

                if any(strcmp(type, {'era5', 'erai'}))
                    tmp(end) = srfc.zs(lo,la,mo); % add surface height
                else
                    tmp(1) = srfc.zs(lo,la,mo); % add surface height
                end

                notnan = find(~isnan(squeeze(tmp))); % only keep nonnan data and redo interpolation

                zgsi(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);

                clear tmp

            end
        end
    end

    zgsi = permute(zgsi, [2 3 1 4]); % reorder to lon x lat x si x mon

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif contains(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='zgsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'zgsi', '-v7.3');
end
function make_pz(type, par)
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 'geopoth';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
        var = 'aps';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ps_orig = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data

    if strcmp(type, 'echam_ml')
        % compute sigma from a and b
        ps_vert = repmat(ps_orig, [1 1 1 size(zg, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        a = permute(repmat(grid.dim3.a, [1 size(ps_orig)]), [2 3 1 4]);
        b = permute(repmat(grid.dim3.b, [1 size(ps_orig)]), [2 3 1 4]);
        plev = a + b.*ps_vert;
    else
        plev = grid.dim3.plev;
    end

    zg = permute(zg, [3 1 2 4]);
    plev = permute(plev, [3 1 2 4]);

    pb=CmdLineProgressBar("Calculating pz..."); % track progress of this loop
    for lo = 1:length(grid.dim3.lon)
        pb.print(lo, length(grid.dim3.lon));
        for la = 1:length(grid.dim3.lat)
            for mo = 1:12
                if strcmp(type, 'echam_ml')
                    pz(:,lo,la,mo) = interp1(zg(:,lo,la,mo), plev(:,lo,la,mo), grid.dim3.z, 'linear', 'extrap');
                else
                    % only keep nonnan data and interpolate
                    notnan = find(~isnan(zg(:,lo,la,mo)));

                    pz(:,lo,la,mo) = interp1(zg(notnan,lo,la,mo), grid.dim3.plev(notnan), grid.dim3.z, 'linear', 'extrap');
                end
            end
        end
    end

    pz = permute(pz, [2 3 1 4]); % reorder to lon x lat x z x mon

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='pz.mat';
    save(sprintf('%s/%s', newdir, filename), 'pz', '-v7.3');
end
function make_psi(type, par)
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif contains(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 'aps';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ps_orig = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 length(grid.dim3.plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 length(grid.dim3.plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 length(grid.dim3.plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam_ml')
        % compute sigma from a and b
        ps_vert = repmat(ps_orig, [1 1 1 length(grid.dim3.a)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        a = permute(repmat(grid.dim3.a, [1 size(ps_orig)]), [2 3 1 4]);
        b = permute(repmat(grid.dim3.b, [1 size(ps_orig)]), [2 3 1 4]);
        pa = a + b.*ps_vert;
    elseif contains(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 length(grid.dim3.plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    sm = nan(size(pa));
    sm(pa < ps_vert) = 1;
    pa_sm = pa.*sm; % filter pa with surface mask

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    pa_sm = permute(pa_sm, [3 1 2 4]);

    pb = CmdLineProgressBar("Calculaing psi...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), pa_sm(:,lo,la,mo), 1e-5*grid.dim3.plev, 'linear');

                % add surface dapa
                if strcmp(type, 'era5') | strcmp(type, 'erai')
                    tmp(end) = srfc.sp(lo,la,mo);
                elseif strcmp(type, 'merra2')
                    tmp = srfc.PS(lo,la,mo);
                elseif strcmp(type, 'gcm')
                    tmp = srfc.ps(lo,la,mo);
                elseif contains(type, 'echam')
                    tmp = srfc.aps(lo,la,mo);
                end

                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(tmp)));
                pa_si(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);

            end
        end
    end

    pa_si = permute(pa_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml','echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='pa_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'pa_si', '-v7.3');
end
function make_thetaeqsi(type, par)
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 'T');
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
        var = 'hur';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        hur_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif contains(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    esat = calc_esat(ta_orig, 0); % compute saturation vapor pressure
    p = permute(repmat(grid.dim3.plev, [1 size(ta_orig,1) size(ta_orig,2) size(ta_orig,4)]), [2 3 1 4]);
    e = esat.*hur_orig/100;
    r = calc_r(p, e, par);
    pd = p - e; % partial pressure of dry air
    clear p e esat;

    % compute eq potential temperature following AMS glossary definition
    thetaeq = ta_orig .* (1e5./pd).^(par.Rd/par.cpd).*(hur_orig/100).^(-r*par.Rv/par.cpd).*exp(par.L*r./(par.cpd*ta_orig));
    clear ta_orig hur_orig;

    % surface eq pot temp
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        tas = srfc.t2m; % surface is the last element in era plev
        ps = srfc.sp;
    elseif strcmp(type, 'merra2')
        tas = srfc.T2M;
        ps = srfc.PS;
    elseif strcmp(type, 'gcm')
        tas = srfc.tas;
        ps = srfc.ps;
        hurs = srfc.hurs;
    elseif contains(type, 'echam')
        tas = srfc.temp2;
    end

    esats = calc_esat(tas, 0); % compute saturation vapor pressure
    es = esats.*hurs/100;
    rs = calc_r(ps, es, par);
    pds = ps - es; % partial pressure of dry air
    clear ps es esats;

    thetaeqs = tas .* (1e5./pds).^(par.Rd/par.cpd).*(hurs/100).^(-rs*par.Rv/par.cpd).*exp(par.L*rs./(par.cpd*tas));
    clear tas hurs;

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    elseif contains(type, 'echam_pl')
        ps_vert = repmat(srfc.aps, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    sm = nan(size(thetaeq));
    sm(pa < 0.9961*ps_vert) = 1;
    thetaeq_sm = thetaeq.*sm; % filter thetaeq with surface mask

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    thetaeq_sm = permute(thetaeq_sm, [3 1 2 4]);

    pb = CmdLineProgressBar("Sorting and interpolating equivalent potential temperature to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), thetaeq_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);

                % add surface data
                tmp(1) = thetaeqs(lo,la,mo);

                % only keep nonnan dathetaeq and redo interpolation
                notnan = find(~isnan(squeeze(tmp)));

                % thetaeq_si.lin(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'linear', nan);
                % thetaeq_si.cub(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'pchip', nan);
                thetaeq_si.spl(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);
                % thetaeq_si.mak(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'makima', nan);

                clear tmp

            end
        end
    end

    % thetaeq_si.lin = permute(thetaeq_si.lin, [2 3 1 4]); % reorder to lon x lat x si x mon
    % thetaeq_si.cub = permute(thetaeq_si.cub, [2 3 1 4]); % reorder to lon x lat x si x mon
    thetaeq_si.spl = permute(thetaeq_si.spl, [2 3 1 4]); % reorder to lon x lat x si x mon
    % thetaeq_si.mak = permute(thetaeq_si.mak, [2 3 1 4]); % reorder to lon x lat x si x mon

    % thetaeq_si = permute(thetaeq_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif conthetaeqins(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='thetaeq_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'thetaeq_si', '-v7.3');
end
