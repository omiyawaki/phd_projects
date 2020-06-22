clc; clear variables; close all;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
par.erai.yr_span = 2000:2012; % spanning years for ERA-Interim
par.erai.yr_text = cellstr(num2str(par.erai.yr_span'))';
par.era5.yr_span = '1979_2019'; % spanning years for ERA5
par.gcm.yr_span = 50; % number of years that I am considering in the GCM climatology
par.era.vars.rad = {'ssr', 'str', 'tsr', 'ttr'}; % radiation variables to read
par.era.vars.pe = {'cp', 'lsp', 'e'}; % radiation variables to read
par.era.vars.stf = {'sshf', 'slhf'}; % surface turbulent flux variables to read
par.era.vars.vert = {'t'}; % 3d variables to read (t = temp)
par.era.vars.srfc = {'sp', 't2m', 'd2m'}; % surface variables to read (sp = surface pressure, t2m = 2m temp, d2m = 2m dew point temp)
par.gcm.vars.rad = {'rsus', 'rsds', 'rlus', 'rlds', 'rsdt', 'rsut', 'rlut'}; % radiation variables to read
par.gcm.vars.pe = {'prc', 'pr', 'evspsbl'}; % radiation variables to read
par.gcm.vars.stf = {'hfss', 'hfls'}; % surface turbulent flux variables to read
par.gcm.vars.vert = {'ta'}; % 3d variables to read (t = temp)
par.gcm.vars.srfc = {'ps', 'tas', 'hurs'}; % surface variables to read (sp = surface pressure, t2m = 2m temp, d2m = 2m dew point temp)
gcm_info
% useful constants
par.cpd = 1005.7; par.Rd = 287; par.L = 2.501e6; par.g = 9.81;

%% call functions
type='era5';
run_func(type, par);
for k=1:length(par.gcm_models); par.model=par.gcm_models{k};
    type='gcm';
    run_func(type, par);
end

%% define functions
function run_func(type, par)
    read_grid(type, par)
    read_rad(type, par)
    read_pe(type, par)
    read_stf(type, par)
    read_srfc(type, par)
end
function read_grid(type, par)
    % read data net SW and LW radiation data downloaded from Era5
    % first read lon and lat vectors since this is different from the Donohoe grid
    if strcmp(type, 'era5')
        grid.dim2.lon = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era5/rad/era5_rad_1979_2019.nc', 'longitude');
        grid.dim2.lat = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era5/rad/era5_rad_1979_2019.nc', 'latitude');
        grid.dim2.plev =  ncread('/project2/tas1/miyawaki/projects/002/data/raw/era5/temp/era5_temp_1979_2019.nc', 'level');
        grid.dim3 = grid.dim2;
        save('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat', 'grid')
    elseif strcmp(type, 'erai')
        grid.dim2.lon = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/interim_rad_2000.nc', 'longitude');
        grid.dim2.lat = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/interim_rad_2000.nc', 'latitude');
        grid.dim2.plev =  ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/temp/interim_temp_2000.nc', 'level');
        grid.dim3 = grid.dim2;
        save('/project2/tas1/miyawaki/projects/002/data/read/era-interim/grid.mat', 'lat_era', 'lon_era', 'plev_era')
    elseif strcmp(type, 'gcm')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.nc', par.model, 'tas', par.model));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.nc', par.model, 'ta', par.model));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=ncread(fullpath.dim2, 'lon');
        grid.dim3.lon=ncread(fullpath.dim3, 'lon');
        grid.dim2.lat=ncread(fullpath.dim2, 'lat');
        grid.dim3.lat=ncread(fullpath.dim3, 'lat');
        grid.dim3.plev=ncread(fullpath.dim3, 'plev');
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    end

    % save grid
end
function read_rad(type, par)
    if strcmp(type, 'erai') | strcmp(type, 'era5')
        rad_vars=par.era.vars.rad;
        for i=1:length(rad_vars)
            % dimensions are (lon x lat x time)
            % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
            rad.(rad_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), rad_vars{i});
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            rad.(rad_vars{i}) = rad.(rad_vars{i})/86400;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/rad.mat', type), 'rad', 'rad_vars');
    elseif strcmp(type, 'gcm')
        rad_vars=par.gcm.vars.rad;
        for i=1:length(par.gcm.vars.rad); var = par.gcm.vars.rad{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.ymonmean.nc', par.model, var, par.model));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            rad.(var)=ncread(fullpath, var);
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='rad.mat';
        save(sprintf('%s/%s', newdir, filename), 'rad', 'rad_vars');
    end
end
function read_pe(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        pe_vars=par.era.vars.pe;
        for i=1:length(pe_vars)
            % dimensions are (lon x lat x time)
            pe.(pe_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/pe/%s_pe_%s.ymonmean.nc', type, type, par.(type).yr_span), pe_vars{i});
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            pe.(pe_vars{i}) = pe.(pe_vars{i})/86400;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/pe.mat', type), 'pe', 'pe_vars');

    elseif strcmp(type, 'gcm')
        pe_vars=par.gcm.vars.pe;
        for i=1:length(par.gcm.vars.pe); var = par.gcm.vars.pe{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.ymonmean.nc', par.model, var, par.model));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            pe.(var)=ncread(fullpath, var);
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='pe.mat';
        save(sprintf('%s/%s', newdir, filename), 'pe', 'pe_vars');
    end
end
function read_stf(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        stf_vars=par.era.vars.stf;
        for i=1:length(stf_vars)
            % dimensions are (lon x lat x time)
            stf.(stf_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/stf/%s_stf_%s.ymonmean.nc', type, type, par.era5.yr_span), stf_vars{i});
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            stf.(stf_vars{i}) = stf.(stf_vars{i})/86400;
        end
        save('/project2/tas1/miyawaki/projects/002/data/read/era5/stf.mat', 'stf', 'stf_vars');

    elseif strcmp(type, 'gcm')
        stf_vars=par.gcm.vars.stf;
        for i=1:length(par.gcm.vars.stf); var = par.gcm.vars.stf{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.ymonmean.nc', par.model, var, par.model));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            stf.(var)=ncread(fullpath, var);
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='stf.mat';
        save(sprintf('%s/%s', newdir, filename), 'stf', 'stf_vars');
    end
end
function read_srfc(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        srfc_vars=par.era.vars.srfc;
        for i=1:length(srfc_vars)
            % dimensions are (lat x time); note that the data is already zonally averaged
            srfc.(srfc_vars{i}) = squeeze(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/srfc/%s_srfc_%s.ymonmean.nc', type, type, par.era5.yr_span), srfc_vars{i}));
        end
        save('/project2/tas1/miyawaki/projects/002/data/read/era5/srfc.mat', 'srfc', 'srfc_vars');

    elseif strcmp(type, 'gcm')
        srfc_vars=par.gcm.vars.srfc;
        for i=1:length(par.gcm.vars.srfc); var = par.gcm.vars.srfc{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.ymonmean.nc', par.model, var, par.model));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            if ~exist(fullpath)
                if strcmp(var, 'hurs')
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.nc', par.model, 'hur', par.model));
                    fullpath=sprintf('%s/%s', file.folder, file.name);
                    hur=ncread(fullpath, 'hur');
                    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/grid.mat', par.model));
                    if ~isequal(grid.dim2.lon, grid.dim3.lon); hur=interp1(grid.dim3.lon, hur, grid.dim2.lon); end; % interpolate to 2D lon if different from 3D
                    hur=permute(hur, [2 1 3 4]); % bring lat to first dim
                    if ~isequal(grid.dim2.lat, grid.dim3.lat); hur=interp1(grid.dim3.lat, hur, grid.dim2.lat); end;
                    hur=permute(hur, [3 2 1 4]); % bring plev to first dim
                    pb=CmdLineProgressBar("Calculating hurs..."); % track progress of this loop
                    for id_lon=1:length(1:length(grid.dim2.lon))
                        pb.print(id_lon, length(1:length(grid.dim3.lon)));
                        for id_lat=1:length(1:length(grid.dim2.lat))
                            for id_time=1:length(1:size(srfc.ps, 3))
                                srfc.hurs(id_lon, id_lat, id_time)=interp1(grid.dim3.plev, hur(:,id_lon,id_lat,id_time), srfc.ps(id_lat, id_lat, id_time));
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
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='srfc.mat';
        save(sprintf('%s/%s', newdir, filename), 'srfc', 'srfc_vars');
    end
end
