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
% read_grid(type, par)
% read_rad(type, par)
% TODO convert all of this below:
% read_pe(type, par)
% read_era5_stf(type, par)
% read_era5_vert(type, par)
% read_era5_sfc(type, par)
for k=1:length(par.gcm_models); par.model=par.gcm_models{k};
    % read_grid('gcm', par)
    % read_rad('gcm', par)
    read_pe('gcm', par)
    % read_gcm_2d(par, model)
    % read_gcm_3d(par, model)
end

%% define functions
function read_grid(type, par)
    % read data net SW and LW radiation data downloaded from Era5
    % first read lon and lat vectors since this is different from the Donohoe grid
    if strcmp(type, 'era5')
        grid.lon = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era5/rad/era5_rad_1979_2019.nc', 'longitude');
        grid.lat = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era5/rad/era5_rad_1979_2019.nc', 'latitude');
        grid.plev =  ncread('/project2/tas1/miyawaki/projects/002/data/raw/era5/temp/era5_temp_1979_2019.nc', 'level');
        save('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat', 'grid')
    elseif strcmp(type, 'erai')
        grid.lon = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/interim_rad_2000.nc', 'longitude');
        grid.lat = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/interim_rad_2000.nc', 'latitude');
        grid.plev =  ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/temp/interim_temp_2000.nc', 'level');
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
            raw.(rad_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.nc', type, type, par.(type).yr_span), rad_vars{i});
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            raw.(rad_vars{i}) = raw.(rad_vars{i})/86400;
            % calculate monthly climatology
            if ~mod(size(raw.(rad_vars{i}),3), 12)
                n_years = size(raw.(rad_vars{i}),3)/12;
            else
                error('Data does not end in a full year. Please make sure the data is available until the end of the year (December).');
            end
            for month = 1:12
                get_months = month + [0:12:(n_years-1)*12];
                rad.(rad_vars{i})(:,:,month) = nanmean(raw.(rad_vars{i})(:,:,get_months),3);
            end
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/rad_climatology.mat', type), 'rad', 'rad_vars');
    elseif strcmp(type, 'gcm')
        rad_vars=par.gcm.vars.rad;
        for i=1:length(par.gcm.vars.rad); var = par.gcm.vars.rad{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.nc', par.model, var, par.model));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            rad.(var)=ncread(fullpath, var);
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='rad_climatology.mat';
        save(sprintf('%s/%s', newdir, filename), 'rad', 'rad_vars');
    end
end
function read_pe(type, par) % TODO finish converting this to be type agnostic
    if strcmp(type, {'era5', 'erai'})
        pe_vars=par.era.vars.pe;
        for i=1:length(pe_vars)
            % dimensions are (lon x lat x time)
            % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
            raw.(pe_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/pe/pe_%s.nc', type, par.(type).yr_span), pe_vars{i});
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            raw.(pe_vars{i}) = raw.(pe_vars{i})/86400;
            % calculate monthly climatology
            if ~mod(size(raw.(pe_vars{i}),3), 12)
                n_years = size(raw.(pe_vars{i}),3)/12;
            else
                error('Data does not end in a full year. Please make sure the data is available until the end of the year (December).');
            end
            for month = 1:12
                get_months = month + [0:12:(n_years-1)*12];
                pe.(pe_vars{i})(:,:,month) = nanmean(raw.(pe_vars{i})(:,:,get_months),3);
            end
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/pe_climatology.mat', type), 'pe', 'pe_vars');

    elseif strcmp(type, 'gcm')
        pe_vars=par.gcm.vars.pe;
        for i=1:length(par.gcm.vars.pe); var = par.gcm.vars.pe{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.nc', par.model, var, par.model));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            pe.(var)=ncread(fullpath, var);
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='pe_climatology.mat';
        save(sprintf('%s/%s', newdir, filename), 'pe', 'pe_vars');
    end
end
function read_era5_stf(stf_vars, par) % TODO
    for i=1:length(stf_vars)
        % dimensions are (lon x lat x time)
        % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
        era5_raw.(stf_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/era5/stf/era5_stf_%s.nc',par.yr_span_era5), stf_vars{i});
        % the data is originally reported as J m^-2 per day, so
        % divide by 86400 s to get the conventional W m^-2 flux
        % over the full day
        era5_raw.(stf_vars{i}) = era5_raw.(stf_vars{i})/86400;
        % calculate monthly climatology
        if ~mod(size(era5_raw.(stf_vars{i}),3), 12)
            n_years = size(era5_raw.(stf_vars{i}),3)/12;
        else
            error('Data does not end in a full year. Please make sure the data is available until the end of the year (December).');
        end
        for month = 1:12
            get_months = month + [0:12:(n_years-1)*12];
            stf.(stf_vars{i})(:,:,month) = nanmean(era5_raw.(stf_vars{i})(:,:,get_months),3);
        end
    end

    save('/project2/tas1/miyawaki/projects/002/data/read/era5/turbfluxes_climatology.mat', 'stf', 'stf_vars');
end
function read_era5_vert(vars_vert, par) % TODO
    for i=1:length(vars_vert)
        % dimensions are (lat x plev x time); note that this data is already zonally-averaged
        era5_raw.(vars_3d{i}) = squeeze(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/era5/temp/era5_tempz_%s.nc', par.yr_span_era5), vars_3d{i}));
    end
    % calculate monthly climatology
    if ~mod(size(era5_raw.(vars_3d{i}),3), 12)
        n_years = size(era5_raw.(vars_3d{i}),3)/12;
    else
        error('Data does not end in a full year. Please make sure the data is available until the end of the year (December).');
    end
    for month = 1:12
        get_months = month + [0:12:(n_years-1)*12];
        vert.(vars_3d{i})(:,:,month) = nanmean(era5_raw.(vars_3d{i})(:,:,get_months),3);
    end

    save('/project2/tas1/miyawaki/projects/002/data/read/era5/temp_climatology.mat', 'vert', 'vars_3d');
end
function read_era5_srfc(vars_srfc, par) % TODO
    for i=1:length(vars_srfc)
        % dimensions are (lat x time); note that the data is already zonally averaged
        era5_raw.(vars_srfc{i}) = squeeze(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/era5/sfc/era5_sfcz_%s.nc', par.yr_span_era5), vars_srfc{i}));
        % calculate monthly climatology
        if ~mod(size(era5_raw.(vars_srfc{i}),2), 12)
            n_years = size(era5_raw.(vars_srfc{i}),2)/12;
        else
            error('Data does not end in a full year. Please make sure the data is available until the end of the year (December).');
        end
        for month = 1:12
            get_months = month + [0:12:(n_years-1)*12];
            srfc.(vars_srfc{i})(:,month) = nanmean(era5_raw.(vars_srfc{i})(:,get_months),2);
        end
    end

    save('/project2/tas1/miyawaki/projects/002/data/read/era5/srfc_climatology.mat', 'srfc', 'vars_srfc');
end
function read_gcm_2d(par, model)
    for i=1:length(par.vars_gcm_2d); var = par.vars_gcm_2d{i};
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.nc', model, var, model));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        if ~exist(fullpath) & strcmp(var, 'hurs')
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.nc', model, 'hur', model));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            hur=ncread(fullpath, 'hur');
            load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/grid.mat', model));
            if ~isequal(grid_2d.lon, grid_3d.lon); hur=interp1(grid_3d.lon, hur, grid_2d.lon); end; % interpolate to 2D lon if different from 3D
            hur=permute(hur, [2 1 3 4]); % bring lat to first dim
            if ~isequal(grid_2d.lat, grid_3d.lat); hur=interp1(grid_3d.lat, hur, grid_2d.lat); end;
            hur=permute(hur, [3 2 1 4]); % bring plev to first dim
            pb=CmdLineProgressBar("Calculating hurs..."); % track progress of this loop
            for id_lon=1:length(1:length(grid_2d.lon))
                pb.print(id_lon, length(1:length(grid_3d.lon)));
                for id_lat=1:length(1:length(grid_2d.lat))
                    for id_time=1:length(1:size(data_2d.ps, 3))
                        data_2d.hurs(id_lon, id_lat, id_time)=interp1(grid_3d.plev, hur(:,id_lon,id_lat,id_time), data_2d.ps(id_lat, id_lat, id_time));
                    end
                end
            end
        else
        end
    end
    newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', model);
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='2d_climatology.mat';
    save(sprintf('%s/%s', newdir, filename), 'data_2d');
end
function read_gcm_3d(vars_gcm_3d, startend_gcm, par)
    gcm_3d.lon = ncread(sprintf('/project2/tas1/CMIP5_piControl/GCM-ESM-LR/%s_Amon_GCM-ESM-LR_piControl_r1i1p1_%s.nc', 'ta', startend_gcm{1}), 'lon');
    gcm_3d.lat = ncread(sprintf('/project2/tas1/CMIP5_piControl/GCM-ESM-LR/%s_Amon_GCM-ESM-LR_piControl_r1i1p1_%s.nc', 'ta', startend_gcm{1}), 'lat');
    gcm_3d.plev = ncread(sprintf('/project2/tas1/CMIP5_piControl/GCM-ESM-LR/%s_Amon_GCM-ESM-LR_piControl_r1i1p1_%s.nc', 'ta', startend_gcm{1}), 'plev');
    for i=1:length(vars_gcm_3d); var = vars_gcm_3d{i};
        for j=1:length(startend_gcm); startend = startend_gcm{j};
            gcm_raw.(var){j} = ncread(sprintf('/project2/tas1/CMIP5_piControl/GCM-ESM-LR/%s_Amon_GCM-ESM-LR_piControl_r1i1p1_%s.nc', var, startend), var);
        end
        % concatenate all years
        gcm_cat.(var) = cat(4, gcm_raw.(var){:});
    end
    % gcm_cat.pa = permute( repmat(gcm_3d.plev, [1 size(gcm_cat.ta, 1) size(gcm_cat.ta, 2) size(gcm_cat.ta, 4)]), [2 3 1 4]);
    % gcm_cat.rho = gcm_cat.pa./(par.Rd*gcm_cat.ta);
    gcm_cat.h = par.cpd*gcm_cat.ta + par.L*gcm_cat.hus + par.g*gcm_cat.zg;
    gcm_cat.vh = gcm_cat.va.*gcm_cat.h;
    gcm_cat.uh = gcm_cat.ua.*gcm_cat.h;
    % take monthly climatology
    for fn={'ta', 'h', 'vh', 'uh'};
        for month=1:12
            index = month+[0:12:12*(par.yr_span_gcm-1)];
            gcm_3d.(fn{1})(:,:,:,month) = nanmean(gcm_cat.(fn{1})(:,:,:,index),4);
        end
    end

    save('/project2/tas1/miyawaki/projects/002/data/read/gcm/3d_climatology.mat', 'gcm_3d');
end
