clc; clear variables; close all;

%% set parameters
par.yr_span = 2000:2012; % spanning years for ERA-Interim
par.yr_span_era5 = '1979_2019'; % spanning years for ERA5
par.yr_text = cellstr(num2str(par.yr_span'))';
rad_vars = {'ssr', 'str', 'tsr', 'ttr'}; % radiation variables to read
stf_vars = {'sshf', 'slhf'}; % surface turbulent flux variables to read
vars_3d = {'t'}; % 3d variables to read
vars_mpi_2d = {'rsdt', 'rsut', 'rsus', 'rsds', 'rlus', 'rlds', 'rlut', 'hfls', 'hfss', 'ps'}; % 2D MPI-ESM-LR variables to read
vars_mpi_3d = {'ta', 'hus', 'zg', 'ua', 'va'}; % 3d MPI-ESM-LR variables
startend_mpi = {'280001-280912', '281001-281912', '282001-282912', '283001-283912', '284001-284912'};
par.yr_span_mpi = 50; % number of years that I am considering in the MPI climatology
% useful constants
par.cpd = 1005.7; par.Rd = 287; par.L = 2.501e6; par.g = 9.81;

%% call functions
% read_era_grid()
% read_era_rad(rad_vars, par)
% read_era_stf(stf_vars, par)
% read_era_3d(vars_3d, par)
% read_era5_grid()
% read_era5_rad(rad_vars, par)
read_era5_stf(stf_vars, par)
% read_era5_3d(vars_3d, par)
% read_mpi_2d(vars_mpi_2d, par)
% read_mpi_3d(vars_mpi_3d, startend_mpi, par)

%% define functions
function read_era_grid()
    % read data net SW and LW radiation data downloaded from ERA-Interim
    % first read lon and lat vectors since this is different from the Donohoe grid
    lon_era = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/interim_rad_2000.nc', 'longitude');
    lat_era = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/interim_rad_2000.nc', 'latitude');
    plev_era =  ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/temp/interim_temp_2000.nc', 'level');

    % save grid
    save('/project2/tas1/miyawaki/projects/002/data/read/era-interim/grid.mat', 'lat_era', 'lon_era', 'plev_era')
end
function read_era_rad(rad_vars, par)
    for i=1:length(rad_vars)
    text.(rad_vars{i}) = strcat(rad_vars{i}, par.yr_text);
    % dimensions are (lon x lat x time)
    % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
    for j=1:length(par.yr_span)
        interim_raw.(text.(rad_vars{i}){j}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/interim_rad_%g.nc',par.yr_span(j)), rad_vars{i});
        return
        % the data is originally reported as J m^-2 per day, so
        % divide by 86400 s to get the conventional W m^-2 flux
        % over the full day
        for month = 1:12
            interim_raw.(text.(rad_vars{i}){j})(:,:,month) = interim_raw.(text.(rad_vars{i}){j})(:,:,month)/86400;
        end
    end
                end
    for month=1:12
    rad.ssr(:,:,month) = nanmean( cat(3, interim_raw.ssr2000(:,:,month), interim_raw.ssr2001(:,:,month), interim_raw.ssr2002(:,:,month), interim_raw.ssr2003(:,:,month), interim_raw.ssr2004(:,:,month), interim_raw.ssr2005(:,:,month), interim_raw.ssr2006(:,:,month), interim_raw.ssr2007(:,:,month), interim_raw.ssr2008(:,:,month), interim_raw.ssr2009(:,:,month), interim_raw.ssr2010(:,:,month), interim_raw.ssr2011(:,:,month), interim_raw.ssr2012(:,:,month)), 3);
    rad.str(:,:,month) = nanmean( cat(3, interim_raw.str2000(:,:,month), interim_raw.str2001(:,:,month), interim_raw.str2002(:,:,month), interim_raw.str2003(:,:,month), interim_raw.str2004(:,:,month), interim_raw.str2005(:,:,month), interim_raw.str2006(:,:,month), interim_raw.str2007(:,:,month), interim_raw.str2008(:,:,month), interim_raw.str2009(:,:,month), interim_raw.str2010(:,:,month), interim_raw.str2011(:,:,month), interim_raw.str2012(:,:,month)), 3);
    rad.tsr(:,:,month) = nanmean( cat(3, interim_raw.tsr2000(:,:,month), interim_raw.tsr2001(:,:,month), interim_raw.tsr2002(:,:,month), interim_raw.tsr2003(:,:,month), interim_raw.tsr2004(:,:,month), interim_raw.tsr2005(:,:,month), interim_raw.tsr2006(:,:,month), interim_raw.tsr2007(:,:,month), interim_raw.tsr2008(:,:,month), interim_raw.tsr2009(:,:,month), interim_raw.tsr2010(:,:,month), interim_raw.tsr2011(:,:,month), interim_raw.tsr2012(:,:,month)), 3);
    rad.ttr(:,:,month) = nanmean( cat(3, interim_raw.ttr2000(:,:,month), interim_raw.ttr2001(:,:,month), interim_raw.ttr2002(:,:,month), interim_raw.ttr2003(:,:,month), interim_raw.ttr2004(:,:,month), interim_raw.ttr2005(:,:,month), interim_raw.ttr2006(:,:,month), interim_raw.ttr2007(:,:,month), interim_raw.ttr2008(:,:,month), interim_raw.ttr2009(:,:,month), interim_raw.ttr2010(:,:,month), interim_raw.ttr2011(:,:,month), interim_raw.ttr2012(:,:,month)), 3);
        end

    save('/project2/tas1/miyawaki/projects/002/data/read/era-interim/radiation_climatology.mat', 'rad', 'rad_vars');
end
function read_era_stf(stf_vars, par)
    for i=1:length(stf_vars)
    text.(stf_vars{i}) = strcat(stf_vars{i}, par.yr_text);
    % dimensions are (lon x lat x time)
    % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
    for j=1:length(par.yr_span)
        interim_raw.(text.(stf_vars{i}){j}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/stf/interim_stf_%g.nc',par.yr_span(j)), stf_vars{i});
        % the data is originally reported as J m^-2 per day, so
        % we add the 00-12 and 12-24 steps and then
        % divide by 86400 s to get the conventional W m^-2 flux
        % over the full day
        for month = 1:12
            interim_raw.(text.(stf_vars{i}){j})(:,:,month) = interim_raw.(text.(stf_vars{i}){j})(:,:,month)/86400;
        end
    end
                end
    for month=1:12
    stf.sshf(:,:,month) = nanmean( cat(3, interim_raw.sshf2000(:,:,month), interim_raw.sshf2001(:,:,month), interim_raw.sshf2002(:,:,month), interim_raw.sshf2003(:,:,month), interim_raw.sshf2004(:,:,month), interim_raw.sshf2005(:,:,month), interim_raw.sshf2006(:,:,month), interim_raw.sshf2007(:,:,month), interim_raw.sshf2008(:,:,month), interim_raw.sshf2009(:,:,month), interim_raw.sshf2010(:,:,month), interim_raw.sshf2011(:,:,month), interim_raw.sshf2012(:,:,month)), 3);
    stf.slhf(:,:,month) = nanmean( cat(3, interim_raw.slhf2000(:,:,month), interim_raw.slhf2001(:,:,month), interim_raw.slhf2002(:,:,month), interim_raw.slhf2003(:,:,month), interim_raw.slhf2004(:,:,month), interim_raw.slhf2005(:,:,month), interim_raw.slhf2006(:,:,month), interim_raw.slhf2007(:,:,month), interim_raw.slhf2008(:,:,month), interim_raw.slhf2009(:,:,month), interim_raw.slhf2010(:,:,month), interim_raw.slhf2011(:,:,month), interim_raw.slhf2012(:,:,month)), 3);
    end

    save('/project2/tas1/miyawaki/projects/002/data/read/era-interim/turbfluxes_climatology.mat', 'stf', 'stf_vars');
end
function read_era_3d(vars_3d, par)
    for i=1:length(vars_3d)
        text.(vars_3d{i}) = strcat(vars_3d{i}, par.yr_text);
        % dimensions are (lon x lat x plev x time)
        for j=1:length(par.yr_span)
            interim_raw.(text.(vars_3d{i}){j}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/temp/interim_temp_%g.nc',par.yr_span(j)), vars_3d{i});
            for month = 1:12
                interim_raw.(text.(vars_3d{i}){j})(:,:,:,month) = interim_raw.(text.(vars_3d{i}){j})(:,:,:,month);
            end
        end
    end
    for month=1:12
        interim_raw.t(:,:,:,month) = nanmean( cat(4, interim_raw.t2000(:,:,:,month), interim_raw.t2001(:,:,:,month), interim_raw.t2002(:,:,:,month), interim_raw.t2003(:,:,:,month), interim_raw.t2004(:,:,:,month), interim_raw.t2005(:,:,:,month), interim_raw.t2006(:,:,:,month), interim_raw.t2007(:,:,:,month), interim_raw.t2008(:,:,:,month), interim_raw.t2009(:,:,:,month), interim_raw.t2010(:,:,:,month), interim_raw.t2011(:,:,:,month), interim_raw.t2012(:,:,:,month)), 4);
    end
    % take zonal average to save space
    vert.t = squeeze(nanmean(interim_raw.t, 1));

    save('/project2/tas1/miyawaki/projects/002/data/read/era-interim/temp_climatology.mat', 'vert', 'vars_3d');
end
function read_era5_grid()
    % read data net SW and LW radiation data downloaded from Era5
    % first read lon and lat vectors since this is different from the Donohoe grid
    lon_era = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era5/rad/era5_rad_1979_2019.nc', 'longitude');
    lat_era = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era5/rad/era5_rad_1979_2019.nc', 'latitude');
    plev_era =  ncread('/project2/tas1/miyawaki/projects/002/data/raw/era5/temp/era5_temp_1979_2019.nc', 'level');

    % save grid
    save('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat', 'lat_era', 'lon_era', 'plev_era')
end
function read_era5_rad(rad_vars, par)
    for i=1:length(rad_vars)
        % dimensions are (lon x lat x time)
        % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
        era5_raw.(rad_vars{i}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/era5/rad/era5_rad_%s.nc',par.yr_span_era5), rad_vars{i});
        % the data is originally reported as J m^-2 per day, so
        % divide by 86400 s to get the conventional W m^-2 flux
        % over the full day
        era5_raw.(rad_vars{i}) = era5_raw.(rad_vars{i})/86400;
        % calculate monthly climatology
        if ~mod(size(era5_raw.(rad_vars{i}),3), 12)
            n_years = size(era5_raw.(rad_vars{i}),3)/12;
        else
            error('Data does not end in a full year. Please make sure the data is available until the end of the year (December).');
        end
        for month = 1:12
            get_months = month + [0:12:(n_years-1)*12];
            rad.(rad_vars{i})(:,:,month) = nanmean(era5_raw.(rad_vars{i})(:,:,get_months),3);
        end
    end

    save('/project2/tas1/miyawaki/projects/002/data/read/era5/radiation_climatology.mat', 'rad', 'rad_vars');
end
function read_era5_stf(stf_vars, par)
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
function read_era5_3d(vars_3d, par)
    for i=1:length(vars_3d)
        text.(vars_3d{i}) = strcat(vars_3d{i}, par.yr_text);
        % dimensions are (lon x lat x plev x time)
        for j=1:length(par.yr_span)
            interim_raw.(text.(vars_3d{i}){j}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/era5-interim/temp/interim_temp_%g.nc',par.yr_span(j)), vars_3d{i});
            for month = 1:12
                interim_raw.(text.(vars_3d{i}){j})(:,:,:,month) = interim_raw.(text.(vars_3d{i}){j})(:,:,:,month);
            end
        end
    end
    for month=1:12
        interim_raw.t(:,:,:,month) = nanmean( cat(4, interim_raw.t2000(:,:,:,month), interim_raw.t2001(:,:,:,month), interim_raw.t2002(:,:,:,month), interim_raw.t2003(:,:,:,month), interim_raw.t2004(:,:,:,month), interim_raw.t2005(:,:,:,month), interim_raw.t2006(:,:,:,month), interim_raw.t2007(:,:,:,month), interim_raw.t2008(:,:,:,month), interim_raw.t2009(:,:,:,month), interim_raw.t2010(:,:,:,month), interim_raw.t2011(:,:,:,month), interim_raw.t2012(:,:,:,month)), 4);
    end
    % take zonal avera5ge to save space
    vert.t = squeeze(nanmean(interim_raw.t, 1));

    save('/project2/tas1/miyawaki/projects/002/data/read/temp_climatology.mat', 'vert', 'vars_3d');
end
function read_mpi_2d(vars_mpi_2d, par)
    for i=1:length(vars_mpi_2d); var = vars_mpi_2d{i};
        mpi_raw.(var) = ncread(sprintf('/project2/tas1/CMIP5_piControl/MPI-ESM-LR/%s_Amon_MPI-ESM-LR_piControl_r1i1p1_280001-284912.nc', var), var);
        for month=1:12
            index = month+[0:12:12*(par.yr_span_mpi-1)];
            mpi_2d.(var)(:,:,month) = nanmean(mpi_raw.(var)(:,:,index),3);
        end
    end

    save('/project2/tas1/miyawaki/projects/002/data/read/mpi/2d_climatology.mat', 'mpi_2d');
end
function read_mpi_3d(vars_mpi_3d, startend_mpi, par)
    mpi_3d.lon = ncread(sprintf('/project2/tas1/CMIP5_piControl/MPI-ESM-LR/%s_Amon_MPI-ESM-LR_piControl_r1i1p1_%s.nc', 'ta', startend_mpi{1}), 'lon');
    mpi_3d.lat = ncread(sprintf('/project2/tas1/CMIP5_piControl/MPI-ESM-LR/%s_Amon_MPI-ESM-LR_piControl_r1i1p1_%s.nc', 'ta', startend_mpi{1}), 'lat');
    mpi_3d.plev = ncread(sprintf('/project2/tas1/CMIP5_piControl/MPI-ESM-LR/%s_Amon_MPI-ESM-LR_piControl_r1i1p1_%s.nc', 'ta', startend_mpi{1}), 'plev');
    for i=1:length(vars_mpi_3d); var = vars_mpi_3d{i};
        for j=1:length(startend_mpi); startend = startend_mpi{j};
            mpi_raw.(var){j} = ncread(sprintf('/project2/tas1/CMIP5_piControl/MPI-ESM-LR/%s_Amon_MPI-ESM-LR_piControl_r1i1p1_%s.nc', var, startend), var);
        end
        % concatenate all years
        mpi_cat.(var) = cat(4, mpi_raw.(var){:});
    end
    % mpi_cat.pa = permute( repmat(mpi_3d.plev, [1 size(mpi_cat.ta, 1) size(mpi_cat.ta, 2) size(mpi_cat.ta, 4)]), [2 3 1 4]);
    % mpi_cat.rho = mpi_cat.pa./(par.Rd*mpi_cat.ta);
    mpi_cat.h = par.cpd*mpi_cat.ta + par.L*mpi_cat.hus + par.g*mpi_cat.zg;
    mpi_cat.vh = mpi_cat.va.*mpi_cat.h;
    mpi_cat.uh = mpi_cat.ua.*mpi_cat.h;
    % take monthly climatology
    for fn={'ta', 'h', 'vh', 'uh'};
        for month=1:12
            index = month+[0:12:12*(par.yr_span_mpi-1)];
            mpi_3d.(fn{1})(:,:,:,month) = nanmean(mpi_cat.(fn{1})(:,:,:,index),4);
        end
    end

    save('/project2/tas1/miyawaki/projects/002/data/read/mpi/3d_climatology.mat', 'mpi_3d');
end
