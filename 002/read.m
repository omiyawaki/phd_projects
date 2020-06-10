clc; clear variables; close all;

%% set parameters
% spanning years
par.yr_span = 2000:2012;
par.yr_text = cellstr(num2str(par.yr_span'))';
% radiation variables to read
rad_vars = {'ssr', 'str', 'tsr', 'ttr'};
% surface turbulent flux variables to read
stf_vars = {'sshf', 'slhf'};
% 3d variables to read
vars_3d = {'t'};

%% call functions
read_grid()
% read_rad(rad_vars, par)
% read_stf(stf_vars, par)
% read_3d(vars_3d, par)

%% define functions
function read_grid()
    % read data net SW and LW radiation data downloaded from ERA-Interim
    % first read lon and lat vectors since this is different from the Donohoe grid
    lon_era = ncread('/project2/tas1/miyawaki/projects/002/data/raw/rad/interim_rad_2000.nc', 'longitude');
    lat_era = ncread('/project2/tas1/miyawaki/projects/002/data/raw/rad/interim_rad_2000.nc', 'latitude');
    plev_era =  ncread('/project2/tas1/miyawaki/projects/002/data/raw/temp/interim_temp_2000.nc', 'level');

    % save grid
    save('/project2/tas1/miyawaki/projects/002/data/read/era_grid.mat', 'lat_era', 'lon_era', 'plev_era')
end
function read_rad(rad_vars, par)
    for i=1:length(rad_vars)
    text.(rad_vars{i}) = strcat(rad_vars{i}, par.yr_text);
    % dimensions are (lon x lat x time)
    % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
    for j=1:length(par.yr_span)
        interim_raw.(text.(rad_vars{i}){j}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/rad/interim_rad_%g.nc',par.yr_span(j)), rad_vars{i});
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

    save('/project2/tas1/miyawaki/projects/002/data/read/radiation_climatology.mat', 'rad', 'rad_vars');
end
function read_stf(stf_vars, par)
    for i=1:length(stf_vars)
    text.(stf_vars{i}) = strcat(stf_vars{i}, par.yr_text);
    % dimensions are (lon x lat x time)
    % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
    for j=1:length(par.yr_span)
        interim_raw.(text.(stf_vars{i}){j}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/stf/interim_stf_%g.nc',par.yr_span(j)), stf_vars{i});
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

    save('/project2/tas1/miyawaki/projects/002/data/read/turbfluxes_climatology.mat', 'stf', 'stf_vars');
end
function read_3d(vars_3d, par)
    for i=1:length(vars_3d)
        text.(vars_3d{i}) = strcat(vars_3d{i}, par.yr_text);
        % dimensions are (lon x lat x plev x time)
        for j=1:length(par.yr_span)
            interim_raw.(text.(vars_3d{i}){j}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/temp/interim_temp_%g.nc',par.yr_span(j)), vars_3d{i});
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

    save('/project2/tas1/miyawaki/projects/002/data/read/temp_climatology.mat', 'vert', 'vars_3d');
end
