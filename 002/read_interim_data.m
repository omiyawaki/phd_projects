clc; clear variables; close all;

% read data net SW and LW radiation data downloaded from ERA-Interim
% first read lon and lat vectors since this is different from the Donohoe grid
interim.lon = ncread('/project2/tas1/miyawaki/projects/002/data/rad/interim_rad_2000.nc', 'longitude');
interim.lat = ncread('/project2/tas1/miyawaki/projects/002/data/rad/interim_rad_2000.nc', 'latitude');

% spanning years
yr_span = 2000:2012;
yr_text = cellstr(num2str(yr_span'))';

% radiation variables to read
rad_vars = {'ssr', 'str', 'tsr', 'ttr'};
for i=1:length(rad_vars)
    text.(rad_vars{i}) = strcat(rad_vars{i}, yr_text);
    % dimensions are (lon x lat x time)
    % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
    for j=1:length(yr_span)
        interim_raw.(text.(rad_vars{i}){j}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/rad/interim_rad_%g.nc',yr_span(j)), rad_vars{i});
        % the data is originally reported as J m^-2 per day, so
        % divide by 86400 s to get the conventional W m^-2 flux
        % over the full day
        for month = 1:12
            interim.(text.(rad_vars{i}){j})(:,:,month) = interim_raw.(text.(rad_vars{i}){j})(:,:,month)/86400;
        end
    end
end
for month=1:12
    interim.ssr(:,:,month) = nanmean( cat(3, interim.ssr2000(:,:,month), interim.ssr2001(:,:,month), interim.ssr2002(:,:,month), interim.ssr2003(:,:,month), interim.ssr2004(:,:,month), interim.ssr2005(:,:,month), interim.ssr2006(:,:,month), interim.ssr2007(:,:,month), interim.ssr2008(:,:,month), interim.ssr2009(:,:,month), interim.ssr2010(:,:,month), interim.ssr2011(:,:,month), interim.ssr2012(:,:,month)), 3);
    interim.str(:,:,month) = nanmean( cat(3, interim.str2000(:,:,month), interim.str2001(:,:,month), interim.str2002(:,:,month), interim.str2003(:,:,month), interim.str2004(:,:,month), interim.str2005(:,:,month), interim.str2006(:,:,month), interim.str2007(:,:,month), interim.str2008(:,:,month), interim.str2009(:,:,month), interim.str2010(:,:,month), interim.str2011(:,:,month), interim.str2012(:,:,month)), 3);
    interim.tsr(:,:,month) = nanmean( cat(3, interim.tsr2000(:,:,month), interim.tsr2001(:,:,month), interim.tsr2002(:,:,month), interim.tsr2003(:,:,month), interim.tsr2004(:,:,month), interim.tsr2005(:,:,month), interim.tsr2006(:,:,month), interim.tsr2007(:,:,month), interim.tsr2008(:,:,month), interim.tsr2009(:,:,month), interim.tsr2010(:,:,month), interim.tsr2011(:,:,month), interim.tsr2012(:,:,month)), 3);
    interim.ttr(:,:,month) = nanmean( cat(3, interim.ttr2000(:,:,month), interim.ttr2001(:,:,month), interim.ttr2002(:,:,month), interim.ttr2003(:,:,month), interim.ttr2004(:,:,month), interim.ttr2005(:,:,month), interim.ttr2006(:,:,month), interim.ttr2007(:,:,month), interim.ttr2008(:,:,month), interim.ttr2009(:,:,month), interim.ttr2010(:,:,month), interim.ttr2011(:,:,month), interim.ttr2012(:,:,month)), 3);
end

% turbulent flux variables to read
tf_vars = {'sshf', 'slhf'};
for i=1:length(tf_vars)
    text.(tf_vars{i}) = strcat(tf_vars{i}, yr_text);
    % dimensions are (lon x lat x time)
    % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
    for j=1:length(yr_span)
        interim_raw.(text.(tf_vars{i}){j}) = ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/tf/interim_tf_%g.nc',yr_span(j)), tf_vars{i});
        % the data is originally reported as J m^-2 per day, so
        % we add the 00-12 and 12-24 steps and then
        % divide by 86400 s to get the conventional W m^-2 flux
        % over the full day
        for month = 1:12
            interim.(text.(tf_vars{i}){j})(:,:,month) = interim_raw.(text.(tf_vars{i}){j})(:,:,month)/86400;
        end
    end
end
for month=1:12
    interim.sshf(:,:,month) = nanmean( cat(3, interim.sshf2000(:,:,month), interim.sshf2001(:,:,month), interim.sshf2002(:,:,month), interim.sshf2003(:,:,month), interim.sshf2004(:,:,month), interim.sshf2005(:,:,month), interim.sshf2006(:,:,month), interim.sshf2007(:,:,month), interim.sshf2008(:,:,month), interim.sshf2009(:,:,month), interim.sshf2010(:,:,month), interim.sshf2011(:,:,month), interim.sshf2012(:,:,month)), 3);
    interim.slhf(:,:,month) = nanmean( cat(3, interim.slhf2000(:,:,month), interim.slhf2001(:,:,month), interim.slhf2002(:,:,month), interim.slhf2003(:,:,month), interim.slhf2004(:,:,month), interim.slhf2005(:,:,month), interim.slhf2006(:,:,month), interim.slhf2007(:,:,month), interim.slhf2008(:,:,month), interim.slhf2009(:,:,month), interim.slhf2010(:,:,month), interim.slhf2011(:,:,month), interim.slhf2012(:,:,month)), 3);
end

save('/project2/tas1/miyawaki/projects/002/data/radiation_climatology.mat', 'interim');
